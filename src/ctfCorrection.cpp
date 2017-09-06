#include "ctfCorrection.h"
#include "exception.h"
#include <sstream>
#include "fftRoutines.h"
#include "volumeIO.h"
#include "defocusFileFormats.h"

CTFCorrection::CTFCorrection(ParameterSetup& aParams)
{
    params=aParams;
    projSet = new ProjectionSetIdentity();
    inputStack = new MRCStack(params.InputStackName(),false, false, true);
    projSet->init(inputStack->getNumberOfProjections());
}

CTFCorrection::~CTFCorrection()
{
    delete inputStack;
}

void CTFCorrection::run()
{
    initVariables();
    computeFrequencyArray();
    checkAstigmatism();
    correctCTF();
}

void CTFCorrection::correctCTF()
{
    std::vector<std::vector<float>> correctedProjections;
    std::vector<float> ctfAmp;
    std::vector<float> ctfFilter;

    correctedProjections.resize(inputStack->getNumberOfProjections());

    for(ProjectionSet::iterator it=projSet->begin(); it!=projSet->end(); it++)
    {
        unsigned int projIndex = it.second();
        correctedProjections[projIndex].resize(inputStack->getProjectionSize());
        inputStack->readProjections(&correctedProjections[projIndex][0],1,projIndex);

        computeCTF(ctfAmp,ctfFilter,defocusFileValues[projIndex][0],defocusFileValues[projIndex][1],defocusFileValues[projIndex][2],defocusFileValues[projIndex][3]);

        if(ctfCorrectionType=="phaseflip")
        {
            if(!correctAstigmatism)
                FFTRoutines::real2DTransform(arraySizeX,arraySizeY,correctedProjections[projIndex],ctfFilter);
            else
                FFTRoutines::complex2DTransform(arraySizeX,arraySizeY,correctedProjections[projIndex],ctfFilter);
        }
        else //multiplication
        {
            if(!correctAstigmatism)
                FFTRoutines::real2DTransform(arraySizeX,arraySizeY,correctedProjections[projIndex],ctfAmp);
            else
                FFTRoutines::complex2DTransform(arraySizeX,arraySizeY,correctedProjections[projIndex],ctfAmp);
        }
    }

    VolumeIO::writeMRCStack(inputStack->getStackHeader(),correctedProjections,params.OutputFilename(),inputStack->getExtraData());

}

/* If astigmatism correction is required it performs check whether the values are available.
 * If they are not, than the warning is written out and the correction is performed without
 * the astigmatism correction.
 */
void CTFCorrection::checkAstigmatism()
{
    if(!correctAstigmatism)
        return;

    float diff = 0.0f;

    for(unsigned int i=0; i<defocusFileValues.size(); i++)
    {
        diff= diff+(defocusFileValues[i][0]-defocusFileValues[i][1]);
    }

    if(diff<0.0001)
    {
        cout << "The values necessary for astigmatism correction are not present in the defocus file!!!" << std::endl;
        cout << "The CTF correction will be performed without the astigmatism correction!!!" << std::endl << std::endl;
        correctAstigmatism=false;
    }
}

void CTFCorrection::initVariables()
{
    defocusFileFormat = params.DefocusFileFormat();
    ctfCorrectionType = params.CtfCorrectionType();

    pixelSize = params.PixelSize();
    amplitude = params.Amplitude();
    cs=params.Cs();
    evk=params.Evk();

    //Convert spherical aberration term from mm to Angstroms
    cs = cs*(1.0e7);

    //Convert from nm to Angstroms
    pixelSize=pixelSize*10.0f;

    //readDefocusFile();
    defocusFileValues.resize(inputStack->getNumberOfProjections());
    DefocusFileFormat* defFile=DefocusFileFormat::createFileFormat(defocusFileFormat);
    defFile->read(params.DefocusFile(),*projSet);
    defFile->getValues(defocusFileValues,*projSet,"microns");

    correctAstigmatism=params.CorrectAstigmatism();
}

void CTFCorrection::computeCTF(std::vector<float>& ctfAmp, std::vector<float>& ctfFilter, float defocus1, float defocus2, float astigmatism, float phaseShift)
{
    // Convert defocii from microns to Angstroms
    defocus1 = defocus1*(1.0e4);
    defocus2 = defocus2*(1.0e4);


    // Calculate electron wavelength
    // Most of the equations here are from Mindell and Grigorieff (2003)

    double h = 6.62606957e-34;
    double c = 299792458;
    float eRest = 511000;
    float v = evk*1000;
    float eCharge = 1.602e-19;
    double lambda = (c*h)/sqrt(((2*eRest*v)+(v*v))*(eCharge*eCharge))*(pow(10.0,10));

    // Calculate astigmatic defocus array

    // Initialize defocus array
    std::vector<std::vector<float> > defocusArray;

    // Calculate center of image
    float centerX = arraySizeX/2.0f-1.0f;
    float centerY = arraySizeY/2.0f-1.0f;

    // For no astigmatism
    if(!correctAstigmatism)
    {
        initMatrixWithValue(defocusArray,defocus1);
    }
    else	// For astigmatism
    {
        // Precalculate some numbers
        float defSum=defocus1+defocus2;
        float defDiff=defocus1-defocus2;

        // Calculate the astigmatic defocus and store in array
        defocusArray.resize(arraySizeX);
        for(size_t x = 0; x < arraySizeX; x++)
        {
            defocusArray[x].resize(arraySizeY);

            for (size_t y = 0; y < arraySizeY; y++)
            {
                float angle = (-atan2((centerY-(y+1)),((x+1)-centerX)))*180/M_PI;
                defocusArray[x][y] = (defSum + defDiff*cosd(2.0f*(angle-astigmatism)))/2.0f;
            }
        }
    }

    // CTF calculation

    // Calculate weighting factors
    float w = sqrt(1.0f-amplitude*amplitude);

    // Calculate phase array

    ctfFilter.resize(arraySizeX*arraySizeY);
    std::fill(ctfFilter.begin(),ctfFilter.end(),1.0f);

    ctfAmp.resize(arraySizeX*arraySizeY);
    for(size_t x = 0; x < arraySizeX; x++)
    {
        for (size_t y = 0; y < arraySizeY; y++)
        {
            float phase = ((M_PI*lambda*(frequencyArray[x][y]*frequencyArray[x][y]))*(defocusArray[x][y] - 0.5*(lambda*lambda)*(frequencyArray[x][y]*frequencyArray[x][y])*cs))+phaseShift;
            // Calculate CTF
            // Has to be transposed here for the case of astigmatism!!!
            ctfAmp[y+x*arraySizeY]= (w*sin(phase)) + (amplitude*cos(phase));
            if(ctfAmp[y+x*arraySizeY] <= 0)
                ctfFilter[y+x*arraySizeY]=-1.0f;
        }
    }
}

void CTFCorrection::initMatrixWithValue(std::vector<std::vector<float> >& matrix, float value)
{
    matrix.resize(arraySizeX);

    for(size_t i=0; i<matrix.size(); i++)
    {
        matrix[i].resize(arraySizeY);
        std::fill(matrix[i].begin(), matrix[i].end(), value);
    }
}

void CTFCorrection::computeFrequencyArray()
{
    size_t resX=inputStack->getResolution().x;
    size_t resY=inputStack->getResolution().y;

    std::vector<float> xArray;
    std::vector<float> yArray;

    for(int i = -floor(resX/2); i<floor(resX/2)+(resX%2); i++ )
        xArray.push_back(i/(resX*pixelSize));

    for(int i = -floor(resY/2); i<floor(resY/2)+(resY%2); i++ )
        yArray.push_back(i/(resY*pixelSize));

    frequencyArray.resize(xArray.size());

    for(size_t x=0; x<xArray.size(); x++)
        for (size_t y=0; y<yArray.size(); y++)
            frequencyArray[x].push_back(sqrt(xArray[x]*xArray[x] + yArray[y]*yArray[y]));

    arraySizeX = xArray.size();
    arraySizeY = yArray.size();
}
