#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "parameterSetup.h"
#include <algorithm>

ParameterSetup::ParameterSetup()
{

}

void separateNameAndValue(std::string& aLine, std::string& aParamName, std::string& aParamValue)
{
    for(unsigned int i = 0; i < aLine.length(); i++)
    {
        if(isspace(aLine[i]))
        {
            //do not want the space neither in name or value
            aParamName = aLine.substr(0, i);
            aParamValue = aLine.substr(i+1, aLine.length());
            return;
        }
    }
    aParamName = aLine;
    aParamValue.clear();
}

std::string removeEmptySpacesFromString(string input)
{
    string output=input;
    output.erase(std::remove(output.begin(), output.end(), ' '), output.end());
    return output;
}

void parseDimensions2D(std::string& aParamValue,int& aX, int& aY, const char separator)
{
    for(unsigned int i = 0; i < aParamValue.length(); i++)
    {
        if(aParamValue[i]==separator)
        {
            //do not want the space neither in name or value
            string tmpX = aParamValue.substr(0, i);
            string tmpY = aParamValue.substr(i+1, aParamValue.length());
            aX = atoi(tmpX.c_str());
            aY = atoi(tmpY.c_str());
            return;
        }
    }
}

void parseFloats(std::string& aParamValue,float& aX, float& aY, const char separator)
{
    for(unsigned int i = 0; i < aParamValue.length(); i++)
    {
        if(aParamValue[i]==separator)
        {
            //do not want the space neither in name or value
            string tmpX = aParamValue.substr(0, i);
            string tmpY = aParamValue.substr(i+1, aParamValue.length());
            aX = atof(tmpX.c_str());
            aY = atof(tmpY.c_str());
            return;
        }
    }
}

void ParameterSetup::parseComFile(const string tiltcom)
{

    std::string line;
    ifstream ifs(tiltcom.c_str());

    if (!ifs)
    {
        cout << "Could not open following file: " << tiltcom.c_str() << endl;
        return;
    }

    while(getline(ifs,line))
    {
        if(line[0]!='#' && line[0]!='$')
        {
            //cout << "Processing following line:" << line << endl;

            string paramName;
            string paramValue;
            separateNameAndValue(line, paramName, paramValue);
            storeValues(paramName,paramValue,' ');
        }
    }
}

bool ParameterSetup::parameterWithoutValue(string paramName)
{
    bool paramWithoutValue = false;

    if(paramName == "-PERPENDICULAR" || paramName == "-AdjustOrigin")
    {
        paramWithoutValue = true;
    }

    return paramWithoutValue;
}

void ParameterSetup::parseCommandLine(std::vector<string> argList)
{
    for(std::vector<string>::iterator it=argList.begin(); it!=argList.end();it++)
    {
        string paramName=*it;
        string paramValue;
        std::vector<string>::iterator nextArgument = it+1;
        if(nextArgument!=argList.end() && !parameterWithoutValue(paramName))
        {
            paramValue=*nextArgument;
            it++;
        }
        else if(!parameterWithoutValue(paramName))
        {
            cout << "Missing value of the following parameter: " << *it << endl;
            return;
        }

        //if(paramName!="tilt.com")
        {
            paramName.erase(0,1); //erase "-" at the beginning of command line arguments
            storeValues(paramName,paramValue,',');
        }
    }
}

void ParameterSetup::storeValues(string paramName, string paramValue, char separator)
{
    if(paramName == "IMAGEBINNED")
    {
        binning = atoi(paramValue.c_str());
    }
    else if(paramName == "OutputFile")
    {
        outputFilename = removeEmptySpacesFromString(paramValue);
    }
    else if(paramName == "InputProjections")
    {
        inputStackName = removeEmptySpacesFromString(paramValue);
    }
    else if(paramName == "TILTFILE")
    {
        tiltAnglesFilename = removeEmptySpacesFromString(paramValue);
    }
    else if(paramName == "THICKNESS")
    {
        depth = atoi(paramValue.c_str());
    }
    else if(paramName == "XAXISTILT")
    {
        xAxisTilt = atof(paramValue.c_str());
        if(xAxisTilt!=0.0)
        {
            cout << "WARNING: novaCTF currently does not support non-zero x-tilt!!!" << endl;
            cout << "\t\t The x-tilt will be set to zero!!!" << endl;
            xAxisTilt = 0.0f;
        }
    }
    else if(paramName == "FULLIMAGE")
    {
        //width and height before binning. Will be divided at the end of parseComFile function.
        int tempW;
        int tempH;
        parseDimensions2D(paramValue, tempW, tempH, separator);
        width=(unsigned int) tempW;
        height=(unsigned int) tempH;
    }
    else if((paramName == "EXCLUDELIST") || (paramName == "EXCLUDELIST2"))
    {
        if(!skipProjectionsList.empty())
        {
            skipProjectionsList+=",";
        }
        skipProjectionsList +=paramValue;
    }
    else if(paramName == "SUBSETSTART")
    {
        parseDimensions2D(paramValue,subsetWidth, subsetHeight ,separator);
    }
    else if (paramName == "SHIFT")
    {
        parseDimensions2D(paramValue,zShiftX,zShiftY,separator);
    }
    else if(paramName == "MODE")
    {
        outputMode = atoi(paramValue.c_str());
    }
    else if (paramName == "RADIAL")
    {
        useRadialFilter = true;
        parseFloats(paramValue,radialCutOff,radialFallOff,separator);
    }
    else if (paramName == "SCALE")
    {
        parseFloats(paramValue,scalingOffset,scaling,separator);
        if(scaling!=1.0 || scalingOffset!=0.0)
            useScaling = true;
    }
    else if (paramName == "OFFSET")
    {
        if(paramValue.find(separator)!=string::npos)
        {
            parseFloats(paramValue,offset.x,offset.y,separator);
        }
        else
        {
            offset.x = atof(paramValue.c_str());
            offset.y = 0.0f;
        }
    }
    else if(paramName == "WIDTH")
    {
        width = atoi(paramValue.c_str());
    }
    else if(paramName == "LOG")
    {
        logarithmizeData = true;
        cout << "WARNING: novaCTF currently does not support logarithmic scaling of densities!!!" << endl;
        cout << "\t\t If you wish to apply this scaling prior the reconstruction use clip function from IMOD (should be applied before the filtering)!" << endl;
    }
    else if(paramName == "Algorithm")
    {
        algorithmType =  removeEmptySpacesFromString(paramValue.c_str());
    }
    else if(paramName == "NumberOfInputStacks")
    {
        numberOfInputStacks =  atoi(paramValue.c_str());
    }
    else if(paramName == "CorrectAstigmatism")
    {
        correctAstigmatism = atoi(paramValue.c_str());
    }
    else if(paramName == "WriteOutDefocusSlices")
    {
        writeOutDefocusSlices = atoi(paramValue.c_str());
    }
    else if(paramName == "DefocusFile")
    {
        defocusFile =  removeEmptySpacesFromString(paramValue.c_str());
    }
    else if(paramName == "param")
    {
        // do nothing, it was already processed
    }
    else if(paramName == "DefocusFileFormat")
    {
        defocusFileFormat =  removeEmptySpacesFromString(paramValue.c_str());
        if(defocusFileFormat!="ctffind4" && defocusFileFormat!="imod")
            cout << "Unknown type of defocus file format: " << defocusFileFormat << endl;
    }
    else if(paramName == "PixelSize")
    {
        pixelSize = atof(paramValue.c_str());
    }
    else if(paramName=="Use3DCTF")
    {
        use3DCTF=atoi(paramValue.c_str());
    }
    else if(paramName == "VolumeThicknessType")
    {
        string thicknessType = removeEmptySpacesFromString(paramValue.c_str());
        if(thicknessType=="maximal")
            volumeThicknessType = novaCTF::VolumeThickness::MAXIMAL;
        else if(thicknessType=="minimal")
            volumeThicknessType = novaCTF::VolumeThickness::MINIMAL;
        else
            cout << "Unknown type of volume thickness: " << thicknessType << endl;
    }
    else if(paramName == "FilterProjections")
    {
        filterProjections=atoi(paramValue.c_str());
    }
    else if(paramName == "DefocusShiftFile")
    {
        defocusShiftFile =  removeEmptySpacesFromString(paramValue.c_str());
        useAdditionalShift = true;
    }
    else if(paramName == "MemoryLimit")
    {
        memoryLimit=atof(paramValue.c_str());
    }
    else if(paramName == "AmplitudeContrast")
    {
        amplitude=atof(paramValue.c_str());
    }
    else if(paramName == "Cs")
    {
        cs=atof(paramValue.c_str());
    }
    else if(paramName == "Volt")
    {
        evk=atof(paramValue.c_str());
    }
    else if(paramName == "CorrectionType")
    {
        ctfCorrectionType=removeEmptySpacesFromString(paramValue.c_str());
        if(ctfCorrectionType!="phaseflip" && ctfCorrectionType!="multiplication")
            cout << "Unknown type of CTF correction: " << ctfCorrectionType << endl;
    }
    else if(paramName == "KeepFilesOpen")
    {
        keepFilesOpen =  atoi(paramValue.c_str());
    }
    else if(paramName == "StackOrientation")
    {
        string orientationType = removeEmptySpacesFromString(paramValue.c_str());
        if(orientationType=="xy")
            stackOrientation =  novaCTF::VolumeRotation::ALONG_XY;
        else if(orientationType=="xz")
            stackOrientation =  novaCTF::VolumeRotation::ALONG_XZ;
        else
            cout << "Unknown type of stack orientation: " << orientationType << endl;
    }
    else if(paramName == "DefocusStep")
    {
        defocusStep = atof(paramValue.c_str());
    }
    else if(paramName == "FakeSIRTiterations")
    {
        fakeSirtIterations = atoi(paramValue.c_str());
        useFakeSirtIterations=true;
    }
    else
    {
        if(paramName!="")
            cout << "Ignoring following parameter: " << paramName << endl;
    }
}

void ParameterSetup::initVariables()
{
    width = 0;
    height = 0;
    depth = 0;          //output volume dimensions AFTER binning
    binning = 1;

    subsetWidth = 0;    //can be also negative according to eTomo specification
    subsetHeight = 0;   //can be also negative according to eTomo specification
    zShiftX = 0;
    zShiftY = 0;

    offset.x = 0.0f;
    offset.y = 0.0f;

    logarithmizeData = false;

    algorithmType = "none";

    defocusFileFormat = "imod";
    writeOutDefocusSlices = 0;
    volumeThicknessType = novaCTF::VolumeThickness::MAXIMAL;
    correctAstigmatism = false;

    useRadialFilter = false;
    radialCutOff = 0.0f;
    radialFallOff = 0.0f;

    use3DCTF = true;

    filterProjections = true;
    useAdditionalShift = false;
    defocusStep = 0.0f;
    numberOfInputStacks = 0;

    scalingOffset = 0.0f;
    scaling = 1.0;

    stackOrientation = novaCTF::VolumeRotation::ALONG_XY;

    keepFilesOpen = true;

    ctfCorrectionType = "phaseflip";
    memoryLimit=0;

    xAxisTilt = 0.0f;
    outputMode = 2;

    useFakeSirtIterations= false;

}

ParameterSetup::ParameterSetup(std::vector<string> argList)
{
    initVariables();

    //check if tilt.com is among the parameters and if so parse it first - the command line parameters have more priority and can later on change any parameters set by tilt.com
    for(std::vector<string>::iterator it=argList.begin(); it!=argList.end();it++)
    {
        if((*it)=="-param")
        {
            std::vector<string>::iterator nextArgument = it+1;
            cout << "Parsing the following parameter file: " << *nextArgument << endl;
            parseComFile((*nextArgument));
        }
    }

    parseCommandLine(argList);

    if(!skipProjectionsList.empty())
    {
        skipProjections = skipProjectionsList;
    }

    if(binning!=0)
    {
        width = width / binning;
        height = height / binning;
        depth = depth / binning;

        subsetWidth = subsetWidth / binning;		// usually zero and not cleared at all how it is computed within imod (not really set via etomo)
        subsetHeight = subsetHeight / binning;

        zShift.x = (float)(zShiftX/binning);
        zShift.y = (float)(zShiftY/binning);		//actually shift in Z direction
    }

    volumeDimensions = novaCTF::Vec3ui(width,height,depth);
    subsetStart=novaCTF::Vec2i(subsetWidth,subsetHeight);
}

ParameterSetup::~ParameterSetup()
{}

novaCTF::Vec3ui ParameterSetup::VolumeDimensions()
{
    return volumeDimensions;
}

string ParameterSetup::InputStackName()
{
    return inputStackName;
}

string ParameterSetup::TiltAnglesFilename()
{
    return tiltAnglesFilename;
}

string ParameterSetup::OutputFilename()
{
    return outputFilename;
}

float ParameterSetup::XAxisTilt()
{
    return xAxisTilt;
}

novaCTF::Vec2f ParameterSetup::Offset()
{
    return offset;
}

string ParameterSetup::SkipProjections()
{
    return skipProjections;
}

novaCTF::Vec2f ParameterSetup::ZShift()
{
    return zShift;
}

novaCTF::Vec2i ParameterSetup::SubsetStart()
{
    return subsetStart;
}

unsigned int ParameterSetup::OutputMode()
{
    return outputMode;
}

bool ParameterSetup::UseRadialFilter()
{
    return useRadialFilter;
}

bool ParameterSetup::LogarithmizeData()
{
    return logarithmizeData;
}

novaCTF::VolumeThickness ParameterSetup::VolumeThicknessType()
{
    return volumeThicknessType;
}

bool ParameterSetup::UseScaling()
{
    return useScaling;
}

float ParameterSetup::RadialCutOff()
{
    return radialCutOff;
}

float ParameterSetup::RadialFallOff()
{
    return radialFallOff;
}

float ParameterSetup::Scaling()
{
    return scaling;
}

float ParameterSetup::PixelSize()
{
    return pixelSize;
}

string ParameterSetup::Algorithm()
{
    return algorithmType;
}

float ParameterSetup::ScalingOffset()
{
    return scalingOffset;
}


bool ParameterSetup::WriteOutDefocusSlices()
{
    return writeOutDefocusSlices;
}


unsigned int ParameterSetup::NumberOfInputStacks()
{
    return numberOfInputStacks;
}

string ParameterSetup::DefocusFile()
{
    return defocusFile;
}

string ParameterSetup::DefocusFileFormat()
{
    return defocusFileFormat;
}


bool ParameterSetup::Use3DCTF()
{
    return use3DCTF;
}


bool ParameterSetup::FilterProjections()
{
    return filterProjections;
}
string ParameterSetup::DefocusShiftFile()
{
    return defocusShiftFile;
}
float ParameterSetup::DefocusStep()
{
    return defocusStep;
}

bool ParameterSetup::UseAdditionalShift()
{
    return useAdditionalShift;
}

novaCTF::VolumeRotation ParameterSetup::StackOrientation()
{
    return stackOrientation;
}

bool ParameterSetup::KeepFilesOpen()
{
    return keepFilesOpen;
}

int ParameterSetup::Binning()
{
    return binning;
}

float ParameterSetup::Amplitude()
{
    return amplitude;
}

float ParameterSetup::Cs()
{
    return cs;
}

float ParameterSetup::Evk()
{
    return evk;
}

std::string ParameterSetup::CtfCorrectionType()
{
    return ctfCorrectionType;
}

float ParameterSetup::MemoryLimit()
{
    return memoryLimit;
}

bool ParameterSetup::CorrectAstigmatism()
{
    return correctAstigmatism;
}

bool ParameterSetup::UseFakeSirtIterations()
{
    return useFakeSirtIterations;
}
int ParameterSetup::FakeSirtIterations()
{
    return fakeSirtIterations;
}
