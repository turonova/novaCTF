#include "ctf3d.h"
#include "common.h"
#include "volumeIO.h"
#include <math.h>
#include <algorithm>
#include <float.h>
#include <sstream>

using namespace novaCTF;
using namespace std;


CTF3d::CTF3d(ParameterSetup& aParams)
{
	params=aParams;
	projSet = new ProjectionSetIdentity();
	if(params.Use3DCTF())
		firstFileName=generateFilename(params.InputStackName(),0);
	else
		firstFileName=params.InputStackName();

	initialStack = new MRCStack(firstFileName,true,true,false);

	projSet->init(initialStack->getNumberOfProjections());
	microscopeGeometry = new Geometry(*initialStack,params.VolumeDimensions(),params.TiltAnglesFilename(),params.XAxisTilt(), params.ZShift());
}

CTF3d::~CTF3d()
{
	for(unsigned int inputStackIndex = 0; inputStackIndex<numberOfStacks; inputStackIndex++)
		delete inputStacks[inputStackIndex];


	delete microscopeGeometry;
	delete initialStack;
	delete projSet;
}

void CTF3d::run()
{
	VolumeIO::writeHeader(firstFileName,params.OutputFilename(),params.VolumeDimensions(),params.OutputMode(),novaCTF::VolumeRotation::ALONG_XZ);

	initVariables();

	prepareCTF();
	loadInputStacks();
	computeProjectionBoundaries();
	setDataSize();
	backproject();

	VolumeIO::convertVolumeValues(params.OutputFilename(), globalMin, globalMax,globalMean/(float)(params.VolumeDimensions().z*params.VolumeDimensions().y*params.VolumeDimensions().x),params.OutputMode());

}

void CTF3d::loadInputStacks()
{

	if(!params.Use3DCTF())
	{
		inputStacks.push_back(new MRCStack(firstFileName,false,false,params.KeepFilesOpen()));
		return;
	}

	for(unsigned int inputStackIndex = 0; inputStackIndex<numberOfStacks; inputStackIndex++)
	{
		inputStacks.push_back(new MRCStack(generateFilename(params.InputStackName(),inputStackIndex),false,false,params.KeepFilesOpen()));
	}
}

void CTF3d::prepareCTF()
{

	sliceDefocusSplits.resize(initialStack->getNumberOfProjections());
	size_t xzSliceSize = params.VolumeDimensions().x*params.VolumeDimensions().z;

	//In case of no CTF just use one slice with 0 to point to the only input stack used
	if(!params.Use3DCTF())
	{
		for(ProjectionSet::iterator it=projSet->begin(); it!=projSet->end(); it++)
		{
			unsigned int projIndex = it.second();
			sliceDefocusSplits[projIndex].resize(xzSliceSize);
			std::fill(sliceDefocusSplits[projIndex].begin(),sliceDefocusSplits[projIndex].end(),0);
		}
		numberOfStacks = 1;
		currentRows.resize(1);
		return;
	}

	if(params.DefocusStep()!=0.0f)
		numberOfStacks = microscopeGeometry->computeNumberOfParts(params.DefocusStep(),params.PixelSize(),params.VolumeThicknessType());
	else if(params.NumberOfInputStacks()!=0)
		numberOfStacks = params.NumberOfInputStacks();
	else
	{
		cout << "Either defocus step size (DefocusStep in nm) or number of input stacks (NumberOfInputStacks) has to be specified!" << endl;
		return;
	}

	currentRows.resize(numberOfStacks);

	for(ProjectionSet::iterator it=projSet->begin(); it!=projSet->end(); it++)
	{
		unsigned int projIndex = it.second();
		microscopeGeometry->setProjectionGeometry(projIndex);
		sliceDefocusSplits[projIndex].reserve(xzSliceSize);
		microscopeGeometry->generateFocusGrid(sliceDefocusSplits[projIndex],numberOfStacks,params.VolumeThicknessType());
	}

	writeOutDefocusSlices(sliceDefocusSplits);
}

void CTF3d::writeOutDefocusSlices(vector<vector<unsigned int> >& sliceDefocusSplits)
{
	if(!params.WriteOutDefocusSlices())
		return;

	vector<float> sliceData;

	for(ProjectionSet::iterator it=projSet->begin(); it!=projSet->end(); it++)
	{
		unsigned int projIndex = it.second();
		for(size_t i = 0; i < sliceDefocusSplits[projIndex].size(); i++)
			sliceData.push_back(sliceDefocusSplits[projIndex][i]);
	}

	VolumeIO::write(sliceData,Vec3ui(params.VolumeDimensions().x,params.VolumeDimensions().z,initialStack->getNumberOfProjections()),params.OutputFilename()+"_defocusSlices.mrc",0,novaCTF::VolumeRotation::ALONG_XY);

}

void CTF3d::initVariables()
{
	projRes.x = (int) initialStack->getResolution().x;
	projRes.y = (int) initialStack->getResolution().y;
	projRes.z = (int) initialStack->getResolution().z;

	nviews=projSet->getSize();

	sliceZ = params.VolumeDimensions().z;
	sliceX = params.VolumeDimensions().x;

	volumeSliceSize = sliceX*sliceZ;

	addition=params.ScalingOffset()*nviews;
	scale=params.Scaling()/nviews;

	initAngles();

	edgeFill = initialStack->getStackMean();
	cout << "Using computed mean value: " << edgeFill << endl;

    globalMin = FLT_MAX;
    globalMax = -FLT_MAX;
    globalMean= 0.0;

}

void CTF3d::initAngles()
{
	// Set up trigonometric tables, then convert angles to radians
	cbet.resize(nviews);
	sbet.resize(nviews);
	angles.resize(nviews);

	for(ProjectionSet::iterator it=projSet->begin(); it!=projSet->end(); it++)
		listOfIndices.push_back(it.second());


	float degreesToRadiansFactor = M_PI/180.0;

	for(unsigned int i=0; i < nviews; i++)
	{
		angles[i] = microscopeGeometry->getAngleInDegrees(listOfIndices[i]);
		float thetanv=angles[i]+params.Offset().x;
		if(thetanv > 180.0)
		thetanv=thetanv-360.0;

		if(thetanv<=-180.0)
		thetanv=thetanv+360.0;

		cbet[i]=cos(thetanv*degreesToRadiansFactor);

		//Keep cosine from going to zero so it can be divided by
		if (fabs(cbet[i]) < 1.e-6)
			cbet[i] = sign(cbet[i])*1.e-6;

		sbet[i]=-sin(thetanv*degreesToRadiansFactor);
		angles[i]=-degreesToRadiansFactor*(angles[i]+params.Offset().x);
	}

}

/* Precomputes projection boundaries for each row in xz slice
 * Assuming zero x-tilt, these values are same for all slices
 */
void CTF3d::computeProjectionBoundaries()
{
	projectionBoundaries.resize(nviews);

	float delxx=params.Offset().y;		// Offset of tilt axis from center of input images. Default is no offset or rotation
	float yoffset = params.ZShift().y; 	// Vertical shift of data in an output slice. In reality z shift

	float xcenin = projRes.x / 2.0 + 0.5 - params.SubsetStart().x; 	//center coordinate of input slice

	float xoffAdj = params.ZShift().x - (projRes.x / 2 + params.SubsetStart().x - projRes.x / 2);
	float xcen = sliceX / 2 + 0.5 + delxx + xoffAdj;	// center x coordinate of output slice
	float ycen = sliceZ / 2 + 0.5 + yoffset;			// center y coordinate of output slice

	for(unsigned int iv=0; iv<nviews; iv++)
	{
		projectionBoundaries[iv].resize(sliceZ);

		// Set view angle
		float cosBeta=cbet[iv];
		float sinBeta=sbet[iv];

		for(unsigned int i=0; i<sliceZ; i++)
		{
			float zz=(i+1-ycen);
	        float zpart=zz*sinBeta+xcenin+delxx;

	        float x = cosBeta;
	        if (fabs(cosBeta) <0.001)
	        	x = 0.001*sign(cosBeta);

	        float projStart=(1.0-zpart) / x + xcen;
	        float projEnd=(projRes.x-zpart) / x + xcen;

	        if (projEnd < projStart)
	        	novaCTF::swap(projStart,projEnd);

	        int projStartPixelIndex=projStart;

	        if(projStartPixelIndex < projStart)
	        	projStartPixelIndex=projStartPixelIndex+1;

	        projStartPixelIndex=max(projStartPixelIndex-1,0);

	        int projEndPixelIndex=projEnd;

	        if(projEndPixelIndex==projEnd)
	        	projEndPixelIndex=projEndPixelIndex-1;

	        projEndPixelIndex=min(projEndPixelIndex-1,(int)sliceX-1);

	        if (projStartPixelIndex <= projEndPixelIndex)
	        {
	            x=projStartPixelIndex-xcen+1;
	            projectionBoundaries[iv][i].position=zpart+x*cosBeta-1;
	            projectionBoundaries[iv][i].projectionImpact=projEndPixelIndex-projStartPixelIndex+1;
	            projectionBoundaries[iv][i].fillWithProjectionValue=true;

	        }
	        else
	        {
	        	projEndPixelIndex = 0;
	        	projectionBoundaries[iv][i].fillWithProjectionValue=false;
			}

	        projectionBoundaries[iv][i].startingIndex=projStartPixelIndex;
	        projectionBoundaries[iv][i].endingIndex=projEndPixelIndex;

		}
	}
}

/* Actual projection of one or more slices
 * Goes over all slices in the batch and over all projections.
 * For the currently processed slice it takes every row and based on precomputed projection boundaries
 * adds to the voxels values from input projections
 */
void CTF3d::projectSlice()
{
	int projectionIndex = 0;
	for(unsigned int sliceID=0; sliceID<slicesToProcessAtOnce; sliceID++)
	{
		size_t volumeSliceOffset=sliceID*volumeSliceSize;

		for(unsigned int iv=0; iv<nviews; iv++)
		{
			int volumeIndex = volumeSliceSize*sliceID;
			bool computeSliceStatistics = false;
			if(iv==nviews-1)
				computeSliceStatistics=true;

			for(unsigned int i=0; i<sliceZ; i++)
			{
				if (projectionBoundaries[iv][i].fillWithProjectionValue)
				{
					int startIndex = volumeIndex;
					int endIndex = volumeIndex + projectionBoundaries[iv][i].startingIndex;

					for(int ind = startIndex; ind < endIndex; ind++)
					{
						currentVolumeSlice[ind] =  currentVolumeSlice[ind] + edgeFill;
						volumeIndex++;
						if(computeSliceStatistics)
						{
							currentVolumeSlice[ind]=currentVolumeSlice[ind]*scale+addition;
							globalMax=max(globalMax,currentVolumeSlice[ind]);
							globalMin=min(globalMin,currentVolumeSlice[ind]);
							globalMean+=currentVolumeSlice[ind];
						}
					}

					computeOneRow(volumeIndex,projectionIndex,projectionBoundaries[iv][i].projectionImpact, projectionBoundaries[iv][i].position,cbet[iv],computeSliceStatistics,iv,volumeSliceOffset);
				}

				int startIndex2 = volumeIndex;
				int endIndex2 = volumeIndex+sliceX-projectionBoundaries[iv][i].endingIndex-1;
				for(int ind = startIndex2; ind<endIndex2; ind++)
				{
					currentVolumeSlice[ind] =  currentVolumeSlice[ind] + edgeFill;
					volumeIndex++;

					if(computeSliceStatistics)
					{
						currentVolumeSlice[ind]=currentVolumeSlice[ind]*scale+addition;
						globalMax=max(globalMax,currentVolumeSlice[ind]);
						globalMin=min(globalMin,currentVolumeSlice[ind]);
						globalMean+=currentVolumeSlice[ind];
					}
				}

			}
			projectionIndex=projectionIndex+projRes.x;
		}
	}
}

/* Computes values for all voxels affected by the currently processed projection
 * It uses simple linear interpolation between two adjacent projection pixel values
 */
void CTF3d::computeOneRow(int& volumeSliceIndex, int projectionIndex, int numberOfVoxels, double xPosition, float cosAngle, bool computeSliceStatistics, unsigned int defocusProjIndex, size_t volumeSliceOffset)
{
	for(int j=0; j<numberOfVoxels; j++)
	{
	  int pixelIndex=xPosition;
	  float fraction=xPosition-pixelIndex;
	  unsigned int defocusIndex = sliceDefocusSplits[defocusProjIndex][volumeSliceIndex-volumeSliceOffset];

	  currentVolumeSlice[volumeSliceIndex]=currentVolumeSlice[volumeSliceIndex]+(1.0f-fraction)*currentRows[defocusIndex][projectionIndex+pixelIndex]+fraction*currentRows[defocusIndex][projectionIndex+pixelIndex+1];

	  if(computeSliceStatistics)
	  {
		  currentVolumeSlice[volumeSliceIndex]=currentVolumeSlice[volumeSliceIndex]*scale+addition;
		  globalMax=max(globalMax,currentVolumeSlice[volumeSliceIndex]);
	  	  globalMin=min(globalMin,currentVolumeSlice[volumeSliceIndex]);
	  	  globalMean+=currentVolumeSlice[volumeSliceIndex];
	  }

	  volumeSliceIndex=volumeSliceIndex+1;
	  xPosition=xPosition+cosAngle;
	}
}

void CTF3d::setDataSize()
{

	float safeBuffer=50.0f; //in MB - just to be safe, the actual program needs around 20MB only (libraries etc.)

	size_t parameterSetupSize=sizeof(ParameterSetup);
	size_t mrcStackSize=sizeof(MRCStack)*numberOfStacks;
	size_t ctfClassSize=sizeof(CTF3d);
	size_t defocusGridSize=sizeof(unsigned int)*volumeSliceSize*nviews;
	size_t boundariesSize=sizeof(projectionBoundaries)+sizeof(voxelProjection)*projectionBoundaries.capacity();

	float baseSum=parameterSetupSize+mrcStackSize+ctfClassSize+defocusGridSize+boundariesSize;
	baseSum=baseSum/1000000.f + safeBuffer;

	float memoryLimit=params.MemoryLimit() - baseSum;

	if(memoryLimit<=0)
	{
		memoryLimit=0;
		reportIndex = 200;
		slicesToProcessAtOnce=1;
		processedDataSize=volumeSliceSize;

		if(params.MemoryLimit()!=0)
			cout << "Memory limit was too low!!!" << std::endl;

		cout << "The number of slices to be processed was set to 1." << std::endl << std::endl;

		return;
	}

	size_t doubleVolumeSliceSize=2*sizeof(float)*volumeSliceSize;				//two times because of writing in load
	size_t loadingCurrentRowsSize=projRes.x*nviews*sizeof(float);
	size_t currentRowsSize=projRes.x*nviews*numberOfStacks*sizeof(float);

	float sumForOneSlice=(doubleVolumeSliceSize+loadingCurrentRowsSize+currentRowsSize)/1000000.f;

	slicesToProcessAtOnce=floor(memoryLimit/sumForOneSlice);

	//To avoid more complications make sure Y dimension can be divided by the number of slices to be processed
	while((projRes.y%slicesToProcessAtOnce)!=0)
	{
		slicesToProcessAtOnce--;
	}

	reportIndex = slicesToProcessAtOnce;
	while(reportIndex<200)
	{
		reportIndex = 2*reportIndex;
	}

	processedDataSize=volumeSliceSize*slicesToProcessAtOnce;
	cout << "Number of slices processed at once was set to:" << slicesToProcessAtOnce << std::endl;

}

// Main loop over slices perpendicular to tilt axis
void CTF3d::backproject()
{
	currentVolumeSlice.resize(processedDataSize);
	currentRows.resize(numberOfStacks);

	for(unsigned int inputStackIndex = 0; inputStackIndex<numberOfStacks; inputStackIndex++)
		currentRows[inputStackIndex].resize(projRes.x*nviews*slicesToProcessAtOnce);

	for(int sliceNumber=0; sliceNumber<projRes.y; sliceNumber+=slicesToProcessAtOnce)
	{
		if(sliceNumber%reportIndex==0)
		{
			cout << "Started processing slice number " << (size_t) sliceNumber << "." << endl;
		}

		for(unsigned int inputStackIndex = 0; inputStackIndex<numberOfStacks; inputStackIndex++)
		{
			inputStacks[inputStackIndex]->readProjections(&currentRows[inputStackIndex][0],slicesToProcessAtOnce,sliceNumber);
		}

		// clear out the slice
		fill(currentVolumeSlice.begin(), currentVolumeSlice.end(),0.0f);
		projectSlice();

		VolumeIO::writeVolumeSliceInFloat(params.OutputFilename(),currentVolumeSlice,processedDataSize,params.OutputMode());
	}
}

string CTF3d::generateFilename(string originalFilename, unsigned int number)
{
	stringstream newFilename;
	newFilename << originalFilename << "_" << number;
	string ret = newFilename.str();
	return ret;
}
