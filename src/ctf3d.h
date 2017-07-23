#pragma once
#include "parameterSetup.h"
#include "mrcStack.h"
#include "projectionSet.h"
#include "microscopeGeometry.h"
#include "common.h"

using namespace novaCTF;

class CTF3d
{
public:

	struct voxelProjection{
		int startingIndex;
		int endingIndex;
		double position;
		int projectionImpact;
		bool fillWithProjectionValue;
	};

	CTF3d(ParameterSetup& aParams);
	~CTF3d();

	void run();
private:

	void writeVolumeSlice();
	void backproject();
	void initVariables();
    void initAngles();

   // void setNeededSlices();
    void computeOneRow(int& volumeSliceIndex, int projectionIndex, int numberOfVoxels, double xPosition, float cosAngle, bool computeSliceStatistics, unsigned int defocusProjIndex, size_t volumeSliceOffset);
    void projectSlice();
    void computeProjectionBoundaries();

    void prepareCTF();
    void loadInputStacks();

    void setDataSize();
    void writeOutDefocusSlices(vector<vector<unsigned int> >& sliceDefocusSplits);

    string generateFilename(string originalFilename, unsigned int number);

    MRCStack* initialStack;

	std::vector<MRCStack*> inputStacks;
	ParameterSetup params;
	ProjectionSet* projSet;
	Geometry* microscopeGeometry;

	vector<vector<unsigned int> > sliceDefocusSplits;
	unsigned int numberOfStacks;
	string firstFileName;

	float globalMin;
	float globalMax;
	float globalMean;

	size_t volumeSliceSize;

	std::vector<std::vector<voxelProjection>> projectionBoundaries;

	unsigned int  slicesToProcessAtOnce;	//number of slices to process at once
	size_t processedDataSize;
	Vec3i projRes;
	unsigned int nviews;    	// number of input views, or number used
	unsigned int sliceX;     	// width of output slice
	unsigned int sliceZ;    	// thickness of output slice made by PROJECT

	vector<float> angles;		// tilt angles
	float edgeFill;				// filling value for FS
	vector<float> cbet;			// cos angles in radians;
	vector<float> sbet;			// sin angles in radians;
	float scale;
	float addition;

	vector<vector<float> > currentRows;		// stores one or more rows from all views that correspond to the current slice
	vector<float> currentVolumeSlice;		// for xz slice (or slices);
	vector<int>	  listOfIndices;			// list of indices of projections after the exlusions have been processed

	unsigned int reportIndex;
};
