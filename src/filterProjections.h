#pragma once
#include "parameterSetup.h"
#include "mrcStack.h"
#include "projectionSet.h"
#include "microscopeGeometry.h"

using namespace novaCTF;

class FilterProjections
{
public:
	FilterProjections(ParameterSetup& aParams);
	~FilterProjections();

	void run();

private:
	void radialWeighting(float cutOff, float fallOff);
	void transform();
	void filter();
	void initParameters();
	void taperEndToStart(unsigned int projNumber);
	void loadTiltAngles();

	MRCStack* inputStack;
	ParameterSetup params;


	Vec3i projRes;
	unsigned int numberOfAngles;    	// number of input views, or number used

	int npad;
	vector<float> angles;				// B: tilt angles
	int paddedResX;
	float scale;

	vector<float> filteredProjections;
	vector<float> currentRows;				// for now stores one row from all views that correspond to the current slice
	vector<float> radialFilter;				// for radial filter storage;
	vector<int>	  listOfIndices;			// list of indices of projections after the exclusion have been processed
};
