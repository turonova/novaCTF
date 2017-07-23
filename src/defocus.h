#pragma once

#include "parameterSetup.h"
#include "projectionSet.h"
#include "microscopeGeometry.h"
#include "defocusFileFormats.h"

using namespace std;
using namespace novaCTF;


class Defocus
{
public:

	Defocus(ParameterSetup& aParams);
	~Defocus();

	void run();

private:

	MRCStack* inputStack;
	ParameterSetup params;
	ProjectionSet* projSet;
	Geometry* microscopeGeometry;

	void calculateDefocus();
	string generateFilename(string originalFilename, unsigned int number);
	void readDefocusFile();
	void writeDefocusFiles();
	void computeCompleteAdditionalShift();

	Vec3ui sliceDimensions;
	unsigned int numberOfInputStacks;
	unsigned int numberOfFileValues;
	string defocusFileFormat;

	vector<vector<float>> defocusFileValues;
	vector<vector<vector<float>>> newDefocusValues;

	float completeDefocusZShift;
	DefocusFileFormat* defFile;
};
