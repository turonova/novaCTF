#pragma once

#include "parameterSetup.h"
#include "projectionSet.h"
#include "microscopeGeometry.h"

using namespace std;
using namespace novaCTF;


class CTFCorrection
{
public:

	CTFCorrection(ParameterSetup& aParams);
	~CTFCorrection();

	void run();

private:

	void initVariables();
	void computeFrequencyArray();
	void correctCTF();
	void initMatrixWithValue(std::vector<std::vector<float> >& matrix, float value);
	void computeCTF(std::vector<float>& ctfAmp, std::vector<float>& ctfFilter, float defocus1, float defocus2, float astigmatism, float phaseShift);
	void checkAstigmatism();

	ParameterSetup params;
	MRCStack* inputStack;
	ProjectionSet* projSet;

	float pixelSize;
	float amplitude;
	float cs;
	float evk;

	size_t arraySizeX;
	size_t arraySizeY;

	string ctfCorrectionType;
	string defocusFileFormat;

	bool correctAstigmatism;

	std::vector<std::vector<float>> defocusFileValues;

	std::vector<std::vector<float> > frequencyArray;

};
