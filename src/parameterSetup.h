#pragma once
#include <vector>
#include "common.h"
#include <string>
using namespace std;

class ParameterSetup
{

public:
	ParameterSetup(std::vector<string> argList);
	ParameterSetup();

	~ParameterSetup();

	novaCTF::Vec3ui VolumeDimensions();
	string InputStackName();
	string TiltAnglesFilename();
	string OutputFilename();
	float XAxisTilt();
	string SkipProjections();
	novaCTF::Vec2f ZShift();
	novaCTF::Vec2i SubsetStart();
	unsigned int OutputMode();

	novaCTF::Vec2f Offset();

	bool UseRadialFilter();
	bool LogarithmizeData();
	bool UseScaling();

	float RadialCutOff();
	float RadialFallOff();

	float Scaling();
	float ScalingOffset();

	string Algorithm();
	unsigned int NumberOfInputStacks();
	novaCTF::VolumeThickness VolumeThicknessType();
	bool WriteOutDefocusSlices();

	string DefocusFile();
	float PixelSize();

	float MaskThreshold();

	string DefocusFileFormat();

	bool Use3DCTF();

	bool FilterProjections();
	string DefocusShiftFile();
	float DefocusStep();

	bool UseAdditionalShift();

	novaCTF::VolumeRotation StackOrientation();

	bool KeepFilesOpen();
	int Binning();

	float Amplitude();
	float Cs();
	float Evk();

	string CtfCorrectionType();
	float MemoryLimit();

	bool CorrectAstigmatism();


private:

	void initVariables();

	void parseComFile(const string tiltcom);
	void parseCommandLine(vector<string> argList);

	void storeValues(string paramName, string paramValue, char separator);
	bool parameterWithoutValue(string paramName);

	novaCTF::Vec3ui volumeDimensions;
	string inputStackName;
	string tiltAnglesFilename;
	string outputFilename;
	float xAxisTilt;
	string skipProjections;
	novaCTF::Vec2f zShift;
	novaCTF::Vec2i subsetStart;
	unsigned int outputMode;

	unsigned int width;
	unsigned int height;
	unsigned int depth; 	//output volume dimensions AFTER binning
	int binning;


	int subsetWidth;		//can be also negative according to eTomo specification
	int subsetHeight;		//can be also negative according to eTomo specification
	int zShiftX;			//can be also negative according to eTomo specification
	int zShiftY;			//can be also negative according to eTomo specification
	std::string skipProjectionsList;

	bool useRadialFilter;
	bool logarithmizeData;
	bool useScaling;

	float radialCutOff;
	float radialFallOff;

	float scaling;
	float scalingOffset;

	novaCTF::Vec2f offset;
	string algorithmType;

	unsigned int numberOfInputStacks;
	novaCTF::VolumeThickness volumeThicknessType;
	bool writeOutDefocusSlices;

	string defocusFile;
	string defocusFileFormat;
	float pixelSize;

	bool use3DCTF;

	bool filterProjections;
	string defocusShiftFile;
	bool useAdditionalShift;
	float defocusStep;

	novaCTF::VolumeRotation stackOrientation;
	bool keepFilesOpen;

	float amplitude;
	float cs;
	float evk;

	float memoryLimit;

	string ctfCorrectionType;

	bool correctAstigmatism;

};
