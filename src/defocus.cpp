#include "defocus.h"
#include "mrcStack.h"
#include "volumeIO.h"
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include "fftRoutines.h"
#include "exception.h"
#include <iomanip>


using namespace novaCTF;
using namespace std;

Defocus::Defocus(ParameterSetup& aParams)
{
	params=aParams;
	projSet = new ProjectionSetIdentity();
	inputStack = new MRCStack(params.InputStackName(),false, false, false);
	projSet->init(inputStack->getNumberOfProjections());
	microscopeGeometry = new Geometry(*inputStack,params.VolumeDimensions(),params.TiltAnglesFilename(),params.XAxisTilt(),params.ZShift());
	defocusFileFormat = params.DefocusFileFormat();
}

Defocus::~Defocus()
{
	delete microscopeGeometry;
	delete projSet;
	delete inputStack;
}

void Defocus::run()
{
	readDefocusFile();
	calculateDefocus();
	writeDefocusFiles();

	cout << "Number of generated files: " << numberOfInputStacks << endl;
}

void Defocus::readDefocusFile()
{
	defocusFileValues.resize(inputStack->getNumberOfProjections());
	defFile=DefocusFileFormat::createFileFormat(defocusFileFormat);
	defFile->read(params.DefocusFile(),*projSet);
	defFile->getValues(defocusFileValues,*projSet,"nanometers");

}

void Defocus::writeDefocusFiles()
{
	for(unsigned int i=0; i<numberOfInputStacks;i++)
	{
		defFile->writeWithShiftedDefocii(newDefocusValues[i],generateFilename(params.DefocusFile(),numberOfInputStacks-i-1),*projSet,"nanometers");
	}
}

void Defocus::computeCompleteAdditionalShift()
{
	if(params.UseAdditionalShift())
	{
		ifstream infile;
		infile.open(params.DefocusShiftFile().c_str());

		if (!infile.good())
		{
			throw ExceptionFileOpen(params.DefocusShiftFile().c_str());
		}

		float additionalZShift;
		infile >> additionalZShift;

		additionalZShift = (float)params.VolumeDimensions().z*0.5f - additionalZShift;

		completeDefocusZShift = params.ZShift().y + additionalZShift;
		infile.close();
	}
	else
		completeDefocusZShift = params.ZShift().y;

	cout << "Using following z-shift (in voxels): " << completeDefocusZShift << endl;
}

void Defocus::calculateDefocus()
{
	if(params.DefocusStep()!=0.0f)
		numberOfInputStacks = microscopeGeometry->computeNumberOfParts(params.DefocusStep(),params.PixelSize(),params.VolumeThicknessType());
	else if(params.NumberOfInputStacks()!=0)
		numberOfInputStacks = params.NumberOfInputStacks();
	else
	{
		cout << "Either defocus step size (DefocusStep in nm) or number of input stacks (NumberOfInputStacks) has to be specified!" << endl;
		return;
	}

	newDefocusValues.resize(numberOfInputStacks);

	for(unsigned int i=0; i<numberOfInputStacks;i++)
	{
		newDefocusValues[i].resize(inputStack->getNumberOfProjections());
	}

	computeCompleteAdditionalShift();

	microscopeGeometry->calculateGridCenters(newDefocusValues,defocusFileValues,numberOfInputStacks,params.VolumeThicknessType(),params.PixelSize(),completeDefocusZShift);
}

string Defocus::generateFilename(string originalFilename, unsigned int number)
{
	stringstream newFilename;
	newFilename << originalFilename << "_" << number;
	return newFilename.str();
}
