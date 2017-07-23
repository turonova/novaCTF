#include "defocusFileFormats.h"
#include "exception.h"
#include <sstream>
#include "fftRoutines.h"
#include "volumeIO.h"
#include <fstream>
#include <iomanip>

DefocusFileFormat* DefocusFileFormat::createFileFormat(string fileFormatName)
{

	if(fileFormatName=="imod")
		return new ImodCTFPlotter;
	else if(fileFormatName=="ctffind4")
		return new CTFFind4;
	else
	{
		cout << "The following defocus file format is unknown: " << fileFormatName << std::endl;
		exit(EXIT_FAILURE);
	}
}

void CTFFind4::skipComments(ifstream& infile)
{
	std::string line;
	unsigned int linesToSkip=0;

	while(getline(infile,line))
	{
		if(line[0]=='#')
			linesToSkip++;
		else
			break;
	}

	infile.seekg(0,infile.beg);

	for(unsigned int i=0; i<linesToSkip;i++)
	{
		string line;
		getline(infile,line);
	}
}

void CTFFind4::getValues(std::vector<std::vector<float>>& values, ProjectionSet& projSet, std::string units)
{
	double conversionFactor=1.0f;

	if(units=="nanometers")
	{
		conversionFactor=10e-2;
	}
	else if(units=="microns")
	{
		conversionFactor=10e-5;
	}



	for(ProjectionSet::iterator it=projSet.begin(); it!=projSet.end(); it++)
	{
		unsigned int projIndex = it.second();

		values[projIndex].push_back(originalValues[projIndex][1]*conversionFactor);	//defocus #1
		values[projIndex].push_back(originalValues[projIndex][2]*conversionFactor);	//defocus #2
		values[projIndex].push_back(originalValues[projIndex][3]);			//azimuth of astigmatism in degrees
		values[projIndex].push_back(originalValues[projIndex][4]);			//phase shift in radians
	}
}

/* Assuming format from CTFFind4 (version 4.0.16):
 * First few lines are comments and are skipped
 * 1 column: tilt number
 * 2 column: defocus #1 in Angstroms
 * 3 column: defocus #2 in Angstroms
 * 4 column: azimuth of astigmatism in degrees
 * 5 column: phase shift in radians
 * 6 column: cross correlation
 * 7 column: spacing (in Angstroms) up to which CTF rings were fit successfully
 */
void CTFFind4::read(string fileName, ProjectionSet& projSet)
{
	ifstream infile;
	infile.open(fileName.c_str());

	if (!infile.good())
	{
		throw ExceptionFileOpen(fileName.c_str());
	}

	skipComments(infile);

	originalValues.resize(projSet.getSize());

	for(ProjectionSet::iterator it=projSet.begin(); it!=projSet.end(); it++)
	{
		unsigned int projIndex = it.second();

		for(unsigned int i=0; i<7; i++)
		{
			float value;
			infile >> value;
			originalValues[projIndex].push_back(value);
		}
	}

	infile.close();
}

void CTFFind4::writeWithShiftedDefocii(std::vector<std::vector<float>>& newValues, string fileName, ProjectionSet& projSet, std::string units)
{
	ofstream file;
	file.open(fileName.c_str());

	if (!file.good())
	{
		throw ExceptionFileOpen(fileName.c_str());
	}

	double conversionFactor = 1.0;
	if(units=="nanometers")
	{
		conversionFactor=10.0;
	}
	else if(units=="microns")
	{
		conversionFactor=1.0e4;
	}

	for(ProjectionSet::iterator it=projSet.begin(); it!=projSet.end(); it++)
	{
		unsigned int projIndex = it.second();

		// Write out the first 1 value that correspond to tilt number
		file << std::setprecision(6) << std::fixed << originalValues[projIndex][0];
		file << "\t";

		// Change both defocus values in the file (it does not matter if we correct for astigmatism or not, the shifted
		// values are always computed for both defocii in the file and written out to avoid confusions
		for(unsigned int j=0; j<2; j++)
		{
			file << std::setprecision(2) << std::fixed << newValues[projIndex][j]*conversionFactor;
			file << "\t";
		}

		// Write out the rest of the original file values
		for(unsigned int i=3; i<originalValues[projIndex].size(); i++)
		{
			file  << std::setprecision(6) << std::fixed << originalValues[projIndex][i];
			file << "\t";
		}

		file << "\n";
	}

	file.close();
}

int ImodCTFPlotter::parseHeader(ifstream& infile)
{
	std::string line;
	getline(infile,line);
	originalHeader=line;

	std::vector<float> headerValues;
	unsigned int stringStart = 0;

	//parse first line and store values in floats
	for(unsigned int i = 0; i < line.length(); i++)
	{
		if(line[i]==' ')
		{
			string value = line.substr(stringStart, i);

			//skip over spaces
			while(line[i]==' ')
				i++;

			stringStart=i;
			headerValues.push_back(atof(value.c_str()));
		}
	}

	// last value on the line
	string value = line.substr(stringStart, line.length());
	headerValues.push_back(atof(value.c_str()));

	unsigned int numberOfValues=headerValues.size();

	hasAstigmatismValues = false;

	//no header in the file -> set pointer to the beginning of the file and return
	if(numberOfValues==5)
	{
		infile.seekg(0,infile.beg);
		numberOfColumns=5;
		containsHeader=false;
		originalHeader="";
		return 0;
	}

	//the file is in the format of version 2 (see ctfphaseflip man page), i.e. the first row contains
	//5 values as in case of no header but also 6th value with version specification that has to be skipped
	// -> set pointer to the beginning of the file and return flag to indicate the necessary skip
	if(numberOfValues==6 && headerValues[numberOfValues-1]==2)
	{
		infile.seekg(0,infile.beg);
		numberOfColumns=5;
		containsHeader=false;
		originalHeader="";
		return -1;
	}

	//the file is in the format of version 3 (see ctfphaseflip man page), i.e. the first row contains header
	//-> keep the pointer at the second line (where the actual values start) and return the first value from
	//the first row - it contains specification of values that were used
	if(numberOfValues==6 && headerValues[numberOfValues-1]==3)
	{
			infile.seekg(0,infile.beg);
			containsHeader=true;
			if(headerValues[0]==16)
			{
				numberOfColumns=5;
			}
			else if(headerValues[0]==4 || headerValues[0]==12 || headerValues[0]==20 || headerValues[0]==28)
			{
				numberOfColumns=6;
			}
			else if(headerValues[0]==1 || headerValues[0]==3 || headerValues[0]==17 || headerValues[0]==19)
			{
				numberOfColumns=7;
				hasAstigmatismValues=true;
			}
			else if( headerValues[0]==5 || headerValues[0]==7 || headerValues[0]==13 || headerValues[0]==15
				  || headerValues[0]==21 || headerValues[0]==23 || headerValues[0]==29 || headerValues[0]==31)
			{
				numberOfColumns=8;
				hasAstigmatismValues=true;
			}
			else
			{
				cout << "Following flag is unknown: " << headerValues[0] << ". See specification of flags on ctfphaseflip man page on IMOD webpage)!" << std::endl;
				exit(EXIT_FAILURE);
			}

			return (int)headerValues[0];
	}

	cout << "The defocus file does not correspond to the specification either of version 2 or 3 (see ctfphaseflip man page on IMOD webpage)!" << std::endl;
	exit(EXIT_FAILURE);
}

float degreesToRadians(float value)
{
	return value*M_PI/180.f;
}

float radiansToDegrees(float value)
{
	return value*180.f/M_PI;
}

/* Assuming format from version 2 or 3 from IMOD version 4.9 (see ctfphaseflip man page)
 * */
void ImodCTFPlotter::read(string fileName, ProjectionSet& projSet)
{
	ifstream infile;
	infile.open(fileName.c_str());

	if (!infile.good())
	{
		throw ExceptionFileOpen(fileName.c_str());
	}

	formatFlag=parseHeader(infile);

	originalValues.resize(projSet.getSize());

	for(ProjectionSet::iterator it=projSet.begin(); it!=projSet.end(); it++)
	{
		unsigned int projIndex = it.second();

		for(unsigned int i=0; i<numberOfColumns; i++)
		{
			float value;
			infile >> value;
			originalValues[projIndex].push_back(value);
		}

		// If the format is from version 2 it contains one more value on the first line that has to be read
		if(formatFlag==-1 && it.first()==0)
		{
			float value;
			infile >> value;
			originalValues[projIndex].push_back(value);
			formatFlag=0;
		}
	}

	infile.close();
}

void ImodCTFPlotter::writeWithShiftedDefocii(std::vector<std::vector<float>>& newValues, string fileName, ProjectionSet& projSet, std::string units)
{
	ofstream file;
	file.open(fileName.c_str());

	if (!file.good())
	{
		throw ExceptionFileOpen(fileName.c_str());
	}

	if(containsHeader)
	{
		file << originalHeader;
		file << "\n";
	}

	unsigned int valuesToChange=1;
	if(hasAstigmatismValues)
		valuesToChange=2;

	double conversionFactor = 1.0;
	if(units=="angstrom")
	{
		conversionFactor=10e-2;
	}
	else if(units=="microns")
	{
		conversionFactor=1.0e4;
	}

	for(ProjectionSet::iterator it=projSet.begin(); it!=projSet.end(); it++)
	{
		unsigned int projIndex = it.second();

		// Write out the first 4 values that correspond to tilts and agles
		for(unsigned int i=0; i<4; i++)
		{
			file << originalValues[projIndex][i];
			file << " ";
		}

		// Change all defocus values in the file (it does not matter if we correct for astigmatism or not, the shifted
		// values are always computed for all defocii in the file and written out to avoid confusions
		for(unsigned int j=0; j<valuesToChange; j++)
		{
			file << round(newValues[projIndex][j]);
			file << " ";
		}

		// Write out the rest of the original file values
		for(unsigned int i=4+valuesToChange; i<originalValues[projIndex].size(); i++)
		{
			file << originalValues[projIndex][i];
			file << " ";
		}

		file << "\n";
	}

	file.close();

}

void ImodCTFPlotter::getValues(std::vector<std::vector<float>>& values, ProjectionSet& projSet, std::string units)
{
	double conversionFactor = 1.0f;

	if(units=="angstroms")
	{
		conversionFactor=10.0;
	}
	else if(units=="microns")
	{
		conversionFactor=10e-4;
	}

	for(ProjectionSet::iterator it=projSet.begin(); it!=projSet.end(); it++)
	{
		unsigned int projIndex = it.second();

		if(formatFlag==0 || formatFlag==16)
		{
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 		// defocus #1
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 		// defocus #2 same as #1 in case of no astigmatism
			values[projIndex].push_back(0.0f);								// no astigmatims
			values[projIndex].push_back(0.0f);								// no phase shift
		}
		else if(formatFlag==1 || formatFlag==17)		// the file has astigmatism values in degrees
		{
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 		// defocus #1
			values[projIndex].push_back(originalValues[projIndex][5]*conversionFactor); 		// defocus #2
			values[projIndex].push_back(originalValues[projIndex][6]);		// astigmatims in degrees
			values[projIndex].push_back(0.0f);								// no phase shift
		}
		else if(formatFlag==3 || formatFlag==19)		// the file has astigmatism values in radians
		{
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 						// defocus #1
			values[projIndex].push_back(originalValues[projIndex][5]*conversionFactor); 						// defocus #2
			values[projIndex].push_back(radiansToDegrees(originalValues[projIndex][6]));	// astigmatims in radians
			values[projIndex].push_back(0.0f);												// no phase shift

		}
		else if(formatFlag==4 || formatFlag==20)		//  the file has phase shifts in degrees
		{
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 						// defocus #1
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 						// defocus #2 same as #1 in case of no astigmatism
			values[projIndex].push_back(0.0f);												// no astigmatims
			values[projIndex].push_back(degreesToRadians(originalValues[projIndex][5]));	//phase shift in degrees
		}
		else if(formatFlag==12 || formatFlag==28)		//  the file has phase shifts in radians
		{
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 			// defocus #1
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 			// defocus #2 same as #1 in case of no astigmatism
			values[projIndex].push_back(0.0f);									// no astigmatims
			values[projIndex].push_back(originalValues[projIndex][5]);			//phase shift in radians
		}
		else if(formatFlag==13 || formatFlag==29)		//  the file has phase shifts in radians and astigmatism in degrees
		{
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 			// defocus #1
			values[projIndex].push_back(originalValues[projIndex][5]*conversionFactor); 			// defocus #2
			values[projIndex].push_back(originalValues[projIndex][6]);			// astigmatims in degrees
			values[projIndex].push_back(originalValues[projIndex][7]);			// phase shift in radians
		}
		else if(formatFlag==5 || formatFlag==21)		//  the file has phase shifts in degrees and astigmatism in degrees
		{
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 						// defocus #1
			values[projIndex].push_back(originalValues[projIndex][5]*conversionFactor); 						// defocus #2
			values[projIndex].push_back(originalValues[projIndex][6]);						// astigmatims in degrees
			values[projIndex].push_back(degreesToRadians(originalValues[projIndex][7]));	// phase shift in degrees
		}
		else if(formatFlag==7 || formatFlag==23)		//  the file has phase shifts in degress and astigmatism in radians
		{
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 						// defocus #1
			values[projIndex].push_back(originalValues[projIndex][5]*conversionFactor); 						// defocus #2
			values[projIndex].push_back(radiansToDegrees(originalValues[projIndex][6]));	// astigmatims in radians
			values[projIndex].push_back(degreesToRadians(originalValues[projIndex][7]));	// phase shift in degrees
		}
		else if(formatFlag==15 || formatFlag==31)		//  the file has phase shifts in radians and astigmatism in radians
		{
			values[projIndex].push_back(originalValues[projIndex][4]*conversionFactor); 						// defocus #1
			values[projIndex].push_back(originalValues[projIndex][5]*conversionFactor); 						// defocus #2
			values[projIndex].push_back(radiansToDegrees(originalValues[projIndex][6]));	// astigmatims in radians
			values[projIndex].push_back(originalValues[projIndex][7]);	// phase shift in radians
		}
		else
		{
			cout << "The defocus file does not correspond to specification either of version 2 or 3 (see ctfphaseflip man page on IMOD webpage)!" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}
