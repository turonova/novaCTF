#pragma once

#include "parameterSetup.h"
#include "projectionSet.h"
#include "microscopeGeometry.h"


using namespace std;
using namespace novaCTF;


class DefocusFileFormat
{
public:

    static DefocusFileFormat* createFileFormat(string fileFormatName);
    virtual void read(string fileName, ProjectionSet& projSet) = 0;
    virtual void getValues(std::vector<std::vector<float>>& values, ProjectionSet& projSet, std::string units) = 0;
    virtual void writeWithShiftedDefocii(std::vector<std::vector<float>>& newValues, string fileName, ProjectionSet& projSet, std::string units) = 0;
};


class CTFFind4: public DefocusFileFormat
{
public:
    void read(string fileName, ProjectionSet& projSet);
    void getValues(std::vector<std::vector<float>>& values, ProjectionSet& projSet, std::string units);
    void writeWithShiftedDefocii(std::vector<std::vector<float>>& newValues, string fileName, ProjectionSet& projSet, std::string units);
private:
    void skipComments(ifstream& infile);

    std::vector<std::vector<float>> originalValues;
};

class ImodCTFPlotter: public DefocusFileFormat
{
public:
    void read(string fileName, ProjectionSet& projSet);
    void getValues(std::vector<std::vector<float>>& values, ProjectionSet& projSet, std::string units);
    void writeWithShiftedDefocii(std::vector<std::vector<float>>& newValues, string fileName, ProjectionSet& projSet, std::string units);
private:

    int parseHeader(ifstream& infile);

    bool hasAstigmatismValues;
    int formatFlag;
    unsigned int numberOfColumns;
    bool containsHeader;
    std::string originalHeader;
    std::vector<std::vector<float>> originalValues;
};

class GCTF: public DefocusFileFormat
{
public:
    void read(string fileName, ProjectionSet& projSet);
    void getValues(std::vector<std::vector<float>>& values, ProjectionSet& projSet, std::string units);
    void writeWithShiftedDefocii(std::vector<std::vector<float>>& newValues, string fileName, ProjectionSet& projSet, std::string units);
private:
    void parseAndSaveHeader(ifstream& infile);

    std::vector<std::vector<string>> originalValues;
    std::vector<std::string> originalHeader;
    unsigned int numberOfColumns;
    std::vector<unsigned int> columnIndices;
    bool containsPhaseShift;
};

