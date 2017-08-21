#pragma once

#include <stdio.h>
#include <stdlib.h>
//#include "volume.h"
#include "common.h"
#include "mrcStack.h"
#include "parameterSetup.h"


    class VolumeIO
    {
    public:

        //static Volume* load(std::string fileName);
       // static Volume* load(MRCStack& dataSource);

       static void write(std::vector<float>& data, novaCTF::Vec3ui resolution, std::string outputVolumeFileName, int mode);
       static void write(std::vector<float>& data, novaCTF::Vec3ui resolution, std::string outputVolumeFileName, int mode, novaCTF::VolumeRotation rotation);

         static int getFileFormatModeForOutputVolume( ParameterSetup& parameters);

        static void writeHeader(std::string inputStackFileName, std::string outputVolumeFileName, novaCTF::Vec3ui volumeDimensions, int mode, novaCTF::VolumeRotation rotation);

        static void writeVolumeSliceInFloat(std::string outputVolumeFileName, std::vector<float> data, size_t sliceSize, int mode);
        static void writeVolumeSliceInFloat(std::string outputVolumeFileName, float& data,  size_t sliceSize, int mode);
        static void convertVolumeValues(std::string outputVolumeFileName, float min, float max, float mean, int mode);

        static void writeMRCStack(MRCHeader header, std::vector<float>& data, string outputName, char* extraData);
        static void writeMRCStack(MRCHeader header, std::vector<std::vector< float > >& data, string outputName, char* extraData);
    protected:

        static void computeMinMaxMean(std::vector<float>& data, size_t dataSize, float& min, float& max, float& mean);
        static void writeProjections(std::ofstream& outfile, std::vector<float>& data, size_t projectionSize, unsigned int numberOfProjections);

        static std::string generateTempName(std::string filename, int mode);

        template <typename voxelType>
        static void convertSliceData(std::string outputVolumeFileName, size_t sliceSize, unsigned int sliceNumber, float min, float max, float newMin, float newMax);

        static void correctHeader(std::string outputVolumeFileName, MRCHeader& header, float min, float max, float mean, int mode);

        template <typename voxelType>
        static void writeData(std::vector<float>& data, novaCTF::Vec3ui resolution, std::ofstream&  outfile, novaCTF::VolumeRotation rotation, float oldMin, float oldMax,float newMin, float newMax);

        template <typename voxelType>
        static void writeDataAlongXZPlane(std::vector<float>& data, novaCTF::Vec3ui resolution, std::ofstream&  outfile, float oldMin, float oldMax, float newMin, float newMax);

        template <typename voxelType>
        static void writeDataAlongXYPlane(std::vector<float>& data, novaCTF::Vec3ui resolution, std::ofstream&  outfile, float oldMin, float oldMax,float newMin, float newMax);

    };

    std::ostream& operator<<( std::ostream& oss, int mode );
