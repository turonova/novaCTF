    #include <fstream>
    #include "microscopeGeometry.h"
    #include "exception.h"
    #include <cmath>
    #include <math.h>
    #include <float.h>
    #include <algorithm>


        Geometry::Geometry(MRCStack& aStack,novaCTF::Vec3ui volumeResolution, string aTiltAnglesFileName, float aXAxisTiltAngle, Vec2f zShift, float additionalTilt)
        {
            GeomHeader* header=aStack.getHeader();

            mXAxisTiltAngle = -aXAxisTiltAngle/180.f*M_PI;

            mDetector = novaCTF::makeVec3f(-0.5,-0.5,-1.f);
            mSource = novaCTF::makeVec3f(-0.5,-0.5,1.f);
            mHorizontalPitch = novaCTF::makeVec3f(1.f,0.f,0.f);
            mVerticalPitch = novaCTF::makeVec3f(0.f,1.f,0.f);

            novaCTF::expand2Vec3f(&mDetector,&header->mDetectorPosition);
            novaCTF::expand2Vec3f(&mSource,&header->mSourcePosition);
            novaCTF::expand2Vec3f(&mHorizontalPitch,&header->mHorizontalPitch);
            novaCTF::expand2Vec3f(&mVerticalPitch,&header->mVerticalPitch);

            //load tilt angles from given file
            mTiltAngles.resize(header->mDepth);
            loadTiltAngles(aTiltAnglesFileName);

            pretilt=additionalTilt;

            //volume
            mVoxelSize.x = 1.0f;
            mVoxelSize.y = 1.0f;
            mVoxelSize.z = 1.0f;

            mPosition = novaCTF::makeVec3f(mVoxelSize.x * (-0.5f*(float)volumeResolution.x), mVoxelSize.y * (-0.5f*(float)volumeResolution.y), mVoxelSize.z * (-0.5f*(float)volumeResolution.z)-zShift.y);

            setVolumeGeometry(volumeResolution);
        }

        Geometry::~Geometry()
        {
        }

        void Geometry::setVolumeGeometry(novaCTF::Vec3ui volumeResolution)
        {

            //Dimensions of the volume
            mSetup.c_volumeDimComplete = novaCTF::makeVec3f((float)volumeResolution.x,(float)volumeResolution.y,(float)volumeResolution.z);

            //Bounding boxes
            Vec3f bBoxComplete = novaCTF::makeVec3f(mVoxelSize.x*(float)volumeResolution.x,mVoxelSize.y*volumeResolution.y,mVoxelSize.z*(float)volumeResolution.z);

            mSetup.c_bBoxMinComplete=mPosition;
            mSetup.c_bBoxMaxComplete=novaCTF::addVec3f(&mPosition, &bBoxComplete);

            //Voxel sizes
            mSetup.c_voxelSize=mVoxelSize;
        }

        novaCTF::Vec3f Geometry::getVolumeBBoxMin()
        {
            return mSetup.c_bBoxMinComplete;
        }

        novaCTF::Vec3f Geometry::getVolumeBBoxMax()
        {
            return mSetup.c_bBoxMaxComplete;
        }

        float Geometry::getAngleInDegrees(unsigned int aProjectionIndex)
        {
            return mTiltAngles[aProjectionIndex];
        }

        float Geometry::getAngleInRadians(unsigned int aProjectionIndex)
        {
            float tiltAngle=mTiltAngles[aProjectionIndex];
            tiltAngle=tiltAngle/180.f*M_PI;

            return tiltAngle;
        }

        unsigned int Geometry::getDefocusID(std::vector<novaCTF::Vec4f>& focusGrid, Vec2f point)
        {
            for(unsigned int i = 0; i<focusGrid.size(); i++)
            {
                float position = ((focusGrid[i].z - focusGrid[i].x) * (point.y - focusGrid[i].y) - (focusGrid[i].w - focusGrid[i].y) * (point.x - focusGrid[i].x));
                if(position<=0.0f)
                    return i;
            }

            return focusGrid.size()-1;
        }

        float Geometry::computeVolumeThickness(VolumeThickness volumeThickness)
        {
            if(volumeThickness==VolumeThickness::MAXIMAL)
            {
                float maxDifferenceInZ = -FLT_MAX;
                for(unsigned int i = 0; i<mTiltAngles.size();i++)
                {
                    float maxTiltAngleInDegrees=mTiltAngles[i]/180.f*M_PI;
                    Vec3f cornerTR = mSetup.c_bBoxMaxComplete;
                    Vec3f cornerBL = mSetup.c_bBoxMinComplete;
                    Vec3f cornerBR = mSetup.c_bBoxMinComplete;
                    Vec3f cornerTL = mSetup.c_bBoxMinComplete;
                    cornerTL.z=cornerTR.z;  // cortneTL.z=-cornerTL.z is incorrect in case the z-shift is not zero!!!
                    cornerBR.x=cornerTR.x;
                    rotate(cornerTR,maxTiltAngleInDegrees);
                    rotate(cornerBL,maxTiltAngleInDegrees);
                    rotate(cornerBR,maxTiltAngleInDegrees);
                    rotate(cornerTL,maxTiltAngleInDegrees);
                    float maxZ=max(cornerTR.z,max(cornerBL.z,max(cornerBR.z,cornerTL.z)));
                    float minZ=min(cornerTR.z,min(cornerBL.z,min(cornerBR.z,cornerTL.z)));
                    maxDifferenceInZ = max(maxDifferenceInZ,fabs(maxZ-minZ));
                }

                return maxDifferenceInZ;
            }
            else
                return mSetup.c_bBoxMaxComplete.z-mSetup.c_bBoxMinComplete.z;
        }

        unsigned int Geometry::computeNumberOfParts(float stripeSize, float pixelSize, VolumeThickness volumeThickness)
        {
            float volumeThicknessInNm = computeVolumeThickness(volumeThickness)*pixelSize;
            unsigned int numberOfStripes = floor(volumeThicknessInNm/stripeSize);

            if(numberOfStripes==0)
                numberOfStripes=1;

            if(numberOfStripes%2==0)	//number of stripes is even
            {
                numberOfStripes--;
            }

            return numberOfStripes;
        }

        void Geometry::generateFocusGrid(std::vector<unsigned int>& sliceDefocusSplit, unsigned int numberOfParts, VolumeThickness volumeThickness)
        {
            std::vector<novaCTF::Vec4f> focusGrid;
            int endBoundary =  (int)numberOfParts/2+1;
            int startBoundary = -(int)numberOfParts/2+1;
            int centerAddCoefficient = 1;

            if((numberOfParts%2)==0)
            {
                endBoundary = (int)numberOfParts/2;
                centerAddCoefficient = 0;
            }

            float stripeSize = computeVolumeThickness(volumeThickness)/numberOfParts;

            novaCTF::Vec2f startPoint = Vec2f((mSetup.c_bBoxMaxComplete.x+mSetup.c_bBoxMinComplete.x)*0.5f-mSetup.c_ray.x*stripeSize*0.5*centerAddCoefficient,(mSetup.c_bBoxMaxComplete.z+mSetup.c_bBoxMinComplete.z)*0.5f-mSetup.c_ray.z*stripeSize*0.5*centerAddCoefficient);
            novaCTF::Vec2f endPoint;
            endPoint.x = startPoint.x-mSetup.c_yPitch.x*mSetup.c_volumeDimComplete.x; //rotated horizontal pitch
            endPoint.y = startPoint.y-mSetup.c_yPitch.z*mSetup.c_volumeDimComplete.x; //rotated vertical pitch


            for(int i = startBoundary; i <= endBoundary; i++)
            {
                focusGrid.push_back(Vec4f(startPoint.x+mSetup.c_ray.x*i*stripeSize,startPoint.y+mSetup.c_ray.z*i*stripeSize,endPoint.x+mSetup.c_ray.x*i*stripeSize,endPoint.y+mSetup.c_ray.z*i*stripeSize));
            }

            for(float zCoord = mSetup.c_bBoxMinComplete.z + 0.5f; zCoord < mSetup.c_bBoxMaxComplete.z; zCoord+=1.0f)
            {
                for(float xCoord = mSetup.c_bBoxMinComplete.x + 0.5f; xCoord < mSetup.c_bBoxMaxComplete.x; xCoord+=1.0f)
                {
                    sliceDefocusSplit.push_back(getDefocusID(focusGrid,Vec2f(xCoord,zCoord)));
                }

            }

        }

        void Geometry::calculateGridCenters(std::vector<std::vector<std::vector<float>>>& newDefocusCenters, std::vector<std::vector<float>>& defocusCenters, unsigned int numberOfParts, VolumeThickness volumeThickness, float pixelSize,float zShift)
        {
            int endBoundary =  (int)numberOfParts/2;
            int startBoundary = -(int)numberOfParts/2;
            int centerAddCoefficient = 0;

            if((numberOfParts%2)==0)
            {
                endBoundary = (int)numberOfParts/2-1;
                centerAddCoefficient = 1;
            }

            float stripeSize = computeVolumeThickness(volumeThickness)/numberOfParts;

            int partIndex = 0;
            for(int i = startBoundary; i <= endBoundary; i++)
            {
                for(size_t projIndex = 0; projIndex < newDefocusCenters[partIndex].size(); projIndex++)
                {
                    newDefocusCenters[partIndex][projIndex].push_back(defocusCenters[projIndex][0]-(zShift*pixelSize)+(i*stripeSize*pixelSize)+(stripeSize*pixelSize*0.5)*centerAddCoefficient);
                    newDefocusCenters[partIndex][projIndex].push_back(defocusCenters[projIndex][1]-(zShift*pixelSize)+(i*stripeSize*pixelSize)+(stripeSize*pixelSize*0.5)*centerAddCoefficient);
                }
                partIndex++;
            }
        }

        void Geometry::setProjectionGeometry(unsigned int aProjectionIndex)
        {

            float tiltAngle=mTiltAngles[aProjectionIndex]-pretilt; //??? sign
            tiltAngle=-tiltAngle/180.f*M_PI;

            Vec3f detector = mDetector;
            Vec3f source = mSource;
            Vec3f horizontalPitch = mHorizontalPitch;
            Vec3f verticalPitch = mVerticalPitch;

            rotate(detector,tiltAngle);
            rotate(source,tiltAngle);
            rotate(horizontalPitch,tiltAngle);
            rotate(verticalPitch,tiltAngle);

            mSetup.c_detektor=detector;
            mSetup.c_source=source;
            mSetup.c_yPitch=horizontalPitch;
            mSetup.c_zPitch=verticalPitch;

            mSetup.c_ray = novaCTF::subVec3f(&detector,&source);
            mSetup.c_ray = novaCTF::normalizeVec3f(&mSetup.c_ray);
        }

        void Geometry::rotate(Vec3f &aCoord, float aTiltAngle)
        {

            if(mXAxisTiltAngle==0.0f)
            {	// there is no rotation around x axis
                aCoord = novaCTF::makeVec3f(aCoord.x * cos(aTiltAngle) - aCoord.z * sin(aTiltAngle), aCoord.y, aCoord.x * sin(aTiltAngle) + aCoord.z * cos(aTiltAngle));
            }
            else
            {
                aCoord = novaCTF::makeVec3f(aCoord.x, aCoord.y * cos(mXAxisTiltAngle) + aCoord.z * sin(mXAxisTiltAngle), - aCoord.y * sin(mXAxisTiltAngle) + aCoord.z * cos(mXAxisTiltAngle));	//rotate the point around x axis

                novaCTF::Vec3f axis = novaCTF::makeVec3f(0.f, cos(mXAxisTiltAngle),-sin(mXAxisTiltAngle));	//y axis around x
                axis=novaCTF::normalizeVec3f(&axis);

                novaCTF::Vec3f tempPoint = aCoord;

                aCoord.x = (cos(aTiltAngle)+(1-cos(aTiltAngle))*axis.x*axis.x)		  * tempPoint.x + ((1-cos(aTiltAngle))*axis.x*axis.y-sin(aTiltAngle)*axis.z) * tempPoint.y + ((1-cos(aTiltAngle))*axis.x*axis.z+sin(aTiltAngle)*axis.y) * tempPoint.z;
                aCoord.y = ((1-cos(aTiltAngle))*axis.x*axis.y+sin(aTiltAngle)*axis.z) * tempPoint.x + (cos(aTiltAngle)+(1-cos(aTiltAngle))*axis.y*axis.y)        * tempPoint.y + ((1-cos(aTiltAngle))*axis.y*axis.z-sin(aTiltAngle)*axis.x) * tempPoint.z;
                aCoord.z = ((1-cos(aTiltAngle))*axis.x*axis.z-sin(aTiltAngle)*axis.y) * tempPoint.x + ((1-cos(aTiltAngle))*axis.y*axis.z+sin(aTiltAngle)*axis.x) * tempPoint.y + (cos(aTiltAngle)+(1-cos(aTiltAngle))*axis.z*axis.z)        * tempPoint.z;
            }

        }


        void Geometry::loadTiltAngles(string aTiltAnglesFileName)
        {
            ifstream infile;
            infile.open(aTiltAnglesFileName.c_str());

            if (!infile.good())
            {
                throw ExceptionFileOpen(aTiltAnglesFileName);
            }

            for(size_t i = 0; i<mTiltAngles.size();i++)
            {
                infile >>mTiltAngles[i];
                //mTiltAngles[i]+=additionalTilt;
            }

            infile.close();
        }

        novaCTF::Vec2f Geometry::projectPointOnDetector(Vec3f point)
        {
            Vec3f pointDet = point - mSetup.c_detektor;
            float x = novaCTF::dotVec3f(&mSetup.c_yPitch, &pointDet);
            float y = novaCTF::dotVec3f(&mSetup.c_zPitch, &pointDet);

            return Vec2f(x,y);
        }

        novaCTF::Vec4f Geometry::voxelProjectionBoundaries()
        {
            Vec2f maxValues(-FLT_MAX,-FLT_MAX);
            Vec2f minValues(FLT_MAX,FLT_MAX);

            for(float x = -0.5f; x <=0.5f; x+=1.0f)
                for(float y = -0.5f; y <=0.5f; y+=1.0f)
                    for(float z = -0.5f; z <=0.5f; z+=1.0f)
                    {
                        Vec2f projectedCoordinates = projectPointOnDetector(Vec3f(x,y,z));
                        minValues = minVec2f(&minValues,&projectedCoordinates);
                        maxValues = maxVec2f(&maxValues,&projectedCoordinates);
                    }

            Vec2f projectedCenter = projectPointOnDetector(0.0f);

            return Vec4f(projectedCenter.x-minValues.x,projectedCenter.y-minValues.y,projectedCenter.x-maxValues.x,projectedCenter.y-maxValues.y);

        }
