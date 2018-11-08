# pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "common.h"
#include "mrcStack.h"
//#include "volume.h"
#include "geometrySetup.h"


/**	Class to store geometry of the stage
*	Handles rotations according to given angle
**/
class Geometry
{

public:

	/*Projection set data dependent*/

	novaCTF::Vec3f mDetector;			//coordinates of left upper corner of the detector (without rotation)
	novaCTF::Vec3f mSource;				//coordinates of left upper corner of the source (without rotation)
	novaCTF::Vec3f mHorizontalPitch;		//pitch in horizontal direction (without rotation)
	novaCTF::Vec3f mVerticalPitch;		//pitch in vertical direction (without rotation)

	std::vector<float> mTiltAngles;			//list of the angles for the projection set
	float mXAxisTiltAngle;		//tilt of the beam

	/* Volume data dependent*/
	novaCTF::Vec3f	mPosition;					//position of the volume - front bottom left corner
    novaCTF::Vec3f	mVoxelSize;					//voxel size

	/* General */
	GeometrySetup mSetup;		//geometry setup

	//Geometry(MRCStack& aProjData,Volume& aVolume, string aTiltAnglesFileName, float aXAxisTiltAngle, Vec2f zShift);
	Geometry(MRCStack& aProjData,novaCTF::Vec3ui volumeResolution, string aTiltAnglesFileName, float aXAxisTiltAngle, Vec2f zShift,float additionalTilt);

	~Geometry();

	void setProjectionGeometry(unsigned int aProjectionIndex);

	float computeRayTraversalLength();
	novaCTF::Vec4f voxelProjectionBoundaries();
	float getAngleInDegrees(unsigned int aProjectionIndex);
	float getAngleInRadians(unsigned int aProjectionIndex);
	void generateFocusGrid(std::vector<unsigned int>& sliceDefocusSplit, unsigned int numberOfParts,VolumeThickness thicknessType);
	float computeVolumeThickness(VolumeThickness volumeThickness);
	void calculateGridCenters(std::vector<std::vector<std::vector<float>>>& newDefocusCenters, std::vector<std::vector<float>>& defocusCenters, unsigned int numberOfParts, VolumeThickness volumeThickness, float pixelSize,float zShift);
	novaCTF::Vec3f getVolumeBBoxMin();
	novaCTF::Vec3f getVolumeBBoxMax();
	unsigned int computeNumberOfParts(float stripeSize, float pixelSize, VolumeThickness volumeThickness);

private:

		void rotate(novaCTF::Vec3f &aCoord, float aTiltAngle);
		void setVolumeGeometry(novaCTF::Vec3ui volumeResolution);
		void loadTiltAngles(string aTiltAnglesFileName);
		unsigned int getDefocusID(std::vector<novaCTF::Vec4f>& focusGrid, Vec2f point);
		novaCTF::Vec2f projectPointOnDetector(Vec3f point);

		float pretilt;

};
