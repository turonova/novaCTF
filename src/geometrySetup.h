#pragma once

#include "common.h"

using namespace novaCTF;

class GeometrySetup
{
public:

	Vec3f c_bBoxMinComplete;
	Vec3f c_bBoxMaxComplete;

	Vec3f c_volumeDimComplete;
	Vec3f c_voxelSize;

	Vec3f c_ray;

	Vec3f c_source;
	Vec3f c_detektor;
	Vec3f c_zPitch;
	Vec3f c_yPitch;

	GeometrySetup(){}

};
