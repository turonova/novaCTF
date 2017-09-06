#pragma once

#include "common.h"

class GeomHeader
{
public:
    unsigned int mWidth;    //width of the volume/projections in the stack
    unsigned int mHeight;   //height of the volume/projections in the stack
    unsigned int mDepth;    //depth of the volume or number of projections in the stack

    novaCTF::Vec3f mSourcePosition;
    novaCTF::Vec3f mDetectorPosition;
    novaCTF::Vec3f mHorizontalPitch;
    novaCTF::Vec3f mVerticalPitch;

    GeomHeader(){}
    ~GeomHeader(){}

}; //end of class geomHeader

