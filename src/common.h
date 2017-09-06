#pragma once

#include <sstream>
#include <cstddef>
#include <float.h>
/**
	Defines some constants, which are used throughout the whole renderer.
*/

#ifndef M_PI
    #define M_PI 3.141592654
#endif

#ifndef M_PI2
    #define M_PI2 6.2831853
#endif

#define PRINTF2

#define __NO_STD_VECTOR
#define __NO_STD_STRING

#define cosd(x) (cos(fmod((x),360) * M_PI / 180))

namespace novaCTF
{
    enum VolumeRotation
    {
        ALONG_XY    = 0, //normal rotation => no flipping
        ALONG_XZ    = 1  //eTomo rotation
    };

    enum VolumeThickness
    {
        MINIMAL     = 0, //corresponds to volume thickness
        MAXIMAL     = 1  //corresponds to volume span in z for the max tilt angle
    };

    struct DataStats
    {
        float mean;
        float max;
        float min;

        DataStats()
        {
            mean = 0.0;
            max = -FLT_MAX;
            min = FLT_MAX;
        }

    };

    struct Vec2i
    {
        int x;
        int y;

        Vec2i()
        {
            x=0;
            y=0;
        }

        Vec2i(int aX, int aY)
        {
            x=aX;
            y=aY;
        }
    };

    struct Vec2f
    {
    public:
      float x;
      float y;

        ~Vec2f()
        {}

        Vec2f()
        {
          x=0.0f;
            y=0.0f;
        }

        Vec2f(float a)
        {
            x=a;
            y=a;
        }

        Vec2f(float aX, float aY)
        {
            x=aX;
            y=aY;
        }

        Vec2f operator+(const Vec2f& other) const;
        Vec2f operator*(float scalar) const;
        Vec2f operator-(const Vec2f& other) const;
    };

    struct Vec2ui
    {
        unsigned int x;
        unsigned int y;
    };

    struct Vec4f
    {
    public:
        float x;
        float y;
        float z;
        float w;

        ~Vec4f()
        {}

        Vec4f()
        {
            x=0.0f;
            y=0.0f;
            z=0.0f;
            w=0.0f;
        }

        Vec4f(float aX, float aY, float aZ,float aW)
        {
            x=aX;
            y=aY;
            z=aZ;
            w=aW;
        }
    };

    struct Vec4i
    {
    public:
        int x;
        int y;
        int z;
        int w;

        Vec4i():x(0),y(0),z(0),w(0)
        {}

        Vec4i(int aX, int aY, int aZ, int aW)
        {
            x=aX;
            y=aY;
            z=aZ;
            w=aW;
        }

        Vec4i operator+(const Vec4i& other) const;
    };

    struct Vec4ui
    {
        unsigned int x;
        unsigned int y;
        unsigned int z;
        unsigned int w;
    };

    struct Vec3i
    {
    public:
        int x;
        int y;
        int z;

        Vec3i():x(0),y(0),z(0)
        {}

        Vec3i(int aX, int aY, int aZ)
        {
            x=aX;
            y=aY;
            z=aZ;
        }

        Vec3i(int value)
        {
            x=value;
            y=value;
            z=value;
        }
    };

    struct Vec3ui
    {
    public:
        unsigned int x;
        unsigned int y;
        unsigned int z;

        Vec3ui():x(0),y(0),z(0)
        {}

        Vec3ui(unsigned int aX, unsigned int aY, unsigned int aZ)
        {
            x=aX;
            y=aY;
            z=aZ;
        }
    };

    struct Vec3t
    {
    public:
        size_t x;
        size_t y;
        size_t z;

        Vec3t():x(0),y(0),z(0)
        {}

        Vec3t(size_t aX, size_t aY, size_t aZ)
        {
            x=aX;
            y=aY;
            z=aZ;
        }

        Vec3t(int aX, int aY, int aZ)
        {
            x=(size_t) aX;
            y=(size_t) aY;
            z=(size_t) aZ;
        }

        Vec3t(unsigned int aX, unsigned int aY, unsigned int aZ)
        {
            x=(size_t) aX;
            y=(size_t) aY;
            z=(size_t) aZ;
        }

        Vec3t(Vec3ui vec)
        {
            x=(size_t) vec.x;
            y=(size_t) vec.y;
            z=(size_t) vec.z;
        }

    };

    struct Vec3f
    {
    public:
        float x;
        float y;
        float z;

        Vec3f():x(0.0f),y(0.0f),z(0.0f)
        {}

        Vec3f(float a)
        {
            x=a;
            y=a;
            z=a;
        }

        Vec3f(float aX, float aY, float aZ)
        {
            x=aX;
            y=aY;
            z=aZ;
        }

        Vec3f operator+(const Vec3f& other) const;
        Vec3f operator*(float scalar) const;
        Vec3f operator-(const Vec3f& other) const;
    };


    Vec3f makeVec3f(float aX, float aY, float aZ);

    Vec3f minVec3f(Vec3f *a, Vec3f *b);
    Vec3f maxVec3f(Vec3f *a, Vec3f *b);
    Vec2f minVec2f(Vec2f *a, Vec2f *b);
    Vec2f maxVec2f(Vec2f *a, Vec2f *b);

    Vec3f addVec3f(Vec3f *a, Vec3f *b);

    Vec3f subVec3f(Vec3f *a, Vec3f *b);

    Vec3f expandVec3f(Vec3f *a, Vec3f *b);

    void expand2Vec3f(Vec3f *a, Vec3f *b);

    Vec3f mulVec3f(Vec3f *a, float b);

    Vec3f divideVec3f(Vec3f *a, Vec3f *b);

    Vec3f normalizeVec3f(Vec3f *a);

    Vec3f inverseVec3f(Vec3f *a);

    Vec3ui makeVec3ui(unsigned int aX, unsigned int aY, unsigned int aZ);

    float convertRange(float value, float oldMin, float oldMax, float newMin, float newMax);

    float clamp(float n, float lower, float upper);

    Vec3f matrixMult(float* matrix, Vec3f vec);

    float dotVec3f(Vec3f *a, Vec3f *b);

    float sign(float val);

    std::string convertSeconds(double seconds_total);

    template<typename T>
    void swap(T& a, T& b)
    {
        T tmpValue = b;
        b = a;
        a = tmpValue;
    }

};
