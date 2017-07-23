#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "common.h"
#include <sstream>

using namespace novaCTF;

float novaCTF::convertRange(float value, float oldMin, float oldMax, float newMin, float newMax)
{
	float oldRange = (oldMax - oldMin);
	float newRange = (newMax - newMin);
	float newValue;

	if (oldRange == 0)
	{
	    newValue = newMin;
	}
	else
	{
		newValue= (((value - oldMin) * newRange) / oldRange) + newMin;
	}

	return newValue;

}

Vec3f novaCTF::minVec3f(Vec3f *a, Vec3f *b)
{
	Vec3f minVec;
	minVec.x = fmin(a->x,b->x);
	minVec.y = fmin(a->y,b->y);
	minVec.z = fmin(a->z,b->z);

	return minVec;
}

Vec3f novaCTF::maxVec3f(Vec3f *a, Vec3f *b)
{
	Vec3f maxVec;
	maxVec.x = fmax(a->x,b->x);
	maxVec.y = fmax(a->y,b->y);
	maxVec.z = fmax(a->z,b->z);

	return maxVec;
}

Vec2f novaCTF::minVec2f(Vec2f *a, Vec2f *b)
{
	Vec2f minVec;
	minVec.x = fmin(a->x,b->x);
	minVec.y = fmin(a->y,b->y);

	return minVec;
}

Vec2f novaCTF::maxVec2f(Vec2f *a, Vec2f *b)
{
	Vec2f maxVec;
	maxVec.x = fmax(a->x,b->x);
	maxVec.y = fmax(a->y,b->y);

	return maxVec;
}

Vec4i novaCTF::Vec4i::operator+(const Vec4i& other) const
{
	return Vec4i(x + other.x, y + other.y, z + other.z, w + other.w);
}

Vec3f novaCTF::Vec3f::operator+(const Vec3f& other) const
{
	return Vec3f(x + other.x, y + other.y, z + other.z);
}

Vec3f novaCTF::Vec3f::operator-(const Vec3f& other) const
{
	return Vec3f(x - other.x, y - other.y, z - other.z);
}

Vec3f novaCTF::Vec3f::operator*(float scalar) const
{
    return Vec3f(x * scalar, y * scalar, z * scalar);
}

Vec2f novaCTF::Vec2f::operator+(const Vec2f& other) const
{
	return Vec2f(x + other.x, y + other.y);
}

Vec2f novaCTF::Vec2f::operator-(const Vec2f& other) const
{
	return Vec2f(x - other.x, y - other.y);
}

Vec2f novaCTF::Vec2f::operator*(float scalar) const
{
    return Vec2f(x * scalar, y * scalar);
}


Vec3ui novaCTF::makeVec3ui(unsigned int aX, unsigned int aY, unsigned int aZ)
{
	Vec3ui ret;
	ret.x=aX;
	ret.y=aY;
	ret.z=aZ;
	return ret;
}

Vec3f novaCTF::makeVec3f(float aX, float aY, float aZ)
{
	Vec3f ret;
	ret.x=aX;
	ret.y=aY;
	ret.z=aZ;
	return ret;
}

Vec3f novaCTF::addVec3f(Vec3f *a, Vec3f *b)
{
	Vec3f ret;
	ret.x=a->x+b->x;
	ret.y=a->y+b->y;
	ret.z=a->z+b->z;
	return ret;
}

Vec3f novaCTF::subVec3f(Vec3f *a, Vec3f *b)
{
	Vec3f ret;
	ret.x=a->x-b->x;
	ret.y=a->y-b->y;
	ret.z=a->z-b->z;
	return ret;
}

//TODO rename

Vec3f novaCTF::expandVec3f(Vec3f *a, Vec3f *b)
{
	Vec3f ret;
	ret.x=a->x*b->x;
	ret.y=a->y*b->y;
	ret.z=a->z*b->z;
	return ret;
}

//TODO - rewrite
void novaCTF::expand2Vec3f(Vec3f *a, Vec3f *b)
{
	a->x=a->x*b->x;
	a->y=a->y*b->y;
	a->z=a->z*b->z;
}

Vec3f novaCTF::normalizeVec3f(Vec3f *a)
{
	Vec3f ret;

	float invDenom = 1.f/sqrt(a->x*a->x + a->y*a->y + a->z*a->z);

	ret.x=a->x*invDenom;
	ret.y=a->y*invDenom;
	ret.z=a->z*invDenom;

	return ret;
}

Vec3f novaCTF::mulVec3f(Vec3f *a, float b)
{
	Vec3f ret;
	ret.x=a->x*b;
	ret.y=a->y*b;
	ret.z=a->z*b;
	return ret;
}

float novaCTF::dotVec3f(Vec3f *a, Vec3f *b)
{
	return a->x*b->x + a->y*b->y + a->z*b->z;
}

Vec3f novaCTF::divideVec3f(Vec3f *a, Vec3f *b)
{
	Vec3f ret;

	if(b->x!=0.f) ret.x=a->x / b->x;
	else ret.x = 0.0f;

	if(b->y!=0.f) ret.y=a->y / b->y;
	else ret.y = 0.0f;

	if(b->z!=0.f) ret.z=a->z / b->z;
	else ret.z = 0.0f;

	return ret;
}

Vec3f novaCTF::inverseVec3f(Vec3f *a)
{
	Vec3f ret;

	if(a->x!=0.f) ret.x=1.0f / a->x;
	else ret.x = 0.0f;
	//else ret.x = 1000000.0f;

	if(a->y!=0.f) ret.y=1.0f / a->y;
	else ret.y = 0.0f;
	//else ret.y = 1000000.0f;

	if(a->z!=0.f) ret.z=1.0f / a->z;
	else ret.z = 0.0f;
	//else ret.z = 1000000.0f;

	return ret;
}

float novaCTF::clamp(float n, float lower, float upper)
{
	return fmax(lower, fmin(n, upper));
}

Vec3f novaCTF::matrixMult(float* matrix, Vec3f vec)
{
	Vec3f ret;

	ret.x = matrix[0] * vec.x + matrix[1] * vec.y + matrix[2] * vec.z +  matrix[3];
	ret.y = matrix[4] * vec.x + matrix[5] * vec.y + matrix[6] * vec.z +  matrix[7];
	ret.z = matrix[8] * vec.x + matrix[9] * vec.y + matrix[10] * vec.z +  matrix[11];

	return ret;
}

float novaCTF::sign(float val)
{
	if (val > 0.0) return 1.0;
	if (val < 0.0) return -1.0;

	return 0.0;
}

std::string novaCTF::convertSeconds(double seconds_total)
{
	unsigned int hours = ((unsigned int)seconds_total)/3600;
	unsigned int minutes = (((int)seconds_total)%3600)/60;
	unsigned int seconds = (((int)seconds_total)%3600)%60;
	//timeinfo = localtime(&rawtime);
	std::stringstream timeDifference;


	if(hours>0)
		timeDifference << hours << "h ";

	if(minutes > 0 || hours > 0)
		timeDifference << minutes << "m ";

	timeDifference << seconds << "s";

	//strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
	//std::string stamp(buffer);

	return timeDifference.str();
}
