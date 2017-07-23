#pragma once
#include <fftw3.h>
#include <vector>
#include "common.h"

class FFTRoutines
{
public:
	static void real1DTransform(size_t size, std::vector<float>& data, std::vector<float>& fft);
	static void real2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<double>& mask, novaCTF::DataStats& dataStats);
	static void real2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& filter);
	static void computePowerSpectrum(std::vector<float>& powerSpectrum, fftw_complex* fftOut, size_t sizeX, size_t nyh, bool logarithmizeData);
	static void complex2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& filter);
	static void many1DTransform(float* data, int sizeX, int sizeY, int direction);
private:
	static double logarithmizeValue(double value, bool logarithmizeData);
	static void maskFFT(fftw_complex* fftOut, std::vector<double>& mask, size_t sizeX, size_t nyh);
	static void filterFFT(fftw_complex* fftOut, std::vector<float>& filter, size_t sizeX, size_t nyh);
	static void normalizeValues(std::vector<float>& normalizedData, std::vector<double>& originalData, size_t dataSize, novaCTF::DataStats& dataStats);
	static void normalize(float *array, float scale, size_t dataSize);

};
