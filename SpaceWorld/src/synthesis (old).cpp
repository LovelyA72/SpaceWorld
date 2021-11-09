#include "world.h"

#include <stdio.h> // for debug
#include <stdlib.h>
#include <float.h>

// spectrum, cepstrum�͖���malloc, free����̂��ʓ|������D
/*
int getOneFrameSegment(double *f0, int tLen, double **specgram, double **aperiodicity, int fftl, double framePeriod, double currentTime, int fs, double defaultF0,
						fft_complex *spectrum, fft_complex *cepstrum, 
						double *response, int xLen);
*/

void getMinimumPhaseSpectrum(double *inputSpec, fft_complex *spectrum, fft_complex *cepstrum, int fftl)
{
	int i;
	double real, imag;
	fft_plan forwardFFT, inverseFFT;
	forwardFFT = fft_plan_dft_1d(fftl, spectrum, cepstrum, FFT_FORWARD, FFT_ESTIMATE);
	inverseFFT = fft_plan_dft_1d(fftl, cepstrum, spectrum, FFT_BACKWARD, FFT_ESTIMATE);

	// �l�����o��
	for(i = 0;i <= fftl/2;i++)	
	{
		spectrum[i][0] = log(inputSpec[i]?inputSpec[i]:1.0e-20)/2.0;
		spectrum[i][1] = 0.0;
	}
	for(;i < fftl;i++)
	{
		spectrum[i][0] = spectrum[fftl-i][0];
		spectrum[i][1] = 0.0;
	}
	fft_execute(forwardFFT);
	for(i = 1;i < fftl/2;i++)
	{
		cepstrum[i][0] = 0.0;
		cepstrum[i][1] = 0.0;
	}
	for(;i < fftl;i++)
	{
		cepstrum[i][0] *= 2.0;
		cepstrum[i][1] *= 2.0;
	}
	fft_execute(inverseFFT);
	for(i = 0;i < fftl;i++)
	{
		real = exp(spectrum[i][0]/(double)fftl)*cos(spectrum[i][1]/(double)fftl);
		imag = exp(spectrum[i][0]/(double)fftl)*sin(spectrum[i][1]/(double)fftl);
		spectrum[i][0] = real;
		spectrum[i][1] = imag;
	}
}

// ���莞���̉������擾����D
void getOneFrameSegment(double *f0, int tLen, double **specgram, double **residualSpecgram, int fftl, double framePeriod, double currentTime, int fs, double defaultF0,
						fft_complex *spectrum, fft_complex *cepstrum, 
						double *response, int xLen)
{
	int i;
	double real, imag, tmp;
	fft_plan	inverseFFT_RP;				// FFT�Z�b�g

	int currentFrame, currentPosition;

	inverseFFT_RP = fft_plan_dft_c2r_1d(fftl, spectrum, response ,  FFT_ESTIMATE);

	currentFrame = (int)(currentTime/(framePeriod/1000.0) + 0.5);	
	currentPosition = (int)(currentTime*(double)fs);

	tmp = currentTime + 1.0/(f0[currentFrame] == 0.0 ? defaultF0 : f0[currentFrame]);

	// �l�����o��
	getMinimumPhaseSpectrum(specgram[currentFrame], spectrum, cepstrum, fftl);

	spectrum[0][0] *= residualSpecgram[currentFrame][0];
	for(i = 1;i < fftl/2;i++)
	{
		real = spectrum[i][0]*residualSpecgram[currentFrame][(i-1)*2+1] - spectrum[i][1]*residualSpecgram[currentFrame][i*2];
		imag = spectrum[i][0]*residualSpecgram[currentFrame][i*2] + spectrum[i][1]*residualSpecgram[currentFrame][(i-1)*2+1];
		spectrum[i][0] = real;
		spectrum[i][1] = imag;
	}
	spectrum[fftl/2][0] *= residualSpecgram[currentFrame][fftl-1];
	spectrum[fftl/2][1] = 0;
	fft_execute(inverseFFT_RP);

	fft_destroy_plan(inverseFFT_RP);
}



void synthesis(double *f0, int tLen, double **specgram, double **residualSpecgram, int fftl, double framePeriod, int fs, 
			   double *synthesisOut, int xLen)
{
	int i,j;
	double *impulseResponse;
	impulseResponse = (double *)malloc(sizeof(double) * fftl);
	fft_complex		*cepstrum, *spectrum;	// �P�v�X�g�����ƃX�y�N�g��
	cepstrum = (fft_complex *)malloc(sizeof(fft_complex) * fftl);
	spectrum = (fft_complex *)malloc(sizeof(fft_complex) * fftl);

	double currentTime = 0.0;
	int currentPosition = 0;//currentTime / framePeriod;
	int currentFrame = 0;
	for(i = 0;;i++)
	{
		for(j = 0;j < fftl;j++) impulseResponse[j] = 0.0; // �z��͖��񏉊���

		getOneFrameSegment(f0, tLen, specgram, residualSpecgram, fftl, framePeriod, currentTime, fs, DEFAULT_F0,
						spectrum, cepstrum, impulseResponse, xLen);

		currentPosition = (int)(currentTime*(double)fs);
//		for(j = 0;j < fftl/2;j++)
		for(j = 0;j < 3*fftl/4;j++)
		{
			if(j+currentPosition >= xLen) break;
			synthesisOut[j+currentPosition] += impulseResponse[j];
		}

		// �X�V
		currentTime += 1.0/(f0[currentFrame] == 0.0 ? DEFAULT_F0 : f0[currentFrame]);
		currentFrame = (int)(currentTime/(framePeriod/1000.0) + 0.5);
		currentPosition = (int)(currentTime*(double)fs);
		if(fftl/8+currentPosition >= xLen || currentFrame >= tLen) break;
	}

	free(cepstrum); free(spectrum);
	free(impulseResponse);
	return;
}
