#include "world.h"
#include "world/dio.h"

#include <stdio.h> // for debug
#include <stdlib.h>

int getSamplesForDIO2(int fs, int x_length, double frame_period) {
	return static_cast<int>(1000.0 * x_length / fs / frame_period) + 1;
}
void starGeneralBody(double *x, int xLen, int fs, double f0, double t, int fftl,
					double * sliceSTAR);

// �T���v�����O���g������K�v��FFT����v�Z
// �{���͍��������̍Œ�F0��K�v�����ǖ���
int getFFTLengthForStar(int fs)
{
	return (int)pow(2.0, 1.0+(int)(log(3.0*fs/FLOOR_F0+1) / log(2.0)));
}

// STAR�ɂ��X�y�N�g�������
void star(double *x, int xLen, int fs, double *timeAxis, double *f0,
		 double **specgram)
{
	int i,j;
	double framePeriod = (timeAxis[1]-timeAxis[0])*1000.0;
	double f0LowLimit = FLOOR_F0; // F0��FLOOR Hz�ȉ��̏ꍇ�͖������Ƃ��Ĉ���
	double currentF0;

	int	fftl = (int)pow(2.0, 1.0+(int)(log(3.0*fs/f0LowLimit+1) / log(2.0)));
	int tLen = getSamplesForDIO2(fs, xLen, framePeriod);

	double *sliceSTAR;
	sliceSTAR = (double *)malloc(sizeof(double) * fftl);

	for(i = 0;i < tLen;i++)
	{
		currentF0 = f0[i] <= FLOOR_F0 ? DEFAULT_F0 : f0[i];
		starGeneralBody(x, xLen, fs, currentF0, timeAxis[i], fftl, sliceSTAR);
		for(j = 0;j <= fftl/2;j++) specgram[i][j] = sliceSTAR[j];
	}

	free(sliceSTAR);
}

void starGeneralBody(double *x, int xLen, int fs, double f0, double t, int fftl,
							   double * sliceSTAR)
{
	int i,j;
	double t0 = 1.0 / f0;

	int *baseIndex, *index; // i�t���̂��܂��� (Matlab�ŎQ��)
	int nFragment = (int)(0.5 + 3.0*(double)fs/f0/2.0);
	baseIndex = (int *)malloc(sizeof(int) * (nFragment*2+1));
	index  = (int *)malloc(sizeof(int) * (nFragment*2+1));

	for(i = -nFragment, j = 0;i <= nFragment;i++, j++) 
		baseIndex[j] = i;
	for(i = 0;i <= nFragment*2;i++) 
		index[i]  = min(xLen, max(1, round(t*(double)fs+1+baseIndex[i]) ) ) - 1;

	double *segment, *window;
	double position, average;
	segment  = (double *)malloc(sizeof(double) * (nFragment*2+1));
	window   = (double *)malloc(sizeof(double) * (nFragment*2+1));
	average  = 0.0;
	for(i = 0;i <= nFragment*2;i++)
	{
		segment[i]  = x[index[i]];
		position  = (double)(baseIndex[i]/(double)fs/(3.0/2.0) ) + 
			(t*(double)fs - (double)(round(t*(double)fs))) / (double)fs;
		window[i]  = 0.5*cos(PI*position*f0) +0.5;
		average  += window[i]*window[i];
	}
	average  = sqrt(average);
	for(i = 0;i <= nFragment*2;i++) window[i]  /= average;

	// �g�`�̃X�y�N�g����v�Z
	double				*waveform;
	double				*powerSpec;
	waveform  = (double *)malloc(sizeof(double) * fftl);
	powerSpec = (double *)malloc(sizeof(double) * fftl);

	fft_plan			forwardFFT;	// FFT�Z�b�g
	fft_complex		*ySpec;		// �X�y�N�g��
	ySpec = (fft_complex *)malloc(sizeof(fft_complex) * fftl);
	forwardFFT = fft_plan_dft_r2c_1d(fftl, waveform, ySpec, FFT_ESTIMATE);

	// �p���[�X�y�N�g���̌v�Z
	for(i = 0;i <= nFragment*2;i++) 
		waveform[i] = segment[i] * window[i];
	for(;i < fftl;i++) 
		waveform[i] = 0.0;
	fft_execute(forwardFFT); // FFT�̎��s
	for(i = 1;i <= fftl/2;i++) 
		powerSpec[i] = ySpec[i][0]*ySpec[i][0] + ySpec[i][1]*ySpec[i][1];
	powerSpec[0] = powerSpec[1];

	free(ySpec);
	fft_destroy_plan(forwardFFT);
	free(segment);
	free(window);
	free(baseIndex); free(index);

	// adroit smoothing
	// ���ʂ�d�����D
	int limit;
	limit = (int)(f0 / (fs/(double)fftl))+1;
	double *dSpectrum, dFrequencyAxis, dShift;
	dSpectrum		= (double *)malloc(sizeof(double) * (fftl+limit*2 + 1) );

	// �v�Z�R�X�g������ł���炷
	dFrequencyAxis = -((double)limit-0.5)*(double)fs/(double)fftl;
	dShift = (double)fs/(double)fftl;

	for(i = 0;i < limit;i++) 
		dSpectrum[i] = powerSpec[limit-i];
	for(j = 0;i < fftl/2+limit;i++,j++) 
		dSpectrum[i] = powerSpec[j];
	for(j=1;i < fftl/2+limit*2+1;i++,j++) 
		dSpectrum[i] = powerSpec[fftl/2-j];

	int tmp = (int)(f0*fftl/(double)fs);
	double *dSegment, *centers;
	dSegment	= (double *)malloc(sizeof(double) * fftl*2);
	centers		= (double *)malloc(sizeof(double) * (fftl/2 + 1) );
	
	dSegment[0]	= log(dSpectrum[0]?dSpectrum[0]:1.0e-20)*(double)fs/(double)fftl;
	for(i = 1;i < fftl/2+limit*2+1;i++)
		dSegment[i] = log(dSpectrum[i]?dSpectrum[i]:1.0e-20)*(double)fs/(double)fftl + dSegment[i-1];

	for(i = 0;i <= fftl/2;i++)
		centers[i] = (double)i / (double)fftl * (double)fs - f0/2.0;

	double *lowLevels, *highLevels;
	lowLevels  = (double *)malloc(sizeof(double) * (fftl/2+1));
	highLevels = (double *)malloc(sizeof(double) * (fftl/2+1));
	interp1Q(dFrequencyAxis, dShift, dSegment, fftl/2+limit*2+1, centers, (fftl/2 + 1), lowLevels);
	for(i = 0;i <= fftl/2;i++)
		centers[i] += f0;
	interp1Q(dFrequencyAxis, dShift, dSegment, fftl/2+limit*2+1, centers, (fftl/2 + 1), highLevels);

	for(i = 0;i <= fftl/2;i++)
		sliceSTAR[i] = exp( (highLevels[i]-lowLevels[i])/f0);

	free(lowLevels); free(highLevels);
	free(dSegment); free(centers);
	free(dSpectrum);
	free(waveform);
	free(powerSpec); 
}
