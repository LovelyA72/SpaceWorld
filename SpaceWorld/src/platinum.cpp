#include "world.h"
#include "world/dio.h"

#include <stdio.h> // for debug
#include <stdlib.h>
#include <math.h>

void getOneFrameResidualSpec(double* x, int xLen, int fs, int positionIndex, double framePeriod, double f0, double* specgram, int fftl, double* pulseLocations, int pCount,
	double* residualSpec)
{
	int i;
	double* tmpWave;
	double T0;
	int index, tmpIndex, wLen;

	double tmp;
	double tmpValue = 100000.0; // safeGuard
	for (i = 0;i < pCount;i++)
	{
		tmp = fabs(pulseLocations[i] - (double)positionIndex * framePeriod);
		if (tmp < tmpValue)
		{
			tmpValue = tmp;
			tmpIndex = i;
		}
		index = 1 + (int)(0.5 + pulseLocations[tmpIndex] * fs);
	}

	T0 = (double)fs / f0;
	wLen = (int)(0.5 + T0 * 2.0);
	tmpWave = (double*)malloc(sizeof(double) * fftl);

	if (wLen + index - (int)(0.5 + T0) >= xLen)
	{
		for (i = 0;i < fftl;i++) tmpWave[i] = 0;
		free(tmpWave);
		return;
	}

	for (i = 0;i < wLen;i++)
	{
		tmpIndex = i + index - (int)(0.5 + T0);
		tmpWave[i] = x[min(xLen - 1, max(0, tmpIndex))] *
			(0.5 - 0.5 * cos(2.0 * PI * (double)(i + 1) / ((double)(wLen + 1))));
	}
	for (;i < fftl;i++)
		tmpWave[i] = 0.0;

	fft_plan			forwardFFT;	// FFTセット
	fft_complex* tmpSpec, * starSpec, * ceps;		// スペクトル
	tmpSpec = (fft_complex*)malloc(sizeof(fft_complex) * fftl);
	starSpec = (fft_complex*)malloc(sizeof(fft_complex) * fftl);
	ceps = (fft_complex*)malloc(sizeof(fft_complex) * fftl);

	//getMinimumPhaseSpectrum(specgram, starSpec, ceps, fftl);

	forwardFFT = fft_plan_dft_r2c_1d(fftl, tmpWave, tmpSpec, FFT_ESTIMATE);
	fft_execute(forwardFFT);

	residualSpec[0] = tmpSpec[0][0] / starSpec[0][0];
	for (i = 0;i < fftl / 2 - 1;i++)
	{
		tmp = starSpec[i + 1][0] * starSpec[i + 1][0] + starSpec[i + 1][1] * starSpec[i + 1][1];
		residualSpec[i * 2 + 1] = (starSpec[i + 1][0] * tmpSpec[i + 1][0] + starSpec[i + 1][1] * tmpSpec[i + 1][1]) / tmp;
		residualSpec[i * 2 + 2] = (-starSpec[i + 1][1] * tmpSpec[i + 1][0] + starSpec[i + 1][0] * tmpSpec[i + 1][1]) / tmp;
	}
	residualSpec[fftl - 1] = tmpSpec[fftl / 2][0] / starSpec[fftl / 2][0];

	free(tmpSpec); free(ceps); free(starSpec);
	free(tmpWave);
}

int getPulseLocations(double* x, int xLen, double* totalPhase, int vuvNum, int* stList, int* edList, int fs, double framePeriod, int* wedgeList, double* pulseLocations)
{
	int i, j;
	int stIndex, edIndex;

	int pCount = 0;
	int numberOfLocation;
	double* tmpPulseLocations, * basePhase;
	tmpPulseLocations = (double*)malloc(sizeof(double) * xLen);
	basePhase = (double*)malloc(sizeof(double) * xLen);

	double tmp;
	for (i = 0;i < vuvNum;i++)
	{
		stIndex = max(0, (int)((double)fs * (stList[i]) * framePeriod / 1000.0));  //島の先頭のサンプル位置
		edIndex = min(xLen - 1, (int)((double)fs * (edList[i] + 1) * framePeriod / 1000.0 + 0.5) - 1);//島の末尾のサンプル位置

		tmp = totalPhase[wedgeList[i]];

		for (j = stIndex;j < edIndex;j++) {
			basePhase[j] = fmod(totalPhase[j + 1] - tmp, 2 * PI) - fmod(totalPhase[j] - tmp, 2 * PI);//フレーム間の位相差？
//			basePhase[j] = fmod(totalPhase[j]-tmp+PI*0.5, 2*PI);  //tn_fnds 各サンプルの位相を補正
		}

		basePhase[0] = 0; numberOfLocation = 0;
		//		numberOfLocation = 0;
		for (j = stIndex;j < edIndex - 1;j++)
		{
			if (fabs(basePhase[j]) > PI / 1.0)
				//			if(basePhase[j+1] < basePhase[j])  //位相が2*PIを超えて0に戻った
			{
				tmpPulseLocations[numberOfLocation++] = (double)j / (double)fs;//ゼロクロス位置の時刻
			}
		}
		for (j = 0;j < numberOfLocation;j++) pulseLocations[pCount++] = tmpPulseLocations[j];//
	}

	free(tmpPulseLocations);
	free(basePhase);
	return pCount;
}

void getWedgeList(double* x, int xLen, int vuvNum, int* stList, int* edList, int fs, double framePeriod, double* f0, int* wedgeList)
{
	int i, j;
	double LowestF0 = 40.0;
	int center, T0;
	double peak;
	int peakIndex = 0;
	double* tmpWav;
	double currentF0;
	tmpWav = (double*)malloc(sizeof(double) * (int)(fs * 2 / LowestF0));

	for (i = 0;i < vuvNum;i++)
	{
		center = (int)((stList[i] + edList[i] + 1) / 2);           //島の中央のフレーム位置
		currentF0 = f0[center] == 0.0 ? DEFAULT_F0 : f0[center];//島の中央のF0 ノイズ領域の場合はデフォルト
		T0 = (int)((fs / currentF0) + 0.5);                //島の中央のF0の１周期のサンプル数
//		peakIndex = (int)(((1+center)*framePeriod*fs/1000.0)+0.5);//島の中央のサンプル位置
		peakIndex = (int)(((center)*framePeriod * fs / 1000.0) + 0.5);//島の中央のサンプル位置
//		for(j = 0;j < T0*2;j++)
		for (j = 0;j < T0 * 2 + 1;j++)
		{
			//			tmpWav[j] = x[peakIndex-T0+j-1];
			tmpWav[j] = x[max(0, min(xLen - 1, peakIndex - T0 + j - 1))];//島の中央の２周期分のデータ
		}
		peak = 0.0;
		peakIndex = 0;
		for (j = 0;j < T0 * 2 + 1;j++)//波形のピークのサンプル位置を検出
		{
			if (fabs(tmpWav[j]) > peak)
			{
				peak = tmpWav[j];
				peakIndex = j;
			}
		}
		//		wedgeList[i] = max(0, min(xLen-1, (int)(0.5 + ((center+1)*framePeriod*fs/1000.0)-T0+peakIndex+1.0) - 1));//島の中央のフレームのピークのサンプル位置
		wedgeList[i] = max(0, min(xLen - 1, (int)(0.5 + ((center)*framePeriod * fs / 1000.0) - T0 + peakIndex + 1.0) - 1));//島の中央のフレームのピークのサンプル位置
	}
	free(tmpWav);
}

// PLATINUM Version 0.0.4. 恐らくこの仕様で確定です．
// Aperiodicity estimation based on PLATINUM

void platinum(double* x, int xLen, int fs, double* timeAxis, double* f0, double** specgram,
	double** residualSpecgram)
{
	int i, j, index;
	double framePeriod = (timeAxis[1] - timeAxis[0]) * 1000.0;

	int	fftl = (int)pow(2.0, 1.0 + (int)(log(3.0 * fs / FLOOR_F0 + 1) / log(2.0)));
	int tLen = GetSamplesForDIO(fs, xLen, framePeriod);

	int vuvNum;
	//	vuvNum = 0;
	vuvNum = 1;	//tn_fuds
	for (i = 1;i < tLen;i++)
	{
		if (f0[i] != 0.0 && f0[i - 1] == 0.0) vuvNum++;	//無声→有声
		if (f0[i] == 0.0 && f0[i - 1] != 0.0) vuvNum++;	//有声→無声  tn_fnds
	}
	//	vuvNum+=vuvNum-1; // 島数の調整 (有声島と無声島)  tn_fnds コメントアウト
	//	if(f0[0] == 0) vuvNum++;  tn_fnds コメントアウト
	//	if(f0[tLen-1] == 0) vuvNum++;  tn_fnds コメントアウト

	int stCount, edCount;
	int* stList, * edList;
	stList = (int*)malloc(sizeof(int) * vuvNum);
	edList = (int*)malloc(sizeof(int) * vuvNum);
	edCount = 0;

	stList[0] = 0;
	stCount = 1;
	index = 1;
	if (f0[0] != 0)
	{
		for (i = 1;i < tLen;i++)
		{
			if (f0[i] == 0 && f0[i - 1] != 0)
			{
				edList[0] = i - 1;
				edCount++;
				stList[1] = i;
				stCount++;
				index = i;

				break;
			}
		}
	}

	edList[vuvNum - 1] = tLen - 1;
	for (i = index;i < tLen;i++)
	{
		if (f0[i] != 0.0 && f0[i - 1] == 0.0)
		{
			edList[edCount++] = i - 1;
			stList[stCount++] = i;
		}
		if (f0[i] == 0.0 && f0[i - 1] != 0.0)
		{
			edList[edCount++] = i - 1;
			stList[stCount++] = i;
		}
	}

	int* wedgeList;
	wedgeList = (int*)malloc(sizeof(int) * vuvNum);
	getWedgeList(x, xLen, vuvNum, stList, edList, fs, framePeriod, f0, wedgeList);

	double* signalTime, * f0interpolatedRaw, * totalPhase;
	double* fixedF0;
	fixedF0 = (double*)malloc(sizeof(double) * tLen);
	signalTime = (double*)malloc(sizeof(double) * xLen);
	f0interpolatedRaw = (double*)malloc(sizeof(double) * xLen);
	totalPhase = (double*)malloc(sizeof(double) * xLen);

	for (i = 0;i < tLen;i++) fixedF0[i] = f0[i] == 0 ? DEFAULT_F0 : f0[i];
	for (i = 0;i < xLen;i++) signalTime[i] = (double)i / (double)fs;
	interp1(timeAxis, fixedF0, tLen, signalTime, xLen, f0interpolatedRaw);
	totalPhase[0] = f0interpolatedRaw[0] * 2 * PI / (double)fs;
	for (i = 1;i < xLen;i++) totalPhase[i] = totalPhase[i - 1] + f0interpolatedRaw[i] * 2 * PI / (double)fs;

	double* pulseLocations;
	pulseLocations = (double*)malloc(sizeof(double) * xLen);
	int pCount;
	pCount = getPulseLocations(x, xLen, totalPhase, vuvNum, stList, edList, fs, framePeriod, wedgeList, pulseLocations);

	double* tmpResidualSpec;
	tmpResidualSpec = (double*)malloc(sizeof(double) * fftl);
	double currentF0;
	for (j = 0;j < fftl;j++) residualSpecgram[0][j] = 0.0;
	for (i = 1;i < tLen;i++)
	{
		currentF0 = f0[i] <= FLOOR_F0 ? DEFAULT_F0 : f0[i];
		getOneFrameResidualSpec(x, xLen, fs, i, framePeriod / 1000.0, currentF0, specgram[i], fftl, pulseLocations, pCount,
			tmpResidualSpec);
		for (j = 0;j < fftl;j++) residualSpecgram[i][j] = tmpResidualSpec[j];
	}

	free(fixedF0);
	free(tmpResidualSpec);
	free(pulseLocations);
	free(totalPhase); free(f0interpolatedRaw); free(signalTime);
	free(wedgeList);
	free(edList); free(stList);
	return;
}