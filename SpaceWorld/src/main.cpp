/*
 * ※このプログラムはGPL (General Public License)でライセンスされています．
 * ライセンスの内容は，copying.txtをご確認ください．
 * 以下、WORLD(0.0.4)のreadme.txtより引用。
 * > 注意すべきは以下の2点です．
 * > ・本プログラムを利用して作ったアプリケーションは，
 * > 　ソースコードも一緒に配布する必要がある．
 * > ・そのアプリケーションもGPLライセンスにしなければならない．
 */
 //
 // world4utau.cpp
 //
 //このソースコードはWORLD0.0.4のtest.cppを元にしています
 //引数
 // 1 入力ファイル（OK）
 // 2 出力ファイル（OK）
 // 3 音階（OK）
 // 4 タイムパーセント→Velocity
 // 5 フラグ　t,g,P 実装
 // 6 オフセット
 // 7 長さ：要求長
 // 8 前半の固定部分 fixed
 // 9 最後の使用しない部分。マイナスの場合offsetからの使用部分 blank
 // 10 ボリューム (OK) volm
 // 11 モジュレーション (OK) modulation
 // 12~ ピッチベンド  pitches


#include <stdio.h>
#include <stdlib.h>

#include <windows.h>

#include "world.h"
#include "wavread.h"

#include <math.h>
#include <float.h>
#include <memory.h>
#include "world/stonemask.h"
#include "world/dio.h"
#include "world/cheaptrick.h"
#include "world/d4c.h"
#include "world/synthesis.h"
#include "getWorldValues77.h"
#include "matlabmyfunctions.h"
#include "audio_io.h"

#include <math.h>
#include <string>
using namespace std;


// 分析シフト量 [msec]
#define FRAMEPERIOD (1000.0*256/44100)
//#define FRAMEPERIOD 5.80498866213152

#pragma comment(lib, "winmm.lib")

int get64(int c)
{
	if (c >= '0' && c <= '9')
	{
		return c - '0' + 52;
	}
	else if (c >= 'A' && c <= 'Z')
	{
		return c - 'A';
	}
	else if (c >= 'a' && c <= 'z')
	{
		return c - 'a' + 26;
	}
	else if (c == '+')
	{
		return 62;
	}
	else if (c == '/')
	{
		return 63;
	}
	else
	{
		return 0;
	}
}

int decpit(char* str, int* dst, int cnt)
{
	int len = 0;
	int i, n = 0;
	int k = 0, num, ii;
	if (str != NULL)
	{
		len = strlen(str);
		for (i = 0; i < len; i += 2)
		{
			if (str[i] == '#')
			{
				i++;
				sscanf(str + i, "%d", &num);
				for (ii = 0; ii < num && k < cnt; ii++) {
					dst[k++] = n;
				}
				while (str[i] != '#' && str[i] != 0) i++;
				i--;
			}
			else
			{
				n = get64(str[i]) * 64 + get64(str[i + 1]);
				if (n > 2047) n -= 4096;
				if (k < cnt) {
					dst[k++] = n;
				}
			}
		}
	}
	return len;
}

double name2freq(char* tname, int trim)
{
	char c;
	int n, m, oct, num;
	//01234567890A
	//C D EF G A B
	c = tname[0];
	if (c >= 'A' && c <= 'G')
	{
		if (c <= 'B')
		{
			n = 9 + (c - 'A') * 2;
		}
		else if (c <= 'E')
		{
			n = (c - 'C') * 2;
		}
		else
		{
			n = 5 + (c - 'F') * 2;
		}

		c = tname[1];

		m = 2;
		if (c == '#') {
			n++;
		}
		else if (c == 'b')
		{
			n--;
		}
		else
		{
			m = 1;
		}

		if (tname[m] == 0)
		{
			return 0;
		}

		sscanf(tname + m, "%d", &oct);

		num = (n + oct * 12 - 21) * 10 + trim;
		//return num;
		//0 で 55Hz
		return 55.0 * pow(2, (double)num / 120.0);
	}
	return 0;
}
//
void makeFilename(const char* filename, const char* ext, char* output)
{
	strcpy(output, filename);
	char* cp = strrchr(output, '.');
	if (cp) *cp = 0;
	strcat(output, ext);
}
//TODO 分析ファイルが揃ってたらwaveは読まなくて良い。（今は必ず読んでる）

double getFreqAvg(double f0[], int tLen)
{
	int i, j;
	double value = 0, r;
	double p[6], q;
	double freq_avg = 0;
	double base_value = 0;
	for (i = 0; i < tLen; i++)
	{
		value = f0[i];
		if (value < 1000.0 && value > 55.0)
		{
			r = 1.0;
			//連続して近い値の場合のウエイトを重くする
			for (j = 0; j <= 5; j++)
			{
				if (i > j) {
					q = f0[i - j - 1] - value;
					p[j] = value / (value + q * q);
				}
				else {
					p[j] = 1 / (1 + value);
				}
				r *= p[j];
			}
			freq_avg += value * r;
			base_value += r;
		}
	}
	if (base_value > 0) freq_avg /= base_value;
	return freq_avg;
}

double getUTAUfrq(char* fname, int dioFrames, double fs, double framePeriod,
	double* f0)
{
	// 入力ファイル名作成
	char frqfname[512];
	strcpy(frqfname, fname);
	strcat(frqfname, ".frq");

	// check *.aiff.frq
	FILE* fp = fopen(frqfname, "rb");
	if (NULL == fp) {
		strcpy(frqfname, fname);
		char* cptr = strrchr(frqfname, '.');
		// '*.wav' => '*_wav'
		if (NULL != cptr) *cptr = '_';
		strcat(frqfname, ".frq");

		// check *_wav.frq
		fp = fopen(frqfname, "rb");
		if (NULL == fp) {
			printf("There is no UTAU frequency file.\n");
			return 0.0;
		}
	}

	// ヘッダ確認
	char tempStr[128];
	fread(tempStr, sizeof(char), 4, fp); // "FREQ"
	if (0 != strncmp(tempStr, "FREQ", 4)) {
		fclose(fp);
		printf("Header Error \"FREQ\".\n");
		return 0.0;
	}

	fread(tempStr, sizeof(char), 4, fp); // "0003"
	if (0 != strncmp(tempStr, "0003", 4)) {
		fclose(fp);
		printf("Header Error version.\n");
		return 0.0;
	}

	int pitchStep;
	fread(&pitchStep, sizeof(int), 1, fp); // デフォールトは256

	double avrFo;
	fread(&avrFo, sizeof(double), 1, fp); // 平均周波数

	fread(tempStr, sizeof(char), 8 + 8, fp); // よく判らない、飛ばしてるだけ

	int utauFrames;
	fread(&utauFrames, sizeof(int), 1, fp); // データ数

	double* utauF0 = new double[utauFrames];
	for (int i = 0; i < utauFrames; i++) {
		fread(&utauF0[i], sizeof(double), 1, fp); // f0
		fread(tempStr, sizeof(double), 1, fp); // 振幅は読み飛ばす
	}
	fclose(fp);
#if RUNTIME_CHECK
	for (int i = 0;i < utauFrames;i++) {
		if (0 == _finite(utauF0[i])) {
			printf("in %s\n", __func__);
			printf("utauF0[%d] error!! inf or NaN !!\n", i);
			getchar();
		}
	}
#endif

	// f0列整形
	// f0が一定範囲外は0にする
	for (int i = 0; i < utauFrames; i++) {
		if ((utauF0[i] < FLOOR_F0) || (CEIL_F0 < utauF0[i]))
			utauF0[i] = 0.0;
	}


	// utauF0列の有声部分を検出する
	int* vIndex = new int[utauFrames];
	int* uIndex = new int[utauFrames];
	int vCnt, uCnt;
	vCnt = uCnt = 0;
	if (FLOOR_F0 <= utauF0[0]) vIndex[vCnt++] = 0;
	// 必ずvIndex[X],uIndex[X]をペアにする
	for (int i = 1; i < utauFrames; i++) {
		// 無声 => 有声（無声島終わり）
		if ((utauF0[i - 1] < FLOOR_F0) && (FLOOR_F0 <= utauF0[i]))
			vIndex[vCnt++] = i;
		// 有声 => 無声（有声島終わり）
		else if ((FLOOR_F0 <= utauF0[i - 1]) && (utauF0[i] < FLOOR_F0))
			uIndex[uCnt++] = i - 1;
	}
	// 無声で終わってるのは構わない
	if (uCnt < vCnt) uIndex[uCnt++] = utauFrames - 1;


	// 出力f0列クリア
	for (int i = 0; i < dioFrames; i++) f0[i] = 0.0;
	double utauFramePeriod, stride;
	utauFramePeriod = Sample2Time(pitchStep, fs); // [msec]
	stride = framePeriod / utauFramePeriod;


	// f0が4個しか連続しない（23msecのはず）ところは、子音とみなして0にする
	// 4個以下の連続なら：　dioF0列を触らないで0のままにする
	// 4個 超 の連続なら：　dioF0列をutauF0列から補間して作成する
	const int minLen = 4;
	for (int i = 0; i < vCnt; i++) {
		// 小さな連続は、子音とみなしてスキップ
		int LenU = uIndex[i] - vIndex[i] + 1; // 間の数ではなくて、サンプル数なので +1
		if (minLen < LenU) {
			// 有声部分は必ず周波数表より狭くしてたが、最近傍の方がよさげ
			double stUtime, edUtime;
			stUtime = Frame2Time(vIndex[i], utauFramePeriod);
			edUtime = Frame2Time(uIndex[i], utauFramePeriod);
			int stD, edD, LenD;
			stD = static_cast<int>(Time2Frame(stUtime, framePeriod) + 0.5);
			edD = static_cast<int>(Time2Frame(edUtime, framePeriod) + 0.5);
			if (fmod(edUtime, framePeriod) <= 0) edD--;
			if (dioFrames <= edD) edD = dioFrames - 1;
			LenD = edD - stD + 1; // 間の数ではなくて、サンプル数なので +1
			double stDtime, edDtime;
			stDtime = Frame2Time(stD, framePeriod);
			edDtime = Frame2Time(edD, framePeriod);
			double offset = Time2Frame(stDtime - stUtime, utauFramePeriod);

			// 区間外はクリップの高速trimスプライン補間
			itrp1Qtrim_clip(offset, &utauF0[vIndex[i]], LenU,
				stride, LenD, &f0[stD]);
		}
	}
	// 念のため、再度f0が一定範囲外は0にする
	for (int i = 0; i < dioFrames; i++) {
		if ((f0[i] < FLOOR_F0) || (CEIL_F0 < f0[i]))
			f0[i] = 0.0;
		else if (f0[i] < avrFo * 0.55)
			f0[i] = f0[i] * 2;
		else if (f0[i] > avrFo * 1.65)
			f0[i] = f0[i] / 2;
	}

	// メモリ解放
	delete[] vIndex; delete[] uIndex; delete[] utauF0;
	return avrFo;
} // getUTAUfrq

double getavgUTAUfrq(char* fname)
{
	// 入力ファイル名作成
	char frqfname[512];
	strcpy(frqfname, fname);
	strcat(frqfname, ".frq");

	// check *.aiff.frq
	FILE* fp = fopen(frqfname, "rb");
	if (NULL == fp) {
		strcpy(frqfname, fname);
		char* cptr = strrchr(frqfname, '.');
		// '*.wav' => '*_wav'
		if (NULL != cptr) *cptr = '_';
		strcat(frqfname, ".frq");

		// check *_wav.frq
		fp = fopen(frqfname, "rb");
		if (NULL == fp) {
			printf("There is no UTAU frequency file.\n");
			printf("Please generate .frq files using resampler.exe or fresamp. \n");
			return 0.0;
		}
	}

	// ヘッダ確認
	char tempStr[128];
	fread(tempStr, sizeof(char), 4, fp); // "FREQ"
	if (0 != strncmp(tempStr, "FREQ", 4)) {
		fclose(fp);
		printf("Header Error \"FREQ\".\n");
		return 0.0;
	}

	fread(tempStr, sizeof(char), 4, fp); // "0003"
	if (0 != strncmp(tempStr, "0003", 4)) {
		fclose(fp);
		printf("Header Error Version Problem.\n");
		return 0.0;
	}
	int pitchStep;
	fread(&pitchStep, sizeof(int), 1, fp); // デフォールトは256

	double avrFo;
	fread(&avrFo, sizeof(double), 1, fp); // 平均周波数

	fclose(fp);
	return avrFo;
}

double** readCTParam(int signalLen, int fs, const char* filename, int tLen, int fftl)
{
	int i, j;
	char fname2[512];
	DWORD elapsedTime;
	unsigned short tn = 0;
	unsigned short us = 0;
	int siglen = 0;
	int rate = 0;
	double** specgram = 0;

	makeFilename(filename, ".ctspec", fname2);

	printf("read .ctspec:");
	elapsedTime = timeGetTime();

	FILE* fp = fopen(fname2, "rb");
	if (fp)
	{
		char st[9];
		fread(st, 1, 8, fp);
		if (strncmp(st, "world-ct", 8) != 0)
		{
			fclose(fp);
			printf(" bad file.\n");
			return 0;
		}
		fread(&siglen, sizeof(int), 1, fp);
		fread(&rate, sizeof(int), 1, fp);
		fread(&tn, sizeof(unsigned short), 1, fp);
		fread(&us, sizeof(unsigned short), 1, fp);
		if (tn == tLen && us == (fftl / 2 + 1) && signalLen == siglen && fs == rate)
		{
			specgram = (double**)malloc(tLen * sizeof(double*));
			if (specgram)
			{
				for (i = 0; i < tLen; i++)
				{
					specgram[i] = (double*)malloc((fftl / 2 + 1) * sizeof(double));
					memset(specgram[i], 0, (fftl / 2 + 1) * sizeof(double));
					if (specgram[i])
					{
						for (j = 0; j <= fftl / 2; j++)
						{
							unsigned short v;
							fread(&v, sizeof(unsigned short), 1, fp);
							//= (unsigned short)(log(specgram[i][j]*(2048.0*2048*2048)+1) * 512.0 + 0.5);
							specgram[i][j] = (exp(v / 1024.0) - 1) * 1.16415321826935E-10;// /(2048.0*2048*2048);
						}
					}
					else
					{
						break;
					}
				}
				if (i < tLen)
				{
					for (j = 0; j < i; j++)
					{
						free(specgram[i]);
					}
					free(specgram);
					specgram = 0;
					fprintf(stderr, " メモリーが確保できません。%d\n", i);
				}
			}
			else
			{
				fprintf(stderr, " メモリーが確保できません。\n");
			}
		}
		else
		{
			tn = 0;
		}
		fclose(fp);
	}
	printf(" %d [msec]\n", timeGetTime() - elapsedTime);
	return specgram;
}
double** getCTParam(double x[], int signalLen, int fs, double t[], double f0[], int tLen, int fftl)
{
	printf("CHEAPTRICK:");
	DWORD elapsedTime = timeGetTime();

	double** specgram = (double**)malloc(sizeof(double*) * tLen);
	if (specgram)
	{
		int i, j;
		for (i = 0;i < tLen;i++)
		{
			specgram[i] = (double*)malloc(sizeof(double) * (fftl / 2 + 1));
			memset(specgram[i], 0, sizeof(double) * (fftl / 2 + 1));
			if (specgram[i])
			{
				memset(specgram[i], 0, sizeof(double) * (fftl / 2 + 1));
			}
			else
			{
				break;
			}
		}
		if (i == tLen)
		{
			CheapTrick(x, signalLen, fs, t, f0, tLen, specgram, fftl);
		}
		else
		{
			for (j = 0; j < i; j++)
			{
				free(specgram[i]);
			}
			free(specgram);
			specgram = 0;
			fprintf(stderr, " メモリーが確保できません。%d\n", i);
		}
	}
	else
	{
		fprintf(stderr, " メモリーが確保できません。\n");
	}
	printf(" %d [msec]\n", timeGetTime() - elapsedTime);
	return specgram;
}
void writeCTParam(int signalLen, int fs, const char* filename, double* specgram[], int tLen, int fftl)
{
	unsigned short tn;
	unsigned short us;
	int i, j;
	char fname2[512];
	makeFilename(filename, ".ctspec", fname2);


	printf("write .ctspec:");
	DWORD elapsedTime = timeGetTime();

	short max = -32767, min = 32767;
	//FILE *ft = fopen("star0.txt", "wt");
	FILE* f1 = fopen(fname2, "wb");
	if (f1)
	{
		fwrite("world-ct", 1, 8, f1);
		tn = (unsigned short)tLen;
		us = (unsigned short)fftl / 2 + 1;
		fwrite(&signalLen, sizeof(int), 1, f1);
		fwrite(&fs, sizeof(int), 1, f1);
		fwrite(&tn, sizeof(unsigned short), 1, f1);
		fwrite(&us, sizeof(unsigned short), 1, f1);
		for (i = 0; i < tLen; i++)
		{
			for (j = 0; j <= fftl / 2; j++)
			{
				int un;
				if ((un = _fpclass(specgram[i][j]) & 0x0087) != 0)
				{
					specgram[i][j] = 0;
#ifdef _DEBUG
					printf("un[%d][%d]=%04x!\n", i, j, un);
#endif
				}
				unsigned short v = (unsigned short)(log(specgram[i][j] * (2048.0 * 2048 * 2048) + 1) * 1024.0 + 0.5);
				fwrite(&v, sizeof(unsigned short), 1, f1);
				//fprintf(ft, "%0.9lf\t", specgram[i][j]*1000000.0);
				if (max < v) max = v;
				if (min > v) min = v;
			}
			//fprintf(ft, "\n");
		}
		fclose(f1);
	}
	//fclose(ft);
	printf(" %d [msec]\n", timeGetTime() - elapsedTime);
	printf("max = %d, min = %d\n", max, min);
}
double** readD4CParam(int signalLen, int fs, const char* filename, int tLen, int fftl)
{
	int i, j;

	DWORD elapsedTime;
	unsigned short tn = 0;
	unsigned short us = 0;
	int siglen = 0;
	int rate = 0;
	double** residualSpecgram = 0;

	char fname3[512];
	makeFilename(filename, ".d4c", fname3);

	printf("read .d4c:");
	elapsedTime = timeGetTime();

	FILE* fp = fopen(fname3, "rb");
	if (fp)
	{
		char b[9];
		fread(b, 1, 8, fp);
		if (strncmp(b, "wrld-d4c", 8) != 0)
		{
			fclose(fp);
			printf(" bad file.\n");
			return 0;
		}
		fread(&siglen, sizeof(int), 1, fp);
		fread(&rate, sizeof(int), 1, fp);
		fread(&tn, sizeof(unsigned short), 1, fp);
		fread(&us, sizeof(unsigned short), 1, fp);
		if (tn == tLen && (fftl / 2 + 1) == us && signalLen == siglen && fs == rate)
		{
			residualSpecgram = (double**)malloc(tLen * sizeof(double*));
			if (residualSpecgram)
			{
				for (i = 0; i < tLen; i++)
				{
					residualSpecgram[i] = (double*)malloc((fftl / 2 + 1) * sizeof(double));
					if (residualSpecgram[i])
					{
						memset(residualSpecgram[i], 0, (fftl / 2 + 1) * sizeof(double));
						for (j = 0; j <= fftl / 2; j++)
						{
							short v;
							fread(&v, sizeof(short), 1, fp);
							residualSpecgram[i][j] = v * 3.90625E-03;// /256.0;
						}
					}
					else
					{
						break;
					}
				}
				if (i < tLen)
				{
					for (j = 0; j < i; j++)
					{
						free(residualSpecgram[i]);
					}
					free(residualSpecgram);
					residualSpecgram = 0;
					fprintf(stderr, " メモリーが確保できません。%d\n", i);
				}
			}
			else
			{
				fprintf(stderr, " メモリーが確保できません。\n");
			}
		}
		fclose(fp);
	}
	printf(" %d [msec]\n", timeGetTime() - elapsedTime);
	return residualSpecgram;
}
double** getD4CParam(double x[], int signalLen, int fs, double t[], double f0[], int tLen, int fftl)
{
	printf("D4C:");
	DWORD elapsedTime = timeGetTime();

	double** residualSpecgram = (double**)malloc(sizeof(double*) * tLen);
	if (residualSpecgram)
	{
		int i, j;
		for (i = 0;i < tLen;i++)
		{
			residualSpecgram[i] = (double*)malloc(sizeof(double) * (fftl / 2 + 1));
			memset(residualSpecgram[i], 0, sizeof(double) * (fftl / 2 + 1));
			if (residualSpecgram[i])
			{
				memset(residualSpecgram[i], 0, sizeof(double) * (fftl / 2 + 1));
			}
			else
			{
				break;
			}
		}
		if (i == tLen)
		{
			D4C(x, signalLen, fs, t, f0, tLen, fftl, residualSpecgram);
		}
		else
		{
			for (j = 0; j < i; j++)
			{
				free(residualSpecgram[i]);
			}
			free(residualSpecgram);
			residualSpecgram = 0;
			fprintf(stderr, " メモリーが確保できません。%d\n", i);
		}
	}
	else
	{
		fprintf(stderr, " メモリーが確保できません。\n");
	}
	printf(" %d [msec]\n", timeGetTime() - elapsedTime);
	return residualSpecgram;
}
void writeD4CParam(int signalLen, int fs, const char* filename, double* residualSpecgram[], int tLen, int fftl)
{
	int i, j;
	DWORD elapsedTime;
	printf("write .d4c:");
	elapsedTime = timeGetTime();

	unsigned short tn = 0;
	unsigned short us = 0;
	short max = -32767, min = 32767;

	char fname3[512];
	makeFilename(filename, ".d4c", fname3);

	FILE* f1 = fopen(fname3, "wb");
	if (f1)
	{
		fwrite("wrld-d4c", 1, 8, f1);
		tn = (unsigned short)tLen;
		us = (unsigned short)(fftl / 2 + 1);
		fwrite(&signalLen, sizeof(int), 1, f1);
		fwrite(&fs, sizeof(int), 1, f1);
		fwrite(&tn, sizeof(unsigned short), 1, f1);
		fwrite(&us, sizeof(unsigned short), 1, f1);
		for (i = 0; i < tLen; i++)
		{
			for (j = 0; j <= fftl / 2; j++)
			{
				int un;
				if ((un = _fpclass(residualSpecgram[i][j]) & 0x0087) != 0)
				{
					residualSpecgram[i][j] = 0;
#ifdef _DEBUG
					printf("unr[%d][%d]=%04x!\n", i, j, un);
#endif
				}
				short v = (short)(residualSpecgram[i][j] * 256.0);
				//v = log(v * (2048.0*2048.0*2048.0) + 1) * 1024.0;
				if (v > 32767) v = 32767;
				else if (v < -32768) v = -32768;//一応飽和計算しておく（したら不味いけど）
				fwrite(&v, sizeof(short), 1, f1);
				if (max < v) max = v;
				if (min > v) min = v;
			}
		}
		fclose(f1);
	}
	printf("max = %d, min = %d\n", max, min);
	printf(" %d [msec]\n", timeGetTime() - elapsedTime);
}
int readDIOParam(const char* filename, double* p_t[], double* p_f0[], int* p_fs, int* p_siglen)//return tLen
{
	char fname1[512];
	int n = 0;
	int len = 0, sps = 0;
	int i;
	int siglen = 0, fs = 0, tLen = 0;
	double* t = 0;
	double* f0 = 0;

	makeFilename(filename, ".dio", fname1);

	DWORD elapsedTime;

	printf("read .dio:");
	elapsedTime = timeGetTime();

	FILE* fp = fopen(fname1, "rb");
	if (fp)
	{
		char d[9];
		fread(d, 8, 1, fp);
		if (strncmp(d, "wrld-dio", 8) != 0)
		{
			fclose(fp);
			printf(" bad file.\n");
			return 0;
		}
		fread(&siglen, sizeof(int), 1, fp);
		fread(&fs, sizeof(int), 1, fp);
		fread(&tLen, sizeof(int), 1, fp);
		if (tLen > 0)
		{
			t = (double*)malloc(tLen * sizeof(double));
			f0 = (double*)malloc(tLen * sizeof(double));
			if (t && f0)
			{
				for (i = 0; i < tLen; i++)
				{
					fread(&(t[i]), sizeof(double), 1, fp);
					fread(&(f0[i]), sizeof(double), 1, fp);
				}
			}
			else
			{
				if (t) free(t);
				if (f0) free(f0);
				t = 0;
				f0 = 0;
				tLen = 0;
				fprintf(stderr, " メモリーが確保できません。%d\n");
			}
		}
		fclose(fp);
	}
	*p_t = t;
	*p_f0 = f0;
	*p_fs = fs;
	*p_siglen = siglen;
	printf(" %d [msec]\n", timeGetTime() - elapsedTime);
	return tLen;
}
int getDIOParam(double x[], int signalLen, int fs, double framePeriod, double* p_t[], double* p_f0[], char* filename)
{
	DWORD elapsedTime;
	printf("DIO:");
	elapsedTime = timeGetTime();
	int tLen = GetSamplesForDIO(fs, signalLen, framePeriod);
	int i, useutaufrq = 0;
	double utauavgf0, avgf0;
	double* t = (double*)malloc(tLen * sizeof(double));
	double* f0 = (double*)malloc(tLen * sizeof(double));
	double* utauf0 = (double*)malloc(tLen * sizeof(double));
	double* refined_f0 = (double*)malloc(tLen * sizeof(double));
	if (t && f0)
	{
		//	    utauavgf0 = getavgUTAUfrq(filename);
		utauavgf0 = getUTAUfrq(filename, tLen, fs, FRAMEPERIOD, utauf0);
		if (utauavgf0 == 0.0)
		{
			fprintf(stderr, "ERROR: No frq file found!\n");
			//If StAkira can implement my suggestion, this -60 error code will be useful
			exit(-60);
		}
		if (utauavgf0 < FLOOR_F0) return EXIT_FAILURE;
		Dio(x, signalLen, fs, t, f0, 0);
		printf("%d", f0[0]);
		avgf0 = getFreqAvg(f0, tLen);
		if (0.95 * utauavgf0 > avgf0 || avgf0 > 1.05 * utauavgf0)
			avgf0 = utauavgf0;
		for (i = 0;i < tLen;i++)
			if (0.95 * avgf0 > f0[i] || f0[i] > 1.05 * avgf0)
				f0[i] = utauf0[i];
		StoneMask(x, signalLen, fs, t, f0, tLen, refined_f0);
		for (int i = 0; i < tLen; ++i)
			f0[i] = refined_f0[i];
	}
	else
	{
		fprintf(stderr, "Memory cannot be secured. \n");
		if (t) free(t);
		if (f0) free(f0);
		if (refined_f0) free(refined_f0);
		if (utauf0) free(utauf0);
		t = 0;
		f0 = 0;
		utauf0 = 0;
		refined_f0 = 0;
		tLen = 0;
	}
	*p_t = t;
	*p_f0 = f0;
	printf(" %d [msec]\n", timeGetTime() - elapsedTime);
	return tLen;
}
int writeDIOParam(int signalLen, int fs, int tLen, const char* filename, double t[], double f0[])
{
	char fname1[512];
	makeFilename(filename, ".dio", fname1);


	printf("write .dio");

	//FILE *ft = fopen("dio0.txt", "wt");
	FILE* f = fopen(fname1, "wb");
	if (f)
	{
		fwrite("wrld-dio", 1, 8, f);
		fwrite(&signalLen, sizeof(int), 1, f);
		fwrite(&fs, sizeof(int), 1, f);
		fwrite(&tLen, sizeof(int), 1, f);
		int i = 0;
		for (i = 0; i < tLen; i++)
		{
			int un;
			int test = f0[0];
			if ((un = _fpclass(f0[i]) & 0x0087) != 0)//NaN,+Inf,-Inf,denormalを除外する
			{
#ifdef _DEBUG
				printf("un[%d]=%04x!\n", i, un);
#endif
				f0[i] = 0;
			}
			fwrite(&(t[i]), sizeof(double), 1, f);
			fwrite(&(f0[i]), sizeof(double), 1, f);
			//fprintf(ft, "%lf\t%lf\n", t[i], f0[i]);
		}
		fclose(f);
	}
	return 0;
}
void freeSpecgram(double** spec, int n)
{
	int i;
	if (spec && n)
	{
		for (i = 0; i < n; i++)
		{
			free(spec[i]);
		}
		free(spec);
	}
}
// スペクトル伸縮
void stretchSpectrum(double** specgram, int oLen, double ratio, int fs, int fftl)
{
	int i, j;
	double w = 0;
	//ratio = atof(argv[4]);
	if (ratio != 1.0)
	{
		double* freqAxis1, * freqAxis2;
		double* spec1, * spec2;
		freqAxis1 = (double*)malloc(sizeof(double) * fftl);
		freqAxis2 = (double*)malloc(sizeof(double) * fftl);
		spec1 = (double*)malloc(sizeof(double) * fftl);
		spec2 = (double*)malloc(sizeof(double) * fftl);

		// 周波数伸縮の前処理
		for (i = 0;i <= fftl / 2;i++)
		{
			freqAxis1[i] = (double)i / (double)fftl * (double)fs * ratio;
			freqAxis2[i] = (double)i / (double)fftl * (double)fs;
		}
		for (i = 0;i < oLen;i++)
		{
			for (j = 0;j <= fftl / 2;j++)
				spec1[j] = log(specgram[i][j]);
			interp1(freqAxis1, spec1, fftl / 2 + 1, freqAxis2, fftl / 2 + 1, spec2);
			for (j = 0;j <= fftl / 2;j++)
				specgram[i][j] = exp(spec2[j]);
			if (ratio < 1.0)
			{
				for (j = int((double)fftl / 2 * ratio);j <= fftl / 2;j++)
				{
					specgram[i][j] = specgram[i][(int)((double)fftl / 2 * ratio) - 1];
				}
			}
		}

		free(spec1); free(spec2);
		free(freqAxis1); free(freqAxis2);
	}
}
void makeHeader(char* header, int samples, int fs, int nbit)
{
	memcpy(header, "RIFF", 4);
	*(long*)(header + 4) = samples * 2 + 36;
	memcpy(header + 8, "WAVE", 4);
	memcpy(header + 12, "fmt ", 4);
	*(long*)(header + 16) = 16;

	*(short*)(header + 20) = 1;
	*(unsigned short*)(header + 22) = 1;
	*(unsigned long*)(header + 24) = fs;
	*(unsigned long*)(header + 28) = fs * nbit / 8;
	*(unsigned short*)(header + 32) = nbit / 8;
	*(unsigned short*)(header + 34) = nbit;
	memcpy(header + 36, "data", 4);
	*(long*)(header + 40) = samples * 2;
}

//tn_fnds
//イコライジング用スペクトル作成
void createWaveSpec(double* x, int xLen, int fftl, int equLen, fft_complex** waveSpecgram)
{
	int i, j;

	double* waveBuff;
	fft_plan			wave_f_fft;				// fftセット
	fft_complex* waveSpec;	// スペクトル
	waveBuff = (double*)malloc(sizeof(double) * fftl);
	waveSpec = (fft_complex*)malloc(sizeof(fft_complex) * fftl);
	wave_f_fft = fft_plan_dft_r2c_1d(fftl, waveBuff, waveSpec, FFT_ESTIMATE);

	int offset;

	for (i = 0;i < equLen;i++)
	{
		offset = i * fftl / 2;
		//データをコピー
		for (j = 0;j < fftl; j++) waveBuff[j] = x[offset + j] *
			(0.5 - 0.5 * cos(2.0 * PI * (double)j / (double)fftl));//窓を掛ける;

//fft実行
		fft_execute(wave_f_fft);

		//スペクトルを格納
		for (j = 0;j < fftl / 2 + 1; j++)
		{
			waveSpecgram[i][j][0] = waveSpec[j][0];
			waveSpecgram[i][j][1] = waveSpec[j][1];
		}
	}

	fft_destroy_plan(wave_f_fft);
	free(waveBuff);
	free(waveSpec);

}

//スペクトルから波形を再構築
void rebuildWave(double* x, int xLen, int fftl, int equLen, fft_complex** waveSpecgram)
{
	int i, j;
	double* waveBuff;
	fft_plan			wave_i_fft;				// fftセット
	fft_complex* waveSpec;	// スペクトル
	waveBuff = (double*)malloc(sizeof(double) * fftl);
	waveSpec = (fft_complex*)malloc(sizeof(fft_complex) * fftl);
	wave_i_fft = fft_plan_dft_c2r_1d(fftl, waveSpec, waveBuff, FFT_ESTIMATE);

	int offset;
	for (i = 0;i < xLen;i++) x[i] = 0;

	for (i = 0;i < equLen;i++)
	{
		offset = i * fftl / 2;

		//スペクトルを格納
		for (j = 0;j < fftl / 2 + 1; j++)
		{
			waveSpec[j][0] = waveSpecgram[i][j][0];
			waveSpec[j][1] = waveSpecgram[i][j][1];
		}


		//fft実行
		fft_execute(wave_i_fft);

		for (j = 0;j < fftl; j++) waveBuff[j] /= fftl;

		//データをコピー
		for (j = 0;j < fftl; j++) x[offset + j] += waveBuff[j];

	}

	fft_destroy_plan(wave_i_fft);
	free(waveBuff);
	free(waveSpec);

}

//Bフラグ（息）を適用する
void breath2(double* f0, int tLen, int fs, double* x, int xLen, fft_complex** waveSpecgram, int equLen, int fftl, int flag_B)
{
	int i, j;

	//ノイズfftの準備
	double* noiseData;
	double* noiseBuff;
	double* noise;
	fft_plan			noise_f_fft;				// fftセット
	fft_plan			noise_i_fft;				// fftセット
	fft_complex* noiseSpec;	// スペクトル

	noiseData = (double*)malloc(sizeof(double) * xLen);
	for (i = 0;i < xLen; i++) noiseData[i] = (double)rand() / (RAND_MAX + 1) - 0.5;
	noise = (double*)malloc(sizeof(double) * xLen);
	for (i = 0;i < xLen; i++) noise[i] = 0.0;
	//	for(i=0;i < xLen; i++) noiseData[i] *= noiseData[i] * (noiseData[i] < 0)? -1 : 1;//ノイズの分布をいじる
	noiseBuff = (double*)malloc(sizeof(double) * fftl);
	noiseSpec = (fft_complex*)malloc(sizeof(fft_complex) * fftl);
	noise_f_fft = fft_plan_dft_r2c_1d(fftl, noiseBuff, noiseSpec, FFT_ESTIMATE);
	noise_i_fft = fft_plan_dft_c2r_1d(fftl, noiseSpec, noiseBuff, FFT_ESTIMATE);

	//wavefftの準備
	fft_complex* waveSpec;	// スペクトル
	waveSpec = (fft_complex*)malloc(sizeof(fft_complex) * fftl);

	int offset;
	double volume;

	int SFreq, MFreq, EFreq;

	SFreq = (int)(fftl * 1500 / fs);//ブレス開始周波数
	MFreq = (int)(fftl * 5000 / fs);//ブレス開始周波数
	EFreq = (int)(fftl * 20000 / fs);//ブレスの周波数帯

	double nowIndex;
	int sIndex, eIndex;
	double nowF0;
	int specs, spece;
	double hs, he;
	int baion;

	for (i = 0; i < equLen; i++)
	{
		offset = i * fftl / 2;
		//データをコピー
		for (j = 0;j < fftl; j++) noiseBuff[j] = noiseData[offset + j] *
			(0.5 - 0.5 * cos(2.0 * PI * (double)j / (double)fftl));//窓を掛ける;

//fft実行
		fft_execute(noise_f_fft);

		//スペクトル包絡（超手抜き）
		for (j = 0;j < fftl / 2 + 1; j++) waveSpec[j][0] = sqrt(waveSpecgram[i][j][0] * waveSpecgram[i][j][0] + waveSpecgram[i][j][1] * waveSpecgram[i][j][1]);
		for (j = 0;j < fftl / 2 + 1; j++) waveSpec[j][0] = log10(waveSpec[j][0] + 0.00000001);//対数化
		for (j = 0;j < fftl / 2 + 1; j++) waveSpec[j][1] = waveSpec[j][0];

		nowIndex = max(0, min(tLen - 1, (double)(offset + fftl / 2) / fs * 1000 / FRAMEPERIOD));
		sIndex = min(tLen - 2, (int)nowIndex);
		eIndex = sIndex + 1;

		nowF0 = (f0[sIndex] == 0 && f0[eIndex] == 0) ? DEFAULT_F0 :
			(f0[sIndex] == 0) ? f0[eIndex] :
			(f0[eIndex] == 0) ? f0[sIndex] :
			(f0[eIndex] - f0[sIndex]) * (nowIndex - sIndex) + f0[sIndex];

		specs = 0;
		hs = 0.0;
		j = 0;
		baion = 1;
		spece = 0;
		for (baion = 1;spece != fftl / 2 + 1;baion++)
		{
			spece = min(fftl / 2 + 1, (int)((double)fftl / fs * nowF0 * baion + 0.5));
			he = waveSpec[spece][1];
			for (j = specs;j < spece;j++)
			{
				waveSpec[j][0] = (he - hs) / (spece - specs) * (j - specs) + hs;
			}
			specs = spece;
			hs = he;
		}

		for (j = 0;j < fftl / 2 + 1; j++) waveSpec[j][0] = pow(10, waveSpec[j][0]);//振幅化

		//ノイズのスペクトルを変形
		for (j = 0;j < SFreq; j++)
		{
			noiseSpec[j][0] = 0.0;
			noiseSpec[j][1] = 0.0;
		}

		for (;j < MFreq; j++)
		{
			volume = waveSpec[j][0] * (0.5 - 0.5 * cos(PI * (j - SFreq) / (double)(MFreq - SFreq)));
			noiseSpec[j][0] *= volume;
			noiseSpec[j][1] *= volume;
		}
		for (;j < EFreq; j++)
		{
			volume = waveSpec[j][0] * (0.5 - 0.5 * cos(PI + PI * (j - MFreq) / (double)(EFreq - MFreq)));
			noiseSpec[j][0] *= volume;
			noiseSpec[j][1] *= volume;
		}

		for (;j < fftl / 2 + 1; j++)
		{
			noiseSpec[j][0] = 0.0;
			noiseSpec[j][1] = 0.0;
		}

		noiseSpec[0][1] = 0.0;
		noiseSpec[fftl / 2][1] = 0.0;

		//逆fft
		fft_execute(noise_i_fft);
		for (j = 0;j < fftl; j++) noiseBuff[j] /= fftl;

		//窓を掛ける
	//	for(j = 0;j < fftl; j++) noiseBuff[j] *= 0.5 - 0.5*cos(2.0*PI*(double)j/(double)fftl);

		//ノイズを加える
		for (j = 0;j < fftl; j++)
		{
			noise[offset + j] += noiseBuff[j] * 0.2;
		}
	}

	//ノイズを合成
	double noiseRatio = max(0, (double)(flag_B - 50) / 50.0);
	double waveRatio = 1 - noiseRatio;
	for (i = 0;i < xLen;i++) x[i] = x[i] * waveRatio + noise[i] * noiseRatio;

	//後処理
	fft_destroy_plan(noise_f_fft);
	fft_destroy_plan(noise_i_fft);
	free(noise);
	free(noiseData);
	free(noiseBuff);
	free(noiseSpec);
	free(waveSpec);
}

//Oフラグ（声の強さ）
void Opening(double* f0, int tLen, int fs, fft_complex** waveSpecgram, int equLen, int fftl, int flag_O)
{
	int i, j;
	double opn = (double)flag_O / 100.0;
	int sFreq = (int)(fftl * 500 / fs);//制御周波数1
	int eFreq = (int)(fftl * 2000 / fs);//制御周波数2
	double sRatio = -10.0;//制御周波数1の振幅倍率デシベル
	double eRatio = 10.0;//制御周波数2の振幅倍率デシベル

	//周波数ごとの音量マップ作成
	double volume;
	double* volumeMap;
	volumeMap = (double*)malloc(sizeof(double) * fftl / 2 + 1);

	volume = pow(10, sRatio * opn / 20);
	for (j = 0;j < sFreq;j++)
	{
		volumeMap[j] = volume;
	}
	for (;j < eFreq;j++)
	{
		volume = pow(10, ((0.5 + 0.5 * cos(PI + PI / (eFreq - sFreq) * (j - sFreq))) * (eRatio - sRatio) + sRatio) * opn / 20);
		volumeMap[j] = volume;
	}
	volume = pow(10, eRatio * opn / 20);
	for (;j < fftl / 2 + 1;j++)
	{
		volumeMap[j] = volume;
	}

	//周波数ごとの音量を変更
	int f0Frame;
	for (i = 0;i < equLen;i++)
	{
		f0Frame = max(0, min(tLen - 1, (int)((double)((i + 1) * fftl / 2) / fs * 1000 / FRAMEPERIOD + 0.5)));
		if (f0[f0Frame] == 0.0) continue;
		for (j = 0;j < fftl / 2 + 1;j++)
		{
			waveSpecgram[i][j][0] *= volumeMap[j];
			waveSpecgram[i][j][1] *= volumeMap[j];
		}
	}

	free(volumeMap);
}

// フィードフォワードコムフィルタ
// UTAUのBREはこれだと思うんだ
// 息成分の追加、追加で悩んでたけど、声成分の削除とは！？　飴屋さん凄ぇ
// y[n]= x[n] - c*x[n-SampleNumOfWaveLength]
// for B flag
void FeedForwardCombFilter(double* waves, int sNum, double fs,
	double* f0, int fNum, double coef)
{
	// 時間遅れの波形メモリ
	// Frame2Sample(2, FRAMEPERIOD, fs) より大きいので大丈夫
	double dWave[MAX_FFT_LENGTH];
	const int margin = 3; // spline でもまあ使える余裕代
	// 時間遅れの波形にローパスフィルタを掛けて少しでもささやきぽくする
	// ごく軽く利かせないといけない
	const double a = 0.05; // [0 .. 0.2 1.0]
	double delayedWave = 0.0;


	// 初期位置を計算する
	fNum = imin(fNum, GetFramesForDIO(fs, sNum, FRAMEPERIOD));
	int Index = fNum - 1;
	int ed = sNum;
	int st = static_cast<int>(Frame2Sample(Index, FRAMEPERIOD, fs));
	int dWLen = ed - st;
	double lastF0 = 300.0; // f0_target を使ってるから必要ないけども
	for (int i = 0; i < fNum; i++) {
		if (FLOOR_F0 <= f0[i]) lastF0 = f0[i];
	}
	double T0Len_f = Frequency2Sample(lastF0, fs);
	int T0Len_i = static_cast<int>(T0Len_f);
	double itrpOffset = margin - (T0Len_f - T0Len_i);
	int itrp_st = st - T0Len_i - margin;
	int itrpLen = margin + dWLen + margin;

	// フィルタリング(backward)
	// 子音に母音f0のフィルタを掛けたいため、後ろから前に作業している
	Index--;
	while (0 <= itrp_st) {
		// 時間遅れ波形の作成
		itrp1Qtrim_clip(itrpOffset, &waves[itrp_st], itrpLen, 1.0, dWLen, dWave);

		// フィードフォワードコムフィルタ with ローパスフィルタ
		for (int i = ed - 1, j = dWLen - 1; st <= i; i--, j--) {
			delayedWave = dWave[j] + a * delayedWave;
			waves[i] -= coef * delayedWave;
		}

		// 位置を計算する
		if (FLOOR_F0 <= f0[Index]) {
			T0Len_f = Frequency2Sample(f0[Index], fs);
			T0Len_i = static_cast<int>(T0Len_f);
			itrpOffset = margin - (T0Len_f - T0Len_i);
		}
		else {
			// 今までの値を使う
		}
		st = static_cast<int>(Frame2Sample(Index, FRAMEPERIOD, fs));
		ed = static_cast<int>(Frame2Sample(Index + 1, FRAMEPERIOD, fs));
		dWLen = ed - st;
		itrp_st = st - T0Len_i - margin;
		itrpLen = margin + dWLen + margin;

		Index--;
	}
}  // FeedForwardCombFilter


int main(int argc, char* argv[])
{
	// メモリリーク検出
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

	int i, j;
	char* thisChar;
	for (i = 1; i < argc; i++) { // argv[0] may be the file name (no guarantee, see Peter M's comment)
		thisChar = argv[i]; // If the parameter is "abc", thisChar = 'a'
		printf("%s ", thisChar);
	}
	printf("\n");
	//*
	if (argc <= 4)
	{
		fprintf(stderr, "error: 引数の数が不正です．\n");
		return 0;
	}
	//*/

	FILE* fp;
	double fsread;
	int fs, nbit = 16, ch;
	//	int fs, nbit = 16;
	int flag_G = 0;
	if (argc > 5)
	{
		flag_G = strchr(argv[5], 'G') != 0;
	}
	int flag_R = 0;
	if (argc > 5)
	{
		flag_R = strchr(argv[5], 'R') != 0;
	}

	double* x = 0, * f0, * t, * y;
	double** specgram = 0;
	double** residualSpecgram = 0;
	int fftl;

	int signalLen;
	int tLen;

	tLen = readDIOParam(argv[1], &t, &f0, &fs, &signalLen);
	if (flag_R != 0)
	{
		tLen = 0;
	}
	if (tLen != 0)
	{
		fftl = getFFTLengthForStar(fs);
		specgram = readCTParam(signalLen, fs, argv[1], tLen, fftl);
		if (specgram)
		{
			residualSpecgram = readD4CParam(signalLen, fs, argv[1], tLen, fftl);
			if (!residualSpecgram)
			{
				tLen = 0;
			}
		}
		else
		{
			tLen = 0;
		}
	}
	if (tLen == 0)
	{
		freeSpecgram(specgram, tLen);
		freeSpecgram(residualSpecgram, tLen);

		//		x = wavread(argv[1], &fs, &nbit, &signalLen);
		x = audio_read(argv[1], &fsread, &nbit, &ch, &signalLen);
		fs = static_cast<int>(fsread);
		if (x == NULL)
		{
			fprintf(stderr, "error: 指定されたファイルは存在しません．\n");
			return 0;
		}
		tLen = getDIOParam(x, signalLen, fs, FRAMEPERIOD, &t, &f0, argv[1]);
		if (tLen != 0)
		{
			writeDIOParam(signalLen, fs, tLen, argv[1], t, f0);
		}
		else
		{
			free(x);
			fprintf(stderr, "error: DIOパラメータ作成失敗．\n");
			return 0;
		}
		fftl = getFFTLengthForStar(fs);
		specgram = getCTParam(x, signalLen, fs, t, f0, tLen, fftl);
		if (!specgram)
		{
			free(x);
			free(t);
			free(f0);
			fprintf(stderr, "error: CHEAPTRICKパラメータ作成失敗．\n");
			return 0;
		}
		else
		{
			writeCTParam(signalLen, fs, argv[1], specgram, tLen, fftl);
		}
		residualSpecgram = getD4CParam(x, signalLen, fs, t, f0, tLen, fftl);
		if (!residualSpecgram)
		{
			free(x);
			free(t);
			free(f0);
			free(specgram);
			fprintf(stderr, "error: D4Cパラメータ作成失敗．\n");
			return 0;
		}
		else
		{
			writeD4CParam(signalLen, fs, argv[1], residualSpecgram, tLen, fftl);
		}
	}

	// 注意！
	// 全ての分析が完了するまで，F0は操作しないこと．

	//引数を取り込む
	double offset = 0;
	double wavelength = (double)signalLen / (double)fs * 1000;//音源の長さをmsecにする
	double length_req = wavelength;//初期値いれとく
	double fixed = 0;
	double blank = 0;
	double* target_freqs = NULL;
	double velocity = 1.0;
	double value = 100;
	if (argc > 4) sscanf(argv[4], "%lf", &value);
	velocity = pow(2, value / 100 - 1.0);
	if (argc > 6) sscanf(argv[6], "%lf", &offset);
	if (argc > 7) sscanf(argv[7], "%lf", &length_req);
	if (argc > 8) sscanf(argv[8], "%lf", &fixed);
	if (argc > 9) sscanf(argv[9], "%lf", &blank);
#ifdef _DEBBUG
	printf("Parameters\n");
	printf("velocity      :%lf\n", velocity);
	printf("offset        :%lf\n", offset);
	printf("request length:%lf\n", length_req);
	printf("fixed         :%lf\n", fixed);
	printf("blank         :%lf\n", blank);
#endif
	//伸縮の概念図
	//  offset    fixed      m2      blank
	//|--------|--------|---------|---------| 原音
	//         |        |          |
	//         |   l1   |    l2     |
	//         |--------|------------|  出力
	// l1  = fixed / velocity
	// l2  = m2    / stretch
	// l1 + l2  = 要求長＝argv[7]
	if (blank < 0)
	{
		blank = wavelength - offset + blank;
		if (blank < 0) blank = 0;
	}
	if (offset + blank >= wavelength)
	{
		fprintf(stderr, "error: パラメータ異常．\n");
		free(x);
		return 0;
	}
	if (offset + blank + fixed >= wavelength)
	{
		fixed = wavelength - offset + blank;
	}
	double l1, l2;
	double m2 = wavelength - offset - fixed - blank;

	l1 = fixed / velocity;
	l2 = length_req - l1;
	if (m2 <= 0 && l2 > 0)
	{
		fprintf(stderr, "error: パラメータ異常2．\n");
		free(x);
		return 0;
	}
	double stretch = m2 / l2;
	if (stretch > 1.0) stretch = 1.0;

	int outSamples = (int)(length_req * 0.001 * fs + 1);
	int oLen = GetSamplesForDIO(fs, outSamples, FRAMEPERIOD);

	printf("File information\n");
	printf("Sampling : %d Hz %d Bit\n", fs, nbit);
	printf("Input:\n");
	printf("Length %d [sample]\n", signalLen);
	printf("Length %f [sec]\n", (double)signalLen / (double)fs);
	printf("Output:\n");
	printf("Length %d [sample]\n", outSamples);
	printf("Length %f [sec]\n", (double)outSamples / (double)fs);


	int flag_t = 0;
	char* cp;
	if (argc > 5 && (cp = strchr(argv[5], 't')) != 0)
	{
		sscanf(cp + 1, "%lf", &flag_t);
	}
	int flag_e = 1; //Take note e flag changes stretching to looping, not the other way around
	if (argc > 5 && (cp = strchr(argv[5], 'e')) != 0)
	{
		flag_e = 0;
	}
	int flag_B = 50;//BRE（息）成分
	if (argc > 5 && (cp = strchr(argv[5], 'B')) != 0)
	{
		sscanf(cp + 1, "%d", &flag_B);
		flag_B = max(0, min(100, flag_B));
	}
	int flag_O = 0;//独自フラグ　声の強さ
	if (argc > 5 && (cp = strchr(argv[5], 'O')) != 0)
	{
		sscanf(cp + 1, "%d", &flag_O);
		flag_O = max(-100, min(100, flag_O));
	}
	double modulation = 100;
	if (argc > 11) sscanf(argv[11], "%lf", &modulation);

	double volume = 1.0;
	if (argc > 10)
	{
		volume = 100;
		sscanf(argv[10], "%lf", &volume);
		volume *= 0.01;
	}

	double target_freq = name2freq(argv[3], flag_t);
	double freq_avg = getFreqAvg(f0, tLen);

#ifdef _DEBUG
	printf("volume        :%lf\n", volume);
	printf("modulation    :%lf\n", modulation);
	printf("target frequency     :%lf\n", target_freq);
	printf("input frequency(avg.):%lf\n", freq_avg);
#endif

	double* f0out = (double*)malloc(oLen * sizeof(double));
	memset(f0out, 0, sizeof(double) * oLen);
	//double *tout = (double*)malloc(oLen * sizeof(double));
	int* pit = NULL;
	double tempo = 120;
	int pLen = oLen;
	int pStep = 256;
	if (argc > 13)
	{
		cp = argv[12];
		sscanf(cp + 1, "%lf", &tempo);
		pStep = (int)(60.0 / 96.0 / tempo * fs + 0.5);
		pLen = outSamples / pStep + 1;
		pit = (int*)malloc((pLen + 1) * sizeof(int));
		memset(pit, 0, (pLen + 1) * sizeof(int));
		decpit(argv[13], pit, pLen);
	}
	else
	{
		pit = (int*)malloc((pLen + 1) * sizeof(int));
		memset(pit, 0, (pLen + 1) * sizeof(int));
	}

	double** specgram_out = (double**)malloc(sizeof(double*) * oLen);
	double** residualSpecgram_out = (double**)malloc(sizeof(double*) * oLen);
	for (i = 0;i < oLen;i++)
	{
		specgram_out[i] = (double*)malloc(sizeof(double) * (fftl / 2 + 1));
		memset(specgram_out[i], 0, sizeof(double) * (fftl / 2 + 1));
		residualSpecgram_out[i] = (double*)malloc(sizeof(double) * (fftl / 2 + 1));
		memset(residualSpecgram_out[i], 0, sizeof(double) * (fftl / 2 + 1));
	}
	//出力f0数列
	double tmo, tmi, laststretch, m2edit;
	//	double ires = 0.0;
	int totalloops, direct = 0;
	int loopcount = 1;
	DWORD elapsedTime = timeGetTime();
	printf("\nTransform\n");
#ifdef _DEBUG
	FILE* fp0 = fopen("time.txt", "wt");
	FILE* fp1 = fopen("dio.txt", "wt");
	FILE* fp2 = fopen("star.txt", "wt");
	FILE* fp3 = fopen("plat.txt", "wt");
#endif

	if (flag_e == 0)
	{
		//totalloops = (int)floor(1.0/stretch);
		m2edit = m2 - FRAMEPERIOD / 2.0;
		totalloops = (int)floor(l2 / m2edit);
		totalloops--;
		if (totalloops > 0)
			laststretch = m2 / (l2 - m2edit * totalloops);
		else
			flag_e = 1;
		//	        printf("Loops: %d \n", totalloops);
	}

	double v, u;
	int n, m;
	if (flag_e == 0)
	{
		for (i = 0; i < oLen; i++)
		{
			tmo = FRAMEPERIOD * i;
			if (tmo < l1)
			{
				tmi = offset + tmo * velocity;
			}
			else if (tmo < l1 + m2edit * totalloops)
			{
				if (direct == 0)
				{
					tmi = offset + fixed + (tmo - l1 - m2edit * (loopcount - 1));
					if (tmo - l1 > m2edit * loopcount)
					{
						loopcount++;
						direct = 1;
					}
					//	            printf("Count: %d \n", loopcount);
				}
				else
				{
					tmi = wavelength - blank - (tmo - l1 - m2edit * (loopcount - 1));
					if (tmo - l1 > m2edit * loopcount)
					{
						loopcount++;
						direct = 0;
					}
					//	            printf("Count: %d \n", loopcount);
				}
			}
			else
			{
				if (direct == 1)
					tmi = offset + fixed + (tmo - l1 - m2edit * totalloops) * laststretch;
				else
					tmi = wavelength - blank - (tmo - l1 - m2edit * totalloops) * laststretch;
			}
			v = tmi / FRAMEPERIOD;
			n = (int)floor(v);
			v -= n;

			double f0i = f0[n];
			if (n < tLen - 1) {
				double f0j = f0[n + 1];
				if (f0i != 0 || f0j != 0)
				{
					if (f0i == 0) f0i = freq_avg;
					if (f0j == 0) f0j = freq_avg;
					f0i = f0i * (1.0 - v) + f0j * v;
				}
			}

			u = tmo * 0.001 * fs / pStep;
			m = (int)floor(u);
			u -= m;
			if (m >= pLen) { m = pLen - 1; v = 0.0; }
			f0out[i] = target_freq * pow(2, (pit[m] * (1.0 - u) + pit[m + 1] * u) / 1200.0);
			f0out[i] *= pow(f0i / freq_avg, modulation * 0.01);

			for (j = 0; j <= fftl / 2; j++)
			{
				if (n < tLen - 1)
				{
					specgram_out[i][j] = specgram[n][j] * (1.0 - v) + specgram[n + 1][j] * v;
				}
				else
				{
					specgram_out[i][j] = specgram[tLen - 1][j];
				}
			}

			int m = n;
			if (v > 0.5) m++;
			//		for (j = 0; j <= fftl; j++)
			for (j = 0; j <= fftl / 2; j++)
			{
				if (m < tLen)//(n < tLen - 1)
				{
					//residualSpecgram_out[i][j] = residualSpecgram[n][j] * (1.0 - v) + residualSpecgram[n + 1][j] * v;
					residualSpecgram_out[i][j] = residualSpecgram[m][j];
				}
				else
				{
					residualSpecgram_out[i][j] = residualSpecgram[tLen - 1][j];
				}
			}
		}
	}
	else
	{
		for (i = 0; i < oLen; i++)
		{
			tmo = FRAMEPERIOD * i;
			if (tmo < l1)
			{
				tmi = offset + tmo * velocity;
			}
			else
			{
				tmi = offset + fixed + (tmo - l1) * stretch;
			}
#ifdef _DEBUG
			fprintf(fp0, "%0.6lf\t%0.6lf\n", tmi, tmo);
#endif
			v = tmi / FRAMEPERIOD;
			n = (int)floor(v);
			v -= n;

			double f0i = f0[n];
			if (n < tLen - 1) {
				double f0j = f0[n + 1];
				if (f0i != 0 || f0j != 0)
				{
					if (f0i == 0) f0i = freq_avg;
					if (f0j == 0) f0j = freq_avg;
					f0i = f0i * (1.0 - v) + f0j * v;
				}
			}

			u = tmo * 0.001 * fs / pStep;
			m = (int)floor(u);
			u -= m;
			if (m >= pLen) { m = pLen - 1; v = 0.0; }
			f0out[i] = target_freq * pow(2, (pit[m] * (1.0 - u) + pit[m + 1] * u) / 1200.0);
			f0out[i] *= pow(f0i / freq_avg, modulation * 0.01);
			/*		printf("pStep: %d\n", pStep);
					printf("tempo: %.6f\n", tempo);
					printf("target_freq: %.6f\n", target_freq);
					printf("freq_avg: %.6f\n", freq_avg);
					printf("f0out: %.6f\n", f0out[i]);z
			*/

#ifdef _DEBUG
			fprintf(fp1, "%lf\n", f0out[i]);
#endif
			for (j = 0; j <= fftl / 2; j++)
			{
				if (n < tLen - 1)
				{
					specgram_out[i][j] = specgram[n][j] * (1.0 - v) + specgram[n + 1][j] * v;
				}
				else
				{
					specgram_out[i][j] = specgram[tLen - 1][j];
				}
#ifdef _DEBUG
				fprintf(fp2, "%lf\t", specgram_out[i][j] * 1000000);
#endif
				/*if (_isnan(specgram_out[i][j]))
				{
					printf("nan!\n");
				}
				else if (specgram_out[i][j] == 0)
				{
					printf("(%d)(%d)zero!\n", i, j);
				}*/
			}
#ifdef _DEBUG
			fprintf(fp2, "\n");
#endif
			int m = n;
			if (v > 0.5) m++;
			//		for (j = 0; j <= fftl; j++)
			for (j = 0; j <= fftl / 2; j++)
			{
				if (m < tLen)//(n < tLen - 1)
				{
					//residualSpecgram_out[i][j] = residualSpecgram[n][j] * (1.0 - v) + residualSpecgram[n + 1][j] * v;
					residualSpecgram_out[i][j] = residualSpecgram[m][j];
				}
				else
				{
					residualSpecgram_out[i][j] = residualSpecgram[tLen - 1][j];
				}
				/*if (_isnan(residualSpecgram_out[i][j]))
				{
					printf("nan!\n");
				}
				else if (residualSpecgram_out[i][j] == 0)
				{
					printf("(%d)(%d)zero!\n", i, j);
				}*/
			}
			/*
			//debug
					for (j = 0; j < fftl/2; j++)
					{
						ires += residualSpecgram_out[i][j];
					}
					printf("iresavg: %lf\n", ires*2/fftl);
					ires = 0.0;
			*/
#ifdef _DEBUG
			for (j = 0; j < fftl / 2; j += 8)
			{
				fprintf(fp3, "%lf\t", residualSpecgram_out[i][j]);
			}
			fprintf(fp3, "\n");
			for (j = 0; j < fftl / 2; j += 8)
			{
				fprintf(fp3, "%lf\t", residualSpecgram_out[i][j + 1]);
			}
			fprintf(fp3, "\n");
#endif
		}
	}

#ifdef _DEBUG
	fclose(fp0);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
#endif
	// スペクトル伸縮
	//stretchSpectrum(double **specgram, double ratio)
	if (argc > 5 && (cp = strchr(argv[5], 'g')) != 0)
	{
		double w = 0;
		double ratioW = 1.0;
		sscanf(cp + 1, "%lf", &w);
		if (w > 100) w = 100;
		if (w < -100) w = -100;
		ratioW = pow(10, -w / 200);
		stretchSpectrum(specgram_out, oLen, ratioW, fs, fftl);
	}
	printf("TRANSFORM: %d [msec]\n", timeGetTime() - elapsedTime);

	// 合成
	y = (double*)malloc(sizeof(double) * outSamples);
	memset(y, 0, sizeof(double) * outSamples);

	printf("\nSynthesis\n");
	elapsedTime = timeGetTime();
	//synthesis(f0out, oLen, specgram_out, residualSpecgram_out, fftl, FRAMEPERIOD, fs, y, outSamples);
//	synthesis(f0out, oLen, specgram_out, residualSpecgram_out, fftl, FRAMEPERIOD, fs, y, outSamples);
	Synthesis(f0out, oLen, specgram_out, residualSpecgram_out, fftl, FRAMEPERIOD, fs, outSamples, y);
	printf("WORLD: %d [msec]\n", timeGetTime() - elapsedTime);

	//tn_fnds flag start
	//イコライジング
	int equfftL = 1024;//イコライザーのfft長
	int equLen = (outSamples / (equfftL / 2)) - 1; //繰り返し回数
	fft_complex** waveSpecgram;  //スペクトル
	waveSpecgram = (fft_complex**)malloc(sizeof(fft_complex*) * equLen);
	for (i = 0;i < equLen;i++) waveSpecgram[i] = (fft_complex*)malloc(sizeof(fft_complex) * (equfftL / 2 + 1));
	//スペクトル作成
	if (flag_B > 50 || flag_O != 0)
	{
		createWaveSpec(y, outSamples, equfftL, equLen, waveSpecgram);
	}

	//声の強さ
	if (flag_O != 0)
	{
		Opening(f0out, oLen, fs, waveSpecgram, equLen, equfftL, flag_O);
	}

	//イコライズ結果を波形に反映
	if (flag_O != 0)
	{
		rebuildWave(y, outSamples, equfftL, equLen, waveSpecgram);
	}

	//ノイズ
	if (flag_B > 50)
	{
		breath2(f0out, oLen, fs, y, outSamples, waveSpecgram, equLen, equfftL, flag_B);
	}
	else if (flag_B < 50)
	{
		FeedForwardCombFilter(y, outSamples, fs, f0out, oLen, flag_B * 0.01 + 0.005);
	}
	//tn_fnds flags end

	// ファイルの書き出し (内容には関係ないよ)
	char header[44];
	short* output;
	double maxAmp;

	output = (short*)malloc(sizeof(short) * outSamples);
	// 振幅の正規化
	maxAmp = 0.0;
#ifdef _DEBUG
	{
		FILE* f = fopen("synthesis.txt", "wt");
		for (i = 0;i < outSamples;i++)
		{
			fprintf(f, "%f\n", y[i]);
		}
		fclose(f);
	}
#endif
	for (i = 0;i < outSamples;i++)
	{
		if (!_isnan(y[i]))
		{
			maxAmp = maxAmp < fabs(y[i]) ? fabs(y[i]) : maxAmp;
		}
	}
	value = 0.86;
	if (argc > 5)
	{
		cp = strchr(argv[5], 'P');
		if (cp)
		{
			sscanf(cp + 1, "%lf", &value);
			if (value < 0) value = 0;
			else if (value > 100) value = 100;
			value *= 0.01;
		}
	}
	double peekcomp = 3 * 32.0 * pow(512.0 / maxAmp, value);
	//double peekcomp = pow( 16384.0 / maxAmp, value);
	for (i = 0;i < outSamples;i++)
	{
		//****** 空白部分でsynthesis()がnanを吐くので対策した*********
		value = _isnan(y[i]) ? 0 : y[i] * peekcomp * volume;
		if (value > 32767.0) value = 32767.0;
		else if (value < -32767.0) value = -32767.0;
		output[i] = (short)value;
	}

	fp = fopen(argv[1], "rb");
	makeHeader(header, outSamples, fs, nbit);

	fp = fopen(argv[2], "wb");
	fwrite(header, sizeof(char), 44, fp);
	fwrite(output, sizeof(short), outSamples, fp);
	fclose(fp);

	printf("complete.\n");

	// メモリの解放
	free(output);
	free(x); free(t); free(f0); free(y);
	for (i = 0;i < tLen;i++)
	{
		free(specgram[i]);
		free(residualSpecgram[i]);
	}
	free(specgram);
	free(residualSpecgram);

	free(f0out);
	for (i = 0;i < oLen;i++)
	{
		free(specgram_out[i]);
		free(residualSpecgram_out[i]);
	}
	free(specgram_out);
	free(residualSpecgram_out);
	free(pit);

	for (i = 0;i < equLen;i++) free(waveSpecgram[i]);
	free(waveSpecgram);

	printf("complete.\n");

	return 0;
}

