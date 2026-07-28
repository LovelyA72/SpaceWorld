/*
 * This program is licensed under the GNU General Public License (GPL).
 * See the LICENSE file in the repository root for the license terms.
 *
 * This source started from test.cpp in WORLD 0.0.4 and was adapted for
 * World4UTAU.
 *
 * Arguments:
 *  1  Input file
 *  2  Output file
 *  3  Target note
 *  4  Velocity
 *  5  Flags
 *  6  Offset
 *  7  Requested length
 *  8  Fixed leading region
 *  9  Unused trailing region, or usable length from offset when negative
 * 10  Volume
 * 11  Modulation
 * 12+ Tempo and pitch bend
 */


#include <stdio.h>
#include <stdlib.h>

#include <windows.h>

#include "world.h"
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


// Analysis frame shift [ms].
#define FRAMEPERIOD (1000.0*256/44100)
//#define FRAMEPERIOD 5.80498866213152

#pragma comment(lib, "winmm.lib")

int get64(int ch)
{
	if (ch >= '0' && ch <= '9')
	{
		return ch - '0' + 52;
	}
	else if (ch >= 'A' && ch <= 'Z')
	{
		return ch - 'A';
	}
	else if (ch >= 'a' && ch <= 'z')
	{
		return ch - 'a' + 26;
	}
	else if (ch == '+')
	{
		return 62;
	}
	else if (ch == '/')
	{
		return 63;
	}
	else
	{
		return 0;
	}
}

int decpit(char* pitch_text, int* pitch_values, int value_count)
{
	int text_len = 0;
	int i, pitch = 0;
	int k = 0, repeat_count, ii;
	if (pitch_text != NULL)
	{
		text_len = static_cast<int>(strlen(pitch_text));
		for (i = 0; i < text_len; i += 2)
		{
			if (pitch_text[i] == '#')
			{
				i++;
				sscanf(pitch_text + i, "%d", &repeat_count);
				for (ii = 0; ii < repeat_count && k < value_count; ii++) {
					pitch_values[k++] = pitch;
				}
				while (pitch_text[i] != '#' && pitch_text[i] != 0) i++;
				i--;
			}
			else
			{
				pitch = get64(pitch_text[i]) * 64 + get64(pitch_text[i + 1]);
				if (pitch > 2047) pitch -= 4096;
				if (k < value_count) {
					pitch_values[k++] = pitch;
				}
			}
		}
	}
	return text_len;
}

double name2freq(char* note, double cent_shift)
{
	char note_ch;
	int semitone, octave_pos, octave, cent_steps;
	//01234567890A
	//C D EF G A B
	note_ch = note[0];
	if (note_ch >= 'A' && note_ch <= 'G')
	{
		if (note_ch <= 'B')
		{
			semitone = 9 + (note_ch - 'A') * 2;
		}
		else if (note_ch <= 'E')
		{
			semitone = (note_ch - 'C') * 2;
		}
		else
		{
			semitone = 5 + (note_ch - 'F') * 2;
		}

		note_ch = note[1];

		octave_pos = 2;
		if (note_ch == '#') {
			semitone++;
		}
		else if (note_ch == 'b')
		{
			semitone--;
		}
		else
		{
			octave_pos = 1;
		}

		if (note[octave_pos] == 0)
		{
			return 0;
		}

		sscanf(note + octave_pos, "%d", &octave);

		cent_steps = (int)((semitone + octave * 12 - 21) * 10 + cent_shift);
		// A1 maps to 55 Hz when cent_steps is zero.
		return 55.0 * pow(2, (double)cent_steps / 120.0);
	}
	return 0;
}

void makeFilename(const char* input_name, const char* extension, char* output_name)
{
	strcpy(output_name, input_name);
	char* dot = strrchr(output_name, '.');
	if (dot) *dot = 0;
	strcat(output_name, extension);
}
// TODO: Avoid reading the waveform when all analysis files are available.

double getFreqAvg(double f0[], int frame_count)
{
	int i, j;
	double current_f0 = 0, weight;
	double weights[6], delta;
	double weighted_sum = 0;
	double weight_sum = 0;
	for (i = 0; i < frame_count; i++)
	{
		current_f0 = f0[i];
		if (current_f0 < 1000.0 && current_f0 > 55.0)
		{
			weight = 1.0;
			// give more weight to values close to the preceding values.
			for (j = 0; j <= 5; j++)
			{
				if (i > j) {
					delta = f0[i - j - 1] - current_f0;
					weights[j] = current_f0 / (current_f0 + delta * delta);
				}
				else {
					weights[j] = 1 / (1 + current_f0);
				}
				weight *= weights[j];
			}
			weighted_sum += current_f0 * weight;
			weight_sum += weight;
		}
	}
	if (weight_sum > 0) weighted_sum /= weight_sum;
	return weighted_sum;
}

int hasUTAUfrq(char* input_name)
{
	char frq_name[512];
	strcpy(frq_name, input_name);
	strcat(frq_name, ".frq");

	FILE* fp = fopen(frq_name, "rb");
	if (NULL == fp) {
		strcpy(frq_name, input_name);
		char* extension = strrchr(frq_name, '.');
		if (NULL != extension) *extension = '_';
		strcat(frq_name, ".frq");
		fp = fopen(frq_name, "rb");
	}
	if (NULL == fp) return 0;
	fclose(fp);
	return 1;
}

// The FREQ0003 layout and amplitude calculation was adapted from OpenUtau
// OpenUtau is Copyright (c) StAkira and licensed under the MIT License
// Thank you StAkira for the awesome OpenUtau! You're the best! <3
int writeUTAUfrq(double wave[], int sample_count, int sample_rate,
	char* input_name)
{
	const int pitch_step = 256;
	double frame_period = Sample2Time(pitch_step, sample_rate);
	int frame_count = GetSamplesForDIO(sample_rate, sample_count, frame_period);
	double* time_axis = (double*)malloc(frame_count * sizeof(double));
	double* f0 = (double*)malloc(frame_count * sizeof(double));
	double* refined_f0 = (double*)malloc(frame_count * sizeof(double));
	if (!time_axis || !f0 || !refined_f0)
	{
		if (time_axis) free(time_axis);
		if (f0) free(f0);
		if (refined_f0) free(refined_f0);
		fprintf(stderr, "Unable to allocate memory for .frq generation.\n");
		return EXIT_FAILURE;
	}

	DioWithFramePeriod(wave, sample_count, sample_rate, frame_period,
		time_axis, f0, 0);
	StoneMask(wave, sample_count, sample_rate, time_axis, f0,
		frame_count, refined_f0);

	double avg_f0 = 0.0;
	int voiced_count = 0;
	for (int i = 0; i < frame_count; i++)
	{
		int invalid_class = _fpclass(refined_f0[i]) & 0x0087;
		if (invalid_class != 0) refined_f0[i] = 0.0;
		if (refined_f0[i] > 0.0)
		{
			avg_f0 += refined_f0[i];
			voiced_count++;
		}
	}
	if (voiced_count > 0) avg_f0 /= voiced_count;

	char frq_name[512];
	strcpy(frq_name, input_name);
	char* extension = strrchr(frq_name, '.');
	if (NULL != extension) *extension = '_';
	strcat(frq_name, ".frq");

	printf("write .frq");
	FILE* fp = fopen(frq_name, "wb");
	if (NULL == fp)
	{
		free(time_axis);
		free(f0);
		free(refined_f0);
		fprintf(stderr, "Unable to write %s.\n", frq_name);
		return EXIT_FAILURE;
	}

	fwrite("FREQ0003", 1, 8, fp);
	fwrite(&pitch_step, sizeof(int), 1, fp);
	fwrite(&avg_f0, sizeof(double), 1, fp);
	int reserved = 0;
	for (int i = 0; i < 4; i++)
		fwrite(&reserved, sizeof(int), 1, fp);
	fwrite(&frame_count, sizeof(int), 1, fp);
	for (int i = 0; i < frame_count; i++)
	{
		fwrite(&refined_f0[i], sizeof(double), 1, fp);
		double amp = 0.0;
		int start = i * pitch_step;
		int end = min(start + pitch_step, sample_count);
		for (int j = start; j < end; j++)
			amp += fabs(wave[j]);
		if (start < end) amp = amp * 32768.0 / (end - start);
		fwrite(&amp, sizeof(double), 1, fp);
	}
	fclose(fp);
	printf("\n");

	free(time_axis);
	free(f0);
	free(refined_f0);
	return 0;
}

double getUTAUfrq(char* input_name, int frame_count, double sample_rate, double frame_period,
	double* f0)
{
	// Build the frequency-file name.
	char frq_name[512];
	strcpy(frq_name, input_name);
	strcat(frq_name, ".frq");

	// check *.aiff.frq
	FILE* fp = fopen(frq_name, "rb");
	if (NULL == fp) {
		strcpy(frq_name, input_name);
		char* extension = strrchr(frq_name, '.');
		// '*.wav' => '*_wav'
		if (NULL != extension) *extension = '_';
		strcat(frq_name, ".frq");

		// check *_wav.frq
		fp = fopen(frq_name, "rb");
		if (NULL == fp) {
			printf("There is no UTAU frequency file.\n");
			return 0.0;
		}
	}

	// Validate the header.
	char header[128];
	fread(header, sizeof(char), 4, fp); // "FREQ"
	if (0 != strncmp(header, "FREQ", 4)) {
		fclose(fp);
		printf("Header Error \"FREQ\".\n");
		return 0.0;
	}

	fread(header, sizeof(char), 4, fp); // "0003"
	if (0 != strncmp(header, "0003", 4)) {
		fclose(fp);
		printf("Header Error version.\n");
		return 0.0;
	}

	int pitch_step;
	fread(&pitch_step, sizeof(int), 1, fp); // Default: 256

	double avg_f0;
	fread(&avg_f0, sizeof(double), 1, fp); // Average frequency

	fread(header, sizeof(char), 8 + 8, fp); // Reserved fields; skip them.

	int utau_frames;
	fread(&utau_frames, sizeof(int), 1, fp); // Number of records

	double* utau_f0 = new double[utau_frames];
	for (int i = 0; i < utau_frames; i++) {
		fread(&utau_f0[i], sizeof(double), 1, fp); // F0
		fread(header, sizeof(double), 1, fp); // Skip amplitude.
	}
	fclose(fp);
#if RUNTIME_CHECK
	for (int i = 0;i < utau_frames;i++) {
		if (0 == _finite(utau_f0[i])) {
			printf("in %s\n", __func__);
			printf("utau_f0[%d] error!! inf or NaN !!\n", i);
			getchar();
		}
	}
#endif

	// Clean the F0 contour. Values outside the allowed range are unvoiced.
	for (int i = 0; i < utau_frames; i++) {
		if ((utau_f0[i] < FLOOR_F0) || (CEIL_F0 < utau_f0[i]))
			utau_f0[i] = 0.0;
	}


	// Find the voiced regions in the UTAU F0 contour.
	int* voiced_start = new int[utau_frames];
	int* voiced_end = new int[utau_frames];
	int voiced_count, end_count;
	voiced_count = end_count = 0;
	if (FLOOR_F0 <= utau_f0[0]) voiced_start[voiced_count++] = 0;
	// Keep voiced_start[X] and voiced_end[X] paired.
	for (int i = 1; i < utau_frames; i++) {
		// Unvoiced to voiced: start of a voiced region.
		if ((utau_f0[i - 1] < FLOOR_F0) && (FLOOR_F0 <= utau_f0[i]))
			voiced_start[voiced_count++] = i;
		// Voiced to unvoiced: end of a voiced region.
		else if ((FLOOR_F0 <= utau_f0[i - 1]) && (utau_f0[i] < FLOOR_F0))
			voiced_end[end_count++] = i - 1;
	}
	// A final voiced region may extend to the last frame.
	if (end_count < voiced_count) voiced_end[end_count++] = utau_frames - 1;


	// Clear the output F0 contour.
	for (int i = 0; i < frame_count; i++) f0[i] = 0.0;
	double utau_period, stride;
	utau_period = Sample2Time(pitch_step, sample_rate); // [ms]
	stride = frame_period / utau_period;


	// Treat voiced runs of four frames or fewer (about 23 ms) as consonants.
	// Longer runs are interpolated into the DIO frame spacing.
	const int min_voiced_frames = 4;
	for (int i = 0; i < voiced_count; i++) {
		int voiced_len = voiced_end[i] - voiced_start[i] + 1;
		if (min_voiced_frames < voiced_len) {
			double src_start_time, src_end_time;
			src_start_time = Frame2Time(voiced_start[i], utau_period);
			src_end_time = Frame2Time(voiced_end[i], utau_period);
			int dst_start, dst_end, dst_len;
			dst_start = static_cast<int>(Time2Frame(src_start_time, frame_period) + 0.5);
			dst_end = static_cast<int>(Time2Frame(src_end_time, frame_period) + 0.5);
			if (fmod(src_end_time, frame_period) <= 0) dst_end--;
			if (frame_count <= dst_end) dst_end = frame_count - 1;
			dst_len = dst_end - dst_start + 1;
			double dst_start_time = Frame2Time(dst_start, frame_period);
			double interp_offset = Time2Frame(dst_start_time - src_start_time, utau_period);

			// Use clipped trim-spline interpolation at the region boundaries.
			itrp1Qtrim_clip(interp_offset, &utau_f0[voiced_start[i]], voiced_len,
				stride, dst_len, &f0[dst_start]);
		}
	}
	// Clamp the interpolated contour to the allowed F0 range again.
	for (int i = 0; i < frame_count; i++) {
		if ((f0[i] < FLOOR_F0) || (CEIL_F0 < f0[i]))
			f0[i] = 0.0;
		else if (f0[i] < avg_f0 * 0.55)
			f0[i] = f0[i] * 2;
		else if (f0[i] > avg_f0 * 1.65)
			f0[i] = f0[i] / 2;
	}

	// Release temporary analysis data.
	delete[] voiced_start; delete[] voiced_end; delete[] utau_f0;
	return avg_f0;
} // getUTAUfrq

double getavgUTAUfrq(char* input_name)
{
	// Build the frequency-file name.
	char frq_name[512];
	strcpy(frq_name, input_name);
	strcat(frq_name, ".frq");

	// check *.aiff.frq
	FILE* fp = fopen(frq_name, "rb");
	if (NULL == fp) {
		strcpy(frq_name, input_name);
		char* extension = strrchr(frq_name, '.');
		// '*.wav' => '*_wav'
		if (NULL != extension) *extension = '_';
		strcat(frq_name, ".frq");

		// check *_wav.frq
		fp = fopen(frq_name, "rb");
		if (NULL == fp) {
			printf("There is no UTAU frequency file.\n");
			return 0.0;
		}
	}

	// Validate the header.
	char header[128];
	fread(header, sizeof(char), 4, fp); // "FREQ"
	if (0 != strncmp(header, "FREQ", 4)) {
		fclose(fp);
		printf("Header Error \"FREQ\".\n");
		return 0.0;
	}

	fread(header, sizeof(char), 4, fp); // "0003"
	if (0 != strncmp(header, "0003", 4)) {
		fclose(fp);
		printf("Header Error Version Problem.\n");
		return 0.0;
	}
	int pitch_step;
	fread(&pitch_step, sizeof(int), 1, fp); // Default: 256

	double avg_f0;
	fread(&avg_f0, sizeof(double), 1, fp); // Average frequency

	fclose(fp);
	return avg_f0;
}

double** readCTParam(int sample_count, int sample_rate, const char* input_name,
	int frame_count, int fft_size)
{
	int i, j;
	char cache_name[512];
	DWORD start_time;
	unsigned short cached_frames = 0;
	unsigned short cached_bins = 0;
	int cached_samples = 0;
	int cached_rate = 0;
	double** spectrogram = 0;

	makeFilename(input_name, ".ctspec", cache_name);

	printf("read .ctspec:");
	start_time = timeGetTime();

	FILE* fp = fopen(cache_name, "rb");
	if (fp)
	{
		char signature[9];
		fread(signature, 1, 8, fp);
		if (strncmp(signature, "world-ct", 8) != 0)
		{
			fclose(fp);
			printf(" bad file.\n");
			return 0;
		}
		fread(&cached_samples, sizeof(int), 1, fp);
		fread(&cached_rate, sizeof(int), 1, fp);
		fread(&cached_frames, sizeof(unsigned short), 1, fp);
		fread(&cached_bins, sizeof(unsigned short), 1, fp);
		if (cached_frames == frame_count && cached_bins == (fft_size / 2 + 1) &&
			sample_count == cached_samples && sample_rate == cached_rate)
		{
			spectrogram = (double**)malloc(frame_count * sizeof(double*));
			if (spectrogram)
			{
				for (i = 0; i < frame_count; i++)
				{
					spectrogram[i] = (double*)malloc((fft_size / 2 + 1) * sizeof(double));
					memset(spectrogram[i], 0, (fft_size / 2 + 1) * sizeof(double));
					if (spectrogram[i])
					{
						for (j = 0; j <= fft_size / 2; j++)
						{
							unsigned short packed_value;
							fread(&packed_value, sizeof(unsigned short), 1, fp);
							spectrogram[i][j] =
								(exp(packed_value / 1024.0) - 1) * 1.16415321826935E-10;
						}
					}
					else
					{
						break;
					}
				}
				if (i < frame_count)
				{
					for (j = 0; j < i; j++)
					{
						free(spectrogram[i]);
					}
					free(spectrogram);
					spectrogram = 0;
					fprintf(stderr, "Unable to allocate memory at frame %d.\n", i);
				}
			}
			else
			{
				fprintf(stderr, "Unable to allocate memory.\n");
			}
		}
		else
		{
			cached_frames = 0;
		}
		fclose(fp);
	}
	printf(" %d [msec]\n", timeGetTime() - start_time);
	return spectrogram;
}
double** getCTParam(double wave[], int sample_count, int sample_rate,
	double time_axis[], double f0[], int frame_count, int fft_size)
{
	printf("CHEAPTRICK:");
	DWORD start_time = timeGetTime();

	double** spectrogram = (double**)malloc(sizeof(double*) * frame_count);
	if (spectrogram)
	{
		int i, j;
		for (i = 0;i < frame_count;i++)
		{
			spectrogram[i] = (double*)malloc(sizeof(double) * (fft_size / 2 + 1));
			memset(spectrogram[i], 0, sizeof(double) * (fft_size / 2 + 1));
			if (spectrogram[i])
			{
				memset(spectrogram[i], 0, sizeof(double) * (fft_size / 2 + 1));
			}
			else
			{
				break;
			}
		}
		if (i == frame_count)
		{
			CheapTrick(wave, sample_count, sample_rate, time_axis, f0,
				frame_count, spectrogram, fft_size);
		}
		else
		{
			for (j = 0; j < i; j++)
			{
				free(spectrogram[i]);
			}
			free(spectrogram);
			spectrogram = 0;
			fprintf(stderr, "Unable to allocate memory at frame %d.\n", i);
		}
	}
	else
	{
		fprintf(stderr, "Unable to allocate memory.\n");
	}
	printf(" %d [msec]\n", timeGetTime() - start_time);
	return spectrogram;
}
void writeCTParam(int sample_count, int sample_rate, const char* input_name,
	double* spectrogram[], int frame_count, int fft_size)
{
	unsigned short stored_frames;
	unsigned short stored_bins;
	int i, j;
	char cache_name[512];
	makeFilename(input_name, ".ctspec", cache_name);


	printf("write .ctspec:");
	DWORD start_time = timeGetTime();

	short max_value = -32767, min_value = 32767;
	//FILE *ft = fopen("star0.txt", "wt");
	FILE* fp = fopen(cache_name, "wb");
	if (fp)
	{
		fwrite("world-ct", 1, 8, fp);
		stored_frames = (unsigned short)frame_count;
		stored_bins = (unsigned short)fft_size / 2 + 1;
		fwrite(&sample_count, sizeof(int), 1, fp);
		fwrite(&sample_rate, sizeof(int), 1, fp);
		fwrite(&stored_frames, sizeof(unsigned short), 1, fp);
		fwrite(&stored_bins, sizeof(unsigned short), 1, fp);
		for (i = 0; i < frame_count; i++)
		{
			for (j = 0; j <= fft_size / 2; j++)
			{
				int invalid_class;
				if ((invalid_class = _fpclass(spectrogram[i][j]) & 0x0087) != 0)
				{
					spectrogram[i][j] = 0;
#ifdef _DEBUG
					printf("invalid[%d][%d]=%04x!\n", i, j, invalid_class);
#endif
				}
				unsigned short packed_value = (unsigned short)
					(log(spectrogram[i][j] * (2048.0 * 2048 * 2048) + 1) * 1024.0 + 0.5);
				fwrite(&packed_value, sizeof(unsigned short), 1, fp);
				if (max_value < packed_value) max_value = packed_value;
				if (min_value > packed_value) min_value = packed_value;
			}
		}
		fclose(fp);
	}
	printf(" %d [msec]\n", timeGetTime() - start_time);
	printf("max = %d, min = %d\n", max_value, min_value);
}
double** readD4CParam(int sample_count, int sample_rate, const char* input_name,
	int frame_count, int fft_size)
{
	int i, j;

	DWORD start_time;
	unsigned short cached_frames = 0;
	unsigned short cached_bins = 0;
	int cached_samples = 0;
	int cached_rate = 0;
	double** aperiodicity = 0;

	char cache_name[512];
	makeFilename(input_name, ".d4c", cache_name);

	printf("read .d4c:");
	start_time = timeGetTime();

	FILE* fp = fopen(cache_name, "rb");
	if (fp)
	{
		char signature[9];
		fread(signature, 1, 8, fp);
		if (strncmp(signature, "wrld-d4c", 8) != 0)
		{
			fclose(fp);
			printf(" bad file.\n");
			return 0;
		}
		fread(&cached_samples, sizeof(int), 1, fp);
		fread(&cached_rate, sizeof(int), 1, fp);
		fread(&cached_frames, sizeof(unsigned short), 1, fp);
		fread(&cached_bins, sizeof(unsigned short), 1, fp);
		if (cached_frames == frame_count && (fft_size / 2 + 1) == cached_bins &&
			sample_count == cached_samples && sample_rate == cached_rate)
		{
			aperiodicity = (double**)malloc(frame_count * sizeof(double*));
			if (aperiodicity)
			{
				for (i = 0; i < frame_count; i++)
				{
					aperiodicity[i] = (double*)malloc((fft_size / 2 + 1) * sizeof(double));
					if (aperiodicity[i])
					{
						memset(aperiodicity[i], 0, (fft_size / 2 + 1) * sizeof(double));
						for (j = 0; j <= fft_size / 2; j++)
						{
							short packed_value;
							fread(&packed_value, sizeof(short), 1, fp);
							aperiodicity[i][j] = packed_value * 3.90625E-03;
						}
					}
					else
					{
						break;
					}
				}
				if (i < frame_count)
				{
					for (j = 0; j < i; j++)
					{
						free(aperiodicity[i]);
					}
					free(aperiodicity);
					aperiodicity = 0;
					fprintf(stderr, "Unable to allocate memory at frame %d.\n", i);
				}
			}
			else
			{
				fprintf(stderr, "Unable to allocate memory.\n");
			}
		}
		fclose(fp);
	}
	printf(" %d [msec]\n", timeGetTime() - start_time);
	return aperiodicity;
}
double** getD4CParam(double wave[], int sample_count, int sample_rate,
	double time_axis[], double f0[], int frame_count, int fft_size)
{
	printf("D4C:");
	DWORD start_time = timeGetTime();

	double** aperiodicity = (double**)malloc(sizeof(double*) * frame_count);
	if (aperiodicity)
	{
		int i, j;
		for (i = 0;i < frame_count;i++)
		{
			aperiodicity[i] = (double*)malloc(sizeof(double) * (fft_size / 2 + 1));
			memset(aperiodicity[i], 0, sizeof(double) * (fft_size / 2 + 1));
			if (aperiodicity[i])
			{
				memset(aperiodicity[i], 0, sizeof(double) * (fft_size / 2 + 1));
			}
			else
			{
				break;
			}
		}
		if (i == frame_count)
		{
			D4C(wave, sample_count, sample_rate, time_axis, f0,
				frame_count, fft_size, aperiodicity);
		}
		else
		{
			for (j = 0; j < i; j++)
			{
				free(aperiodicity[i]);
			}
			free(aperiodicity);
			aperiodicity = 0;
			fprintf(stderr, "Unable to allocate memory at frame %d.\n", i);
		}
	}
	else
	{
		fprintf(stderr, "Unable to allocate memory.\n");
	}
	printf(" %d [msec]\n", timeGetTime() - start_time);
	return aperiodicity;
}
void writeD4CParam(int sample_count, int sample_rate, const char* input_name,
	double* aperiodicity[], int frame_count, int fft_size)
{
	int i, j;
	DWORD start_time;
	printf("write .d4c:");
	start_time = timeGetTime();

	unsigned short stored_frames = 0;
	unsigned short stored_bins = 0;
	short max_value = -32767, min_value = 32767;

	char cache_name[512];
	makeFilename(input_name, ".d4c", cache_name);

	FILE* fp = fopen(cache_name, "wb");
	if (fp)
	{
		fwrite("wrld-d4c", 1, 8, fp);
		stored_frames = (unsigned short)frame_count;
		stored_bins = (unsigned short)(fft_size / 2 + 1);
		fwrite(&sample_count, sizeof(int), 1, fp);
		fwrite(&sample_rate, sizeof(int), 1, fp);
		fwrite(&stored_frames, sizeof(unsigned short), 1, fp);
		fwrite(&stored_bins, sizeof(unsigned short), 1, fp);
		for (i = 0; i < frame_count; i++)
		{
			for (j = 0; j <= fft_size / 2; j++)
			{
				int invalid_class;
				if ((invalid_class = _fpclass(aperiodicity[i][j]) & 0x0087) != 0)
				{
					aperiodicity[i][j] = 0;
#ifdef _DEBUG
					printf("invalid[%d][%d]=%04x!\n", i, j, invalid_class);
#endif
				}
				short packed_value = (short)(aperiodicity[i][j] * 256.0);
				// Saturate defensively. Valid analysis data should already fit.
				if (packed_value > 32767) packed_value = 32767;
				else if (packed_value < -32768) packed_value = -32768;
				fwrite(&packed_value, sizeof(short), 1, fp);
				if (max_value < packed_value) max_value = packed_value;
				if (min_value > packed_value) min_value = packed_value;
			}
		}
		fclose(fp);
	}
	printf("max = %d, min = %d\n", max_value, min_value);
	printf(" %d [msec]\n", timeGetTime() - start_time);
}
int readDIOParam(const char* input_name, double* out_time_axis[],
	double* out_f0[], int* out_sample_rate, int* out_sample_count)
{
	char cache_name[512];
	int i;
	int sample_count = 0, sample_rate = 0, frame_count = 0;
	double* time_axis = 0;
	double* f0 = 0;

	makeFilename(input_name, ".dio", cache_name);

	DWORD start_time;

	printf("read .dio:");
	start_time = timeGetTime();

	FILE* fp = fopen(cache_name, "rb");
	if (fp)
	{
		char signature[9];
		fread(signature, 8, 1, fp);
		if (strncmp(signature, "wrld-dio", 8) != 0)
		{
			fclose(fp);
			printf(" bad file.\n");
			return 0;
		}
		fread(&sample_count, sizeof(int), 1, fp);
		fread(&sample_rate, sizeof(int), 1, fp);
		fread(&frame_count, sizeof(int), 1, fp);
		if (frame_count > 0)
		{
			time_axis = (double*)malloc(frame_count * sizeof(double));
			f0 = (double*)malloc(frame_count * sizeof(double));
			if (time_axis && f0)
			{
				for (i = 0; i < frame_count; i++)
				{
					fread(&(time_axis[i]), sizeof(double), 1, fp);
					fread(&(f0[i]), sizeof(double), 1, fp);
				}
			}
			else
			{
				if (time_axis) free(time_axis);
				if (f0) free(f0);
				time_axis = 0;
				f0 = 0;
				frame_count = 0;
				fprintf(stderr, "Unable to allocate memory.\n");
			}
		}
		fclose(fp);
	}
	*out_time_axis = time_axis;
	*out_f0 = f0;
	*out_sample_rate = sample_rate;
	*out_sample_count = sample_count;
	printf(" %d [msec]\n", timeGetTime() - start_time);
	return frame_count;
}
int getDIOParam(double wave[], int sample_count, int sample_rate,
	double frame_period, double* out_time_axis[], double* out_f0[],
	char* input_name)
{
	DWORD start_time;
	printf("DIO:");
	start_time = timeGetTime();
	int frame_count = GetSamplesForDIO(sample_rate, sample_count, frame_period);
	int i;
	double utau_avg_f0, avg_f0;
	double* time_axis = (double*)malloc(frame_count * sizeof(double));
	double* f0 = (double*)malloc(frame_count * sizeof(double));
	double* utau_f0 = (double*)malloc(frame_count * sizeof(double));
	double* refined_f0 = (double*)malloc(frame_count * sizeof(double));
	if (time_axis && f0 && utau_f0 && refined_f0)
	{
		if (!hasUTAUfrq(input_name))
		{
			printf("There is no UTAU frequency file. Generating one now.\n");
			if (writeUTAUfrq(wave, sample_count, sample_rate, input_name) != 0)
			{
				fprintf(stderr, "Error: Failed to create a .frq file.\n");
				exit(-60);
			}
		}
		utau_avg_f0 = getUTAUfrq(input_name, frame_count, sample_rate,
			FRAMEPERIOD, utau_f0);
		if (utau_avg_f0 == 0.0)
		{
			fprintf(stderr, "Error: No .frq file was found.\n");
			exit(-60);
		}
		if (utau_avg_f0 < FLOOR_F0) return EXIT_FAILURE;
		Dio(wave, sample_count, sample_rate, time_axis, f0, 0);
		printf("%f", f0[0]);
		avg_f0 = getFreqAvg(f0, frame_count);
		if (0.95 * utau_avg_f0 > avg_f0 || avg_f0 > 1.05 * utau_avg_f0)
			avg_f0 = utau_avg_f0;
		for (i = 0;i < frame_count;i++)
			if (0.95 * avg_f0 > f0[i] || f0[i] > 1.05 * avg_f0)
				f0[i] = utau_f0[i];
		StoneMask(wave, sample_count, sample_rate, time_axis, f0,
			frame_count, refined_f0);
		for (int i = 0; i < frame_count; ++i)
			f0[i] = refined_f0[i];
	}
	else
	{
		fprintf(stderr, "Unable to allocate memory.\n");
		if (time_axis) free(time_axis);
		if (f0) free(f0);
		if (refined_f0) free(refined_f0);
		if (utau_f0) free(utau_f0);
		time_axis = 0;
		f0 = 0;
		utau_f0 = 0;
		refined_f0 = 0;
		frame_count = 0;
	}
	*out_time_axis = time_axis;
	*out_f0 = f0;
	printf(" %d [msec]\n", timeGetTime() - start_time);
	return frame_count;
}
int writeDIOParam(int sample_count, int sample_rate, int frame_count,
	const char* input_name, double time_axis[], double f0[])
{
	char cache_name[512];
	makeFilename(input_name, ".dio", cache_name);


	printf("write .dio");

	FILE* fp = fopen(cache_name, "wb");
	if (fp)
	{
		fwrite("wrld-dio", 1, 8, fp);
		fwrite(&sample_count, sizeof(int), 1, fp);
		fwrite(&sample_rate, sizeof(int), 1, fp);
		fwrite(&frame_count, sizeof(int), 1, fp);
		int i = 0;
		for (i = 0; i < frame_count; i++)
		{
			int invalid_class;
			// Replace NaN, infinity, and denormal values with unvoiced F0.
			if ((invalid_class = _fpclass(f0[i]) & 0x0087) != 0)
			{
#ifdef _DEBUG
				printf("invalid[%d]=%04x!\n", i, invalid_class);
#endif
				f0[i] = 0;
			}
			fwrite(&(time_axis[i]), sizeof(double), 1, fp);
			fwrite(&(f0[i]), sizeof(double), 1, fp);
		}
		fclose(fp);
	}
	return 0;
}
void freeSpecgram(double** spectrogram, int frame_count)
{
	int i;
	if (spectrogram && frame_count)
	{
		for (i = 0; i < frame_count; i++)
		{
			free(spectrogram[i]);
		}
		free(spectrogram);
	}
}
// Stretch the spectral envelope along the frequency axis.
void stretchSpectrum(double** spectrogram, int frame_count, double ratio,
	int sample_rate, int fft_size)
{
	int i, j;
	if (ratio != 1.0)
	{
		double* source_freq, * target_freq;
		double* source_spec, * target_spec;
		source_freq = (double*)malloc(sizeof(double) * fft_size);
		target_freq = (double*)malloc(sizeof(double) * fft_size);
		source_spec = (double*)malloc(sizeof(double) * fft_size);
		target_spec = (double*)malloc(sizeof(double) * fft_size);

		// Prepare the source and target frequency axes.
		for (i = 0;i <= fft_size / 2;i++)
		{
			source_freq[i] = (double)i / (double)fft_size * (double)sample_rate * ratio;
			target_freq[i] = (double)i / (double)fft_size * (double)sample_rate;
		}
		for (i = 0;i < frame_count;i++)
		{
			for (j = 0;j <= fft_size / 2;j++)
				source_spec[j] = log(spectrogram[i][j]);
			interp1(source_freq, source_spec, fft_size / 2 + 1,
				target_freq, fft_size / 2 + 1, target_spec);
			for (j = 0;j <= fft_size / 2;j++)
				spectrogram[i][j] = exp(target_spec[j]);
			if (ratio < 1.0)
			{
				for (j = int((double)fft_size / 2 * ratio);j <= fft_size / 2;j++)
				{
					spectrogram[i][j] =
						spectrogram[i][(int)((double)fft_size / 2 * ratio) - 1];
				}
			}
		}

		free(source_spec); free(target_spec);
		free(source_freq); free(target_freq);
	}
}
void makeHeader(char* header, int sample_count, int sample_rate, int bit_depth)
{
	memcpy(header, "RIFF", 4);
	*(long*)(header + 4) = sample_count * 2 + 36;
	memcpy(header + 8, "WAVE", 4);
	memcpy(header + 12, "fmt ", 4);
	*(long*)(header + 16) = 16;

	*(short*)(header + 20) = 1;
	*(unsigned short*)(header + 22) = 1;
	*(unsigned long*)(header + 24) = sample_rate;
	*(unsigned long*)(header + 28) = sample_rate * bit_depth / 8;
	*(unsigned short*)(header + 32) = bit_depth / 8;
	*(unsigned short*)(header + 34) = bit_depth;
	memcpy(header + 36, "data", 4);
	*(long*)(header + 40) = sample_count * 2;
}

// Build spectra for waveform post-processing.
void createWaveSpec(double* wave, int sample_count, int fft_size,
	int spectrum_count, fft_complex** wave_spectrum)
{
	int i, j;

	double* wave_buffer;
	fft_plan forward_fft;
	fft_complex* spectrum;
	wave_buffer = (double*)malloc(sizeof(double) * fft_size);
	spectrum = (fft_complex*)malloc(sizeof(fft_complex) * fft_size);
	forward_fft = fft_plan_dft_r2c_1d(fft_size, wave_buffer, spectrum, FFT_ESTIMATE);

	int offset;

	for (i = 0;i < spectrum_count;i++)
	{
		offset = i * fft_size / 2;
		// Copy the frame and apply a Hann window.
		for (j = 0;j < fft_size; j++) wave_buffer[j] = wave[offset + j] *
			(0.5 - 0.5 * cos(2.0 * PI * (double)j / (double)fft_size));

		// Run the forward FFT.
		fft_execute(forward_fft);

		// Store the spectrum.
		for (j = 0;j < fft_size / 2 + 1; j++)
		{
			wave_spectrum[i][j][0] = spectrum[j][0];
			wave_spectrum[i][j][1] = spectrum[j][1];
		}
	}

	fft_destroy_plan(forward_fft);
	free(wave_buffer);
	free(spectrum);

}

// Rebuild a waveform from overlapping spectra.
void rebuildWave(double* wave, int sample_count, int fft_size,
	int spectrum_count, fft_complex** wave_spectrum)
{
	int i, j;
	double* wave_buffer;
	fft_plan inverse_fft;
	fft_complex* spectrum;
	wave_buffer = (double*)malloc(sizeof(double) * fft_size);
	spectrum = (fft_complex*)malloc(sizeof(fft_complex) * fft_size);
	inverse_fft = fft_plan_dft_c2r_1d(fft_size, spectrum, wave_buffer, FFT_ESTIMATE);

	int offset;
	for (i = 0;i < sample_count;i++) wave[i] = 0;

	for (i = 0;i < spectrum_count;i++)
	{
		offset = i * fft_size / 2;

		// Load one spectrum.
		for (j = 0;j < fft_size / 2 + 1; j++)
		{
			spectrum[j][0] = wave_spectrum[i][j][0];
			spectrum[j][1] = wave_spectrum[i][j][1];
		}

		// Run the inverse FFT.
		fft_execute(inverse_fft);

		for (j = 0;j < fft_size; j++) wave_buffer[j] /= fft_size;

		// Add the frame to the output using 50 percent overlap.
		for (j = 0;j < fft_size; j++) wave[offset + j] += wave_buffer[j];

	}

	fft_destroy_plan(inverse_fft);
	free(wave_buffer);
	free(spectrum);

}

// Apply the B (breathiness) flag.
void breath2(double* f0, int frame_count, int sample_rate, double* wave,
	int sample_count, fft_complex** wave_spectrum, int spectrum_count,
	int fft_size, int breath)
{
	int i, j;

	// Prepare the noise FFT.
	double* noise_data;
	double* noise_buffer;
	double* noise_wave;
	fft_plan noise_forward_fft;
	fft_plan noise_inverse_fft;
	fft_complex* noise_spectrum;

	noise_data = (double*)malloc(sizeof(double) * sample_count);
	for (i = 0;i < sample_count; i++)
		noise_data[i] = (double)rand() / (RAND_MAX + 1) - 0.5;
	noise_wave = (double*)malloc(sizeof(double) * sample_count);
	for (i = 0;i < sample_count; i++) noise_wave[i] = 0.0;
	noise_buffer = (double*)malloc(sizeof(double) * fft_size);
	noise_spectrum = (fft_complex*)malloc(sizeof(fft_complex) * fft_size);
	noise_forward_fft = fft_plan_dft_r2c_1d(
		fft_size, noise_buffer, noise_spectrum, FFT_ESTIMATE);
	noise_inverse_fft = fft_plan_dft_c2r_1d(
		fft_size, noise_spectrum, noise_buffer, FFT_ESTIMATE);

	// Prepare storage for the waveform envelope.
	fft_complex* envelope;
	envelope = (fft_complex*)malloc(sizeof(fft_complex) * fft_size);

	int offset;
	double gain;

	int start_bin, mid_bin, end_bin;

	start_bin = (int)(fft_size * 1500 / sample_rate);
	mid_bin = (int)(fft_size * 5000 / sample_rate);
	end_bin = (int)(fft_size * 20000 / sample_rate);

	double f0_pos;
	int f0_start, f0_end;
	double current_f0;
	int bin_start, bin_end;
	double amp_start, amp_end;
	int harmonic;

	for (i = 0; i < spectrum_count; i++)
	{
		offset = i * fft_size / 2;
		// Copy the noise frame and apply a Hann window.
		for (j = 0;j < fft_size; j++) noise_buffer[j] = noise_data[offset + j] *
			(0.5 - 0.5 * cos(2.0 * PI * (double)j / (double)fft_size));

		// Run the forward FFT.
		fft_execute(noise_forward_fft);

		// Estimate a simple log-amplitude envelope.
		for (j = 0;j < fft_size / 2 + 1; j++)
			envelope[j][0] = sqrt(
				wave_spectrum[i][j][0] * wave_spectrum[i][j][0] +
				wave_spectrum[i][j][1] * wave_spectrum[i][j][1]);
		for (j = 0;j < fft_size / 2 + 1; j++)
			envelope[j][0] = log10(envelope[j][0] + 0.00000001);
		for (j = 0;j < fft_size / 2 + 1; j++)
			envelope[j][1] = envelope[j][0];

		f0_pos = max(0, min(frame_count - 1,
			(double)(offset + fft_size / 2) / sample_rate * 1000 / FRAMEPERIOD));
		f0_start = min(frame_count - 2, (int)f0_pos);
		f0_end = f0_start + 1;

		current_f0 = (f0[f0_start] == 0 && f0[f0_end] == 0) ? DEFAULT_F0 :
			(f0[f0_start] == 0) ? f0[f0_end] :
			(f0[f0_end] == 0) ? f0[f0_start] :
			(f0[f0_end] - f0[f0_start]) * (f0_pos - f0_start) + f0[f0_start];

		bin_start = 0;
		amp_start = 0.0;
		j = 0;
		harmonic = 1;
		bin_end = 0;
		for (harmonic = 1;bin_end != fft_size / 2 + 1;harmonic++)
		{
			bin_end = min(fft_size / 2 + 1,
				(int)((double)fft_size / sample_rate * current_f0 * harmonic + 0.5));
			amp_end = envelope[bin_end][1];
			for (j = bin_start;j < bin_end;j++)
			{
				envelope[j][0] = (amp_end - amp_start) /
					(bin_end - bin_start) * (j - bin_start) + amp_start;
			}
			bin_start = bin_end;
			amp_start = amp_end;
		}

		for (j = 0;j < fft_size / 2 + 1; j++)
			envelope[j][0] = pow(10, envelope[j][0]);

		// Shape the noise spectrum with the waveform envelope.
		for (j = 0;j < start_bin; j++)
		{
			noise_spectrum[j][0] = 0.0;
			noise_spectrum[j][1] = 0.0;
		}

		for (;j < mid_bin; j++)
		{
			gain = envelope[j][0] *
				(0.5 - 0.5 * cos(PI * (j - start_bin) / (double)(mid_bin - start_bin)));
			noise_spectrum[j][0] *= gain;
			noise_spectrum[j][1] *= gain;
		}
		for (;j < end_bin; j++)
		{
			gain = envelope[j][0] *
				(0.5 - 0.5 * cos(PI + PI * (j - mid_bin) / (double)(end_bin - mid_bin)));
			noise_spectrum[j][0] *= gain;
			noise_spectrum[j][1] *= gain;
		}

		for (;j < fft_size / 2 + 1; j++)
		{
			noise_spectrum[j][0] = 0.0;
			noise_spectrum[j][1] = 0.0;
		}

		noise_spectrum[0][1] = 0.0;
		noise_spectrum[fft_size / 2][1] = 0.0;

		// Run the inverse FFT.
		fft_execute(noise_inverse_fft);
		for (j = 0;j < fft_size; j++) noise_buffer[j] /= fft_size;

		// Add the noise using 50 percent overlap.
		for (j = 0;j < fft_size; j++)
		{
			noise_wave[offset + j] += noise_buffer[j] * 0.2;
		}
	}

	// Mix the shaped noise into the waveform.
	double noise_ratio = max(0, (double)(breath - 50) / 50.0);
	double wave_ratio = 1 - noise_ratio;
	for (i = 0;i < sample_count;i++)
		wave[i] = wave[i] * wave_ratio + noise_wave[i] * noise_ratio;

	// Release temporary FFT data.
	fft_destroy_plan(noise_forward_fft);
	fft_destroy_plan(noise_inverse_fft);
	free(noise_wave);
	free(noise_data);
	free(noise_buffer);
	free(noise_spectrum);
	free(envelope);
}

// Apply the O (vocal opening/strength) flag.
void Opening(double* f0, int frame_count, int sample_rate,
	fft_complex** wave_spectrum, int spectrum_count, int fft_size, int opening)
{
	int i, j;
	double open_ratio = (double)opening / 100.0;
	int low_bin = (int)(fft_size * 500 / sample_rate);
	int high_bin = (int)(fft_size * 2000 / sample_rate);
	double low_gain_db = -10.0;
	double high_gain_db = 10.0;

	// Build the gain curve for each frequency bin.
	double gain;
	double* gain_map;
	gain_map = (double*)malloc(sizeof(double) * (fft_size / 2 + 1));

	gain = pow(10, low_gain_db * open_ratio / 20);
	for (j = 0;j < low_bin;j++)
	{
		gain_map[j] = gain;
	}
	for (;j < high_bin;j++)
	{
		gain = pow(10, ((0.5 + 0.5 *
			cos(PI + PI / (high_bin - low_bin) * (j - low_bin))) *
			(high_gain_db - low_gain_db) + low_gain_db) * open_ratio / 20);
		gain_map[j] = gain;
	}
	gain = pow(10, high_gain_db * open_ratio / 20);
	for (;j < fft_size / 2 + 1;j++)
	{
		gain_map[j] = gain;
	}

	// Apply the curve to voiced spectra.
	int f0_frame;
	for (i = 0;i < spectrum_count;i++)
	{
		f0_frame = max(0, min(frame_count - 1,
			(int)((double)((i + 1) * fft_size / 2) /
				sample_rate * 1000 / FRAMEPERIOD + 0.5)));
		if (f0[f0_frame] == 0.0) continue;
		for (j = 0;j < fft_size / 2 + 1;j++)
		{
			wave_spectrum[i][j][0] *= gain_map[j];
			wave_spectrum[i][j][1] *= gain_map[j];
		}
	}

	free(gain_map);
}

// Feed-forward comb filter used for the lower half of the B flag.
// y[n] = x[n] - coefficient * x[n - samples_per_period]
// y[n]= x[n] - c*x[n-SampleNumOfWaveLength]
void FeedForwardCombFilter(double* wave, int sample_count, double sample_rate,
	double* f0, int frame_count, double coefficient)
{
	// Temporary storage for the delayed waveform.
	// This is larger than Frame2Sample(2, FRAMEPERIOD, sample_rate).
	double delayed_frame[MAX_FFT_LENGTH];
	const int margin = 3; // Enough extra samples for spline interpolation.
	// Apply a light low-pass filter to make the delayed signal more breath-like.
	const double lowpass = 0.05; // [0 .. 0.2 1.0]
	double delayed_sample = 0.0;


	// Calculate the initial position.
	frame_count = imin(frame_count,
		GetFramesForDIO(sample_rate, sample_count, FRAMEPERIOD));
	int frame_index = frame_count - 1;
	int end_sample = sample_count;
	int start_sample =
		static_cast<int>(Frame2Sample(frame_index, FRAMEPERIOD, sample_rate));
	int block_len = end_sample - start_sample;
	double last_f0 = 300.0;
	for (int i = 0; i < frame_count; i++) {
		if (FLOOR_F0 <= f0[i]) last_f0 = f0[i];
	}
	double period_samples = Frequency2Sample(last_f0, sample_rate);
	int period_floor = static_cast<int>(period_samples);
	double interp_offset = margin - (period_samples - period_floor);
	int interp_start = start_sample - period_floor - margin;
	int interp_len = margin + block_len + margin;

	// Work backward so consonants can use the following vowel's F0.
	frame_index--;
	while (0 <= interp_start) {
		// Build the delayed waveform.
		itrp1Qtrim_clip(interp_offset, &wave[interp_start], interp_len,
			1.0, block_len, delayed_frame);

		// Apply the feed-forward comb filter and light low-pass filter.
		for (int i = end_sample - 1, j = block_len - 1;
			start_sample <= i; i--, j--) {
			delayed_sample = delayed_frame[j] + lowpass * delayed_sample;
			wave[i] -= coefficient * delayed_sample;
		}

		// Calculate the next position.
		if (FLOOR_F0 <= f0[frame_index]) {
			period_samples = Frequency2Sample(f0[frame_index], sample_rate);
			period_floor = static_cast<int>(period_samples);
			interp_offset = margin - (period_samples - period_floor);
		}
		else {
			// Keep the most recent voiced period.
		}
		start_sample =
			static_cast<int>(Frame2Sample(frame_index, FRAMEPERIOD, sample_rate));
		end_sample =
			static_cast<int>(Frame2Sample(frame_index + 1, FRAMEPERIOD, sample_rate));
		block_len = end_sample - start_sample;
		interp_start = start_sample - period_floor - margin;
		interp_len = margin + block_len + margin;

		frame_index--;
	}
}  // FeedForwardCombFilter


int main(int argc, char* argv[])
{
	// Enable this while debugging memory leaks.
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

	int i, j;
	if (argc <= 4)
	{
		fprintf(stderr, "Error: Invalid number of arguments.\n");
		return 0;
	}

	FILE* file;
	double read_sample_rate;
	int sample_rate, bit_depth = 16, channels;
	int legacy_g_flag = 0;
	if (argc > 5)
	{
		legacy_g_flag = strchr(argv[5], 'G') != 0;
	}
	int rebuild_cache = 0;
	if (argc > 5)
	{
		rebuild_cache = strchr(argv[5], 'R') != 0;
	}

	double* input_wave = 0, * f0, * time_axis, * output_wave;
	double** spectrogram = 0;
	double** aperiodicity = 0;
	int fft_size;

	int input_samples;
	int input_frames;

	input_frames = readDIOParam(argv[1], &time_axis, &f0,
		&sample_rate, &input_samples);
	if (rebuild_cache != 0)
	{
		input_frames = 0;
	}
	if (input_frames != 0)
	{
		fft_size = getFFTLengthForStar(sample_rate);
		spectrogram = readCTParam(input_samples, sample_rate, argv[1],
			input_frames, fft_size);
		if (spectrogram)
		{
			aperiodicity = readD4CParam(input_samples, sample_rate, argv[1],
				input_frames, fft_size);
			if (!aperiodicity)
			{
				input_frames = 0;
			}
		}
		else
		{
			input_frames = 0;
		}
	}
	if (input_frames == 0)
	{
		freeSpecgram(spectrogram, input_frames);
		freeSpecgram(aperiodicity, input_frames);

		input_wave = audio_read(argv[1], &read_sample_rate, &bit_depth,
			&channels, &input_samples);
		sample_rate = static_cast<int>(read_sample_rate);
		if (input_wave == NULL)
		{
			fprintf(stderr, "Error: The input file does not exist or is unsupported.\n");
			return 0;
		}
		input_frames = getDIOParam(input_wave, input_samples, sample_rate,
			FRAMEPERIOD, &time_axis, &f0, argv[1]);
		if (input_frames != 0)
		{
			writeDIOParam(input_samples, sample_rate, input_frames,
				argv[1], time_axis, f0);
		}
		else
		{
			free(input_wave);
			fprintf(stderr, "Error: Failed to create DIO parameters.\n");
			return 0;
		}
		fft_size = getFFTLengthForStar(sample_rate);
		spectrogram = getCTParam(input_wave, input_samples, sample_rate,
			time_axis, f0, input_frames, fft_size);
		if (!spectrogram)
		{
			free(input_wave);
			free(time_axis);
			free(f0);
			fprintf(stderr, "Error: Failed to create CheapTrick parameters.\n");
			return 0;
		}
		else
		{
			writeCTParam(input_samples, sample_rate, argv[1],
				spectrogram, input_frames, fft_size);
		}
		aperiodicity = getD4CParam(input_wave, input_samples, sample_rate,
			time_axis, f0, input_frames, fft_size);
		if (!aperiodicity)
		{
			free(input_wave);
			free(time_axis);
			free(f0);
			free(spectrogram);
			fprintf(stderr, "Error: Failed to create D4C parameters.\n");
			return 0;
		}
		else
		{
			writeD4CParam(input_samples, sample_rate, argv[1],
				aperiodicity, input_frames, fft_size);
		}
	}

	// Do not modify F0 until every analysis stage has finished.

	// Read the timing arguments.
	double offset_ms = 0;
	double input_ms = (double)input_samples / (double)sample_rate * 1000;
	double output_ms = input_ms;
	double fixed_ms = 0;
	double cutoff_ms = 0;
	double velocity_ratio = 1.0;
	double velocity_value = 100;
	if (argc > 4) sscanf(argv[4], "%lf", &velocity_value);
	velocity_ratio = pow(2, velocity_value / 100 - 1.0);
	if (argc > 6) sscanf(argv[6], "%lf", &offset_ms);
	if (argc > 7) sscanf(argv[7], "%lf", &output_ms);
	if (argc > 8) sscanf(argv[8], "%lf", &fixed_ms);
	if (argc > 9) sscanf(argv[9], "%lf", &cutoff_ms);
#ifdef _DEBBUG
	printf("Parameters\n");
	printf("velocity      :%lf\n", velocity_ratio);
	printf("offset        :%lf\n", offset_ms);
	printf("request length:%lf\n", output_ms);
	printf("fixed         :%lf\n", fixed_ms);
	printf("cutoff        :%lf\n", cutoff_ms);
#endif
	// Time mapping:
	//   offset    fixed          body       cutoff
	//|--------|---------------|------------|-------| input
	//         |               |              |
	//         |   fixed_out   |    body_out     |
	//         |---------------|--------------------|  output
	// fixed_out = fixed / velocity
	// body_out = body / stretch
	// fixed_out + body_out = requested output length
	if (cutoff_ms < 0)
	{
		cutoff_ms = input_ms - offset_ms + cutoff_ms;
		if (cutoff_ms < 0) cutoff_ms = 0;
	}
	if (offset_ms + cutoff_ms >= input_ms)
	{
		fprintf(stderr, "Error: Invalid offset or cutoff.\n");
		free(input_wave);
		return 0;
	}
	if (offset_ms + cutoff_ms + fixed_ms >= input_ms)
	{
		fixed_ms = input_ms - offset_ms + cutoff_ms;
	}
	double fixed_out_ms, body_out_ms;
	double body_in_ms = input_ms - offset_ms - fixed_ms - cutoff_ms;

	fixed_out_ms = fixed_ms / velocity_ratio;
	body_out_ms = output_ms - fixed_out_ms;
	if (body_in_ms <= 0 && body_out_ms > 0)
	{
		fprintf(stderr, "Error: The stretchable input region is empty.\n");
		free(input_wave);
		return 0;
	}
	double stretch_ratio = body_in_ms / body_out_ms;
	if (stretch_ratio > 1.0) stretch_ratio = 1.0;

	int output_samples = (int)(output_ms * 0.001 * sample_rate + 1);
	int output_frames =
		GetSamplesForDIO(sample_rate, output_samples, FRAMEPERIOD);

	printf("File information\n");
	printf("Sampling : %d Hz %d Bit\n", sample_rate, bit_depth);
	printf("Input:\n");
	printf("Length %d [sample]\n", input_samples);
	printf("Length %f [sec]\n", (double)input_samples / (double)sample_rate);
	printf("Output:\n");
	printf("Length %d [sample]\n", output_samples);
	printf("Length %f [sec]\n", (double)output_samples / (double)sample_rate);


	double cent_shift = 0;
	char* flag_pos;
	if (argc > 5 && (flag_pos = strchr(argv[5], 't')) != 0)
	{
		sscanf(flag_pos + 1, "%lf", &cent_shift);
	}
	int use_stretch = 1; // The e flag selects looping instead of stretching.
	if (argc > 5 && (flag_pos = strchr(argv[5], 'e')) != 0)
	{
		use_stretch = 0;
	}
	int breath = 50;
	if (argc > 5 && (flag_pos = strchr(argv[5], 'B')) != 0)
	{
		sscanf(flag_pos + 1, "%d", &breath);
		breath = max(0, min(100, breath));
	}
	int opening = 0;
	if (argc > 5 && (flag_pos = strchr(argv[5], 'O')) != 0)
	{
		sscanf(flag_pos + 1, "%d", &opening);
		opening = max(-100, min(100, opening));
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

	double target_f0 = name2freq(argv[3], cent_shift);
	double avg_f0 = getFreqAvg(f0, input_frames);

#ifdef _DEBUG
	printf("volume        :%lf\n", volume);
	printf("modulation    :%lf\n", modulation);
	printf("target frequency     :%lf\n", target_f0);
	printf("input frequency(avg.):%lf\n", avg_f0);
#endif

	double* output_f0 = (double*)malloc(output_frames * sizeof(double));
	memset(output_f0, 0, sizeof(double) * output_frames);
	int* pitch_cents = NULL;
	double tempo = 120;
	int pitch_count = output_frames;
	int pitch_step = 256;
	if (argc > 13)
	{
		flag_pos = argv[12];
		sscanf(flag_pos + 1, "%lf", &tempo);
		pitch_step = (int)(60.0 / 96.0 / tempo * sample_rate + 0.5);
		pitch_count = output_samples / pitch_step + 1;
		pitch_cents = (int*)malloc((pitch_count + 1) * sizeof(int));
		memset(pitch_cents, 0, (pitch_count + 1) * sizeof(int));
		decpit(argv[13], pitch_cents, pitch_count);
	}
	else
	{
		pitch_cents = (int*)malloc((pitch_count + 1) * sizeof(int));
		memset(pitch_cents, 0, (pitch_count + 1) * sizeof(int));
	}

	double** output_spectrogram =
		(double**)malloc(sizeof(double*) * output_frames);
	double** output_aperiodicity =
		(double**)malloc(sizeof(double*) * output_frames);
	for (i = 0;i < output_frames;i++)
	{
		output_spectrogram[i] =
			(double*)malloc(sizeof(double) * (fft_size / 2 + 1));
		memset(output_spectrogram[i], 0,
			sizeof(double) * (fft_size / 2 + 1));
		output_aperiodicity[i] =
			(double*)malloc(sizeof(double) * (fft_size / 2 + 1));
		memset(output_aperiodicity[i], 0,
			sizeof(double) * (fft_size / 2 + 1));
	}
	// Output-frame mapping state.
	double output_time, input_time, tail_stretch, loop_span;
	int loop_total, reverse = 0;
	int loop_count = 1;
	DWORD start_time = timeGetTime();
	printf("\nTransform\n");
#ifdef _DEBUG
	FILE* time_log = fopen("time.txt", "wt");
	FILE* f0_log = fopen("dio.txt", "wt");
	FILE* spectrum_log = fopen("star.txt", "wt");
	FILE* aperiodicity_log = fopen("plat.txt", "wt");
#endif

	if (use_stretch == 0)
	{
		loop_span = body_in_ms - FRAMEPERIOD / 2.0;
		loop_total = (int)floor(body_out_ms / loop_span);
		loop_total--;
		if (loop_total > 0)
			tail_stretch =
				body_in_ms / (body_out_ms - loop_span * loop_total);
		else
			use_stretch = 1;
	}

	double frame_frac, pitch_frac;
	int n, m;
	if (use_stretch == 0)
	{
		for (i = 0; i < output_frames; i++)
		{
			output_time = FRAMEPERIOD * i;
			if (output_time < fixed_out_ms)
			{
				input_time = offset_ms + output_time * velocity_ratio;
			}
			else if (output_time < fixed_out_ms + loop_span * loop_total)
			{
				if (reverse == 0)
				{
					input_time = offset_ms + fixed_ms +
						(output_time - fixed_out_ms - loop_span * (loop_count - 1));
					if (output_time - fixed_out_ms > loop_span * loop_count)
					{
						loop_count++;
						reverse = 1;
					}
				}
				else
				{
					input_time = input_ms - cutoff_ms -
						(output_time - fixed_out_ms - loop_span * (loop_count - 1));
					if (output_time - fixed_out_ms > loop_span * loop_count)
					{
						loop_count++;
						reverse = 0;
					}
				}
			}
			else
			{
				if (reverse == 1)
					input_time = offset_ms + fixed_ms +
						(output_time - fixed_out_ms - loop_span * loop_total) *
						tail_stretch;
				else
					input_time = input_ms - cutoff_ms -
						(output_time - fixed_out_ms - loop_span * loop_total) *
						tail_stretch;
			}
			frame_frac = input_time / FRAMEPERIOD;
			n = (int)floor(frame_frac);
			frame_frac -= n;

			double input_f0 = f0[n];
			if (n < input_frames - 1) {
				double next_f0 = f0[n + 1];
				if (input_f0 != 0 || next_f0 != 0)
				{
					if (input_f0 == 0) input_f0 = avg_f0;
					if (next_f0 == 0) next_f0 = avg_f0;
					input_f0 = input_f0 * (1.0 - frame_frac) +
						next_f0 * frame_frac;
				}
			}

			pitch_frac = output_time * 0.001 * sample_rate / pitch_step;
			m = (int)floor(pitch_frac);
			pitch_frac -= m;
			if (m >= pitch_count) {
				m = pitch_count - 1;
				frame_frac = 0.0;
			}
			output_f0[i] = target_f0 * pow(2,
				(pitch_cents[m] * (1.0 - pitch_frac) +
					pitch_cents[m + 1] * pitch_frac) / 1200.0);
			output_f0[i] *= pow(input_f0 / avg_f0, modulation * 0.01);

			for (j = 0; j <= fft_size / 2; j++)
			{
				if (n < input_frames - 1)
				{
					output_spectrogram[i][j] =
						spectrogram[n][j] * (1.0 - frame_frac) +
						spectrogram[n + 1][j] * frame_frac;
				}
				else
				{
					output_spectrogram[i][j] =
						spectrogram[input_frames - 1][j];
				}
			}

			int nearest_frame = n;
			if (frame_frac > 0.5) nearest_frame++;
			for (j = 0; j <= fft_size / 2; j++)
			{
				if (nearest_frame < input_frames)
				{
					output_aperiodicity[i][j] =
						aperiodicity[nearest_frame][j];
				}
				else
				{
					output_aperiodicity[i][j] =
						aperiodicity[input_frames - 1][j];
				}
			}
		}
	}
	else
	{
		for (i = 0; i < output_frames; i++)
		{
			output_time = FRAMEPERIOD * i;
			if (output_time < fixed_out_ms)
			{
				input_time = offset_ms + output_time * velocity_ratio;
			}
			else
			{
				input_time = offset_ms + fixed_ms +
					(output_time - fixed_out_ms) * stretch_ratio;
			}
#ifdef _DEBUG
			fprintf(time_log, "%0.6lf\t%0.6lf\n", input_time, output_time);
#endif
			frame_frac = input_time / FRAMEPERIOD;
			n = (int)floor(frame_frac);
			frame_frac -= n;

			double input_f0 = f0[n];
			if (n < input_frames - 1) {
				double next_f0 = f0[n + 1];
				if (input_f0 != 0 || next_f0 != 0)
				{
					if (input_f0 == 0) input_f0 = avg_f0;
					if (next_f0 == 0) next_f0 = avg_f0;
					input_f0 = input_f0 * (1.0 - frame_frac) +
						next_f0 * frame_frac;
				}
			}

			pitch_frac = output_time * 0.001 * sample_rate / pitch_step;
			m = (int)floor(pitch_frac);
			pitch_frac -= m;
			if (m >= pitch_count) {
				m = pitch_count - 1;
				frame_frac = 0.0;
			}
			output_f0[i] = target_f0 * pow(2,
				(pitch_cents[m] * (1.0 - pitch_frac) +
					pitch_cents[m + 1] * pitch_frac) / 1200.0);
			output_f0[i] *= pow(input_f0 / avg_f0, modulation * 0.01);

#ifdef _DEBUG
			fprintf(f0_log, "%lf\n", output_f0[i]);
#endif
			for (j = 0; j <= fft_size / 2; j++)
			{
				if (n < input_frames - 1)
				{
					output_spectrogram[i][j] =
						spectrogram[n][j] * (1.0 - frame_frac) +
						spectrogram[n + 1][j] * frame_frac;
				}
				else
				{
					output_spectrogram[i][j] =
						spectrogram[input_frames - 1][j];
				}
#ifdef _DEBUG
				fprintf(spectrum_log, "%lf\t",
					output_spectrogram[i][j] * 1000000);
#endif
			}
#ifdef _DEBUG
			fprintf(spectrum_log, "\n");
#endif
			int nearest_frame = n;
			if (frame_frac > 0.5) nearest_frame++;
			for (j = 0; j <= fft_size / 2; j++)
			{
				if (nearest_frame < input_frames)
				{
					output_aperiodicity[i][j] =
						aperiodicity[nearest_frame][j];
				}
				else
				{
					output_aperiodicity[i][j] =
						aperiodicity[input_frames - 1][j];
				}
			}
#ifdef _DEBUG
			for (j = 0; j < fft_size / 2; j += 8)
			{
				fprintf(aperiodicity_log, "%lf\t", output_aperiodicity[i][j]);
			}
			fprintf(aperiodicity_log, "\n");
			for (j = 0; j < fft_size / 2; j += 8)
			{
				fprintf(aperiodicity_log, "%lf\t",
					output_aperiodicity[i][j + 1]);
			}
			fprintf(aperiodicity_log, "\n");
#endif
		}
	}

#ifdef _DEBUG
	fclose(time_log);
	fclose(f0_log);
	fclose(spectrum_log);
	fclose(aperiodicity_log);
#endif
	// Apply spectral stretching for the g flag.
	if (argc > 5 && (flag_pos = strchr(argv[5], 'g')) != 0)
	{
		double gender = 0;
		double spectrum_ratio = 1.0;
		sscanf(flag_pos + 1, "%lf", &gender);
		if (gender > 100) gender = 100;
		if (gender < -100) gender = -100;
		spectrum_ratio = pow(10, -gender / 200);
		stretchSpectrum(output_spectrogram, output_frames,
			spectrum_ratio, sample_rate, fft_size);
	}
	printf("TRANSFORM: %d [msec]\n", timeGetTime() - start_time);

	// Synthesize the transformed WORLD parameters.
	output_wave = (double*)malloc(sizeof(double) * output_samples);
	memset(output_wave, 0, sizeof(double) * output_samples);

	printf("\nSynthesis\n");
	start_time = timeGetTime();
	Synthesis(output_f0, output_frames, output_spectrogram,
		output_aperiodicity, fft_size, FRAMEPERIOD, sample_rate,
		output_samples, output_wave);
	printf("WORLD: %d [msec]\n", timeGetTime() - start_time);

	// Prepare spectra for the B and O post-processing flags.
	int effect_fft_size = 1024;
	int effect_frames = (output_samples / (effect_fft_size / 2)) - 1;
	fft_complex** effect_spectrum;
	effect_spectrum =
		(fft_complex**)malloc(sizeof(fft_complex*) * effect_frames);
	for (i = 0;i < effect_frames;i++)
		effect_spectrum[i] = (fft_complex*)malloc(
			sizeof(fft_complex) * (effect_fft_size / 2 + 1));

	if (breath > 50 || opening != 0)
	{
		createWaveSpec(output_wave, output_samples, effect_fft_size,
			effect_frames, effect_spectrum);
	}

	// Apply vocal opening/strength.
	if (opening != 0)
	{
		Opening(output_f0, output_frames, sample_rate, effect_spectrum,
			effect_frames, effect_fft_size, opening);
	}

	// Rebuild the waveform after spectral changes.
	if (opening != 0)
	{
		rebuildWave(output_wave, output_samples, effect_fft_size,
			effect_frames, effect_spectrum);
	}

	// Apply breath/noise processing.
	if (breath > 50)
	{
		breath2(output_f0, output_frames, sample_rate, output_wave,
			output_samples, effect_spectrum, effect_frames,
			effect_fft_size, breath);
	}
	else if (breath < 50)
	{
		FeedForwardCombFilter(output_wave, output_samples, sample_rate,
			output_f0, output_frames, breath * 0.01 + 0.005);
	}

	// Convert the waveform to PCM and write the output file.
	char header[44];
	short* output_pcm;
	double max_amp;

	output_pcm = (short*)malloc(sizeof(short) * output_samples);
	// Find the peak amplitude for normalization.
	max_amp = 0.0;
#ifdef _DEBUG
	{
		FILE* synthesis_log = fopen("synthesis.txt", "wt");
		for (i = 0;i < output_samples;i++)
		{
			fprintf(synthesis_log, "%f\n", output_wave[i]);
		}
		fclose(synthesis_log);
	}
#endif
	for (i = 0;i < output_samples;i++)
	{
		if (!_isnan(output_wave[i]))
		{
			max_amp = max_amp < fabs(output_wave[i]) ?
				fabs(output_wave[i]) : max_amp;
		}
	}
	double peak_strength = 0.86;
	if (argc > 5)
	{
		flag_pos = strchr(argv[5], 'P');
		if (flag_pos)
		{
			sscanf(flag_pos + 1, "%lf", &peak_strength);
			if (peak_strength < 0) peak_strength = 0;
			else if (peak_strength > 100) peak_strength = 100;
			peak_strength *= 0.01;
		}
	}
	double peak_gain =
		3 * 32.0 * pow(512.0 / max_amp, peak_strength);
	for (i = 0;i < output_samples;i++)
	{
		// WORLD can return NaN in silent regions; write silence instead.
		double sample_value = _isnan(output_wave[i]) ?
			0 : output_wave[i] * peak_gain * volume;
		if (sample_value > 32767.0) sample_value = 32767.0;
		else if (sample_value < -32767.0) sample_value = -32767.0;
		output_pcm[i] = (short)sample_value;
	}

	makeHeader(header, output_samples, sample_rate, bit_depth);

	file = fopen(argv[2], "wb");
	fwrite(header, sizeof(char), 44, file);
	fwrite(output_pcm, sizeof(short), output_samples, file);
	fclose(file);

	printf("complete.\n");

	// Release all run data.
	free(output_pcm);
	free(input_wave); free(time_axis); free(f0); free(output_wave);
	for (i = 0;i < input_frames;i++)
	{
		free(spectrogram[i]);
		free(aperiodicity[i]);
	}
	free(spectrogram);
	free(aperiodicity);

	free(output_f0);
	for (i = 0;i < output_frames;i++)
	{
		free(output_spectrogram[i]);
		free(output_aperiodicity[i]);
	}
	free(output_spectrogram);
	free(output_aperiodicity);
	free(pitch_cents);

	for (i = 0;i < effect_frames;i++) free(effect_spectrum[i]);
	free(effect_spectrum);

	printf("complete.\n");

	return 0;
}
