#include <stdio.h> // for debug
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
 #define _USE_MATH_DEFINES
 #define __func__ __FUNCTION__
#endif
#ifdef __BORLANDC__
 #define __func__ __FUNC__
#endif
#include <math.h>


#include "world.h"
#include "getWorldValues77.h"


namespace {

        // Nothing! Any other functions.

} // end of namespace


// ���O�Ƀ��[�U���������m�ۂł���悤�ɁAF0�O�Ղ̗v�f���𓾂�
// fs                   : Sampling frequency [Hz]
// sNum                 : Number of input signal [sample].
// framePeriod          : Frame shift [msec]
int GetFramesForDIO(double fs, double sNum, double framePeriod) {
  return static_cast<int>(Sample2Frame(sNum, framePeriod, fs));
} // GetFramesForDIO


// You get array datas of time
// framePeriod          : Frame shift [msec]
// timeAxis             : Time axis. [msec]
// fNum                 : Number of frames [frame]
void GetTimeAxisForDIO(double framePeriod, double *timeAxis, int fNum)
{
  for (int i = 0; i < fNum; i++) timeAxis[i] = framePeriod * i;
} // GetTimeAxisForDIO


// �T���v�����O���g������K�v��FFT�����v�Z
// �{���͍��������̍Œ�F0���K�v�����ǖ���
// MAX_FFT_LENGTH�łЂ�������悤�ɂ���
// We can get the length of FFT.
// Essentially, it also depends on the lowest F0 of the input signal.
int GetFFTSizeForStar(double fs)
{
  int fftl = GetSuitableFFTSize(static_cast<int>(3.0 *
                                Frequency2Sample(FLOOR_F0, fs) + 1.0));
  if (MAX_FFT_LENGTH < fftl) fftl = 0;
  return fftl;
} // GetFFTSizeForStar


