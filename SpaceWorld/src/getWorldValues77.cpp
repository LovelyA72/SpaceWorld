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


// Calculate the DIO frame count so callers can allocate F0 and time-axis buffers.
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


// Calculate the FFT size required for the current sampling frequency.
// The required size also depends on the lowest expected F0.
// Return zero when the required size exceeds MAX_FFT_LENGTH.
int GetFFTSizeForStar(double fs)
{
  int fftl = GetSuitableFFTSize(static_cast<int>(3.0 *
                                Frequency2Sample(FLOOR_F0, fs) + 1.0));
  if (MAX_FFT_LENGTH < fftl) fftl = 0;
  return fftl;
} // GetFFTSizeForStar


