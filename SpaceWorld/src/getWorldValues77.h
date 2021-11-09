//-----------------------------------------------------------------------------
// WORLD/XIII (WORLD slash excijie) Library  1.8 karat WORLD Library
//      Tuneup and original functions
//
#ifndef WORLD_GET_WORLD_VALUES_77_H_
#define WORLD_GET_WORLD_VALUES_77_H_


// for User set up and/or user convenience
//-----------------------------------------------------------------------------
// GetFramesForDIO() calculates the number of frames required for Dio().
// Input:
//   fs             : Sampling frequency [Hz]
//   sNum           : Number of input signal [Sample].
//   framePeriod    : Frame shift [msec]
// Output:
//   The number of frames required to store the results of Dio()
//-----------------------------------------------------------------------------
int GetFramesForDIO(double fs, double sNum, double framePeriod);

//-----------------------------------------------------------------------------
// GetTimeAxisForDIO() calculates the time axis required for Dio().
// Input:
//   framePeriod    : Frame shift [msec]
//   timeAxis       : Time axis. [msec]
//   fNum           : Number of frames [frame]
// Output:
//   timeAxis       : Frame time of Dio [msec]
//-----------------------------------------------------------------------------
void GetTimeAxisForDIO(double framePeriod, double *timeAxis, int fNum);

//-----------------------------------------------------------------------------
// GetFFTSizeForStar() calculates the FFT size based on the sampling frequency
// and the lower limit of f0 (It is defined in world.h).
// Input:
//   fs             : Sampling Frequency [Hz]
// Output:
//   FFT size [sample]
// Caution:
//   return Max MAX_FFT_LENGTH. 0 for error.
//-----------------------------------------------------------------------------
int GetFFTSizeForStar(double fs);

//-----------------------------------------------------------------------------
//      inline functions
//
//              �}�N�����͏d�����ǂ�A���x�ቺ�͊����Ȃ��̂ŁA�[��������悤
//              ������̕����\�[�X�̌��ʂ���ǂ��Ȃ�Ǝv���܂�
//
//-----------------------------------------------------------------------------
// These four functions are simple max() and min() function
// for "int" and "double" type.
//-----------------------------------------------------------------------------
inline int imax(int x, int y) {
  if (x <= y) return y;
  else        return x;
}

inline int imin(int x, int y) {
  if (x <= y) return x;
  else        return y;
}

inline double fmax(double x, double y) {
  if (x <= y) return y;
  else        return x;
}

inline double fmin(double x, double y) {
  if (x <= y) return x;
  else        return y;
}

//-----------------------------------------------------------------------------
// ���̃v���O�����́A�P�T���v���͎��Ԃ��L���Ȃ��Ƃ݂Ȃ��Ă��܂�
//
// this program considered, one sample does not occupy one time interval
//

//-----------------------------------------------------------------------------
// Convert number of samples to length of the sample.
// Input:
//   sample         : number of samples [sample]
// Output:
//   length of the sample [distance]
// Caution:
//   rounding is left to the discretion of user.
//-----------------------------------------------------------------------------
inline double SampleNum2SampleLen(double sample) {
  if      (1.0 < sample)
    return sample - 1.0;
  else if (sample < -1.0)
    return sample + 1.0;
  else // -1.0 <= sample <= 1.0
    return 0.0;
} // Frame2Sample

//-----------------------------------------------------------------------------
// �ȉ��̃v���O�����́A�P�T���v����1��Ԃ̎��Ԃ��L����Ƃ݂Ȃ��Ă��܂�
//
// following program considered, one sample occupies one time interval
//

//-----------------------------------------------------------------------------
// Convert number of frames to number of samples.
// Input:
//   frame          : number of frames [frame]
//   framePeriod    : Frame shift [msec]
//   fs             : Sampling Frequency [Hz]
// Output:
//   number of samples, or point of samples. [sample]
// Caution:
//   rounding is left to the discretion of user.
//-----------------------------------------------------------------------------
inline double Frame2Sample(double frame, double framePeriod, double fs) {
  return frame * (0.001 * framePeriod * fs);
} // Frame2Sample

//-----------------------------------------------------------------------------
// Convert length of time to number of samples.
// Input:
//   time           : Length of time [msec]
//   fs             : Sampling Frequency [Hz]
// Output:
//   number of samples, or point of samples. [sample]
// Caution:
//   rounding is left to the discretion of user.
//-----------------------------------------------------------------------------
inline double Time2Sample(double time, double fs) {
  return time * 0.001 * fs;
} // Time2Sample

//-----------------------------------------------------------------------------
// Convert one cycle of frequency to number of samples.
// Input:
//   frequency      : Frequency [Hz]
//   fs             : Sampling Frequency [Hz]
// Output:
//   number of samples, or point of samples. [sample]
// Caution:
//   rounding is left to the discretion of user.
//-----------------------------------------------------------------------------
inline double Frequency2Sample(double frequency, double fs) {
  return fs / frequency;
} // Frequency2Sample

//-----------------------------------------------------------------------------
// Convert number of samples to length of time/framePeriod.
// Input:
//   sample         : number of samples [sample]
//   fs             : Sampling Frequency [Hz]
// Output:
//   length of time/framePeriod [msec]
// Caution:
//   Because, N samples occupy [0 .. N).
//   So X [msec] means to occupy [0 .. X).
//-----------------------------------------------------------------------------
inline double Sample2Time(double sample, double fs) {
  return 1000.0 * sample / fs;
} // Sample2Time

//-----------------------------------------------------------------------------
// Convert number of frames to length of time.
// Input:
//   frame          : number of frames [frame]
//   framePeriod    : Frame shift [msec]
// Output:
//   length of time [msec]
// Caution:
//   Because, N frames occupy [0 .. N).
//   So X [msec] means to occupy [0 .. X).
//-----------------------------------------------------------------------------
inline double Frame2Time(double frame, double framePeriod) {
  return frame * framePeriod;
} // Frame2Time

//-----------------------------------------------------------------------------
// Convert one cycle of frequency to length of time/framePeriod.
// Input:
//   frequency      : Frequency [Hz]
// Output:
//   length of time/framePeriod [msec]
//-----------------------------------------------------------------------------
inline double Frequency2Time(double frequency) {
  return 1000.0 / frequency;
} // Frequency2Time

//-----------------------------------------------------------------------------
// Convert number of samples to frequency.
// Input:
//   sample         : number of samples [sample]
//   fs             : Sampling Frequency [Hz]
// Output:
//   Frequency [Hz]
//-----------------------------------------------------------------------------
inline double Sample2Frequency(double sample, double fs) {
  return fs / sample;
} // Sample2Frequency

//-----------------------------------------------------------------------------
// Convert number of samples to number of frames.
// Input:
//   sample         : number of samples [sample]
//   framePeriod    : Frame shift [msec]
//   fs             : Sampling Frequency [Hz]
// Output:
//   number of frames [frame]
// Caution:
//   rounding is left to the discretion of user.
// Note:
//   Because, N samples occupy [0 .. N).
//   So X [frame] means to occupy [0 .. X).
//-----------------------------------------------------------------------------
inline double Sample2Frame(double sample, double framePeriod, double fs) {
  if       (1.0 < sample)
    return (sample - 1.0) / (0.001 * framePeriod * fs) + 1.0;
  else if  (sample < -1.0)
    return (sample + 1.0) / (0.001 * framePeriod * fs) - 1.0;
  else // -1.0 <= sample <= 1.0
    return 0.0;
} // Sample2Frame

//-----------------------------------------------------------------------------
// Convert length of time to number of frames.
// Input:
//   time           : Length of time [msec]
//   framePeriod    : Frame shift [msec]
// Output:
//   number of frames [frame]
// Caution:
//   rounding is left to the discretion of user.
//-----------------------------------------------------------------------------
inline double Time2Frame(double time, double framePeriod) {
  return time / framePeriod;
} // Time2Frame

//-----------------------------------------------------------------------------
// Convert one cycle of frequency to number of frames.
// Input:
//   frequency      : Frequency [Hz]
//   framePeriod    : Frame shift [msec]
// Output:
//   number of frames [frame]
// Caution:
//   rounding is left to the discretion of user.
//-----------------------------------------------------------------------------
inline double Frequency2Frame(double frequency, double framePeriod) {
  return 1000.0 / (frequency * framePeriod);
} // Frequency2Frame

//-----------------------------------------------------------------------------
// Convert cent (1 octave = 1,200cent) to ratio.
// Input:
//   cent           ; number of cent [cent]
// Output:
//   ratio of note [-]
// Note:
//   ratio is (0 .. +inf.).
//-----------------------------------------------------------------------------
inline double Cent2Ratio(double cent) {
  return pow(2.0, cent / 1200.0);
} // Cent2Ratio

//-----------------------------------------------------------------------------
// Convert dB of amplitude to ratio (ratio 1.0 = 0 dB).
// Input:
//   db             ; number of dB [dB]
// Output:
//   ratio of amplitude [-]
// Caution:
//   dB is not power, it is amplitude.
// Note:
//   ratio is (0 .. +inf.).
//-----------------------------------------------------------------------------
inline double dB2Ratio(double db) {
  return pow(10.0, db / 20.0);
} // dB2Ratio

//-----------------------------------------------------------------------------
// Convert ratio to cent (1 octave = 1,200cent).
// Input:
//   ratio          : ratio of note [-]
// Output:
//   number of cent [cent]
// Caution:
//   ratio must be (0 .. +inf.).
// Note:
//   ratio = 1             0 [cent]
//           2          1200 [cent]
//           1/2       -1200 [cent]
//-----------------------------------------------------------------------------
inline double Ratio2Cent(double ratio) {
  return 1200.0 * (log(ratio) / 1);
} // Ratio2Cent   return 1200.0 * (log(ratio) / M_LN2);

//-----------------------------------------------------------------------------
// Convert ratio to dB of amplitude (ratio 1.0 = 0 dB).
// Input:
//   ratio          : ratio of amplitude [-]
// Output:
//   dB of amplitude [dB]
// Caution:
//   ratio must be (0 .. +inf.).
//   dB is not power, it is amplitude.
// Note:
//   ratio = 1          0 [dB]
//           2          6 [dB]
//           1/2       -6 [dB]
//           10        20 [dB]
//           1/10     -20 [dB]
//-----------------------------------------------------------------------------
inline double Ratio2dB(double ratio) {
  return 20.0 * log10(ratio);
} // Ratio2dB

#endif  // WORLD_GET_WORLD_VALUES_77_H_
