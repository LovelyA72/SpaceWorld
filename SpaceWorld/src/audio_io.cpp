// I do not know well and/or compatibility with MATLAB.
// So programs are turned into the most original.
//               /2012/12/06/   Custom.Maid
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

#include "audio_io.h"
//#include "matlabfunctions.h"

#if (defined (__WIN32__) || defined (_WIN32)) && !defined (__MINGW32__)
  #include <windows.h>
  #ifdef _MSC_VER          // Borland C �ł͕s�v
    // Warning 4996 is hidden
    #pragma warning( disable : 4996 )
  #endif
#endif


namespace {

/****************************************************
 *
 * �@����CPU�����g���^�r�b�O�G���f�B�A���ł��邩�ǂ����Ɋ֌W�Ȃ�
 * �v�Z������@��v�������̂ŁA�������g�����Ƃɂ����B
 *
 * �@���Ȃ݂�CPU�̓��g���^�r�b�O�G���f�B�A���ȊO�ɂ�F�X�ȕ��т�����
 * �炵���̂ł����A����Ƃ�������тɑΉ����Ă���͂��B��
 *
 */

/****************************************************
 *
 * �@aiff�̏ꍇ�̓^�O�`�����N�̃T�C�Y�ɊԈ�����l�������Ă鎖������
 * �̂ł͖����A�`�����N�͕K��2�o�C�g�ŃA���C������Ă�炵��
 * ������̃\�[�X�R�[�h�Q�Ƃ����Ă��������܂���
 *      http://code.google.com/p/bitspersampleconv2/
 *      http://src.gnu-darwin.org/
 *
 * �����Ă���Șb�����̂ŁAwave�ɂ�t���Ă���܂�
 *      http://www.kk.iij4u.or.jp/~kondo/wave/
 *
 * ���������Ƃ�
 *
 */

bool CheckLittleEndian(void)
{
  unsigned long checkNum = 1;
  unsigned char *result = reinterpret_cast<unsigned char *>(&checkNum);

  if(*result) { // Little Endian
    return true;
  } else { // Big Endian
    return false;
  }
} // CheckLittleEndian

void long2memAsBigEndian(long num, unsigned char *array)
{
  for (int i = 3; 0 <= i; i--) {
    array[i] = num & 0xff;
    num /= 256;
  }
} // long2memAsBigEndian

void long2memAsLittleEndian(long num, unsigned char *array)
{
  for (int i = 0; i < 4; i++) {
    array[i] = num & 0xff;
    num /= 256;
  }
} // long2memAsLittleEndian

unsigned long BigEndian2ulong(unsigned char *array)
{
  unsigned long num = 0;
  for (int i = 0; i < 4; i++) {
    num = num * 256 + array[i];
  }
  return num;
} // BigEndian2ulong

unsigned long LittleEndian2ulong(unsigned char *array)
{
  unsigned long num = 0;
  for (int i = 3; 0 <= i; i--) {
    num = num * 256 + array[i];
  }
  return num;
} // LittleEndian2ulong

// lazy conversion routine
// if exponent is out of range, return NaN
// fraction that is out of range, simply truncates
// extendedFloat80bit should be arranged in the order of big endian
float extendedFloat80bit2float(unsigned char *extendedFloat80bit)
{
  // Of course, both "long" and "float" should be 32bit long.
  int sign, exponent, fraction;
  unsigned char signSet, exponentSet;
  unsigned char floatSet[4];
  unsigned long floatPntNum;
  float *floatPnt = reinterpret_cast<float *>(&floatPntNum);

  if (0x80 & extendedFloat80bit[0]) {
    sign = -1;
    signSet = 0x80;
  } else {
    sign = 1;
    signSet = 0x00;
  }

  exponent = (0x7f & extendedFloat80bit[0]) * 256 +
                     extendedFloat80bit[1] - 0x3fff;
  if ((exponent < -0x7f) || (0x7f <exponent)) {
    exponentSet = 0xff; // set for NaN
  } else {
    exponentSet = exponent + 0x7f;
  }

  // msb of extendedFloat80bit[2] should be 1. and do not use this one.
  fraction = ((0x7f & extendedFloat80bit[2]) * 256 +
                      extendedFloat80bit[3]) * 256 +
                      extendedFloat80bit[4];

  floatSet[0] = signSet | (exponentSet >> 1);
  floatSet[1] = (exponentSet << 7) | (0x7f & extendedFloat80bit[2]);
  floatSet[2] = extendedFloat80bit[3];
  floatSet[3] = extendedFloat80bit[4];
  floatPntNum = BigEndian2ulong(floatSet);

  return *floatPnt;
} // extendedFloat80bit2float


} // end of namespace


// only signed 8,16,24bit PCM.  mono and stereo
double *aiff_read(char *filename, double *fs, int *Nbit, int *Ch, int *sNum)
{
  FILE *fp;
  char getString[4];
  unsigned char getLongNumber[4];
  unsigned char getExtendedNumber[10];
  long chunkLength, offset, blockSize;
  long Nbyte, iTmp;
  double normalizeNum;

  fp = fopen(filename, "rb");
  if (NULL == fp) {
    printf("File not found. %s\n", filename);
    return NULL;
  }

  // Header
  fread(getString, sizeof(char), 4, fp);
  if (0 != strncmp("FORM", getString, 4)) {
    fclose(fp);
    printf("FORM error.\n");
    return NULL;
  }
  fseek(fp, 4, SEEK_CUR); // skip chunk size
  fread(getString, sizeof(char), 4, fp);
  if (0 != strncmp("AIFF", getString, 4)) {
    fclose(fp);
    printf("AIFF error.\n");
    return NULL;
  }

  // seek COMM chunk
  fread(getString, sizeof(char), 4, fp); // "COMM"?
  while (0 != strncmp("COMM", getString, 4)) {
    fread(getLongNumber, sizeof(char), 4, fp); // chunk size
    chunkLength = BigEndian2ulong(getLongNumber);
    if (1 == chunkLength % 2) chunkLength++;
/*
    printf("\tsearch COMM chunk: skip %c%c%c%c chunk! %dBytes\n",
           getString[0], getString[1], getString[2], getString[3],
           chunkLength);
 */
    fseek(fp, chunkLength, SEEK_CUR); // skip chunk
    // read next chunk
    fread(getString, sizeof(char), 4, fp); // "COMM"?
    if (feof(fp)) {
      fclose(fp);
      printf("COMM error.\n");
      return NULL;
    }
  }
  fread(getLongNumber, sizeof(char), 4, fp); // chunk size
  chunkLength = BigEndian2ulong(getLongNumber);
  if (18 != chunkLength) { // chunk size shouled be : 00 00 00 12
    fclose(fp);
    printf("COMM chunk size error.\n");
    return NULL;
  }

  // Channels
  fread(getLongNumber, sizeof(char), 2, fp);
  *Ch = getLongNumber[1]; // ignore value of the other byte
  if (2 < *Ch) { // multi channel
    fclose(fp);
    printf("Cannot read %dchannel file\n", *Ch);
    return NULL;
  }
  fseek(fp, 4, SEEK_CUR); // skip sample frames

  // Quantization
  fread(getLongNumber, sizeof(char), 2, fp);
  *Nbit = getLongNumber[1]; // ignore value of the other byte
  if ((24 < *Nbit) || (0 != *Nbit % 8)) { // not 8,16,24bit
    fclose(fp);
    printf("Cannot read %dbit PCM file\n", *Nbit);
    return NULL;
  }
  Nbyte = *Nbit / 8;
  normalizeNum = pow(2.0, *Nbit - 1);

  // Sampling frequency
  fread(getExtendedNumber, sizeof(char), 10, fp);
  *fs = extendedFloat80bit2float(getExtendedNumber);
  if (*fs <= 0.0) {
    fclose(fp);
    printf("Can not calc Sampling frequency %f\n",
           extendedFloat80bit2float(getExtendedNumber));
    return NULL;
  }

  // seek SSND chunk
  fread(getString, sizeof(char), 4, fp); // "SSND"?
  while (0 != strncmp("SSND", getString, 4)) {
    fread(getLongNumber, sizeof(char), 4, fp); // chunk size
    chunkLength = BigEndian2ulong(getLongNumber);
    if (1 == chunkLength % 2) chunkLength++;
/*
    printf("\tsearch SSND chunk: skip %c%c%c%c chunk! %dBytes\n",
           getString[0], getString[1], getString[2], getString[3],
           chunkLength);
 */
    fseek(fp, chunkLength, SEEK_CUR); // skip chunk
    // read next chunk
    fread(getString, sizeof(char), 4, fp); // "SSND"?
    if (feof(fp)) {
      fclose(fp);
      printf("SSND error.\n");
      return NULL;
    }
  }

  // The number of sample.
  fread(getLongNumber, sizeof(char), 4, fp); // data size
  *sNum = BigEndian2ulong(getLongNumber);
  if (1 == *sNum % 2) *sNum++;
  *sNum -= 8; // delete of offset/blockSize
  *sNum /= Nbyte;
  *sNum /= *Ch;

  fread(getLongNumber, sizeof(char), 4, fp); // offset
  offset = BigEndian2ulong(getLongNumber);
  if (0 != offset) {
    fclose(fp);
    printf("SSND offset should be zero error.\n");
    return NULL;
  }

  fread(getLongNumber, sizeof(char), 4, fp); // blockSize
  blockSize = BigEndian2ulong(getLongNumber);
  if (0 != blockSize) {
    fclose(fp);
    printf("SSND blockSize should be zero error.\n");
    return NULL;
  }

  double *waveForm = new double[*sNum];
  if (waveForm == NULL) return NULL;

  for (int i = 0; i < *sNum; i++) {
    // for mono or L chanel
    fread(&getLongNumber[4-Nbyte], sizeof(char),
          Nbyte, fp);
    if (0x80 & getLongNumber[4-Nbyte]) { // check sign
      for (int j = 0; j < 4-Nbyte; j++) getLongNumber[j] = 0xff;
    } else {
      for (int j = 0; j < 4-Nbyte; j++) getLongNumber[j] = 0x00;
    }
    iTmp = BigEndian2ulong(getLongNumber);
    waveForm[i] = iTmp / normalizeNum;

    if (1 == *Ch) continue;
    // for R chanel
    fread(&getLongNumber[4-Nbyte], sizeof(char),
          Nbyte, fp);
    if (0x80 & getLongNumber[4-Nbyte]) { // check sign
      for (int j = 0; j < 4-Nbyte; j++) getLongNumber[j] = 0xff;
    } else {
      for (int j = 0; j < 4-Nbyte; j++) getLongNumber[j] = 0x00;
    }
    iTmp = BigEndian2ulong(getLongNumber);
    // calc. average
    waveForm[i] = (waveForm[i] + iTmp / normalizeNum) / 2.0;
  }

  fclose(fp);
  return waveForm;
} // aiff_read

// only unsigned8, signed 16,24bit PCM and 32bit float.  mono and stereo
double *wave_read(char *filename, double *fs, int *Nbit, int *Ch, int *sNum)
{
  FILE *fp;
  char getString[4];
  unsigned char getLongNumber[4];
  unsigned long floatPntNum;
  float *floatPnt = (float *)&floatPntNum;
  long chunkLength, extdataLength;
  long Nbyte, iTmp;
  double normalizeNum;
  bool flag_float = false;

  fp = fopen(filename, "rb");
  if (NULL == fp) {
    printf("File not found. %s\n", filename);
    return NULL;
  }

  // Header
  fread(getString, sizeof(char), 4, fp);
  if (0 != strncmp("RIFF", getString, 4)) {
    fclose(fp);
    printf("RIFF error.\n");
    return NULL;
  }
  fseek(fp, 4, SEEK_CUR); // skip file size
  fread(getString, sizeof(char), 4, fp);
  if (0 != strncmp("WAVE", getString, 4)) {
    fclose(fp);
    printf("WAVE error.\n");
    return NULL;
  }
  fread(getString, sizeof(char), 4, fp);
  if (0 != strncmp("fmt ", getString, 4)) {
    fclose(fp);
    printf("fmt error.\n");
    return NULL;
  }
  fread(getLongNumber, sizeof(char), 4, fp);
  chunkLength = LittleEndian2ulong(getLongNumber);
  extdataLength = chunkLength - 16; // normaly chunk size: 10 00 00 00

  // Data Format
  fread(getString, sizeof(char), 2, fp);
  if ( !((1 == getString[0]) && (0 == getString[1])) ) { // PCM: 01 00
    if ( !((3 == getString[0]) && (0 == getString[1])) ) { // IEEE float: 03 00
      fclose(fp);
      printf("Format ID error.\n");
      return NULL;
    } else {
      flag_float = true;
    }
  }

  // Channels
  fread(getLongNumber, sizeof(char), 2, fp);
  *Ch = getLongNumber[0]; // ignore value of the other byte
  if (2 < *Ch) { // multi channel
    fclose(fp);
    printf("Cannot read %dchannel file\n", *Ch);
    return NULL;
  }

  // Sampling frequency
  fread(getLongNumber, sizeof(char), 4, fp);
  *fs = LittleEndian2ulong(getLongNumber);
  fseek(fp, 6, SEEK_CUR); // skip bytePerSec, blockAlign

  // Quantization
  fread(getLongNumber, sizeof(char), 2, fp);
  *Nbit = getLongNumber[0]; // ignore value of the other byte
  if (flag_float) {
    if (32 != *Nbit) { // not 32bit float data
      fclose(fp);
      printf("Cannot read %dbit float file\n", *Nbit);
      return NULL;
    }
  } else if ((0 != *Nbit % 8) || (24 < *Nbit)) { // not 8,16,24bit PCM data
    fclose(fp);
    printf("Cannot read %dbit PCM file\n", *Nbit);
    return NULL;
  }
  Nbyte = *Nbit / 8;
  normalizeNum = pow(2.0, *Nbit - 1);

  // skip extra data
  fseek(fp, extdataLength, SEEK_CUR);

  // check data chunk
  fread(getString, sizeof(char), 4, fp); // "data"?
  while (0 != strncmp("data", getString, 4)) {
    fread(getLongNumber, sizeof(char), 4, fp); // chunk size
    chunkLength = LittleEndian2ulong(getLongNumber);
    if (1 == chunkLength % 2) chunkLength++;
/*
    printf("\tsearch data chunk: skip %c%c%c%c chunk! %dBytes\n",
           getString[0], getString[1], getString[2], getString[3],
           chunkLength);
 */
    fseek(fp, chunkLength, SEEK_CUR); // skip chunk
    // read next chunk
    fread(getString, sizeof(char), 4, fp); // "data"?
    if (feof(fp)) {
      fclose(fp);
      printf("data error.\n");
      return NULL;
    }
  }

  // The number of sample.
  fread(getLongNumber, sizeof(char), 4, fp); // data size
  *sNum = LittleEndian2ulong(getLongNumber);
  if (1 == *sNum % 2) *sNum++;
  *sNum /= Nbyte;
  *sNum /= *Ch;

  double *waveForm = new double[*sNum];
  if (waveForm == NULL) return NULL;

  for (int i = 0; i < *sNum; i++) {
    fread(getLongNumber, sizeof(char), Nbyte, fp);
    if (flag_float) { // float data
      floatPntNum = LittleEndian2ulong(getLongNumber);
      waveForm[i] = *floatPnt;
    } else { // PCM data
      if (8 == *Nbit) { // unsigned char
        iTmp = getLongNumber[0];
        iTmp -= 0x80;
      } else {
        if (0x80 & getLongNumber[Nbyte-1]) { // check sign
          for (int j = Nbyte; j < 4; j++) getLongNumber[j] = 0xff;
        } else {
          for (int j = Nbyte; j < 4; j++) getLongNumber[j] = 0x00;
        }
        iTmp = LittleEndian2ulong(getLongNumber);
      }
      waveForm[i] = iTmp / normalizeNum;
  }

    if (1 == *Ch) continue;
    // for R chanel
    fread(getLongNumber, sizeof(char), Nbyte, fp);
    if (flag_float) { // float data
      floatPntNum = LittleEndian2ulong(getLongNumber);
      // calc. average
      waveForm[i] = (waveForm[i] + *floatPnt) / 2.0;
    } else { // PCM data
      if (8 == *Nbit) { // unsigned char
        iTmp = getLongNumber[0];
        iTmp -= 0x80;
      } else {
        if (0x80 & getLongNumber[Nbyte-1]) { // check sign
          for (int j = Nbyte; j < 4; j++) getLongNumber[j] = 0xff;
        } else {
          for (int j = Nbyte; j < 4; j++) getLongNumber[j] = 0x00;
        }
        iTmp = LittleEndian2ulong(getLongNumber);
      }
      // calc. average
      waveForm[i] = (waveForm[i] + iTmp / normalizeNum) / 2.0;
    }
  }

  fclose(fp);
  return waveForm;
} // wave_read

// overall read program
double *audio_read(char *filename, double *fs, int *Nbit, int *Ch, int *sNum)
{
  FILE *fp;
  char getString[4];
  double *waveForm;

  fp = fopen(filename, "rb");
  if (NULL == fp) {
    printf("File not found. %s\n", filename);
    return NULL;
  }

  // Header Check
  fread(getString, sizeof(char), 4, fp);
  fclose(fp);
  if        (0 == strncmp("FORM", getString, 4)) {
    waveForm = aiff_read(filename, fs, Nbit, Ch, sNum);
  } else if (0 == strncmp("RIFF", getString, 4)) {
    waveForm = wave_read(filename, fs, Nbit, Ch, sNum);
  } else {
    printf("Not supported file.\n");
    return NULL;
  }

  return waveForm;
} // audio_read
