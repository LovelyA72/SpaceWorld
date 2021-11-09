// I do not know well and/or compatibility with MATLAB.
// So programs are turned into the most original.
//               /2012/12/06/   Custom.Maid
#ifndef WORLD_AUDIO_IO_H_
#define WORLD_AUDIO_IO_H_

// subcontract programs of audio_read
// only signed 8,16,24bit PCM.  mono and stereo
double *aiff_read(char *filename, double *fs, int *Nbit, int *Ch, int *sNum);
// only unsigned8, signed 16,24bit PCM and 32bit float.  mono and stereo
double *wave_read(char *filename, double *fs, int *Nbit, int *Ch, int *sNum);

// overall read program
double *audio_read(char *filename, double *fs, int *Nbit, int *Ch, int *sNum);

// write wave. only 16bit monaural
void wave_write(char *filename, double fs, int Nbit, double *x, int sNum);

#endif  // WORLD_AUDIO_IO_H_
