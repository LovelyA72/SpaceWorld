// FFT-domain convolution helpers.
#include "../world.h"

void fast_fftfilt(const double *x, int x_length, const double *h, int h_length,
    int fft_size, const ForwardRealFFT *forward_real_fft,
    const InverseRealFFT *inverse_real_fft, double *y) {
  fft_complex *x_spectrum = new fft_complex[fft_size];

  for (int i = 0; i < x_length; ++i)
    forward_real_fft->waveform[i] = x[i] / fft_size;
  for (int i = x_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i) {
    x_spectrum[i][0] = forward_real_fft->spectrum[i][0];
    x_spectrum[i][1] = forward_real_fft->spectrum[i][1];
  }

  for (int i = 0; i < h_length; ++i)
    forward_real_fft->waveform[i] = h[i] / fft_size;
  for (int i = h_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);

  for (int i = 0; i <= fft_size / 2; ++i) {
    inverse_real_fft->spectrum[i][0] =
      x_spectrum[i][0] * forward_real_fft->spectrum[i][0] -
      x_spectrum[i][1] * forward_real_fft->spectrum[i][1];
    inverse_real_fft->spectrum[i][1] =
      x_spectrum[i][0] * forward_real_fft->spectrum[i][1] +
      x_spectrum[i][1] * forward_real_fft->spectrum[i][0];
  }
  fft_execute(inverse_real_fft->inverse_fft);

  for (int i = 0; i < fft_size; ++i)
    y[i] = inverse_real_fft->waveform[i];

  delete[] x_spectrum;
}
