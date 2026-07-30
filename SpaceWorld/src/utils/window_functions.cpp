// Window and sample-array utility functions.
#include "window_functions.h"
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>

void Hamming_window(double *w, int N) {
  for (int i = 0; i < N; i++)
    w[i] = 0.54 - 0.46 * cos(2.0 * M_PI * (i + (N % 2 ? 0.5 : 1.0)) /
        (N % 2 ? N : N + 1));
}

void Hanning_window(double *w, int N) {
  for (int i = 0; i < N; i++)
    w[i] = 0.5 - 0.5 * cos(2.0 * M_PI * (i + (N % 2 ? 0.5 : 1.0)) /
        (N % 2 ? N : N + 1));
}

void Nuttall_window(double *w, int N) {
  for (int i = 0; i < N; i++) {
    double tmp = M_PI * (i + (N % 2 ? 0.5 : 1.0)) / (N % 2 ? N : N + 1);
    w[i] = 0.355768 - 0.487396 * cos(2.0 * tmp) +
        0.144232 * cos(4.0 * tmp) - 0.012604 * cos(6.0 * tmp);
  }
}

void sinc(double *t, int tNum, double *c) {
  for (int i = 0; i < tNum; i++) {
    if (0.0 == t[i]) {
      c[i] = 1.0;
    } else {
      double tmp = M_PI * t[i];
      c[i] = sin(tmp) / tmp;
    }
  }
}

void combsort_ascend(double *x, int sNum) {
  int gap = sNum;
  bool swapped = false;
  do {
    gap = static_cast<int>(gap / 1.247330950103979);
    if (gap < 1) gap = 1;
    swapped = false;
    for (int i = 0; gap + i < sNum; i++) {
      if (x[gap + i] < x[i]) {
        double tmp = x[i]; x[i] = x[gap + i]; x[gap + i] = tmp;
        swapped = true;
      }
    }
  } while (swapped || (1 < gap));
}

void combsort_descend(double *x, int sNum) {
  int gap = sNum;
  bool swapped = false;
  do {
    gap = static_cast<int>(gap / 1.247330950103979);
    if (gap < 1) gap = 1;
    swapped = false;
    for (int i = 0; gap + i < sNum; i++) {
      if (x[i] < x[gap + i]) {
        double tmp = x[i]; x[i] = x[gap + i]; x[gap + i] = tmp;
        swapped = true;
      }
    }
  } while (swapped || (1 < gap));
}

double median(double *x, int sNum) {
  combsort_ascend(x, sNum);
  return sNum % 2 == 0 ? (x[sNum / 2] + x[sNum / 2 - 1]) * 0.5 : x[sNum / 2];
}
