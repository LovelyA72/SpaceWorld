#ifndef SPACEWORLD_UTILS_WINDOW_FUNCTIONS_H_
#define SPACEWORLD_UTILS_WINDOW_FUNCTIONS_H_

void Hamming_window(double *w, int N);
void Hanning_window(double *w, int N);
void Nuttall_window(double *w, int N);
void sinc(double *t, int tNum, double *c);
void combsort_ascend(double *x, int sNum);
void combsort_descend(double *x, int sNum);
double median(double *x, int sNum);

#endif  // SPACEWORLD_UTILS_WINDOW_FUNCTIONS_H_
