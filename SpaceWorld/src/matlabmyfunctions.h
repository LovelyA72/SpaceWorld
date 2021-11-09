//------------------------------------------------------------------------------------
// Matlab functions
//      These are my own functions that I made
#ifndef WORLD_MATLAB_MY_FUNCTIONS_H_
#define WORLD_MATLAB_MY_FUNCTIONS_H_


// ハミング窓
void Hamming_window(double *w, int N);

// ハニング窓
void Hanning_window(double *w, int N);

// ナットール窓
void Nuttall_window(double *w, int N);

// 配列tを元に配列cを作るsinc関数
void sinc(double *t, int tNum, double *c);

// Comb Sort in ascending order
// you must read http://en.wikipedia.org/wiki/Combsort
void combsort_ascend(double *x, int sNum);

// Comb Sort in descending order
void combsort_descend(double *x, int sNum);

// Gets the median
double median(double *x, int sNum);

//------------------------------------------------------------------------------------
// interpolation

// interp1 uses linear interpolation(interp1 default).
// This function make no extrapolation(clip and/or nearest).
void interp1_clip(double *x, double *y, int inNum,
                  double *xo, int outNum, double *yo);

// interp1 uses 2point cubic spline interpolation
//      original idea.
// This function make no extrapolation(clip and/or nearest).
void interp12PC_clip(double *x, double *y, int inNum,
                     double *xo, int outNum, double *yo);

// interp1 uses (T-Spline/trim spline) cubic spline interpolation.
// This function is not interp1 pchip/cubic.
//      original idea.
// This function make no extrapolation(clip and/or nearest).
void interp1trim_clip(double *x, double *y, int inNum,
                      double *xo, int outNum, double *yo);

// interp1 uses (RunOut/Not-a-Knot) cubic spline interpolation.
// interp1 spline default.
//      http://www.vector.co.jp/soft/data/prog/se002453.html
// This function make no extrapolation(clip and/or nearest).
void interp1spline_clip(double *x, double *y, int inNum,
                        double *xo, int outNum, double *yo);

// interp1 uses natural cubic spline interpolation.
// interp1 natural.
//      http://www.vector.co.jp/soft/data/prog/se002453.html
// This function make no extrapolation(clip and/or nearest).
void interp1natural_clip(double *x, double *y, int inNum,
                         double *xo, int outNum, double *yo);

// interp1 uses (Catmull-Rom/cardinal spline a=-1/2) cubic spline interpolation.
// This function is not interp1 pchip/cubic.
//      http://yehar.com/blog/wp-content/uploads/2009/08/deip.pdf
// This function make no extrapolation(clip and/or nearest).
void interp1Catmull_Rom_clip(double *x, double *y, int inNum,
                             double *xo, int outNum, double *yo);

// サンプリング間隔を等間隔に限定し、高速に動作する直線補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
void itrp1Q_clip(double offset, double *y, int inNum, double stride,
                 int outNum, double *yo);

// サンプリング間隔を等間隔に限定し、高速に動作する2PointCubic補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
void itrp1Q2PC_clip(double offset, double *y, int inNum, double stride,
                    int outNum, double *yo);

// サンプリング間隔を等間隔に限定し、高速に動作するトリムスプライン補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
//      両端の区間はクリップに近づくように設定している
void itrp1Qtrim_clip(double offset, double *y, int inNum, double stride,
                            int outNum, double *yo);

// サンプリング間隔を等間隔に限定し、高速に動作するスプライン補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
void itrp1Qspline_clip(double offset, double *y, int inNum, double stride,
                       int outNum, double *yo);

// サンプリング間隔を等間隔に限定し、高速に動作する自然スプライン補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
void itrp1Qnatural_clip(double offset, double *y, int inNum, double stride,
                        int outNum, double *yo);

// サンプリング間隔を等間隔に限定し、高速に動作するCatmull-Romスプライン補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
void itrp1QCatmull_Rom_clip(double offset, double *y, int inNum, double stride,
                            int outNum, double *yo);

#endif  // WORLD_MATLAB_MY_FUNCTIONS_H_
