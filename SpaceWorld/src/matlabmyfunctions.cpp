// Matlab風味の関数たち　FFT不要なもののみ
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
#include "matlabmyfunctions.h"

namespace {


// 単位区間における三次補間関数の係数を計算する
void calc_coeffOfCubic(double dy1, double vy1, double vy2,
                       double *a, double *b, double *c)
{
  *a = vy1 + vy2 - 2.0 * dy1;
  *b = 3.0 * dy1 - 2.0 * vy1 - vy2;
  *c = vy1;
}  // calc_coeffOfCubic


void make_2PC_table(double alpha, double *y, int inNum,
                    double *a, double *b, double *c)
{
/***********************************************
 *      2点三次曲線補間
 *      2PontCubicSplineの補間用 区間関数
 *      要は三次エルミートスプライン(cubic hermit spline)の一種
 *      use fixed start/end tangent cubic hermit spline
 *
 *      f(x) =  2(1+alph)*x**3 - 3(1+alpha)*x**2 + alpha*x + 1  (0 <= x <= 1)
 *              0                                               (1 <= x)
 *              f(-x)                                           (x <= 0)
 *
 *      alphaはx=0, 1での傾き       (-1 <= alpha <= 0)
 *               alpha= 0の時はコサイン補間を近似
 *                      サンプル点を通る二次のB-Splineにも近い
 *               alpha= -1/2が標準的、直線補間代わりにグッドと思う
 *               alpha= -1の時は直線補間に等しい
 *
 *               (-inf. < a < -1) でも動くが、
 *                      曲線の振動が強くなるし、
 *                      サンプル点では鋭く尖がるしで補間には使えない
 *                      CGで模様を描く時に使えるかもくらいか
 *
 *      エルミートさんのエの字も知らないで考え付いた自分えらい！
 *      と自画自賛しておくｗ
 *
 *      Y1(t) = -2(1+alpha)*dy*t**3 + 3(1+alpha)*dy*t**2 - alpha*dy*t + d
 *              d = y1
 *              t = (x - x1)/(x2 - x1)  (0 <= t <= 1)
 *              dy = y2 - y1
 *
 */
  for (int i = 0; i < inNum-1; i++) {
    // 差分の計算をする
    double dy = y[i+1] - y[i];
    calc_coeffOfCubic(dy, -alpha * dy, -alpha * dy, &a[i], &b[i], &c[i]);
  }
}  // make_2PC_table


void make_trim_table(double *x, double *y, int inNum,
                     double *a, double *b, double *c, double *h)
{
/************************************************************
 *      make (T-Spline / trim spline)
 *              cubic spline coefficient table
 *
 *      Catmull-Rom の傾き、接線は必ず中点が作る角に接する
 *      エッジの所がその為に少し緩やかになるのでオーバーシュートが減る
 *      でも上下のピークも少し緩やかにするので少し太ましくなる
 *
 *      このスプラインはエッジや上下のピークはより緩やかに、
 *      直線に近い部分はより急峻なスプライン
 *      いずれにせよ MATLAB の pchip とはだいぶ違う
 *      エッジの所がほぼ傾きゼロになるのでオーバーシュートがほとんど無くなる
 *      直線に近い所は傾き、接線が中点が作る角を貫通する
 *      結果 Catmull-Rom よりふらふら揺れてるようにも見える
 *      見た目は Catmull-Rom より滑らかそうだったりもするんだけども
 *      オリジナルなんだから maid spline でも善かったんだけどｗ
 *
 *      (x0,y0), (x1,y1), (x2,y2), (x3,y3)      (x1 <= x <= x2)
 *      Y1(t) = a*t**3 + b*t**2 + c*t + d       (d = y1)
 *               t  =  (x - x1)/(x2 - x1)       (0 <= t <= 1)
 *              vy1 = (y2 - y1)/(x2 - x1)
 *      Y'0(1)  = Y'1(0) = vy1*dy1
 *      vy1 = 0                                         (vy1 = 0, or vy0 = 0)
 *          = (vy1+vy0)*(4|vy1|*|vy0|/(|vy1|+|vy0|)**2) (same sign vy1,vy0)
 *          = (vy1+vy0)*(2(|vyL|-|vyS|)*|vyS|/|vyL|**2) (opposite sign vy1,vy0)
 *      Y'1(1)  = Y'2(0) = vy2*dy1
 *      vy2 = 0                                         (vy2 = 0, or vy1 = 0)
 *          = (vy2+vy1)*(4|vy2|*|vy1|/(|vy2|+|vy1|)**2) (same sign vy2,vy1)
 *          = (vy1+vy0)*(2(|vyL|-|vyS|)*|vyS|/|vyL|**2) (opposite sign vy1,vy0)
 *
 *      同符号の場合：調和平均と相加平均の比を掛ける            (0 <= q <= 1)
 *      異符号の場合：きつい傾きと緩い傾きの差 と 緩い傾きとの
 *                    調和平均と相加平均の比を1/2して掛ける     (0 <= q <= 1/2)
 *      ※、当初、異符号の場合も(0 <= q <= 1)としてたが、色々 Catmull-Rom に
 *        負けてたので、調整の結果このように決定した
 *                      trim ratio = 1 としている状態
 *
 *      端の傾きはゼロが善いのか、中間値が善いのか？
 *      取り敢えず Catmull-Rom に合わせる事にした
 *      Y'0(0) = 1/2 * vy0：中間値
 *              つまりは vy1 = 1/2 * vy0  vy(-1) = 7*vy0
 *
 */
  double *dy = new double [inNum];
  double *vy = new double [inNum];
  double *vyAdd = new double [inNum];
  double *vyMul = new double [inNum];

  // 正規化の為にサンプル間隔と傾きの計算をする
  for (int i = 0; i < inNum-1; i++) {
    h[i] = x[i+1] - x[i];
    dy[i] = y[i+1] - y[i];
    vy[i] = dy[i] / h[i];
  }
    // 端点の設定 Y'(inNum-1)(0) = 1/2 * vy[inNum-2]
    h[inNum-1] = h[inNum-1-1];
    dy[inNum-1] = dy[inNum-1-1] / 7.0;
    vy[inNum-1] = dy[inNum-1] / h[inNum-1];

    // 端点の設定 Y'0(0) = 1/2 * vy[0]
    vyAdd[0] = vy[0] + (vy[0] / 7.0);
    vyMul[0] = vy[0] * (vy[0] / 7.0);
  for (int i = 1; i < inNum; i++) {
    vyAdd[i] = vy[i] + vy[i-1];
    vyMul[i] = vy[i] * vy[i-1];
  }

  double vy1, vy2;
  // 最初特別
      // 単純にこれでいいんだが、きちんと意識するために計算する
      // vy1 = vy[0] / 2.0;
    if        (0.0 == vyMul[0]) { // zero either
      vy1 = 0.0;
    } else if (0.0 <  vyMul[0]) { // of course same sign
      vy1 = 4.0 * vyMul[0] / vyAdd[0];
    } else { // it's impossible. inf. or Nan ?!
      vy1 = 0.0;
    }
    if        (0.0 == vyMul[1]) { // zero either
      vy2 = 0.0;
    } else if (0.0 <  vyMul[1]) { // same sign
      vy2 = 4.0 * vyMul[1] / vyAdd[1];
    } else if (0.0 < vy[1] * vyAdd[1]) { // vy[i] < vy[i+1]
      vy2 = -2.0 * vyAdd[1] * vyAdd[1] * vy[0] / (vy[1] * vy[1]);
    } else { // vy[i+1] <= vy[i]
      vy2 = -2.0 * vyAdd[1] * vyAdd[1] * vy[1] / (vy[0] * vy[0]);
    }
    calc_coeffOfCubic(dy[0], vy1 * h[0], vy2 * h[0], // 傾きを正規化している
                       &a[0], &b[0], &c[0]);

  // ループ
  for (int i = 1; i < inNum-1; i++) {
    if        (0.0 == vyMul[i]) { // zero either
      vy1 = 0.0;
    } else if (0.0 <  vyMul[i]) { // same sign
      vy1 = 4.0 * vyMul[i] / vyAdd[i];
    } else if (0.0 < vy[i] * vyAdd[i]) { // vy[i-1] < vy[i]
      vy1 = -2.0 * vyAdd[i] * vyAdd[i] * vy[i-1] / (vy[i] * vy[i]);
    } else { // vy[i] <= vy[i-1]
      vy1 = -2.0 * vyAdd[i] * vyAdd[i] * vy[i] / (vy[i-1] * vy[i-1]);
    }
    if        (0.0 == vyMul[i+1]) { // zero either
      vy2 = 0.0;
    } else if (0.0 <  vyMul[i+1]) { // same sign
      vy2 = 4.0 * vyMul[i+1] / vyAdd[i+1];
    } else if (0.0 < vy[i+1] * vyAdd[i+1]) { // vy[i] < vy[i+1]
      vy2 = -2.0 * vyAdd[i+1] * vyAdd[i+1] * vy[i] / (vy[i+1] * vy[i+1]);
    } else { // vy[i+1] <= vy[i]
      vy2 = -2.0 * vyAdd[i+1] * vyAdd[i+1] * vy[i+1] / (vy[i] * vy[i]);
    }
    calc_coeffOfCubic(dy[i], vy1 * h[i], vy2 * h[i], // 傾きを正規化している
                       &a[i], &b[i], &c[i]);
  }

  delete[] vyMul; delete[] vyAdd; delete[] vy; delete[] dy;
}  // make_trim_table


void make_fast_trim_table(double *y, int inNum,
                          double *a, double *b, double *c)
{
  double *dy = new double [inNum];
  double *vyAdd = new double [inNum];
  double *vyMul = new double [inNum];

  // 差分、かつ、傾きの計算をする
  for (int i = 0; i < inNum-1; i++) {
    dy[i] = y[i+1] - y[i];
  }
    // 端点の設定 Y'(inNum-1)(0) = 1/2 * vy[inNum-2]
    dy[inNum-1] = dy[inNum-1-1] / 7.0;

    // 端点の設定 Y'0(0) = 1/2 * vy[0]
    vyAdd[0] = dy[0] + (dy[0] / 7.0);
    vyMul[0] = dy[0] * (dy[0] / 7.0);
  for (int i = 1; i < inNum; i++) {
    vyAdd[i] = dy[i] + dy[i-1];
    vyMul[i] = dy[i] * dy[i-1];
  }

  double vy1, vy2;
  // 最初特別
      // 単純にこれでいいんだが、きちんと意識するために計算する
      // vy1 = dy[0] / 2.0;
    if        (0.0 == vyMul[0]) { // zero either
      vy1 = 0.0;
    } else if (0.0 <  vyMul[0]) { // of course same sign
      vy1 = 4.0 * vyMul[0] / vyAdd[0];
    } else { // it's impossible. inf. or Nan ?!
      vy1 = 0.0;
    }
    if        (0.0 == vyMul[1]) { // zero either
      vy2 = 0.0;
    } else if (0.0 <  vyMul[1]) { // same sign
      vy2 = 4.0 * vyMul[1] / vyAdd[1];
    } else if (0.0 < dy[1] * vyAdd[1]) { // dy[i] < dy[i+1]
      vy2 = -2.0 * vyAdd[1] * vyAdd[1] * dy[0] / (dy[1] * dy[1]);
    } else { // dy[i+1] <= dy[i]
      vy2 = -2.0 * vyAdd[1] * vyAdd[1] * dy[1] / (dy[0] * dy[0]);
    }
    calc_coeffOfCubic(dy[0], vy1, vy2, &a[0], &b[0], &c[0]);

  // ループ
  for (int i = 1; i < inNum-1; i++) {
    if        (0.0 == vyMul[i]) { // zero either
      vy1 = 0.0;
    } else if (0.0 <  vyMul[i]) { // same sign
      vy1 = 4.0 * vyMul[i] / vyAdd[i];
    } else if (0.0 < dy[i] * vyAdd[i]) { // dy[i-1] < dy[i]
      vy1 = -2.0 * vyAdd[i] * vyAdd[i] * dy[i-1] / (dy[i] * dy[i]);
    } else { // dy[i] <= dy[i-1]
      vy1 = -2.0 * vyAdd[i] * vyAdd[i] * dy[i] / (dy[i-1] * dy[i-1]);
    }
    if        (0.0 == vyMul[i+1]) { // zero either
      vy2 = 0.0;
    } else if (0.0 <  vyMul[i+1]) { // same sign
      vy2 = 4.0 * vyMul[i+1] / vyAdd[i+1];
    } else if (0.0 < dy[i+1] * vyAdd[i+1]) { // dy[i] < dy[i+1]
      vy2 = -2.0 * vyAdd[i+1] * vyAdd[i+1] * dy[i] / (dy[i+1] * dy[i+1]);
    } else { // dy[i+1] <= dy[i]
      vy2 = -2.0 * vyAdd[i+1] * vyAdd[i+1] * dy[i+1] / (dy[i] * dy[i]);
    }
    calc_coeffOfCubic(dy[i], vy1, vy2, &a[i], &b[i], &c[i]);
  }

  delete[] vyMul; delete[] vyAdd; delete[] dy;
}  // make_fast_trim_table


void make_spline_table(double *x, double *y, int inNum,
                       double *a, double *b, double *c, double *h)
{
/************************************************************
 *      make Not-A-Knot cubic spline coefficient table
 *
 *      Y1(t) = a*t**3 + b*t**2 + c*t + d       (d = y1)
 *              t = (x - x1)/(x2 - x1)          (0 <= t <= 1)
 *
 *      Not-A-Knot cubic spline
 *      Y'''0(1) = Y'''1(0),  Y'''inNum-1-1(1) =Y'''inNum-1(0)
 *
 */
  // 三重対角行列と差分のメモリ確保
  double *dl = new double [inNum];
  double *dd = new double [inNum];
  double *du = new double [inNum];
  double *dy = new double [inNum];


  // 三重対角行列の作成
  // 正規化の為にサンプル間隔、差分、傾きの計算をする
  for (int i = 0; i < inNum-1; i++) {
    h[i] = x[i+1] - x[i];
    dy[i] = y[i+1] - y[i];
    c[i] = dy[i] / h[i];
  }

  for (int i = 0; i < inNum-3; i++) {
    dl[i] = h[i+1];
    du[i] = h[i+1];
  }
  dl[inNum-4] -= (h[inNum-2] * h[inNum-2] / h[inNum-3]);

  for (int i = 0; i < inNum-2; i++) {
    dd[i] = 2.0 * (h[i] + h[i+1]);
    b[i]  = c[i+1] - c[i];
  }
  dd[0] += (h[0] + h[0] * h[0] / h[1]);
  du[0] -= (h[0] * h[0] / h[1]);
  dd[inNum-3] += (h[inNum-2] + h[inNum-2] * h[inNum-2] / h[inNum-3]);

  // 三重対角行列による連立方程式の解法計算する
  for (int i = 0; i < inNum-2-1; i++) {
    du[i] /= dd[i];
    dd[i+1] -= dl[i] * du[i];
  }

  b[0] /= dd[0];
  for (int i = 1; i < inNum-2; i++)
    b[i] = (b[i] - dl[i-1] * b[i-1]) / dd[i];

  for (int i = inNum-2-2; 0 <= i; i--)
    b[i] -= b[i+1] * du[i];

  // bを整える：両端の計算と正規化
  for (int i = inNum-3; 0 <= i; i--) {
    b[i+1] = b[i];
  }
  b[0] = (1.0 + h[0] / h[1]) * b[1] -
         h[0] / h[1] * b[2];
  b[inNum-1] = (1.0 + h[inNum-2] / h[inNum-3]) * b[inNum-2] -
              h[inNum-2] / h[inNum-3] * b[inNum-3];
  for (int i = 0; i < inNum; i++) {
    b[i] *= h[i] * h[i];
  }

  // 係数を計算する
  for (int i = 0; i < inNum-1; i++) {
    a[i] = (b[i+1] - b[i]);
    c[i] = dy[i] - (b[i+1] + 2.0 * b[i]);
  }
  for (int i = 0; i < inNum; i++) {
    b[i] *= 3.0;
  }

  delete [] dy; delete [] du; delete [] dd; delete [] dl;
}  // make_spline_table


void make_fast_spline_table(double *y, int inNum,
                             double *a, double *b, double *c)
{
  // 三重対角行列と差分のメモリ確保
  double *dl = new double [inNum];
  double *dd = new double [inNum];
  double *du = new double [inNum];
  double *dy = new double [inNum];


  // 三重対角行列の作成
  // 差分、かつ、傾きの計算をする
  for (int i = 0; i < inNum-1; i++) {
    dy[i] = y[i+1] - y[i];
    c[i] = dy[i];
  }

  for (int i = 0; i < inNum-3; i++) {
    dl[i] = 1.0;
    du[i] = 1.0;
  }
  dl[inNum-4] -= 1.0;

  for (int i = 0; i < inNum-2; i++) {
    dd[i] = 4.0;
    b[i]  = c[i+1] - c[i];
  }
  dd[0] += 2.0;
  du[0] -= 1.0;
  dd[inNum-3] += 2.0;

  // 三重対角行列による連立方程式の解法計算
  for (int i = 0; i < inNum-2-1; i++) {
    du[i] /= dd[i];
    dd[i+1] -= dl[i] * du[i];
  }

  b[0] /= dd[0];
  for (int i = 1; i < inNum-2; i++)
    b[i] = (b[i] - dl[i-1] * b[i-1]) / dd[i];

  for (int i = inNum-2-2; 0 <= i; i--)
    b[i] -= b[i+1] * du[i];

  // bを整える：両端の計算
  for (int i = inNum-3; 0 <= i; i--) {
    b[i+1] = b[i];
  }
  b[0] = 2.0 * b[1] - b[2];
  b[inNum-1] = 2.0 * b[inNum-2] - b[inNum-3];

  // 係数を計算する
  for (int i = 0; i < inNum-1; i++) {
    a[i] = (b[i+1] - b[i]);
    c[i] = dy[i] - (b[i+1] + 2.0 * b[i]);
  }
  for (int i = 0; i < inNum-1; i++) {
    b[i] *= 3.0;
  }

  delete [] dy; delete [] du; delete [] dd; delete [] dl;
}  // make_fast_spline_table


void make_natural_table(double *x, double *y, int inNum,
                        double *a, double *b, double *c, double *h)
{
/************************************************************
 *      make naural cubic spline coefficient table
 *
 *      Y1(t) = a*t**3 + b*t**2 + c*t + d       (d = y1)
 *              t = (x - x1)/(x2 - x1)          (0 <= t <= 1)
 *
 *      natural cubic spline
 *      Y''0(0) = 2*b[0] = Y''inNum-1(0) = 2*b[inNum-1] = 0
 *
 */
  double *dy = new double [inNum];

  // 正規化の為にサンプル間隔、差分、傾きの計算をする
  for (int i = 0; i < inNum-1; i++) {
    h[i] = x[i+1] - x[i];
    dy[i] = y[i+1] - y[i];
    c[i+1] = dy[i] / h[i];
  }

  b[0] = 0; b[inNum-1] = 0; // 境界条件：両端点での y''(x) / 6
  b[1] = c[2] - c[1] - h[0] * b[0];
  c[1] = 2.0 * (x[2] - x[0]);
  for (int i = 1; i < inNum-1-1; i++) {
    double temp = h[i] / c[i];
    b[i+1] = c[i+2] - c[i+1] - b[i] * temp;
    c[i+1] = 2.0 * (x[i+2] - x[i]) - h[i] * temp;
  }

  b[inNum-1-1] -= h[inNum-1-1] * b[inNum-1];
  for (int i = inNum-1-1; 0 < i; i--)
    b[i] = (b[i] - h[i] * b[i+1]) / c[i];
  // bを正規化する
  for (int i = 0; i < inNum; i++) {
    b[i] *= h[i] * h[i];
  }

  // 係数を計算する
  for (int i = 0; i < inNum-1; i++) {
    a[i] = (b[i+1] - b[i]);
    c[i] = dy[i] - (b[i+1] + 2.0 * b[i]);
  }
  for (int i = 0; i < inNum; i++) {
    b[i] *= 3.0;
  }

  delete[] dy;
}  // make_natural_table


void make_fast_natural_table(double *y, int inNum,
                             double *a, double *b, double *c)
{
  double *dy = new double [inNum];

  // 差分、かつ、傾きの計算をする
  for (int i = 0; i < inNum-1; i++) {
    dy[i] = y[i+1] - y[i];
    c[i+1] = dy[i];
  }

  b[0] = 0; b[inNum-1] = 0; // 境界条件：両端点での y''(x) / 6
  b[1] = c[2] - c[1] - b[0];
  c[1] = 4.0;
  for (int i = 1; i < inNum-1-1; i++) {
    b[i+1] = c[i+2] - c[i+1] - b[i] / c[i];
    c[i+1] = 4.0 - 1.0 / c[i]; // only set 15/4
  }

  b[inNum-1-1] -= b[inNum-1];
  for (int i = inNum-1-1; 0 < i; i--)
    b[i] = (b[i] - b[i+1]) / c[i];

  // 係数を計算する
  for (int i = 0; i < inNum-1; i++) {
    a[i] = (b[i+1] - b[i]);
    c[i] = dy[i] - (b[i+1] + 2.0 * b[i]);
  }
  for (int i = 0; i < inNum; i++) {
    b[i] *= 3.0;
  }

  delete[] dy;
}  // make_fast_natural_table


void make_hermite_table(double *x, double *y, int inNum,
                        double *a, double *b, double *c, double *h)
{
/************************************************************
 *      make (Catmull-Rom/cardinal spline a=-1/2)
 *              cubic spline coefficient table
 *
 *      (x0,y0), (x1,y1), (x2,y2), (x3,y3)      (x1 <= x <= x2)
 *      Y1(t) = a*t**3 + b*t**2 + c*t + d       (d = y1)
 *              t = (x - x1)/(x2 - x1)          (0 <= t <= 1)
 *      Y'0(1) = Y'1(0) = (y2 - y0)/(x2 - x0) = 1/2(vy1 + vy0)
 *      Y'1(1) = Y'2(0) = (y3 - y1)/(x3 - x1) = 1/2(vy2 + vy1)
 *
 *      端の傾きは中間値を使用
 *      Y'0(0) = (y1 - y0)/(x1 - x0)：直線補間に近づける
 *      Y'0(0) = 1/2 * (y1 - y0)/(x1 - x0)：中間値
 *              つまりは vy(0) = (y1 - y0)/(x1 - x0),  vy(-1) = 0
 *      Y'0(0) = 0：クリッピングに近づけるなら
 *
 */
  double *dy = new double [inNum];
  double *vy = new double [inNum];

  // 正規化の為にサンプル間隔と傾きの計算をする
  for (int i = 0; i < inNum-1; i++) {
    h[i] = x[i+1] - x[i];
    dy[i] = y[i+1] - y[i];
  }
    h[inNum-1] = 1.0;
    dy[inNum-1] = 0.0;

    vy[0] = 0.5 * (dy[0] / h[0]);
  for (int i = 1; i < inNum; i++) {
    vy[i] = 0.5 * (dy[i] / h[i] + dy[i-1] / h[i-1]);
  }

  // ループ
  for (int i = 0; i < inNum-1; i++) {
    calc_coeffOfCubic(dy[i], vy[i] * h[i], vy[i+1] * h[i], // 傾きを正規化している
                      &a[i], &b[i], &c[i]);
  }

  delete[] vy; delete[] dy;
}  // make_hermite_table


void make_fast_hermite_table(double *y, int inNum,
                             double *a, double *b, double *c)
{
  double *dy = new double [inNum];
  double *vy = new double [inNum];

  // 差分、かつ、傾きの計算をする
  for (int i = 0; i < inNum-1; i++) {
    dy[i] = y[i+1] - y[i];
  }
    dy[inNum-1] = 0.0;

    vy[0] = 0.5 * dy[0];
  for (int i = 1; i < inNum; i++) {
    vy[i] = 0.5 * (dy[i] + dy[i-1]);
  }

  // ループ
  for (int i = 0; i < inNum-1; i++) {
    calc_coeffOfCubic(dy[i], vy[i], vy[i+1], &a[i], &b[i], &c[i]);
  }

  delete[] vy; delete[] dy;
}  // make_fast_hermite_table


void calc_cubic_intrp_clip(double *x, double *y, int inNum,
                           double *a, double *b, double *c, double *h,
                           double *xo, int outNum, double *yo)
{
// ----------------------------------------------------
//              出力値の計算
//
  double u;
  int Index, oPnt;

  Index = 0;
  oPnt = 0;
  // 区間最小値以下: (-inf. .. 0]
  while (oPnt < outNum) {
    if (x[Index] < xo[oPnt]) break;
    yo[oPnt++] = y[Index];
  }

  // 補間区間: (0 .. inNum-1]
  while (oPnt < outNum) {
    while (x[Index] < xo[oPnt]) {
      Index++;
      if (inNum-1 < Index) goto ExitWhile;
    }
    // normalize interpolation point
    u = (xo[oPnt] - x[Index-1]) / h[Index-1];
    // 補間の計算
    yo[oPnt++] = y[Index-1] +
                 u * (c[Index-1] + u * (b[Index-1] + u * a[Index-1]));
  }
  ExitWhile:

  // 区間最大値以上: (inNum-1 .. +inf.)
  Index = inNum - 1;
  while (oPnt < outNum) {
    yo[oPnt++] = y[Index];
  }
}  // calc_cubic_intrp_clip


void calc_fast_cubic_intrp_clip(double offset, double *y, int inNum,
                                double stride, double *a, double *b, double *c,
                                int outNum, double *yo)
{
// ----------------------------------------------------
//              出力値の計算
//
  int Index, st, ed;
  double u, currentPosition;

  // 区間最小値以下: (-inf. .. 0)
  Index = 0;
  st = 0;
  currentPosition = offset;
  if (0.0 <= currentPosition) {
    ed = st;
  } else {
    ed = static_cast<int>((0.0 - currentPosition) / stride) + 1;
    if (0 <= offset + stride * ed) ed--;
    if (outNum < ed) ed = outNum;
  }
  for (int i = st; i < ed; i++) {
    yo[i] = y[Index];
  }

  // 補間区間: [0 .. inNum-1)
  st = ed;
  if (st < outNum) {
    currentPosition = offset + stride * ed;
    if (inNum-1 <= currentPosition) {
      ed = st;
    } else {
      ed = st + static_cast<int>((inNum-1 - currentPosition) / stride) + 1;
      if (inNum-1 <= offset + stride * ed) ed--;
      if (outNum < ed) ed = outNum;
    }
  }
  for (int i = st; i < ed; i++) {
    // it should be 0 <= currentPosition
    Index = static_cast<int>(currentPosition);
    // 補間の計算
    u = currentPosition - Index;
    yo[i] = y[Index] +
            u * (c[Index] + u * (b[Index] + u * a[Index]));
    // 位置更新
    currentPosition += stride;
  }

  // 区間最大値以上: [inNum-1 .. +inf.)
  Index = inNum - 1;
  st = ed; ed = outNum;
  for (int i = st; i < ed; i++) {
    yo[i] = y[Index];
  }
}  // calc_fast_cubic_intrp_clip

}  // namespace




// ハミング窓
void Hamming_window(double *w, int N)
{
  if (0 == N % 2) { /* Nが偶数のとき */
    for (int i = 0; i < N; i++)
      w[i] = 0.54 - 0.46 * cos(2.0 * M_PI * (i + 1.0) / (N+1));
  } else { /* Nが奇数のとき */
    for (int i = 0; i < N; i++)
      w[i] = 0.54 - 0.46 * cos(2.0 * M_PI * (i + 0.5) / N);
  }
} // Hamming_window

// ハニング窓
void Hanning_window(double *w, int N)
{
  if (0 == N % 2) { /* Nが偶数のとき */
    for (int i = 0; i < N; i++)
      w[i] = 0.5 - 0.5 * cos(2.0 * M_PI * (i + 1.0) / (N+1));
  } else { /* Nが奇数のとき */
    for (int i = 0; i < N; i++)
      w[i] = 0.5 - 0.5 * cos(2.0 * M_PI * (i + 0.5) / N);
  }
} // Hanning_window

// ナットール窓
void Nuttall_window(double *w, int N)
{
  if (0 == N % 2) { /* Nが偶数のとき */
    for (int i = 0; i < N; i++) {
      double tmp = (M_PI * (i+1)) / (N+1);
      w[i] = 0.355768 - 0.487396 * cos(2.0 * tmp)
                      + 0.144232 * cos(4.0 * tmp)
                      - 0.012604 * cos(6.0 * tmp);
    }
  } else { /* Nが奇数のとき */
    for (int i = 0; i < N; i++) {
      double tmp = (M_PI * (i+0.5)) / N;
      w[i] = 0.355768 - 0.487396 * cos(2.0 * tmp)
                      + 0.144232 * cos(4.0 * tmp)
                      - 0.012604 * cos(6.0 * tmp);
    }
  }
}  // Nuttall_window

// 配列渡しのsinc関数
void sinc(double *t, int tNum, double *c)
{
  for (int i = 0; i < tNum; i++) {
    if (0.0 == t[i]) {
      c[i] = 1.0;
    } else {
      double tmp = M_PI * t[i];
      c[i] = sin(tmp) / tmp;
    }
  }
} // sinc

// 昇順のコムソート qsort使うのが良いかもなんだが、使いこなせないのだよな
void combsort_ascend(double *x, int sNum)
{
  int gap = sNum;
  bool swapped = false;
  do {
    // With this factor, Combsort11 is not required.
    // you must read http://en.wikipedia.org/wiki/Combsort
    gap = static_cast<int>(gap / 1.247330950103979);
    if (gap < 1) gap = 1;
    swapped = false;

    for (int i = 0; gap + i < sNum; i++) {
      if (x[gap+i] < x[i]) { // 小さければ前に持ってくる
        double tmp = x[i];
        x[i] = x[gap+i];
        x[gap+i] = tmp;
        swapped = true;
      }
    }
  } while (swapped || (1 < gap));
} // combsort_ascend

// 降順のコムソート
void combsort_descend(double *x, int sNum)
{
  int gap = sNum;
  bool swapped = false;
  do {
    gap = static_cast<int>(gap / 1.247330950103979);
    if (gap < 1) gap = 1;
    swapped = false;

    for (int i = 0; gap + i < sNum; i++) {
      if (x[i] < x[gap + i]) { // 大きければ前に持ってくる
        double tmp = x[i];
        x[i] = x[gap + i];
        x[gap + i] = tmp;
        swapped = true;
      }
    }
  } while (swapped || (1 < gap));
} // combsort_ascend

// 中央値を返す
double median(double *x, int sNum)
{
  // ソートする
  combsort_ascend(x, sNum);

  // 中央値を計算する
  double result;
  if (0 == sNum % 2) { // 偶数
    sNum /= 2;
    result = (x[sNum] + x[sNum-1]) * 0.5; // 中央の二つの平均
  } else { // 奇数
    sNum /= 2;
    result = x[sNum]; // 中央、ぴったりど真ん中
  }
  return result;
} // median


/**********************************************************
 *
 *  解説のような愚痴
 *
 *  ・linea と 2Point Cubic Spline の違い
 *      WORKD においては、ありとあらゆる補間を linea から 2PC に
 *      変更して音の違いが判る、善い方に変わる、我田引水、自画自賛ｗ
 *
 *      UTAU のピッチ曲線の補間にはオーバーシュートやリンギングは
 *      悪影響を及ぼすと判ったので、自力開発した
 *
 *      感覚的なランキング：
 *      Catmull-Rom >= 2PC >> linea >>> trim spline >= cubic spline
 *
 *      Catmull-Rom が圧倒的と言えないのは、やはりオーバーシュートが
 *      原因なんかね？！
 *      であれば、時々見る尖がったり、角ばったピッチ曲線だと
 *      ランキングが入れ替わる可能性は結構高い
 *
 *  ・Not-a-Knot Spline と Natural Spline の違い
 *      ほとんど両端の3点間の区間をどう補間するかだけの違い
 *      所詮オーバーシュート、リンギングは起きる
 *
 *      声波形のリサンプリングには Catmull-Rom Spline の方が
 *      だいぶ善いです
 *
 *      苦労していろいろ実装した甲斐がない！ orz
 *
 *  ・Catmull-Rom Spline と trim spline, MATLAB pchip Spline の違い
 *      MathWorksの pchip のドキュメントはトンチンカンだが、
 *      こちらは素晴らしい説明だ
 *              www.math.iit.edu/~fass/Notes350_Ch3Print.pdf
 *
 *      しかし、皆オーバーシュートやリンギングでは苦労してるのだな
 *      そして、解決には涙ぐましい努力が必要なのだな。やれやれだ
 *
 *      Catmull-Rom はオーバーシュートこそあるが、割と小さいし
 *      リンギングはしない      速いし、安定していてとても善いと思う
 *
 *      前後の点から傾きを決めるアイデアも cubic spline を弄ってる間に、
 *      一応自力で思いついたんだが、世の中賢い人はたくさん居るので、
 *      とっくの昔にとても有名になってたｗ
 *
 *      とは言え、苦労して pchip や秋間先生の Akima Spline を
 *      実装しても、 WORLD では使いどころが無いような気がしたが、
 *      音声補間の為に trim spline を新たに開発した
 *      結果的に Catmull-Rom 並みの万能選手になったと思う
 *
 *      trim spline は Catmull-Rom よりもオーバーシュートは少なく
 *      直線部分はよりゆらゆら揺れるようにしてある
 *      補間グラフを目で見る限りでは、むしろこっちの方が Catmull-Rom より
 *      滑らかにも思えたりする  親の欲目かもしらんが
 *
 *      オーバーシュートが無い方が善いのであれば、 pchip が善いのだろうが
 *      オーバーシュートがゼロの2PCより Catmull-Rom の方が善い場合もあり
 *      やはり適材適所を確認しないといけない
 *
 **********************************************************
 */
/***********************************************
 *      MATLABのinterp1は基本的に補外はしないNaNを返すのだそうな
 *      なんか昔見た日本語の解説とは動きが違うような？
 *      バージョンで違うのかもですね
 *      わざわざ外挿してたけど、不要と言えば不要
 *      とは言え、クリップするのは使い勝手がいいので、これは残す
 *
 */
// interp1 uses linear interpolation(interp1 default).
// This function make no extrapolation(clip and/or nearest).
void interp1_clip(double *x, double *y, int inNum,
                  double *xo, int outNum, double *yo)
{
// ----------------------------------------------------
// 直線補間をカリカリにスピードアップする意味はあまりない
// これでも充分速いと思うのですよ
//
  const double alpha = -1; // 直線補間

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];
  double *h = new double[inNum];

  make_2PC_table(alpha, y, inNum, a, b, c);

  // 正規化の為にサンプル間隔の計算をする
  for (int i = 0; i < inNum-1; i++) {
    h[i] =  x[i+1] - x[i];
  }

  calc_cubic_intrp_clip(x, y, inNum, a, b, c, h, xo, outNum, yo);

  delete[] h; delete[] c; delete[] b; delete[] a;
}  // interp1_clip


// 最低2ポイントから使える曲線ぽい補間
// interp1 uses 2point cubic spline interpolation
//      original idea.
// This function make no extrapolation(clip and/or nearest).
void interp12PC_clip(double *x, double *y, int inNum,
                     double *xo, int outNum, double *yo)
{
  const double alpha = -1/2; // 標準値

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];
  double *h = new double[inNum];

  make_2PC_table(alpha, y, inNum, a, b, c);

  // 正規化の為にサンプル間隔の計算をする
  for (int i = 0; i < inNum-1; i++) {
    h[i] =  x[i+1] - x[i];
  }

  calc_cubic_intrp_clip(x, y, inNum, a, b, c, h, xo, outNum, yo);

  delete[] h; delete[] c; delete[] b; delete[] a;
}  // interp12PC_clip


// interp1 uses (T-Spline/trim spline) cubic spline interpolation.
// This function is not interp1 pchip/cubic.
//      original idea.
// This function make no extrapolation(clip and/or nearest).
void interp1trim_clip(double *x, double *y, int inNum,
                      double *xo, int outNum, double *yo)
{
  // スプラインの計算には最低4点必要なので、それ以下なら2PC補間を使う
  if (inNum < 4) {
    interp12PC_clip(x, y, inNum, xo, outNum, yo);
    return;
  }

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];
  double *h = new double[inNum];

  make_trim_table(x, y, inNum, a, b, c, h);

  calc_cubic_intrp_clip(x, y, inNum, a, b, c, h, xo, outNum, yo);

  delete[] h; delete[] c; delete[] b; delete[] a;
}  // interp1trim_clip


// interp1 uses (Cubic RunOut/Not-a-Knot) cubic spline interpolation.
// interp1 spline default.
//      http://www.pcs.cnu.edu/~bbradie/cinterpolation.html
// This function make no extrapolation(clip and/or nearest).
void interp1spline_clip(double *x, double *y, int inNum,
                        double *xo, int outNum, double *yo)
{
  // スプラインの計算には最低4点必要なので、それ以下なら2PC補間を使う
  if (inNum < 4) {
    interp12PC_clip(x, y, inNum, xo, outNum, yo);
    return;
  }

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];
  double *h = new double[inNum];

  make_spline_table(x, y, inNum, a, b, c, h);

  calc_cubic_intrp_clip(x, y, inNum, a, b, c, h, xo, outNum, yo);

  delete[] h; delete[] c; delete[] b; delete[] a;
}  // interp1spline_clip


// interp1 uses natural cubic spline interpolation.
// interp1 natural.
//      http://www.vector.co.jp/soft/data/prog/se002453.html
// This function make no extrapolation(clip and/or nearest).
void interp1natural_clip(double *x, double *y, int inNum,
                         double *xo, int outNum, double *yo)
{
  // スプラインの計算には最低4点必要なので、それ以下なら2PC補間を使う
  if (inNum < 4) {
    interp12PC_clip(x, y, inNum, xo, outNum, yo);
    return;
  }

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];
  double *h = new double[inNum];

  make_natural_table(x, y, inNum, a, b, c, h);

  calc_cubic_intrp_clip(x, y, inNum, a, b, c, h, xo, outNum, yo);

  delete[] h; delete[] c; delete[] b; delete[] a;
}  // interp1natural_clip


/***********************************************
 *      Catmull-Romスプライン補間
 *
 */
// interp1 uses (Catmull-Rom/cardinal spline a=-1/2) cubic spline interpolation.
// This function is not interp1 pchip/cubic.
//      http://yehar.com/blog/wp-content/uploads/2009/08/deip.pdf
// This function make no extrapolation(clip and/or nearest).
void interp1Catmull_Rom_clip(double *x, double *y, int inNum,
                             double *xo, int outNum, double *yo)
{
  // スプラインの計算には最低4点必要なので、それ以下なら2PC補間を使う
  if (inNum < 4) {
    interp12PC_clip(x, y, inNum, xo, outNum, yo);
    return;
  }

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];
  double *h = new double[inNum];

  make_hermite_table(x, y, inNum, a, b, c, h);

  calc_cubic_intrp_clip(x, y, inNum, a, b, c, h, xo, outNum, yo);

  delete[] h; delete[] c; delete[] b; delete[] a;
}  // interp1Catmull_Rom_clip


// サンプリング間隔を等間隔に限定し、高速に動作する直線補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
void itrp1Q_clip(double offset, double *y, int inNum, double stride,
                 int outNum, double *yo)
{
// ----------------------------------------------------
// 直線補間をカリカリにスピードアップする意味はあまりないと思うます
// これでも充分速いと思うのですよ
//
  const double alpha = -1; // 直線補間

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];

  make_2PC_table(alpha, y, inNum, a, b, c);

  calc_fast_cubic_intrp_clip(offset, y, inNum, stride, a, b, c, outNum, yo);

  delete[] c; delete[] b; delete[] a;
}  // itrp1Q_clip


// サンプリング間隔を等間隔に限定し、高速に動作する2PointCubic補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
void itrp1Q2PC_clip(double offset, double *y, int inNum, double stride,
                    int outNum, double *yo)
{
  const double alpha = -1/2; // 標準値

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];

  make_2PC_table(alpha, y, inNum, a, b, c);

  calc_fast_cubic_intrp_clip(offset, y, inNum, stride, a, b, c, outNum, yo);

  delete[] c; delete[] b; delete[] a;
}  // itrp1Q2PC_clip


// サンプリング間隔を等間隔に限定し、高速に動作するトリムスプライン補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
//      両端の区間はクリップに近づくように設定している
void itrp1Qtrim_clip(double offset, double *y, int inNum, double stride,
                            int outNum, double *yo)
{
  // スプラインの計算には最低4点必要なので、それ以下なら2PC補間を使う
  if (inNum < 4) {
    itrp1Q2PC_clip(offset, y, inNum, stride, outNum, yo);
    return;
  }


  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];

  make_fast_trim_table(y, inNum, a, b, c);

  calc_fast_cubic_intrp_clip(offset, y, inNum, stride, a, b, c, outNum, yo);

  delete[] c; delete[] b; delete[] a;
}  // itrp1Qtrim_clip


// サンプリング間隔を等間隔に限定し、高速に動作するスプライン補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
void itrp1Qspline_clip(double offset, double *y, int inNum, double stride,
                       int outNum, double *yo)
{
  // スプラインの計算には最低4点必要なので、それ以下なら2PC補間を使う
  if (inNum < 4) {
    itrp1Q2PC_clip(offset, y, inNum, stride, outNum, yo);
    return;
  }

// ----------------------------------------------------
//              make cubic spline coefficient table
//
  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];

  make_fast_spline_table(y, inNum, a, b, c);

  calc_fast_cubic_intrp_clip(offset, y, inNum, stride, a, b, c, outNum, yo);

  delete[] c; delete[] b; delete[] a;
}  // itrp1Qspline_clip


// サンプリング間隔を等間隔に限定し、高速に動作する自然スプライン補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
void itrp1Qnatural_clip(double offset, double *y, int inNum, double stride,
                        int outNum, double *yo)
{
  // スプラインの計算には最低4点必要なので、それ以下なら2PC補間を使う
  if (inNum < 4) {
    itrp1Q2PC_clip(offset, y, inNum, stride, outNum, yo);
    return;
  }

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];

  make_fast_natural_table(y, inNum, a, b, c);

  calc_fast_cubic_intrp_clip(offset, y, inNum, stride, a, b, c, outNum, yo);

  delete[] c; delete[] b; delete[] a;
}  // itrp1Qnatural_clip


// サンプリング間隔を等間隔に限定し、高速に動作するCatmull-Romスプライン補間、区間外はクリップ
// 元データのサンプリング間隔を 1 に限定している
//      両端の区間はクリップに近づくように設定している
void itrp1QCatmull_Rom_clip(double offset, double *y, int inNum, double stride,
                            int outNum, double *yo)
{
  // スプラインの計算には最低4点必要なので、それ以下なら2PC補間を使う
  if (inNum < 4) {
    itrp1Q2PC_clip(offset, y, inNum, stride, outNum, yo);
    return;
  }

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];

  make_fast_hermite_table(y, inNum, a, b, c);

  calc_fast_cubic_intrp_clip(offset, y, inNum, stride, a, b, c, outNum, yo);

  delete[] c; delete[] b; delete[] a;
}  // itrp1QCatmull_Rom_clip


