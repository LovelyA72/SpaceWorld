// Matlab�����̊֐������@FFT�s�v�Ȃ��̂̂�
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


// �P�ʋ�Ԃɂ�����O����Ԋ֐��̌W�����v�Z����
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
 *      2�_�O���Ȑ����
 *      2PontCubicSpline�̕�ԗp ��Ԋ֐�
 *      �v�͎O���G���~�[�g�X�v���C��(cubic hermit spline)�̈��
 *      use fixed start/end tangent cubic hermit spline
 *
 *      f(x) =  2(1+alph)*x**3 - 3(1+alpha)*x**2 + alpha*x + 1  (0 <= x <= 1)
 *              0                                               (1 <= x)
 *              f(-x)                                           (x <= 0)
 *
 *      alpha��x=0, 1�ł̌X��       (-1 <= alpha <= 0)
 *               alpha= 0�̎��̓R�T�C����Ԃ��ߎ�
 *                      �T���v���_��ʂ�񎟂�B-Spline�ɂ��߂�
 *               alpha= -1/2���W���I�A������ԑ���ɃO�b�h�Ǝv��
 *               alpha= -1�̎��͒�����Ԃɓ�����
 *
 *               (-inf. < a < -1) �ł��������A
 *                      �Ȑ��̐U���������Ȃ邵�A
 *                      �T���v���_�ł͉s���낪�邵�ŕ�Ԃɂ͎g���Ȃ�
 *                      CG�Ŗ͗l��`�����Ɏg���邩�����炢��
 *
 *      �G���~�[�g����̃G�̎����m��Ȃ��ōl���t�����������炢�I
 *      �Ǝ��掩�^���Ă�����
 *
 *      Y1(t) = -2(1+alpha)*dy*t**3 + 3(1+alpha)*dy*t**2 - alpha*dy*t + d
 *              d = y1
 *              t = (x - x1)/(x2 - x1)  (0 <= t <= 1)
 *              dy = y2 - y1
 *
 */
  for (int i = 0; i < inNum-1; i++) {
    // �����̌v�Z������
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
 *      Catmull-Rom �̌X���A�ڐ��͕K�����_�����p�ɐڂ���
 *      �G�b�W�̏������ׂ̈ɏ����ɂ₩�ɂȂ�̂ŃI�[�o�[�V���[�g������
 *      �ł��㉺�̃s�[�N�������ɂ₩�ɂ���̂ŏ������܂����Ȃ�
 *
 *      ���̃X�v���C���̓G�b�W��㉺�̃s�[�N�͂��ɂ₩�ɁA
 *      �����ɋ߂������͂��}�s�ȃX�v���C��
 *      ������ɂ��� MATLAB �� pchip �Ƃ͂����ԈႤ
 *      �G�b�W�̏����قڌX���[���ɂȂ�̂ŃI�[�o�[�V���[�g���قƂ�ǖ����Ȃ�
 *      �����ɋ߂����͌X���A�ڐ������_�����p���ђʂ���
 *      ���� Catmull-Rom ���ӂ�ӂ�h��Ă�悤�ɂ�������
 *      �����ڂ� Catmull-Rom ��芊�炩�����������������񂾂��ǂ�
 *      �I���W�i���Ȃ񂾂��� maid spline �ł��P�������񂾂��ǂ�
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
 *      �������̏ꍇ�F���a���ςƑ������ς̔���|����            (0 <= q <= 1)
 *      �ٕ����̏ꍇ�F�����X���Ɗɂ��X���̍� �� �ɂ��X���Ƃ�
 *                    ���a���ςƑ������ς̔��1/2���Ċ|����     (0 <= q <= 1/2)
 *      ���A�����A�ٕ����̏ꍇ��(0 <= q <= 1)�Ƃ��Ă����A�F�X Catmull-Rom ��
 *        �����Ă��̂ŁA�����̌��ʂ��̂悤�Ɍ��肵��
 *                      trim ratio = 1 �Ƃ��Ă�����
 *
 *      �[�̌X���̓[�����P���̂��A���Ԓl���P���̂��H
 *      ��芸���� Catmull-Rom �ɍ��킹�鎖�ɂ���
 *      Y'0(0) = 1/2 * vy0�F���Ԓl
 *              �܂�� vy1 = 1/2 * vy0  vy(-1) = 7*vy0
 *
 */
  double *dy = new double [inNum];
  double *vy = new double [inNum];
  double *vyAdd = new double [inNum];
  double *vyMul = new double [inNum];

  // ���K���ׂ̈ɃT���v���Ԋu�ƌX���̌v�Z������
  for (int i = 0; i < inNum-1; i++) {
    h[i] = x[i+1] - x[i];
    dy[i] = y[i+1] - y[i];
    vy[i] = dy[i] / h[i];
  }
    // �[�_�̐ݒ� Y'(inNum-1)(0) = 1/2 * vy[inNum-2]
    h[inNum-1] = h[inNum-1-1];
    dy[inNum-1] = dy[inNum-1-1] / 7.0;
    vy[inNum-1] = dy[inNum-1] / h[inNum-1];

    // �[�_�̐ݒ� Y'0(0) = 1/2 * vy[0]
    vyAdd[0] = vy[0] + (vy[0] / 7.0);
    vyMul[0] = vy[0] * (vy[0] / 7.0);
  for (int i = 1; i < inNum; i++) {
    vyAdd[i] = vy[i] + vy[i-1];
    vyMul[i] = vy[i] * vy[i-1];
  }

  double vy1, vy2;
  // �ŏ�����
      // �P���ɂ���ł����񂾂��A������ƈӎ����邽�߂Ɍv�Z����
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
    calc_coeffOfCubic(dy[0], vy1 * h[0], vy2 * h[0], // �X���𐳋K�����Ă���
                       &a[0], &b[0], &c[0]);

  // ���[�v
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
    calc_coeffOfCubic(dy[i], vy1 * h[i], vy2 * h[i], // �X���𐳋K�����Ă���
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

  // �����A���A�X���̌v�Z������
  for (int i = 0; i < inNum-1; i++) {
    dy[i] = y[i+1] - y[i];
  }
    // �[�_�̐ݒ� Y'(inNum-1)(0) = 1/2 * vy[inNum-2]
    dy[inNum-1] = dy[inNum-1-1] / 7.0;

    // �[�_�̐ݒ� Y'0(0) = 1/2 * vy[0]
    vyAdd[0] = dy[0] + (dy[0] / 7.0);
    vyMul[0] = dy[0] * (dy[0] / 7.0);
  for (int i = 1; i < inNum; i++) {
    vyAdd[i] = dy[i] + dy[i-1];
    vyMul[i] = dy[i] * dy[i-1];
  }

  double vy1, vy2;
  // �ŏ�����
      // �P���ɂ���ł����񂾂��A������ƈӎ����邽�߂Ɍv�Z����
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

  // ���[�v
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
  // �O�d�Ίp�s��ƍ����̃������m��
  double *dl = new double [inNum];
  double *dd = new double [inNum];
  double *du = new double [inNum];
  double *dy = new double [inNum];


  // �O�d�Ίp�s��̍쐬
  // ���K���ׂ̈ɃT���v���Ԋu�A�����A�X���̌v�Z������
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

  // �O�d�Ίp�s��ɂ��A���������̉�@�v�Z����
  for (int i = 0; i < inNum-2-1; i++) {
    du[i] /= dd[i];
    dd[i+1] -= dl[i] * du[i];
  }

  b[0] /= dd[0];
  for (int i = 1; i < inNum-2; i++)
    b[i] = (b[i] - dl[i-1] * b[i-1]) / dd[i];

  for (int i = inNum-2-2; 0 <= i; i--)
    b[i] -= b[i+1] * du[i];

  // b�𐮂���F���[�̌v�Z�Ɛ��K��
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

  // �W�����v�Z����
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
  // �O�d�Ίp�s��ƍ����̃������m��
  double *dl = new double [inNum];
  double *dd = new double [inNum];
  double *du = new double [inNum];
  double *dy = new double [inNum];


  // �O�d�Ίp�s��̍쐬
  // �����A���A�X���̌v�Z������
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

  // �O�d�Ίp�s��ɂ��A���������̉�@�v�Z
  for (int i = 0; i < inNum-2-1; i++) {
    du[i] /= dd[i];
    dd[i+1] -= dl[i] * du[i];
  }

  b[0] /= dd[0];
  for (int i = 1; i < inNum-2; i++)
    b[i] = (b[i] - dl[i-1] * b[i-1]) / dd[i];

  for (int i = inNum-2-2; 0 <= i; i--)
    b[i] -= b[i+1] * du[i];

  // b�𐮂���F���[�̌v�Z
  for (int i = inNum-3; 0 <= i; i--) {
    b[i+1] = b[i];
  }
  b[0] = 2.0 * b[1] - b[2];
  b[inNum-1] = 2.0 * b[inNum-2] - b[inNum-3];

  // �W�����v�Z����
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

  // ���K���ׂ̈ɃT���v���Ԋu�A�����A�X���̌v�Z������
  for (int i = 0; i < inNum-1; i++) {
    h[i] = x[i+1] - x[i];
    dy[i] = y[i+1] - y[i];
    c[i+1] = dy[i] / h[i];
  }

  b[0] = 0; b[inNum-1] = 0; // ���E�����F���[�_�ł� y''(x) / 6
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
  // b�𐳋K������
  for (int i = 0; i < inNum; i++) {
    b[i] *= h[i] * h[i];
  }

  // �W�����v�Z����
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

  // �����A���A�X���̌v�Z������
  for (int i = 0; i < inNum-1; i++) {
    dy[i] = y[i+1] - y[i];
    c[i+1] = dy[i];
  }

  b[0] = 0; b[inNum-1] = 0; // ���E�����F���[�_�ł� y''(x) / 6
  b[1] = c[2] - c[1] - b[0];
  c[1] = 4.0;
  for (int i = 1; i < inNum-1-1; i++) {
    b[i+1] = c[i+2] - c[i+1] - b[i] / c[i];
    c[i+1] = 4.0 - 1.0 / c[i]; // only set 15/4
  }

  b[inNum-1-1] -= b[inNum-1];
  for (int i = inNum-1-1; 0 < i; i--)
    b[i] = (b[i] - b[i+1]) / c[i];

  // �W�����v�Z����
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
 *      �[�̌X���͒��Ԓl���g�p
 *      Y'0(0) = (y1 - y0)/(x1 - x0)�F������Ԃɋ߂Â���
 *      Y'0(0) = 1/2 * (y1 - y0)/(x1 - x0)�F���Ԓl
 *              �܂�� vy(0) = (y1 - y0)/(x1 - x0),  vy(-1) = 0
 *      Y'0(0) = 0�F�N���b�s���O�ɋ߂Â���Ȃ�
 *
 */
  double *dy = new double [inNum];
  double *vy = new double [inNum];

  // ���K���ׂ̈ɃT���v���Ԋu�ƌX���̌v�Z������
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

  // ���[�v
  for (int i = 0; i < inNum-1; i++) {
    calc_coeffOfCubic(dy[i], vy[i] * h[i], vy[i+1] * h[i], // �X���𐳋K�����Ă���
                      &a[i], &b[i], &c[i]);
  }

  delete[] vy; delete[] dy;
}  // make_hermite_table


void make_fast_hermite_table(double *y, int inNum,
                             double *a, double *b, double *c)
{
  double *dy = new double [inNum];
  double *vy = new double [inNum];

  // �����A���A�X���̌v�Z������
  for (int i = 0; i < inNum-1; i++) {
    dy[i] = y[i+1] - y[i];
  }
    dy[inNum-1] = 0.0;

    vy[0] = 0.5 * dy[0];
  for (int i = 1; i < inNum; i++) {
    vy[i] = 0.5 * (dy[i] + dy[i-1]);
  }

  // ���[�v
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
//              �o�͒l�̌v�Z
//
  double u;
  int Index, oPnt;

  Index = 0;
  oPnt = 0;
  // ��ԍŏ��l�ȉ�: (-inf. .. 0]
  while (oPnt < outNum) {
    if (x[Index] < xo[oPnt]) break;
    yo[oPnt++] = y[Index];
  }

  // ��ԋ��: (0 .. inNum-1]
  while (oPnt < outNum) {
    while (x[Index] < xo[oPnt]) {
      Index++;
      if (inNum-1 < Index) goto ExitWhile;
    }
    // normalize interpolation point
    u = (xo[oPnt] - x[Index-1]) / h[Index-1];
    // ��Ԃ̌v�Z
    yo[oPnt++] = y[Index-1] +
                 u * (c[Index-1] + u * (b[Index-1] + u * a[Index-1]));
  }
  ExitWhile:

  // ��ԍő�l�ȏ�: (inNum-1 .. +inf.)
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
//              �o�͒l�̌v�Z
//
  int Index, st, ed;
  double u, currentPosition;

  // ��ԍŏ��l�ȉ�: (-inf. .. 0)
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

  // ��ԋ��: [0 .. inNum-1)
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
    // ��Ԃ̌v�Z
    u = currentPosition - Index;
    yo[i] = y[Index] +
            u * (c[Index] + u * (b[Index] + u * a[Index]));
    // �ʒu�X�V
    currentPosition += stride;
  }

  // ��ԍő�l�ȏ�: [inNum-1 .. +inf.)
  Index = inNum - 1;
  st = ed; ed = outNum;
  for (int i = st; i < ed; i++) {
    yo[i] = y[Index];
  }
}  // calc_fast_cubic_intrp_clip

}  // namespace




// �n�~���O��
void Hamming_window(double *w, int N)
{
  if (0 == N % 2) { /* N�������̂Ƃ� */
    for (int i = 0; i < N; i++)
      w[i] = 0.54 - 0.46 * cos(2.0 * M_PI * (i + 1.0) / (N+1));
  } else { /* N����̂Ƃ� */
    for (int i = 0; i < N; i++)
      w[i] = 0.54 - 0.46 * cos(2.0 * M_PI * (i + 0.5) / N);
  }
} // Hamming_window

// �n�j���O��
void Hanning_window(double *w, int N)
{
  if (0 == N % 2) { /* N�������̂Ƃ� */
    for (int i = 0; i < N; i++)
      w[i] = 0.5 - 0.5 * cos(2.0 * M_PI * (i + 1.0) / (N+1));
  } else { /* N����̂Ƃ� */
    for (int i = 0; i < N; i++)
      w[i] = 0.5 - 0.5 * cos(2.0 * M_PI * (i + 0.5) / N);
  }
} // Hanning_window

// �i�b�g�[����
void Nuttall_window(double *w, int N)
{
  if (0 == N % 2) { /* N�������̂Ƃ� */
    for (int i = 0; i < N; i++) {
      double tmp = (M_PI * (i+1)) / (N+1);
      w[i] = 0.355768 - 0.487396 * cos(2.0 * tmp)
                      + 0.144232 * cos(4.0 * tmp)
                      - 0.012604 * cos(6.0 * tmp);
    }
  } else { /* N����̂Ƃ� */
    for (int i = 0; i < N; i++) {
      double tmp = (M_PI * (i+0.5)) / N;
      w[i] = 0.355768 - 0.487396 * cos(2.0 * tmp)
                      + 0.144232 * cos(4.0 * tmp)
                      - 0.012604 * cos(6.0 * tmp);
    }
  }
}  // Nuttall_window

// �z��n����sinc�֐�
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

// �����̃R���\�[�g qsort�g���̂��ǂ������Ȃ񂾂��A�g�����Ȃ��Ȃ��̂����
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
      if (x[gap+i] < x[i]) { // ��������ΑO�Ɏ����Ă���
        double tmp = x[i];
        x[i] = x[gap+i];
        x[gap+i] = tmp;
        swapped = true;
      }
    }
  } while (swapped || (1 < gap));
} // combsort_ascend

// �~���̃R���\�[�g
void combsort_descend(double *x, int sNum)
{
  int gap = sNum;
  bool swapped = false;
  do {
    gap = static_cast<int>(gap / 1.247330950103979);
    if (gap < 1) gap = 1;
    swapped = false;

    for (int i = 0; gap + i < sNum; i++) {
      if (x[i] < x[gap + i]) { // �傫����ΑO�Ɏ����Ă���
        double tmp = x[i];
        x[i] = x[gap + i];
        x[gap + i] = tmp;
        swapped = true;
      }
    }
  } while (swapped || (1 < gap));
} // combsort_ascend

// �����l��Ԃ�
double median(double *x, int sNum)
{
  // �\�[�g����
  combsort_ascend(x, sNum);

  // �����l���v�Z����
  double result;
  if (0 == sNum % 2) { // ����
    sNum /= 2;
    result = (x[sNum] + x[sNum-1]) * 0.5; // �����̓�̕���
  } else { // �
    sNum /= 2;
    result = x[sNum]; // �����A�҂�����ǐ^��
  }
  return result;
} // median


/**********************************************************
 *
 *  ����̂悤�ȋ�s
 *
 *  �Elinea �� 2Point Cubic Spline �̈Ⴂ
 *      WORKD �ɂ����ẮA����Ƃ������Ԃ� linea ���� 2PC ��
 *      �ύX���ĉ��̈Ⴂ������A�P�����ɕς��A��c�����A���掩�^��
 *
 *      UTAU �̃s�b�`�Ȑ��̕�Ԃɂ̓I�[�o�[�V���[�g�⃊���M���O��
 *      ���e�����y�ڂ��Ɣ������̂ŁA���͊J������
 *
 *      ���o�I�ȃ����L���O�F
 *      Catmull-Rom >= 2PC >> linea >>> trim spline >= cubic spline
 *
 *      Catmull-Rom �����|�I�ƌ����Ȃ��̂́A��͂�I�[�o�[�V���[�g��
 *      �����Ȃ񂩂ˁH�I
 *      �ł���΁A���X����낪������A�p�΂����s�b�`�Ȑ�����
 *      �����L���O������ւ��\���͌��\����
 *
 *  �ENot-a-Knot Spline �� Natural Spline �̈Ⴂ
 *      �قƂ�Ǘ��[��3�_�Ԃ̋�Ԃ��ǂ���Ԃ��邩�����̈Ⴂ
 *      ���F�I�[�o�[�V���[�g�A�����M���O�͋N����
 *
 *      ���g�`�̃��T���v�����O�ɂ� Catmull-Rom Spline �̕���
 *      �����ԑP���ł�
 *
 *      ��J���Ă��낢����������b�オ�Ȃ��I orz
 *
 *  �ECatmull-Rom Spline �� trim spline, MATLAB pchip Spline �̈Ⴂ
 *      MathWorks�� pchip �̃h�L�������g�̓g���`���J�������A
 *      ������͑f���炵��������
 *              www.math.iit.edu/~fass/Notes350_Ch3Print.pdf
 *
 *      �������A�F�I�[�o�[�V���[�g�⃊���M���O�ł͋�J���Ă�̂���
 *      �����āA�����ɂ͗܂��܂����w�͂��K�v�Ȃ̂��ȁB����ꂾ
 *
 *      Catmull-Rom �̓I�[�o�[�V���[�g�������邪�A���Ə�������
 *      �����M���O�͂��Ȃ�      �������A���肵�Ă��ĂƂĂ��P���Ǝv��
 *
 *      �O��̓_����X�������߂�A�C�f�A�� cubic spline ��M���Ă�ԂɁA
 *      �ꉞ���͂Ŏv�������񂾂��A���̒������l�͂������񋏂�̂ŁA
 *      �Ƃ����̐̂ɂƂĂ��L���ɂȂ��Ă���
 *
 *      �Ƃ͌����A��J���� pchip ��H�Ԑ搶�� Akima Spline ��
 *      �������Ă��A WORLD �ł͎g���ǂ��낪�����悤�ȋC���������A
 *      ������Ԃׂ̈� trim spline ��V���ɊJ������
 *      ���ʓI�� Catmull-Rom ���݂̖��\�I��ɂȂ����Ǝv��
 *
 *      trim spline �� Catmull-Rom �����I�[�o�[�V���[�g�͏��Ȃ�
 *      ���������͂������h���悤�ɂ��Ă���
 *      ��ԃO���t��ڂŌ������ł́A�ނ��낱�����̕��� Catmull-Rom ���
 *      ���炩�ɂ��v�����肷��  �e�̗~�ڂ��������
 *
 *      �I�[�o�[�V���[�g�����������P���̂ł���΁A pchip ���P���̂��낤��
 *      �I�[�o�[�V���[�g���[����2PC��� Catmull-Rom �̕����P���ꍇ������
 *      ��͂�K�ޓK�����m�F���Ȃ��Ƃ����Ȃ�
 *
 **********************************************************
 */
/***********************************************
 *      MATLAB��interp1�͊�{�I�ɕ�O�͂��Ȃ�NaN��Ԃ��̂�������
 *      �Ȃ񂩐̌������{��̉���Ƃ͓������Ⴄ�悤�ȁH
 *      �o�[�W�����ňႤ�̂����ł���
 *      �킴�킴�O�}���Ă����ǁA�s�v�ƌ����Εs�v
 *      �Ƃ͌����A�N���b�v����͎̂g�����肪�����̂ŁA����͎c��
 *
 */
// interp1 uses linear interpolation(interp1 default).
// This function make no extrapolation(clip and/or nearest).
void interp1_clip(double *x, double *y, int inNum,
                  double *xo, int outNum, double *yo)
{
// ----------------------------------------------------
// ������Ԃ��J���J���ɃX�s�[�h�A�b�v����Ӗ��͂��܂�Ȃ�
// ����ł��[�������Ǝv���̂ł���
//
  const double alpha = -1; // �������

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];
  double *h = new double[inNum];

  make_2PC_table(alpha, y, inNum, a, b, c);

  // ���K���ׂ̈ɃT���v���Ԋu�̌v�Z������
  for (int i = 0; i < inNum-1; i++) {
    h[i] =  x[i+1] - x[i];
  }

  calc_cubic_intrp_clip(x, y, inNum, a, b, c, h, xo, outNum, yo);

  delete[] h; delete[] c; delete[] b; delete[] a;
}  // interp1_clip


// �Œ�2�|�C���g����g����Ȑ��ۂ����
// interp1 uses 2point cubic spline interpolation
//      original idea.
// This function make no extrapolation(clip and/or nearest).
void interp12PC_clip(double *x, double *y, int inNum,
                     double *xo, int outNum, double *yo)
{
  const double alpha = -1/2; // �W���l

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];
  double *h = new double[inNum];

  make_2PC_table(alpha, y, inNum, a, b, c);

  // ���K���ׂ̈ɃT���v���Ԋu�̌v�Z������
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
  // �X�v���C���̌v�Z�ɂ͍Œ�4�_�K�v�Ȃ̂ŁA����ȉ��Ȃ�2PC��Ԃ��g��
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
  // �X�v���C���̌v�Z�ɂ͍Œ�4�_�K�v�Ȃ̂ŁA����ȉ��Ȃ�2PC��Ԃ��g��
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
  // �X�v���C���̌v�Z�ɂ͍Œ�4�_�K�v�Ȃ̂ŁA����ȉ��Ȃ�2PC��Ԃ��g��
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
 *      Catmull-Rom�X�v���C�����
 *
 */
// interp1 uses (Catmull-Rom/cardinal spline a=-1/2) cubic spline interpolation.
// This function is not interp1 pchip/cubic.
//      http://yehar.com/blog/wp-content/uploads/2009/08/deip.pdf
// This function make no extrapolation(clip and/or nearest).
void interp1Catmull_Rom_clip(double *x, double *y, int inNum,
                             double *xo, int outNum, double *yo)
{
  // �X�v���C���̌v�Z�ɂ͍Œ�4�_�K�v�Ȃ̂ŁA����ȉ��Ȃ�2PC��Ԃ��g��
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


// �T���v�����O�Ԋu�𓙊Ԋu�Ɍ��肵�A�����ɓ��삷�钼����ԁA��ԊO�̓N���b�v
// ���f�[�^�̃T���v�����O�Ԋu�� 1 �Ɍ��肵�Ă���
void itrp1Q_clip(double offset, double *y, int inNum, double stride,
                 int outNum, double *yo)
{
// ----------------------------------------------------
// ������Ԃ��J���J���ɃX�s�[�h�A�b�v����Ӗ��͂��܂�Ȃ��Ǝv���܂�
// ����ł��[�������Ǝv���̂ł���
//
  const double alpha = -1; // �������

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];

  make_2PC_table(alpha, y, inNum, a, b, c);

  calc_fast_cubic_intrp_clip(offset, y, inNum, stride, a, b, c, outNum, yo);

  delete[] c; delete[] b; delete[] a;
}  // itrp1Q_clip


// �T���v�����O�Ԋu�𓙊Ԋu�Ɍ��肵�A�����ɓ��삷��2PointCubic��ԁA��ԊO�̓N���b�v
// ���f�[�^�̃T���v�����O�Ԋu�� 1 �Ɍ��肵�Ă���
void itrp1Q2PC_clip(double offset, double *y, int inNum, double stride,
                    int outNum, double *yo)
{
  const double alpha = -1/2; // �W���l

  double *a = new double[inNum];
  double *b = new double[inNum];
  double *c = new double[inNum];

  make_2PC_table(alpha, y, inNum, a, b, c);

  calc_fast_cubic_intrp_clip(offset, y, inNum, stride, a, b, c, outNum, yo);

  delete[] c; delete[] b; delete[] a;
}  // itrp1Q2PC_clip


// �T���v�����O�Ԋu�𓙊Ԋu�Ɍ��肵�A�����ɓ��삷��g�����X�v���C����ԁA��ԊO�̓N���b�v
// ���f�[�^�̃T���v�����O�Ԋu�� 1 �Ɍ��肵�Ă���
//      ���[�̋�Ԃ̓N���b�v�ɋ߂Â��悤�ɐݒ肵�Ă���
void itrp1Qtrim_clip(double offset, double *y, int inNum, double stride,
                            int outNum, double *yo)
{
  // �X�v���C���̌v�Z�ɂ͍Œ�4�_�K�v�Ȃ̂ŁA����ȉ��Ȃ�2PC��Ԃ��g��
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


// �T���v�����O�Ԋu�𓙊Ԋu�Ɍ��肵�A�����ɓ��삷��X�v���C����ԁA��ԊO�̓N���b�v
// ���f�[�^�̃T���v�����O�Ԋu�� 1 �Ɍ��肵�Ă���
void itrp1Qspline_clip(double offset, double *y, int inNum, double stride,
                       int outNum, double *yo)
{
  // �X�v���C���̌v�Z�ɂ͍Œ�4�_�K�v�Ȃ̂ŁA����ȉ��Ȃ�2PC��Ԃ��g��
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


// �T���v�����O�Ԋu�𓙊Ԋu�Ɍ��肵�A�����ɓ��삷�鎩�R�X�v���C����ԁA��ԊO�̓N���b�v
// ���f�[�^�̃T���v�����O�Ԋu�� 1 �Ɍ��肵�Ă���
void itrp1Qnatural_clip(double offset, double *y, int inNum, double stride,
                        int outNum, double *yo)
{
  // �X�v���C���̌v�Z�ɂ͍Œ�4�_�K�v�Ȃ̂ŁA����ȉ��Ȃ�2PC��Ԃ��g��
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


// �T���v�����O�Ԋu�𓙊Ԋu�Ɍ��肵�A�����ɓ��삷��Catmull-Rom�X�v���C����ԁA��ԊO�̓N���b�v
// ���f�[�^�̃T���v�����O�Ԋu�� 1 �Ɍ��肵�Ă���
//      ���[�̋�Ԃ̓N���b�v�ɋ߂Â��悤�ɐݒ肵�Ă���
void itrp1QCatmull_Rom_clip(double offset, double *y, int inNum, double stride,
                            int outNum, double *yo)
{
  // �X�v���C���̌v�Z�ɂ͍Œ�4�_�K�v�Ȃ̂ŁA����ȉ��Ȃ�2PC��Ԃ��g��
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


