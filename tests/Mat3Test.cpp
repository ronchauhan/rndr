//------------------------------- Mat3Test.cpp -------------------------------//
//
// This file is a part of the rndr project, under the Apache License v2.0
// See rndr/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0
//
//----------------------------------------------------------------------------//

#include "rndr/Mat3.h"
#include "rndr/Vec3.h"
#include <cassert>
#include <iostream>

int main() {
  using Mat3f = rndr::Mat3f;
  using Vec3f = rndr::Vec3f;

  Vec3f X(1); // (1, 1, 1)
  Mat3f I;
  assert(I * X == X);

  Mat3f M1({1, 3, 2}, {-1, 0, 2}, {3, 1, -1});
  Mat3f M2;
  Mat3f M3;
  Mat3f M4;

  // Operations over self.
  assert(M1.getDeterminant() == 11);
  assert(M1.getTranspose().getTranspose() == M1);
  assert(M1.getInverse().getInverse() == M1);
  assert(M1.getInverse() * M1 == M1 * M1.getInverse());
  assert(M1.getInverse() * M1 == I);

  assert(+M1 == M1);

  M2 = -M1;
  assert(M1 + M2 == Mat3f(0));

  // Matrix-vector multiplication (column-major form).
  Vec3f X_ = M1 * X;
  assert(M1.getInverse() * X_ == X);

  // Operations with constants.
  M3 = I * 3;
  assert(M3 / 3 == I);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      assert(M3[i][j] == I[i][j] * 3);
    }
  }

  M3 = I;
  M3 *= 3;
  M3 /= 3;
  assert(M3 == I);

  // Operations with other matrices.
  M2 = Mat3f(1.5, 2, 3, 4.5, 5, 6, 7.5, 8, 9);

  M3 = M1 + M2;
  M4 = M1 - M2;
  assert(M3 != M4);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      assert(M3[i][j] == M1[i][j] + M2[i][j]);
      assert(M4[i][j] == M1[i][j] - M2[i][j]);
    }
  }

  M1 += M2;
  assert(M1 == M3);
  M1 -= 2 * M2;
  assert(M1 == M4);
}
