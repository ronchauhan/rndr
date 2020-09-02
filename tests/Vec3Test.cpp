//------------------------------- Vec3Test.cpp -------------------------------//
//
// This file is a part of the rndr project, under the Apache License v2.0
// See rndr/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0
//
//----------------------------------------------------------------------------//

#include "rndr/Vec3.h"
#include "rndr/Mat3.h"
#include <cassert>
#include <iostream>

int main() {
  using Vec3f = rndr::Vec3f;
  using Mat3f = rndr::Mat3f;

  Vec3f V1(2.0f, 4.0f, 5.0f);
  Vec3f V2;
  Vec3f V3;

  // Operations over self.
  V3 = +V1;
  assert(V3 == V1);

  V3 = -V1;
  assert(V3 + V1 == Vec3f(0));

  // Operations with scalar values.
  V2 = V1 * 3;
  assert(V2 == Vec3f(3 * 2.0, 3 * 4.0, 3 * 5.0));

  V2 = 3 * V1;
  assert(V2 == V1 * 3);

  V2 = V1;
  V2 *= 3;
  assert(V2 == V1 * 3);

  V2 = V1 / 3;
  assert(V2 == Vec3f(2.0 / 3, 4.0 / 3, 5.0 / 3));

  V2 = V1;
  V2 /= 3;
  assert(V2 == V1 / 3);

  // Operations with other vectors.
  V2 = Vec3f(1.3f, 2.5f, 4.6f);
  V3 = V1 + V2;
  assert(V3 == Vec3f(2.0 + 1.3, 4 + 2.5, 5 + 4.6));

  V3 = V1;
  V3 += V2;
  assert(V3 == V1 + V2);

  V3 = V1 - V2;
  assert(V3 == Vec3f(2.0 - 1.3, 4 - 2.5, 5 - 4.6));

  V3 = V1;
  V3 -= V2;
  assert(V3 == V1 - V2);

  V3 = V1.getNormalizedVector();
  assert(V3.getLength() == 1);

  // Components with small values
  V1 = Vec3f(0.5f, 0.3f, 0.09f);
  V1.normalizeThis();
  assert(1 - V1.getLength() < 0.0000001);

  // Matrix-vector multiplication (row-major form).
  V1 = Vec3f(1, 2, 1);
  Mat3f M1({1, 3, 2}, {-1, 0, 2}, {3, 1, -1});
  assert(M1 * V1 == V1 * M1.getTranspose());

  // Dot and cross products.
  V1 = Vec3f(2.0f, 4.0f, 5.0f);
  V2 = Vec3f(1.3f, 2.5f, 4.6f);

  float P1 = V1.dot(V2);
  assert(P1 == 2.0f * 1.3f + 4.0f * 2.5f + 5.0f * 4.6f);

  float P2 = V2.dot(V1);
  assert(std::abs(P1 - P2) <= 0.000001);

  V3 = V1.cross(V2);
  assert(std::abs(V3.dot(V1) - 0) < 0.000001 &&
         std::abs(V3.dot(V2) - 0) < 0.000001);
  assert(V1.cross(V2) == -V2.cross(V1));
}
