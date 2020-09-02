//--------------------------------- Mat3.h -----------------------------------//
//
// This file is a part of the rndr project, under the Apache License v2.0
// See rndr/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0
//
//----------------------------------------------------------------------------//

#ifndef RNDR_MAT3_H
#define RNDR_MAT3_H

#include <cassert>
#include <cmath>
#include <initializer_list>
#include <iomanip>
#include <iostream>

namespace rndr {

template <typename T> class Vec3;

// Simple implementation of a 3x3 matrix.
template <typename T> class Mat3 {
  T elems[3][3] = {{T(1), T(0), T(0)}, {T(0), T(1), T(0)}, {T(0), T(0), T(1)}};

public:
  using Row = std::initializer_list<T>;

  Mat3<T>(){};
  Mat3<T>(const T e1, const T e2, const T e3, const T e4, const T e5,
          const T e6, const T e7, const T e8, const T e9);
  Mat3<T>(const Row &R1, const Row &R2, const Row &R3);
  Mat3<T>(const T arg);
  Mat3<T>(const Mat3<T> &other);

  T getDeterminant() const;

  Mat3<T> getTranspose() const;
  void transposeThis();

  Mat3<T> getInverse() const;
  void inverseThis();

  // Operations over self.
  Mat3<T> operator+() const { return *this; }
  Mat3<T> operator-() const;

  T *operator[](unsigned int i) { return elems[i]; }
  const T *operator[](unsigned int i) const { return elems[i]; }

  // Operations with constants.
  Mat3<T> operator*(const T &value) const;
  Mat3<T> operator/(const T &value) const;

  Mat3<T> operator*=(const T &value);
  Mat3<T> operator/=(const T &value);

  // Operations with other matrices.
  Mat3<T> operator+(const Mat3<T> &other) const;
  Mat3<T> operator-(const Mat3<T> &other) const;
  Mat3<T> operator*(const Mat3<T> &other) const;

  Mat3<T> operator+=(const Mat3<T> &other);
  Mat3<T> operator-=(const Mat3<T> &other);
  Mat3<T> operator*=(const Mat3<T> &other);

  Mat3<T> operator=(const Mat3<T> &other);
  bool operator==(const Mat3<T> &other) const;
  bool operator!=(const Mat3<T> &other) const;

  // Matrix-vector multiplication (column-major form).
  Vec3<T> operator*(const Vec3<T> &V) const;

  friend Mat3<T> operator*(T value, const Mat3<T> &M) { return M * value; }

  friend std::ostream &operator<<(std::ostream &OS, const Mat3<T> &M) {
    OS << std::setprecision(10);
    OS << M.elems[0][0] << '\t' << M.elems[0][1] << "\t" << M.elems[0][2]
       << '\n';
    OS << M.elems[1][0] << '\t' << M.elems[1][1] << "\t" << M.elems[1][2]
       << '\n';
    OS << M.elems[2][0] << '\t' << M.elems[2][1] << "\t" << M.elems[2][2];
    return OS;
  }

  friend std::istream &operator>>(std::istream &IS, Mat3<T> &M) {
    IS >> M.elems[0][0] >> M.elems[0][1] >> M.elems[0][2];
    IS >> M.elems[1][0] >> M.elems[1][1] >> M.elems[1][2];
    IS >> M.elems[2][0] >> M.elems[2][1] >> M.elems[2][2];
    return IS;
  }
};

template <typename T>
Mat3<T>::Mat3(const T e1, const T e2, const T e3, const T e4, const T e5,
              const T e6, const T e7, const T e8, const T e9) {
  elems[0][0] = e1;
  elems[0][1] = e2;
  elems[0][2] = e3;
  elems[1][0] = e4;
  elems[1][1] = e5;
  elems[1][2] = e6;
  elems[2][0] = e7;
  elems[2][1] = e8;
  elems[2][2] = e9;
}

template <typename T>
Mat3<T>::Mat3(const Row &R1, const Row &R2, const Row &R3) {
  // std::initializer_list doesn't provide a subscript operator, so we use
  // the iterator instead.
  for (int i = 0; i < 3; ++i) {
    elems[0][i] = *(R1.begin() + i);
    elems[1][i] = *(R2.begin() + i);
    elems[2][i] = *(R3.begin() + i);
  }
}

template <typename T> Mat3<T>::Mat3(const T arg) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      elems[i][j] = arg;
    }
  }
}

template <typename T> Mat3<T>::Mat3(const Mat3<T> &other) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      elems[i][j] = other[i][j];
    }
  }
}

template <typename T> Mat3<T> Mat3<T>::operator-() const {
  return Mat3<T>(-elems[0][0], -elems[0][1], -elems[0][2], -elems[1][0],
                 -elems[1][1], -elems[1][2], -elems[2][0], -elems[2][1],
                 -elems[2][2]);
}
template <typename T> T Mat3<T>::getDeterminant() const {
  return elems[0][0] * (elems[1][1] * elems[2][2] - elems[2][1] * elems[1][2]) -
         elems[0][1] * (elems[1][0] * elems[2][2] - elems[2][0] * elems[1][2]) +
         elems[0][2] * (elems[1][0] * elems[2][1] - elems[2][0] * elems[1][1]);
}

template <typename T> Mat3<T> Mat3<T>::getTranspose() const {
  return Mat3<T>(elems[0][0], elems[1][0], elems[2][0], elems[0][1],
                 elems[1][1], elems[2][1], elems[0][2], elems[1][2],
                 elems[2][2]);
}

template <typename T> void Mat3<T>::transposeThis() {
  std::swap(elems[0][1], elems[1][0]);
  std::swap(elems[2][0], elems[0][2]);
  std::swap(elems[1][2], elems[2][1]);
}

template <typename T> Mat3<T> Mat3<T>::getInverse() const {
  T det = getDeterminant();
  assert(std::abs(det - T(0)) > T(0.00001));

  T Cof_00 = elems[1][1] * elems[2][2] - elems[2][1] * elems[1][2];
  T Cof_01 = -1 * (elems[1][0] * elems[2][2] - elems[2][0] * elems[1][2]);
  T Cof_02 = elems[1][0] * elems[2][1] - elems[2][0] * elems[1][1];

  T Cof_10 = -1 * (elems[0][1] * elems[2][2] - elems[2][1] * elems[0][2]);
  T Cof_11 = elems[0][0] * elems[2][2] - elems[2][0] * elems[0][2];
  T Cof_12 = -1 * (elems[0][0] * elems[2][1] - elems[2][0] * elems[0][1]);

  T Cof_20 = elems[0][1] * elems[1][2] - elems[1][1] * elems[0][2];
  T Cof_21 = -1 * (elems[0][0] * elems[1][2] - elems[1][0] * elems[0][2]);
  T Cof_22 = elems[0][0] * elems[1][1] - elems[1][0] * elems[0][1];

  return (1 / det) * Mat3<T>({Cof_00, Cof_10, Cof_20}, {Cof_01, Cof_11, Cof_21},
                             {Cof_02, Cof_12, Cof_22});
}

template <typename T> void Mat3<T>::inverseThis() { *this = getInverse(); }

template <typename T> Mat3<T> Mat3<T>::operator*(const T &value) const {
  return Mat3<T>(elems[0][0] * value, elems[0][1] * value, elems[0][2] * value,
                 elems[1][0] * value, elems[1][1] * value, elems[1][2] * value,
                 elems[2][0] * value, elems[2][1] * value, elems[2][2] * value);
}

template <typename T> Mat3<T> Mat3<T>::operator/(const T &value) const {
  assert(value != T(0) && "Value must be non-zero");
  T dividedValue = T(1) / value;
  return Mat3<T>(elems[0][0] * dividedValue, elems[0][1] * dividedValue,
                 elems[0][2] * dividedValue, elems[1][0] * dividedValue,
                 elems[1][1] * dividedValue, elems[1][2] * dividedValue,
                 elems[2][0] * dividedValue, elems[2][1] * dividedValue,
                 elems[2][2] * dividedValue);
}

template <typename T> Mat3<T> Mat3<T>::operator*=(const T &value) {
  elems[0][0] *= value;
  elems[0][1] *= value;
  elems[0][2] *= value;
  elems[1][0] *= value;
  elems[1][1] *= value;
  elems[1][2] *= value;
  elems[2][0] *= value;
  elems[2][1] *= value;
  elems[2][2] *= value;
  return Mat3<T>(*this);
}

template <typename T> Mat3<T> Mat3<T>::operator/=(const T &value) {
  assert(value != T(0) && "Value must be non-zero");
  T dividedValue = T(1) / value;
  elems[0][0] *= dividedValue;
  elems[0][1] *= dividedValue;
  elems[0][2] *= dividedValue;
  elems[1][0] *= dividedValue;
  elems[1][1] *= dividedValue;
  elems[1][2] *= dividedValue;
  elems[2][0] *= dividedValue;
  elems[2][1] *= dividedValue;
  elems[2][2] *= dividedValue;
  return Mat3<T>(*this);
}

template <typename T> Mat3<T> Mat3<T>::operator+(const Mat3<T> &other) const {
  return Mat3<T>(elems[0][0] + other[0][0], elems[0][1] + other[0][1],
                 elems[0][2] + other[0][2], elems[1][0] + other[1][0],
                 elems[1][1] + other[1][1], elems[1][2] + other[1][2],
                 elems[2][0] + other[2][0], elems[2][1] + other[2][1],
                 elems[2][2] + other[2][2]);
}

template <typename T> Mat3<T> Mat3<T>::operator-(const Mat3<T> &other) const {
  return Mat3<T>(elems[0][0] - other[0][0], elems[0][1] - other[0][1],
                 elems[0][2] - other[0][2], elems[1][0] - other[1][0],
                 elems[1][1] - other[1][1], elems[1][2] - other[1][2],
                 elems[2][0] - other[2][0], elems[2][1] - other[2][1],
                 elems[2][2] - other[2][2]);
}

template <typename T> Mat3<T> Mat3<T>::operator*(const Mat3<T> &other) const {
  return Mat3<T>(other[0][0] * elems[0][0] + other[1][0] * elems[0][1] +
                     other[2][0] * elems[0][2],
                 other[0][1] * elems[0][0] + other[1][1] * elems[0][1] +
                     other[2][1] * elems[0][2],
                 other[0][2] * elems[0][0] + other[1][2] * elems[0][1] +
                     other[2][2] * elems[0][2],

                 other[0][0] * elems[1][0] + other[1][0] * elems[1][1] +
                     other[2][0] * elems[1][2],
                 other[0][1] * elems[1][0] + other[1][1] * elems[1][1] +
                     other[2][1] * elems[1][2],
                 other[0][2] * elems[1][0] + other[1][2] * elems[1][1] +
                     other[2][2] * elems[1][2],

                 other[0][0] * elems[2][0] + other[1][0] * elems[2][1] +
                     other[2][0] * elems[2][2],
                 other[0][1] * elems[2][0] + other[1][1] * elems[2][1] +
                     other[2][1] * elems[2][2],
                 other[0][2] * elems[2][0] + other[1][2] * elems[2][1] +
                     other[2][2] * elems[2][2]);
}

template <typename T> Mat3<T> Mat3<T>::operator+=(const Mat3<T> &other) {
  elems[0][0] += other[0][0];
  elems[0][1] += other[0][1];
  elems[0][2] += other[0][2];
  elems[1][0] += other[1][0];
  elems[1][1] += other[1][1];
  elems[1][2] += other[1][2];
  elems[2][0] += other[2][0];
  elems[2][1] += other[2][1];
  elems[2][2] += other[2][2];
  return Mat3<T>(*this);
}

template <typename T> Mat3<T> Mat3<T>::operator-=(const Mat3<T> &other) {
  elems[0][0] -= other[0][0];
  elems[0][1] -= other[0][1];
  elems[0][2] -= other[0][2];
  elems[1][0] -= other[1][0];
  elems[1][1] -= other[1][1];
  elems[1][2] -= other[1][2];
  elems[2][0] -= other[2][0];
  elems[2][1] -= other[2][1];
  elems[2][2] -= other[2][2];
  return Mat3<T>(*this);
}

template <typename T> Mat3<T> Mat3<T>::operator*=(const Mat3<T> &other) {
  elems[0][0] = other[0][0] * elems[0][0] + other[1][0] * elems[0][1] +
                other[2][0] * elems[0][2];
  elems[0][1] = other[0][1] * elems[0][0] + other[1][1] * elems[0][1] +
                other[2][1] * elems[0][2];
  elems[0][2] = other[0][2] * elems[0][0] + other[1][2] * elems[0][1] +
                other[2][2] * elems[0][2];
  elems[1][0] = other[0][0] * elems[1][0] + other[1][0] * elems[1][1] +
                other[2][0] * elems[1][2];
  elems[1][1] = other[0][1] * elems[1][0] + other[1][1] * elems[1][1] +
                other[2][1] * elems[1][2];
  elems[1][2] = other[0][2] * elems[1][0] + other[1][2] * elems[1][1] +
                other[2][2] * elems[1][2];
  elems[2][0] = other[0][0] * elems[2][0] + other[1][0] * elems[2][1] +
                other[2][0] * elems[2][2];
  elems[2][1] = other[0][1] * elems[2][0] + other[1][1] * elems[2][1] +
                other[2][1] * elems[2][2];
  elems[2][2] = other[0][2] * elems[2][0] + other[1][2] * elems[2][1] +
                other[2][2] * elems[2][2];
  return Mat3<T>(*this);
}

template <typename T> Mat3<T> Mat3<T>::operator=(const Mat3<T> &other) {
  if (this != &other) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        elems[i][j] = other[i][j];
      }
    }
  }

  return Mat3<T>(elems[0][0], elems[0][1], elems[0][2], elems[1][0],
                 elems[1][1], elems[1][2], elems[2][0], elems[2][1],
                 elems[2][2]);
}

template <typename T> bool Mat3<T>::operator==(const Mat3<T> &other) const {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (std::abs(elems[i][j] - other[i][j]) > T(0.000001))
        return false;
    }
  }
  return true;
}

template <typename T> bool Mat3<T>::operator!=(const Mat3<T> &other) const {
  return !(*this == other);
}

template <typename T> Vec3<T> Mat3<T>::operator*(const Vec3<T> &V) const {
  return Vec3<T>(
      elems[0][0] * V.getX() + elems[0][1] * V.getY() + elems[0][2] * V.getZ(),
      elems[1][0] * V.getX() + elems[1][1] * V.getY() + elems[1][2] * V.getZ(),
      elems[2][0] * V.getX() + elems[2][1] * V.getY() + elems[2][2] * V.getZ());
}

using Mat3f = Mat3<float>;

} // namespace rndr

#endif
