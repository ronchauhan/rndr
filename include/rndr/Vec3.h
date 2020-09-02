//--------------------------------- Vec3.h -----------------------------------//
//
// This file is a part of the rndr project, under the Apache License v2.0
// See rndr/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0
//
//----------------------------------------------------------------------------//

#ifndef RNDR_VEC3_H
#define RNDR_VEC3_H

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace rndr {

template <typename T> class Mat3;

template <typename T> class Vec3 {
  T x = (T(0)), y = (T(0)), z = (T(0));

public:
  Vec3<T>() {}
  Vec3<T>(const T argX, T argY, T argZ) : x(argX), y(argY), z(argZ) {}
  Vec3<T>(const T arg) : x(arg), y(arg), z(arg){};
  Vec3<T>(const Vec3<T> &other)
      : x(other.getX()), y(other.getY()), z(other.getZ()) {}

  T getX() const { return x; }
  T getY() const { return y; }
  T getZ() const { return z; }

  void setX(T argX) { x = argX; }
  void setY(T argY) { y = argY; }
  void setZ(T argZ) { z = argZ; }

  T getSquaredLength() const { return x * x + y * y + z * z; }
  T getLength() const { return sqrt(getSquaredLength()); };

  Vec3<T> getNormalizedVector() const;
  void normalizeThis();

  Vec3<T> operator+() const { return *this; }
  Vec3<T> operator-() const { return Vec3<T>(-x, -y, -z); }

  // Operations with scalar values.
  Vec3<T> operator*(const T &value) const;
  Vec3<T> operator/(const T &value) const;

  Vec3<T> operator*=(const T &value);
  Vec3<T> operator/=(const T &value);

  // Operations with other vectors.
  Vec3<T> operator+(const Vec3<T> &other) const;
  Vec3<T> operator-(const Vec3<T> &other) const;

  Vec3<T> operator+=(const Vec3<T> &other);
  Vec3<T> operator-=(const Vec3<T> &other);

  Vec3<T> operator=(const Vec3<T> &other);
  bool operator==(const Vec3<T> &other) const;
  bool operator!=(const Vec3<T> &other) const;

  // Matrix-vector multiplication (row-major form).
  Vec3<T> operator*(const Mat3<T> &M) const;

  T dot(const Vec3<T> &other) const;
  Vec3<T> cross(const Vec3<T> &other) const;

  friend Vec3<T> operator*(T value, const Vec3<T> &vec) { return vec * value; }

  friend std::ostream &operator<<(std::ostream &OS, const Vec3<T> &V) {
    OS << std::setprecision(10) << V.x << " " << V.y << " " << V.z;
    return OS;
  }

  friend std::istream &operator>>(std::istream &IS, Vec3<T> &V) {
    IS >> V.x >> V.y >> V.z;
    return IS;
  }
};

template <typename T> Vec3<T> Vec3<T>::operator*(const T &value) const {
  return Vec3<T>(x * value, y * value, z * value);
}

template <typename T> Vec3<T> Vec3<T>::operator/(const T &value) const {
  assert(value != T(0) && "Value must be non-zero");
  T dividedVal = T(1) / value;
  return Vec3<T>(x * dividedVal, y * dividedVal, z * dividedVal);
}

template <typename T> Vec3<T> Vec3<T>::operator*=(const T &value) {
  x *= value;
  y *= value;
  z *= value;
  return Vec3<T>(x, y, z);
}

template <typename T> Vec3<T> Vec3<T>::operator/=(const T &value) {
  assert(value != T(0) && "Value must be non-zero");
  T dividedVal = T(1) / value;
  x *= dividedVal;
  y *= dividedVal;
  z *= dividedVal;
  return Vec3<T>(x, y, z);
}

template <typename T> Vec3<T> Vec3<T>::operator+(const Vec3<T> &other) const {
  return Vec3<T>(x + other.getX(), y + other.getY(), z + other.getZ());
}

template <typename T> Vec3<T> Vec3<T>::operator-(const Vec3<T> &other) const {
  return Vec3<T>(x - other.getX(), y - other.getY(), z - other.getZ());
}

template <typename T> Vec3<T> Vec3<T>::operator+=(const Vec3<T> &other) {
  x += other.getX();
  y += other.getY();
  z += other.getZ();
  return Vec3<T>(x, y, z);
}

template <typename T> Vec3<T> Vec3<T>::operator-=(const Vec3<T> &other) {
  x -= other.getX();
  y -= other.getY();
  z -= other.getZ();
  return Vec3<T>(x, y, z);
}

template <typename T> Vec3<T> Vec3<T>::operator=(const Vec3<T> &other) {
  if (this != &other) {
    x = other.getX();
    y = other.getY();
    z = other.getZ();
  }
  return Vec3<T>(x, y, z);
}

template <typename T> bool Vec3<T>::operator==(const Vec3<T> &other) const {
  return std::abs(x - other.getX()) <= T(0.000001) &&
         std::abs(y - other.getY()) <= T(0.000001) &&
         std::abs(z - other.getZ()) <= T(0.000001);
}

template <typename T> bool Vec3<T>::operator!=(const Vec3<T> &other) const {
  return !(*this == other);
}

template <typename T> Vec3<T> Vec3<T>::operator*(const Mat3<T> &M) const {
  return Vec3<T>(M[0][0] * x + M[1][0] * y + M[2][0] * z,
                 M[0][1] * x + M[1][1] * y + M[2][1] * z,
                 M[0][2] * x + M[1][2] * y + M[2][2] * z);
}

template <typename T> T Vec3<T>::dot(const Vec3<T> &other) const {
  return x * other.getX() + y * other.getY() + z * other.getZ();
}

template <typename T> Vec3<T> Vec3<T>::cross(const Vec3<T> &other) const {
  return Vec3<T>(y * other.getZ() - z * other.getY(),
                 -(x * other.getZ() - z * other.getX()),
                 x * other.getY() - y * other.getX());
}

template <typename T> Vec3<T> Vec3<T>::getNormalizedVector() const {
  T length = getLength();
  assert(length != T(0));
  T dividedLength = T(1) / length;
  return Vec3<T>(x * dividedLength, y * dividedLength, z * dividedLength);
}

template <typename T> void Vec3<T>::normalizeThis() {
  T length = getLength();
  assert(length != T(0));
  T dividedLength = T(1) / length;
  x *= dividedLength;
  y *= dividedLength;
  z *= dividedLength;
}

using Vec3f = Vec3<float>;

} // namespace rndr

#endif
