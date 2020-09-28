//---------------------------------- Ray.h -----------------------------------//
//
// This file is a part of the rndr project, under the Apache License v2.0
// See rndr/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0
//
//----------------------------------------------------------------------------//

#ifndef RNDR_RAY_H
#define RNDR_RAY_H

#include "rndr/Vec3.h"

namespace rndr {

class Ray {
  // pointAt(t) = origin + t * direction
  // Direction is normalized.

  Vec3f origin;
  Vec3f direction;

public:
  Ray(const Vec3f &start, const Vec3f &end)
      : origin(start), direction((end - start).getNormalizedVector()) {}

  Ray(const Ray &other)
      : origin(other.getOrigin()), direction(other.getDirection()) {}

  const Vec3f &getOrigin() const { return origin; }
  const Vec3f &getDirection() const { return direction; }

  Vec3f getPointAtParam(const float t) const { return origin + t * direction; }
};

} // namespace rndr

#endif
