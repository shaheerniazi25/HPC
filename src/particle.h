#pragma once
#include "vec3.h"

// One charged particle. We use normalized units: q/m carries the sign.
// Position x and velocity v (= v_{n-1/2} before step, v_{n+1/2} after).
struct Particle {
  Vec3 x;
  Vec3 v;
  double q_over_m = -1.0;  // electron-like by default
};
