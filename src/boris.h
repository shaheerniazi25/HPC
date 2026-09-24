#pragma once
#include "particle.h"
#include "vec3.h"

// Standard non-relativistic Boris push (Boris 1970).
// Advances v by dt in given E,B and then x by the new v.
// E,B are evaluated at x_n by the caller.
void boris_step(Particle& p, Vec3 E, Vec3 B, double dt);
