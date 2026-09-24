#pragma once
#include "vec3.h"

enum class FieldCase {
  Cyclotron,  // E = 0, uniform Bz: circular gyro-orbits
  ExB,        // uniform Ex + Bz: gyro-motion + E x B drift
};

struct FieldParams {
  double B0 = 1.0;  // Tesla (normalized)
  double E0 = 1.0;  // V/m (normalized)
};

Vec3 get_electric_field(Vec3 x, double t, FieldCase c, const FieldParams& p);
Vec3 get_magnetic_field(Vec3 x, double t, FieldCase c, const FieldParams& p);

// Theory helpers used for HUD + validation.
double cyclotron_frequency(double q_over_m, double B0);  // |omega_c|
double larmor_radius(double v_perp, double q_over_m, double B0);
Vec3 exb_drift(Vec3 E, Vec3 B);
