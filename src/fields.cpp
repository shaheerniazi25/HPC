#include "fields.h"
#include <cmath>

Vec3 get_electric_field(Vec3 /*x*/, double /*t*/, FieldCase c, const FieldParams& p) {
  if (c == FieldCase::ExB) return {p.E0, 0.0, 0.0};
  return {0.0, 0.0, 0.0};  // Cyclotron: no E
}

Vec3 get_magnetic_field(Vec3 /*x*/, double /*t*/, FieldCase /*c*/, const FieldParams& p) {
  return {0.0, 0.0, p.B0};  // uniform Bz in both demos
}

double cyclotron_frequency(double q_over_m, double B0) {
  return std::abs(q_over_m * B0);
}

double larmor_radius(double v_perp, double q_over_m, double B0) {
  return v_perp / cyclotron_frequency(q_over_m, B0);
}

Vec3 exb_drift(Vec3 E, Vec3 B) {
  double b2 = norm2(B);
  if (b2 == 0.0) return {0.0, 0.0, 0.0};
  return cross(E, B) / b2;
}
