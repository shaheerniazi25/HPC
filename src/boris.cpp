#include "boris.h"

void boris_step(Particle& p, Vec3 E, Vec3 B, double dt) {
  const double qm = p.q_over_m;
  const double half_dt = 0.5 * dt;

  // 1. Half electric acceleration: v^- = v_{n-1/2} + (q/m) E dt/2
  Vec3 v_minus = p.v + E * (qm * half_dt);

  // 2. Magnetic rotation: v' = v^- + v^- x t, with t = (q/m) B dt/2
  Vec3 t = B * (qm * half_dt);
  Vec3 v_prime = v_minus + cross(v_minus, t);

  // 3. Second half of rotation: v^+ = v^- + v' x s, s = 2t/(1+t^2)
  Vec3 s = t * (2.0 / (1.0 + norm2(t)));
  Vec3 v_plus = v_minus + cross(v_prime, s);

  // 4. Second half electric acceleration
  Vec3 v_new = v_plus + E * (qm * half_dt);

  // 5. Drift: x_{n+1} = x_n + v_{n+1/2} dt
  p.v = v_new;
  p.x = p.x + p.v * dt;
}
