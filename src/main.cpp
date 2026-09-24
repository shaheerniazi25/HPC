#include <cmath>
#include <cstdio>
#include <deque>
#include <random>
#include <string>
#include <vector>

#include <SFML/Graphics.hpp>

#include "boris.h"
#include "fields.h"
#include "particle.h"

// ---------- setup ----------

static std::vector<Particle> make_cloud(int n, unsigned seed = 42) {
  std::vector<Particle> ps;
  ps.reserve(n);
  // Quad demo: fixed corner starts so bounces are visible immediately.
  if (n == 4) {
    const double xs[4] = {-5.0, 5.0, -5.0, 5.0};
    const double ys[4] = {4.0, 4.0, -4.0, -4.0};
    const double phs[4] = {0.0, 1.7, 3.1, 4.5};
    for (int i = 0; i < 4; ++i) {
      Particle p;
      p.x = {xs[i], ys[i], 0.0};
      p.v = {2.0 * std::cos(phs[i]), 2.0 * std::sin(phs[i]), 0.0};
      p.q_over_m = -1.0;
      ps.push_back(p);
    }
    return ps;
  }
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> pos_x(-8.0, 8.0);
  std::uniform_real_distribution<double> pos_y(-6.0, 6.0);
  std::uniform_real_distribution<double> phase(0.0, 2.0 * M_PI);
  std::uniform_real_distribution<double> speed(1.5, 2.5);

  for (int i = 0; i < n; ++i) {
    double ph = phase(rng);
    double v = speed(rng);
    Particle p;
    p.x = {pos_x(rng), pos_y(rng), 0.0};
    p.v = {v * std::cos(ph), v * std::sin(ph), 0.2 * (v - 2.0)};
    p.q_over_m = -1.0;  // electron-like, normalized units
    ps.push_back(p);
  }
  return ps;
}

// Headless validation: returns 0 on pass. Prints numbers for learning.
static int run_tests() {
  FieldParams fp;
  fp.B0 = 1.0;
  fp.E0 = 1.0;
  const double dt = 0.02;
  bool ok = true;

  // Test 1: pure B -> speed |v| must be conserved (rotation only).
  {
    Particle p{{0, 0, 0}, {2.0, 0.0, 0.0}, -1.0};
    double v0 = norm(p.v);
    double t = 0.0;
    for (int i = 0; i < 5000; ++i) {
      Vec3 E = get_electric_field(p.x, t, FieldCase::Cyclotron, fp);
      Vec3 B = get_magnetic_field(p.x, t, FieldCase::Cyclotron, fp);
      boris_step(p, E, B, dt);
      t += dt;
    }
    double v1 = norm(p.v);
    double err = std::abs(v1 - v0) / v0;
    std::printf("[cyclotron] |v0|=%.6f |v1|=%.6f rel_err=%.3e (want <1e-6)\n", v0, v1, err);
    if (err > 1e-6) ok = false;

    // Gyro-radius check: particle starting at origin with v=(2,0,0),
    // B=+z, q/m=-1 gyrates with r = v/|q|B = 2.
    double r_expect = larmor_radius(2.0, p.q_over_m, fp.B0);
    std::printf("[cyclotron] expected Larmor radius r=%.4f\n", r_expect);
  }

  // Test 2: E x B drift -> mean <vy> should equal -E0/B0.
  {
    Particle p{{0, 0, 0}, {2.0, 0.0, 0.0}, -1.0};
    double t = 0.0;
    double sum_vy = 0.0;
    const int n = 20000;  // ~64 gyroperiods
    for (int i = 0; i < n; ++i) {
      Vec3 E = get_electric_field(p.x, t, FieldCase::ExB, fp);
      Vec3 B = get_magnetic_field(p.x, t, FieldCase::ExB, fp);
      boris_step(p, E, B, dt);
      t += dt;
      sum_vy += p.v.y;
      // Keep particle near origin for clean mean (drift is uniform so this is exact).
      if (std::abs(p.x.y) > 20.0) p.x.y = 0.0;
      if (std::abs(p.x.x) > 20.0) p.x.x = 0.0;
    }
    Vec3 E0v{fp.E0, 0, 0}, B0v{0, 0, fp.B0};
    Vec3 vd = exb_drift(E0v, B0v);  // expect (0,-1,0)
    double mean_vy = sum_vy / n;
    std::printf("[exb] theory vd=(%.3f,%.3f,%.3f) measured <vy>=%.4f\n", vd.x, vd.y, vd.z,
                mean_vy);
    if (std::abs(mean_vy - vd.y) > 0.05) ok = false;
  }

  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}

// Fixed domain box (half-extents). Particles reflect off these walls in quad mode.
constexpr double BOX_X = 11.0;
constexpr double BOX_Y = 9.0;

// Specular reflection: clamp position, flip normal velocity. Preserves |v|.
static void reflect_walls(Particle& p) {
  if (p.x.x > BOX_X) {
    p.x.x = BOX_X;
    p.v.x = -p.v.x;
  }
  if (p.x.x < -BOX_X) {
    p.x.x = -BOX_X;
    p.v.x = -p.v.x;
  }
  if (p.x.y > BOX_Y) {
    p.x.y = BOX_Y;
    p.v.y = -p.v.y;
  }
  if (p.x.y < -BOX_Y) {
    p.x.y = -BOX_Y;
    p.v.y = -p.v.y;
  }
}

// Color particle by speed: slow=blue, fast=red.
static sf::Color speed_color(double v, double vmin = 0.0, double vmax = 3.5) {
  double f = (v - vmin) / (vmax - vmin);
  f = std::max(0.0, std::min(1.0, f));
  return sf::Color(static_cast<sf::Uint8>(f * 255), 60,
                   static_cast<sf::Uint8>((1.0 - f) * 255));
}

int main(int argc, char** argv) {
  int n_particles = 4;
  FieldCase field_case = FieldCase::Cyclotron;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--test") return run_tests();
    if (a == "--exb") field_case = FieldCase::ExB;
    if (a == "--cyclotron") field_case = FieldCase::Cyclotron;
    if (a == "--particles" && i + 1 < argc) n_particles = std::stoi(argv[++i]);
  }

  const double dt = 0.05;      // resolves gyromotion: omega_c*dt = 0.05
  const int substeps = 4;      // physics steps per rendered frame
  const FieldParams fp{1.0, 1.0};

  std::vector<Particle> ps = make_cloud(n_particles);
  const bool quad = (n_particles <= 4);  // big-circle + box rendering path
  double sim_time = 0.0;
  bool paused = false;

  // Trails for quad mode: one per particle; otherwise first 5.
  const int n_trails = quad ? n_particles : std::min(5, n_particles);
  const size_t trail_len = 150;
  std::vector<std::deque<Vec3>> trails(n_trails);

  sf::RenderWindow window(sf::VideoMode(1000, 800), "Boris box (4 particles, SFML)");
  window.setFramerateLimit(60);

  sf::Font font;
  bool have_font = font.loadFromFile("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf");
  sf::Text hud;
  if (have_font) {
    hud.setFont(font);
    hud.setCharacterSize(15);
    hud.setFillColor(sf::Color::White);
    hud.setPosition(10, 10);
  }

  // World view: half-width in meters; height follows aspect ratio.
  const double view_half_w = 12.0;

  auto world_to_screen = [&](Vec3 w, unsigned W, unsigned H) {
    double half_h = view_half_w * H / W;
    float sx = static_cast<float>((w.x / view_half_w * 0.5 + 0.5) * W);
    float sy = static_cast<float>((-w.y / half_h * 0.5 + 0.5) * H);
    return sf::Vector2f(sx, sy);
  };

  while (window.isOpen()) {
    sf::Event ev;
    while (window.pollEvent(ev)) {
      if (ev.type == sf::Event::Closed) window.close();
      if (ev.type == sf::Event::KeyPressed) {
        if (ev.key.code == sf::Keyboard::Escape) window.close();
        if (ev.key.code == sf::Keyboard::Space) paused = !paused;
        if (ev.key.code == sf::Keyboard::Num1) {
          field_case = FieldCase::Cyclotron;
          ps = make_cloud(n_particles);
          trails = std::vector<std::deque<Vec3>>(n_trails);
          sim_time = 0.0;
        }
        if (ev.key.code == sf::Keyboard::Num2) {
          field_case = FieldCase::ExB;
          ps = make_cloud(n_particles);
          trails = std::vector<std::deque<Vec3>>(n_trails);
          sim_time = 0.0;
        }
        if (ev.key.code == sf::Keyboard::R) {
          ps = make_cloud(n_particles);
          trails = std::vector<std::deque<Vec3>>(n_trails);
          sim_time = 0.0;
        }
      }
    }

    if (!paused) {
      for (int s = 0; s < substeps; ++s) {
        for (auto& p : ps) {
          Vec3 E = get_electric_field(p.x, sim_time, field_case, fp);
          Vec3 B = get_magnetic_field(p.x, sim_time, field_case, fp);
          boris_step(p, E, B, dt);
          if (quad) {
            reflect_walls(p);  // bounce inside fixed box
          } else {
            // Wrap around so the cloud stays on screen (uniform fields: exact).
            if (p.x.x > BOX_X) p.x.x = -BOX_X;
            if (p.x.x < -BOX_X) p.x.x = BOX_X;
            if (p.x.y > BOX_Y) p.x.y = -BOX_Y;
            if (p.x.y < -BOX_Y) p.x.y = BOX_Y;
          }
        }
        sim_time += dt;
      }
      for (int k = 0; k < n_trails; ++k) {
        trails[k].push_back(ps[k].x);
        if (trails[k].size() > trail_len) trails[k].pop_front();
      }
    }

    window.clear(sf::Color(10, 10, 25));

    // Fixed domain box outline (matches BOX_X/BOX_Y physics walls).
    if (quad) {
      sf::Vector2f tl = world_to_screen({-BOX_X, BOX_Y, 0}, window.getSize().x,
                                        window.getSize().y);
      sf::Vector2f br = world_to_screen({BOX_X, -BOX_Y, 0}, window.getSize().x,
                                        window.getSize().y);
      sf::RectangleShape box(sf::Vector2f(br.x - tl.x, br.y - tl.y));
      box.setPosition(tl);
      box.setFillColor(sf::Color::Transparent);
      box.setOutlineColor(sf::Color::White);
      box.setOutlineThickness(2.0f);
      window.draw(box);
    }

    if (quad) {
      // 4 large particles: distinct colors, radius ~10px.
      const sf::Color quad_colors[4] = {sf::Color::Red, sf::Color::Green,
                                        sf::Color::Cyan, sf::Color::Yellow};
      for (int i = 0; i < n_particles; ++i) {
        sf::CircleShape c(10.0f);
        c.setOrigin(10.0f, 10.0f);
        c.setPosition(world_to_screen(ps[i].x, window.getSize().x, window.getSize().y));
        c.setFillColor(quad_colors[i % 4]);
        c.setOutlineColor(sf::Color::White);
        c.setOutlineThickness(1.5f);
        window.draw(c);
      }
    } else {
      // Particles as points.
      sf::VertexArray pts(sf::Points, ps.size());
      for (size_t i = 0; i < ps.size(); ++i) {
        pts[i].position = world_to_screen(ps[i].x, window.getSize().x, window.getSize().y);
        pts[i].color = speed_color(norm(ps[i].v));
      }
      window.draw(pts);
    }

    // Trails as line strips.
    const sf::Color quad_trail[4] = {
        sf::Color(255, 100, 100, 200), sf::Color(100, 255, 100, 200),
        sf::Color(100, 255, 255, 200), sf::Color(255, 255, 100, 200)};
    for (int k = 0; k < n_trails; ++k) {
      if (trails[k].size() < 2) continue;
      sf::VertexArray line(sf::LineStrip, trails[k].size());
      size_t j = 0;
      for (Vec3 w : trails[k]) {
        line[j].position = world_to_screen(w, window.getSize().x, window.getSize().y);
        line[j].color = quad ? quad_trail[k % 4] : sf::Color(255, 255, 100, 160);
        ++j;
      }
      window.draw(line);
    }

    if (have_font) {
      Vec3 E0v{fp.E0, 0, 0}, B0v{0, 0, fp.B0};
      Vec3 vd = exb_drift(E0v, B0v);
      char buf[512];
      std::snprintf(buf, sizeof(buf),
                    "%s  t=%.1f  N=%d  dt=%.3f\n"
                    "omega_c=%.2f  r_L(v=2)=%.2f  vd=(%.2f,%.2f)\n"
                    "[1] cyclotron  [2] ExB drift  [Space] pause  [R] reset",
                    field_case == FieldCase::Cyclotron ? "Cyclotron (E=0,Bz)" : "ExB drift",
                    sim_time, n_particles, dt, cyclotron_frequency(-1.0, fp.B0),
                    larmor_radius(2.0, -1.0, fp.B0), vd.x, vd.y);
      hud.setString(buf);
      window.draw(hud);
    }

    window.display();
  }
  return 0;
}
