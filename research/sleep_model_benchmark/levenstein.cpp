// Independent integration of Levenstein et al. (2019), Eqs. 1–2 and Supplementary Note 1.
// OU convention and default theta/sigma follow author dependencies/OUNoise.m.
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>

struct State { double r, a; };

double sigmoid(double v) { return 1.0 / (1.0 + std::exp(-v)); }

State derivative(State y, double input, double w, double b) {
    return {sigmoid(w*y.r - b*y.a + input - 5.0) - y.r,
            (sigmoid(15.0 * (y.r - 0.5)) - y.a) / 25.0};
}

State step(State y, State direction, double dt) {
    return {y.r + dt*direction.r, y.a + dt*direction.a};
}

int main(int argc, char **argv) {
    if (argc != 9 && argc != 12) {
        std::cerr << "usage: levenstein seed duration_AU steps_per_AU I w b sigma theta [drive_gate recurrence_gate period_AU]\n";
        return 2;
    }
    const auto seed = std::stoul(argv[1]);
    const int duration = std::stoi(argv[2]);
    const int substeps = std::stoi(argv[3]);
    const double input = std::stod(argv[4]), w = std::stod(argv[5]), b = std::stod(argv[6]);
    const double sigma = std::stod(argv[7]), theta = std::stod(argv[8]);
    if (duration <= 0 || substeps < 1 || !std::isfinite(input+w+b+sigma+theta)
            || sigma < 0 || theta <= 0) return 2;
    const double drive_gate = argc == 12 ? std::stod(argv[9]) : 0;
    const double recurrence_gate = argc == 12 ? std::stod(argv[10]) : 0;
    const double period = argc == 12 ? std::stod(argv[11]) : 500;
    if (period <= 125 || drive_gate < 0 || recurrence_gate < 0) return 2;
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> gaussian;
    const double dt = 1.0 / substeps;
    const double diffusion = sigma*std::sqrt(2*theta*dt);
    State y{0.5,0.5};
    double noise = sigma*gaussian(rng);
    std::cout << std::setprecision(12);
    for (int t=0; t<duration; ++t) {
        std::cout << y.r << '\n';
        for (int s=0; s<substeps; ++s) {
            const double next = noise*(1-theta*dt) + diffusion*gaussian(rng);
            const bool closed = std::fmod(t + s*dt, period) < 62.5;
            const double current = input - (closed ? drive_gate : 0);
            const double gain = w - (closed ? recurrence_gate : 0);
            const double middle = current + (noise+next)/2;
            const State a = derivative(y, current+noise, gain, b);
            const State c = derivative(step(y,a,dt/2), middle, gain, b);
            const State d = derivative(step(y,c,dt/2), middle, gain, b);
            const State e = derivative(step(y,d,dt), current+next, gain, b);
            y.r += dt*(a.r + 2*c.r + 2*d.r + e.r)/6;
            y.a += dt*(a.a + 2*c.a + 2*d.a + e.a)/6;
            noise = next;
        }
    }
}
