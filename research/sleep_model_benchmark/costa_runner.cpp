// Minimal CLI for the pinned Costa et al. C++ model, matching the unstimulated MEX loop.
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "Cortical_Column.h"
#include "ODE.h"
#include "Thalamic_Column.h"

extern const int onset = 20;
extern const int res = 10000;
extern const int red = 100;
extern const double dt = 1000.0/res;
extern const double h = std::sqrt(dt);

// Match the friend declaration. The original Data_Storage.h defines a
// reference overload without friend access in this archived revision.
void get_data(int index, Cortical_Column &cortex, Thalamic_Column &thalamus,
              std::vector<double*> output) {
    output[0][index]=cortex.Vp[0];
    output[1][index]=thalamus.Vt[0];
}

int main(int argc, char **argv) {
    if (argc!=3) return 2;
    const auto seed=std::stoul(argv[1]);
    const int duration=std::stoi(argv[2]);
    if (duration<1 || duration>3600) return 2;
    std::array<double,3> cortex_params;
    std::array<double,2> thalamus_params;
    std::array<double,4> connections;
    for (double &d:cortex_params) if (!(std::cin>>d)) return 2;
    for (double &d:thalamus_params) if (!(std::cin>>d)) return 2;
    for (double &d:connections) if (!(std::cin>>d)) return 2;
    std::srand(seed);
    Cortical_Column cortex(cortex_params.data(),connections.data());
    Thalamic_Column thalamus(thalamus_params.data(),connections.data());
    cortex.get_Thalamus(thalamus);
    thalamus.get_Cortex(cortex);
    std::array<double,2> sample;
    std::vector<double*> output{&sample[0],&sample[1]};
    std::cout << std::setprecision(12);
    for (int t=0; t<(duration+onset)*res; ++t) {
        ODE(cortex,thalamus);
        if (t>=onset*res && t%red==0) {
            get_data(0,cortex,thalamus,output);
            std::cout<<sample[0]<<' '<<sample[1]<<'\n';
        }
    }
}
