// Independent discrete branching process based on González et al., Methods.
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>

constexpr int units = 50;

int main(int argc, char **argv) {
    if (argc!=8) return 2;
    const auto seed=std::stoul(argv[1]);
    const int seconds=std::stoi(argv[2]), width=std::stoi(argv[5]), period=std::stoi(argv[6]);
    const double sigma=std::stod(argv[3]);
    const char gate=argv[4][0];
    const double gated_fraction=std::stod(argv[7]);
    if (seconds<1 || seconds>1010 || sigma<0 || width<0 || period<1 || period<width ||
        gated_fraction<0 || gated_fraction>1 || std::string("NDP").find(gate)==std::string::npos) return 2;
    std::mt19937_64 graph(seed), random(seed+200000);
    std::uniform_real_distribution<double> uniform(0,1);
    std::array<uint64_t,units> edges{};
    int connections=0;
    for (int i=0;i<units;++i)
        for (int j=i+1;j<units;++j)
            if (uniform(graph)<.03) {
                edges[i]|=UINT64_C(1)<<j;
                edges[j]|=UINT64_C(1)<<i;
                ++connections;
            }
    const double degree=2.0*connections/units;
    if (degree<sigma || degree<=0) return 3;
    const double p=sigma/degree;
    std::array<double,units+1> propagation{};
    for (int i=0;i<=units;++i) propagation[i]=1-std::pow(1-p,i);
    uint64_t spiking=0,refractory=0;
    const double noise=1-std::exp(-.014);
    const uint64_t gated=(UINT64_C(1)<<static_cast<int>(gated_fraction*units))-1;
    for (int t=0;t<seconds*1000;++t) {
        const bool closed = t%period<width;
        uint64_t next=0;
        for (int i=0;i<units;++i) {
            const uint64_t bit=UINT64_C(1)<<i;
            const double external=uniform(random), internal=uniform(random);
            if ((spiking|refractory)&bit) continue;
            const bool blocked = closed && (bit&gated);
            const uint64_t sources = closed && gate=='P' ? spiking&~gated : spiking;
            if ((!(gate=='D'&&blocked) && external<noise) ||
                internal<propagation[__builtin_popcountll(sources&edges[i])]) next|=bit;
        }
        refractory=spiking;
        spiking=next;
        std::cout.write(reinterpret_cast<const char*>(&spiking),sizeof(spiking));
    }
}
