#include "ceinms2/RienerActivation.h"
using namespace ceinms;
#include <iostream>
#include <fstream>
using namespace std;

int main() {

    ofstream outF("out.csv");
    RienerActivation activation;
    double time{ 0. }, dt{ 0.0001 };
    for (int k{ 0 }; k < 10; ++k) {
        for (int i{ 0 }; i < 1000; ++i) {
            if (i == 1)
                activation.i_.pulse = 1.;
            else
                activation.i_.pulse = 0.;
            activation.integrate(dt);
            time += dt;
            outF << time << ", " << activation.i_.pulse << ", " << activation.s_.beta[0] << 
                ", " << activation.s_.gamma[0] << endl;
        }
    }
    return 0;
}

