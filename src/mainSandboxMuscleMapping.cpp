#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "ceinms2/Types.h"
#include "ceinms2/NMSmodel.h"
#include "ceinms2/ElectromechanicalDelay.h"
#include "ceinms2/ExponentialActivation.h"

using namespace std;
using namespace Eigen;
using namespace ceinms;

int testMatrix() {
    MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    std::cout << "Here is the matrix m:\n" << m << std::endl;
    VectorXd v(2);
    v(0) = 4;
    v(1) = v(0) - 1;
    std::cout << "Here is the vector v:\n" << v << std::endl;
    return 0;
}


int main() {
    size_t N = 10, M = 10;
    using EMGMapping =  ceinms::MultiInputMultiOutput<Excitation, Excitation>;
    EMGMapping emgGenerator(N, M); 
    
    auto f{ [N,M](const vector<Excitation> &in) {
        MatrixXd m{ MatrixXd::Random(M, N) };
        VectorXd v(N);
        for (int i(0); i < N; ++i) {
            v[i] = in[i];
        }
        MatrixXd x = m * v;
        std::vector<Excitation> out(x.data(), x.data() + x.size());
        return out;
    }};
    emgGenerator.setName("emgGenerator");
    emgGenerator.setFunction(f);
    emgGenerator.setInput(vector<Excitation>(N, 1.));
    emgGenerator.evaluate(0.01);
    for (auto &e : emgGenerator.getOutput())
        cout << e << endl;

    ElectromechanicalDelay delay({ 0.05 });
    delay.setName("mtu1");
    ExponentialActivation act;
    act.setName("mtu1");
    NMSmodel<EMGMapping, ElectromechanicalDelay, ExponentialActivation> model;
    model.addComponent(emgGenerator);
    model.addComponent(delay);
    model.addComponent(act);
    model.connect<EMGMapping, ElectromechanicalDelay>({ "emgGenerator", 0 }, "mtu1");
    model.connect<ElectromechanicalDelay, ExponentialActivation>();



    return 0;

}







