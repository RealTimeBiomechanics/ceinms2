/*This is a simple test from the official boost documentation to assess whether boost libraries have
 * been installed correctly
 */
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/differentiation/finite_difference.hpp>

#include <iostream>
#include <vector>

template<typename T>
T fourth_power(T const &x) {
    T x4 = x * x;// retval in operator*() uses x4's memory via NRVO.
    x4 *= x4;// No copies of x4 are made within operator*=() even when squaring.
    return x4;// x4 uses y's memory in main() via NRVO.
}

int main() {
    using namespace boost::math::differentiation;

    constexpr unsigned Order = 5;// Highest order derivative to be calculated.
    auto const x = make_fvar<double, Order>(2.0);// Find derivatives at x=2.
    auto const y = fourth_power(x);
    const std::vector<double> expected{ 16., 32., 48., 48., 24., 0. };
    bool failed = false;
    for (unsigned i = 0; i <= Order; ++i) {
        std::cout << "y.derivative(" << i << ") = " << y.derivative(i) << std::endl;
        failed |= (y.derivative(i) != expected.at(i)); 
    }
    return failed;
}

/*
Output:
y.derivative(0) = 16
y.derivative(1) = 32
y.derivative(2) = 48
y.derivative(3) = 48
y.derivative(4) = 24
y.derivative(5) = 0
*/
