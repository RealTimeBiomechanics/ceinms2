/* -------------------------------------------------------------------------- *
 * CEINMS is a standalone toolbox for neuromusculoskeletal modelling and      *
 * simulation. CEINMS can also be used as a plugin for OpenSim either         *
 * through the OpenSim GUI or API. See https://simtk.org/home/ceinms and the  *
 * NOTICE file for more information. CEINMS development was coordinated       *
 * through Griffith University and supported by the Australian National       *
 * Health and Medical Research Council (NHMRC), the US National Institutes of *
 * Health (NIH), and the European Union Framework Programme 7 (EU FP7). Also  *
 * see the PROJECTS file for more information about the funding projects.     *
 *                                                                            *
 * Copyright (c) 2010-2015 Griffith University and the Contributors           *
 *                                                                            *
 * CEINMS Contributors: C. Pizzolato, M. Reggiani, M. Sartori,                *
 *                      E. Ceseracciu, and D.G. Lloyd                         *
 *                                                                            *
 * Author(s): C. Pizzolato, M. Reggiani                                       *
 *                                                                            *
 * CEINMS is licensed under the Apache License, Version 2.0 (the "License").  *
 * You may not use this file except in compliance with the License. You may   *
 * obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.*
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#ifndef ceinms_Curve_h
#define ceinms_Curve_h

#include <vector>
#include <cmath>
#include <limits>
#include "CircularVector.h"

namespace ceinms {

template<typename T>
bool nearly_equal(T a, T b) {
    return std::nextafter(a, std::numeric_limits<T>::lowest()) <= b
           && std::nextafter(a, std::numeric_limits<T>::max()) >= b;
}

namespace CurveMode {
    enum Mode { Online,
        Offline };
    enum Interpolation { Cubic,
        Linear };
}// namespace CurveMode

template<int v>
struct Int2Type {
    enum { value = v };
};

template<CurveMode::Mode mode, typename T, typename U>
struct Select {
    typedef T Result;
};

template<typename T, typename U>
struct Select<CurveMode::Online, T, U> {
    typedef U Result;
};

template<CurveMode::Mode mode, CurveMode::Interpolation T, size_t N>
class Curve;

template<CurveMode::Mode mode, CurveMode::Interpolation T, size_t N>
std::ostream &operator<<(std::ostream &output, const Curve<mode, T, N> &c);

template<CurveMode::Mode mode, CurveMode::Interpolation T = CurveMode::Cubic, size_t N = 15>
class Curve {
  public:
    typedef typename Select<mode, std::vector<double>, CircularVector<double, N>>::Result VectorType;
    Curve();
    // compute coefficients
    Curve(const std::vector<double> &x, const std::vector<double> &y);
    Curve(const Curve &orig);
    virtual ~Curve() {}
    Curve &operator=(const Curve &orig);
    void reset();
    // add a new points and compute again coefficients
    bool addPoint(double x, double y);
    bool addPointNoUpdate(double x, double y);
    void resetPointsWith(const std::vector<double> &x, const std::vector<double> &y);
    void refresh();
    // remove last point of the Curve_c and compute again coefficients
    void removeLastPoint();
    void removeLastPointNoUpdate();//remove last point without computing the coefficients again
    // interpolation
    double getValue(double xValue) const;
    double getFirstDerivative(double xValue) const;
    double getSecondDerivative(double xValue) const;
    double getMinX() const {
        if (!x_.empty())
            return x_.front();
        return 0;
    }
    double getMaxX() const {
        if (!x_.empty())
            return x_.back();
        return 0;
    }

    VectorType getXNodes() const { return x_; }
    VectorType getYNodes() const { return y_; }

    friend std::ostream &operator<<<>(std::ostream &output, const Curve &c);
    bool empty() { return b_.empty(); }
    std::size_t getNoElements() const { return x_.size(); }

  private:
    void computeCoefficients(Int2Type<CurveMode::Cubic>);
    void computeCoefficients(Int2Type<CurveMode::Linear>);

    double getValue(double xValue, size_t abscissaPoint, Int2Type<CurveMode::Cubic>) const;
    double getValue(double xValue, size_t abscissaPoint, Int2Type<CurveMode::Linear>) const;

    double getFirstDerivative(double xValue, size_t abscissaPoint, Int2Type<CurveMode::Cubic>) const;
    double getFirstDerivative(double xValue, size_t abscissaPoint, Int2Type<CurveMode::Linear>) const;

    size_t getAbscissaPoint(double xValue) const;
    size_t getAbscissaPoint(double xValue, Int2Type<CurveMode::Online>) const;
    size_t getAbscissaPoint(double xValue, Int2Type<CurveMode::Offline>) const;


    VectorType x_;
    VectorType y_;
    std::vector<double> b_;
    std::vector<double> c_;
    std::vector<double> d_;
};
}// namespace ceinms

#include "Curve.cpp"
#endif
