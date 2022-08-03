#ifndef __MATH_COMMON_H__
#define __MATH_COMMON_H__

#include <complex>
#include <memory>
#include <functional>
#include <vector>
#include <cfloat>

namespace maths {
    class IMatrix;
    class IVector;
    class IGMRESSolver;
    class IEigenSolver;

    typedef std::complex<double> complex;
    typedef std::function<complex(int, int)> FuncOfRowCol;
    typedef std::shared_ptr<IMatrix> Matrix;
    typedef std::shared_ptr<IVector> Vector;
    typedef std::shared_ptr<IGMRESSolver> GMRESSolver;
    typedef std::shared_ptr<IEigenSolver> EigenSolver;

    // constants
    const double c = 137.036;                       // speed of light in atomic units
    const double Pi = 3.14159265358979323846;       // pi
    const double LnmToEnergy = 45.5633346744;       // factor to convert from wavelength (nm) to energy (au)
    const double Infinity = DBL_MAX;                // big number
    
    enum DimIndex {
        X = 0,
        Y,
        Z,

        NUM
    };

    // generates a list of all elements in an interval
    template <typename T>
    std::vector<T> range(T start, T end, T step = 1) {
        std::vector<T> r; r.reserve((end - start)/step);
        for (T i = start; i < end; i+=step)
            r.push_back(i);
        return r;
    }
}

#endif