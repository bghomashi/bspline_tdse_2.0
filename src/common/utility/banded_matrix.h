#ifndef __BANDED_MATRIX_H__
#define __BANDED_MATRIX_H__

#include <vector>
#include <complex>

template <typename T>
class TBandedMatrix {
    std::vector<T> _elements;
    int _numBands;
    int _N;
public:
    TBandedMatrix();
    TBandedMatrix(int N, int numBands);
    T& operator() (int row, int col);
    const T& operator() (int row, int col) const;
    TBandedMatrix& operator*= (const T& o);
    void Scale(const T& o);
    void Print() const;
};

typedef TBandedMatrix<std::complex<double>> BandedMatrix;
typedef TBandedMatrix<double> BandedMatrixD;

#endif