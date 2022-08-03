#include "banded_matrix.h"
#include <iostream>

template <typename T>
TBandedMatrix<T>::TBandedMatrix() 
    : _N(0), _numBands(0) {}

template <typename T>
TBandedMatrix<T>::TBandedMatrix(int N, int numBands) 
    : _N(N), _numBands(numBands), _elements(N*numBands) {

    }

template <typename T>
T& TBandedMatrix<T>::operator() (int row, int col) {
    int bandOrder = (_numBands - 1) / 2;
    int index = _numBands*row + (col - row + bandOrder);
    return _elements[index];
}

template <typename T>
const T& TBandedMatrix<T>::operator() (int row, int col) const {
    int bandOrder = (_numBands - 1) / 2;
    int index = _numBands*row + (col - row + bandOrder);
    return _elements[index];
}
template <typename T>
TBandedMatrix<T>& TBandedMatrix<T>::operator*= (const T& o) {
    for (auto& e : _elements)
        e *= o;
    return *this;
}
template <typename T>
void TBandedMatrix<T>::Print() const {
    int bandOrder = (_numBands - 1) / 2;

    for (int i = 0; i < _N; i++) {
        for (int j = 0; j < _N; j++) {
            int diff = std::abs(i-j);
            if (diff <= bandOrder)
                std::cout << (*this)(i,j) << " ";
            else
                std::cout << "0.0" << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template <typename T>
void TBandedMatrix<T>::Scale(const T& o) {
    for (auto& e : _elements)
        e *= o;
}

template class TBandedMatrix<std::complex<double>>;
template class TBandedMatrix<double>;