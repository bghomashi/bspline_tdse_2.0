#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <vector>
#include "common/maths/math_common.h"

// Establishes the interface of any Vector implementation to 
// be used by the TDSE/TISE methods and/or to compute any
// observables.

namespace maths {
    class IVector {
    protected:
        int _len;
    public:
        // retrieve the length of this vector
        inline int Length() const {
            return _len;
        }

        // get the value at some index
        virtual complex Get(int index) const = 0;
        // set the value at some index
        virtual void Set(int index, complex value) = 0;
        // scale all elements in the vector
        virtual void Scale(complex a) = 0;
        // make this vector have the same structure (length) 
        // as another vector. Does not copy values
        virtual void Duplicate(const Vector& o) = 0;
        // only copy all values. Does not create structure
        virtual void Copy(const Vector& o) = 0;
        // Zero all elements of the vector
        virtual void Zero() = 0;
        // Combine a set of vectors into one vector
        virtual void Concatenate(const std::vector<Vector>& vecs) = 0;
        // copy all elements of this vector into a list 
        virtual void CopyTo(std::vector<complex>& values) = 0;
        // apply a transformation (function) to all the elements of 
        // the vector and write them into out.
        virtual void Transform(Vector& out, std::function<std::vector<complex>(const std::vector<complex>&)> f) = 0;

        // retrieve a portion of this vector 
        virtual Vector GetSubVector(int start, int end) = 0;
        // release the portion of the vector retrieved with GetSubVector
        virtual void RestoreSubVector(Vector sub) = 0; 

        // multiply (vector) this vector by (scalar) a then add (vector) X
        // results replaces this
        virtual void AYPX(complex a, Vector X) = 0;
        // multiply (vector) X by (scalar) a then add (vector) this
        // results replaces this
        virtual void AXPY(complex a, Vector X) = 0;
        // that the scalar product of two vectors. result is stored in
        // (scalar) value.
        virtual void Dot(const Vector b, complex& value) const = 0;
        // get the max element of the vector
        virtual double Max() const = 0;
        // get the max element of the vector
        virtual double Min() const = 0;
        // get the max element of the vector
        virtual void Abs() = 0;


        // basically only for PETSc...
        virtual void AssembleBegin() {};
        virtual void AssembleEnd() {};
    };
}

#endif