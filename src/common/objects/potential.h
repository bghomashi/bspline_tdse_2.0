#ifndef __POTENTIAL_H__
#define __POTENTIAL_H__

#include "common/utility/json.hpp"
#include "common/maths/matrix.h"
#include <memory>

#include "common/utility/factory.h"

namespace tdse {
    struct Symmetry {
        enum {
            Central = 3,
            Axial   = 1,
            None = 0,
        };
    };

    class Potential : public utl::Factory<Potential> {
    protected:
        bool _isCentral;
        bool _isAxial;
    public:
        typedef std::shared_ptr<Potential> Ptr_t;
    
        Potential() : _isCentral(true), _isAxial(true) {};


        inline bool isCentral() const {
            return _isCentral;
        }
        inline bool isAxial() const {
            return _isAxial;
        }

        virtual double operator() (double x, double y, double z) const = 0;
        virtual void FillMatrix(maths::Matrix m, int N, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0}) = 0;

        // these are only supported in the axial case right now...
        virtual void FillMatrixGradX(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0}) = 0;
        virtual void FillMatrixGradY(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0}) = 0;
        virtual void FillMatrixGradZ(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0}) = 0;
    };
}


#endif