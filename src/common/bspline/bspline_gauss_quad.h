#ifndef __BSPLINE_GAUSS_QUAD_H__
#define __BSPLINE_GAUSS_QUAD_H__

// class:       bspline::GaussQuadrature
// Description: Encapsulates the quadrature integration of a function over 
//              an interval.

#include "common/bspline/bspline_common.h"
#include <vector>
#include <functional>

namespace bspline {
	class GaussQuadrature {											// how many points

		std::vector<double> _points, _weights;
	public:
		bool Initialize();
		maths::complex Integrate(maths::complex xmin, maths::complex xmax, std::function<maths::complex(maths::complex)> f) const;
	};
}

#endif