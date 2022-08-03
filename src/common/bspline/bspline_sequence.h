#ifndef __BSPLINE_SEQUENCE_H__
#define __BSPLINE_SEQUENCE_H__

// class:       Sequence
// description: The sequence class defines an interface for node sequences 
//              that the main BSpline class will use to compute the knots.
//              
//              A derive class must implement the GetGrid function which 
//              returns the grid. 
//              Additionally, the derived classes need to be added to the:
//              static create/validate functions.

#include <vector>
#include <memory>
#include "common/utility/json.hpp"
#include "common/utility/factory.h"

namespace bspline {
	struct Sequence : utl::Factory<Sequence> {
		typedef std::shared_ptr<Sequence> Ptr_t;

		virtual std::vector<double> GetGrid(double xmin, double xmax, int nodes) = 0;
	};
}

#endif