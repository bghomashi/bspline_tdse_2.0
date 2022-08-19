#ifndef __BSPLINE_H__
#define __BSPLINE_H__

// struct:		bspline::Basis
// description:	Complete a utility-wrapper struct. Maintains no state and only
//				defined the static interface. The allows a single *constant* 
//				(except for Load) global instance of BSpline to be accessed from 
//				anywhere. This should also be more efficent in the case of more 
//				than one dimension using the same basis; it does not need to be 
//				instantiated more than once.
//				** this interface should be preferred
// usage:		bspline::Basis::Load(input);
//				bspline::Basis::Integrate(...);


// class:		bspline::BSpline
// description:	This class describes a BSpline basis on a set of nodes. 
//				Methods are provide to retrieve the value of a single bspline
//				and its n'th derivative at a point and on an interval.
//				Methods are also provided to evaluate a function (and its derivative)
//				expanded in the BSpline basis.
//				Methods are provided to evaluate integrals of functions over the 
//				BSpline basis. 

#include "common/bspline/bspline_common.h"
#include "common/bspline/bspline_gauss_quad.h"
#include "common/bspline/bspline_sequence.h"
#include "common/bspline/bspline_ecs.h"
#include "common/utility/json.hpp"
#include <functional>
#include <vector>

namespace bspline {
	struct Basis {
		static bool Validate(const nlohmann::json& basis);
		static void Load(const nlohmann::json& basis, bool ecs_off = false);

		// get the bs'th Bspline (derivative order dn)
		static maths::complex bspline(maths::complex x, int bs, int dn = 0);		// x - includes ecs
		static maths::complex bspline(double x, int bs, int dn = 0);				// x - does not include ecs
		static std::vector<maths::complex> GetBSpline(const std::vector<double>& x, int bs, int dn = 0);

		// evaluate function expanded in Bspline basis at x - dn is the order of the derivative
		// x-before ecs rotation
		static maths::complex FunctionEvaluate(double x, const std::vector<maths::complex>& fc, int dn = 0);
		static std::vector<maths::complex> FunctionEvaluate(const std::vector<double>& x, const std::vector<maths::complex>& fc, int dn = 0);

		// integrate a function over two Bsplines over the whole grid
		// the parameter passed to 'f' is (complex) x
		static maths::complex Integrate(int bs1, int bs2, int dn1 = 0, int dn2 = 0);
		static maths::complex Integrate(int bs1, int bs2, std::function<maths::complex(maths::complex)> f, int dn1 = 0, int dn2 = 0);
		
		// integrate a function over restricted bounds
		static maths::complex Integrate(double xmin, double xmax, int bs1, int bs2, std::function<maths::complex(maths::complex)> f, int dn1 = 0, int dn2 = 0);

		static const std::vector<double>& GetGrid();
		static int GetNumBSplines();
		static int GetOrder();
		static int GetNumNodes();
		static double GetXmin();
		static double GetXmax();
		static double GetECS_R0();
		static double GetECS_Theta();
	};


    class BSpline {
		// ---- structures to hold the polynomial coefficients ----
		// in a table-like structure 
		struct BSplineCoeff {
			std::vector<std::vector<maths::complex>> d;

			BSplineCoeff(int max_order);
			maths::complex operator() (int interval, int degree) const;
			maths::complex& operator() (int interval, int degree);
		};
		struct CoeffMatrix {
			std::vector<std::vector<BSplineCoeff>> d;
			CoeffMatrix(size_t order, size_t num_knots);
			BSplineCoeff& operator() (size_t order, size_t bspline);
		};


		// ----- member variables -----
        GaussQuadrature _glQuad;					// integration quadrature
        ECS _ecs;									// ecs grid calculator

		std::vector<double> _grid;                  // real locations of the nodes
		std::vector<maths::complex> _ecs_grid, _knots;     //maths::complex nodes, knots
		
		int _nodes;                                 // number of nodes
		int _order;                                 // order of bsplines
		int _numBSplines;                           // number of bsplines
		bool _skipFirst, _skipLast;                 // assert 0 boundary by excluding first/last bspline

		std::vector<BSplineCoeff> _bsCoeffs;        // polynomial coefficients for each interval
		std::vector<double> _partialFactorial, _factorial;// factorial storage
        


		// ----- internal initialization functions ------
		void InitializeFactorials();
		void InitializeBSplines();
		void InitializeBSpline(const std::vector<maths::complex>& knots, BSplineCoeff& out);
		void InitializeBSplineOfOrder(const std::vector<maths::complex>& knots, int order, CoeffMatrix& coeff);
		void PropagateCoeffients(const std::vector<maths::complex>& knots, int tbs, int order, CoeffMatrix& coeff);
    public:
		int Initialize(int order, int nodes, double xmin, double xmax, Sequence::Ptr_t seq, const ECS& ecs = ECS{0.9, maths::Pi/4.});
		void setSkipFirst(bool flag = true);
		void setSkipLast(bool flag = true);
		bool Validate(const nlohmann::json& basis) const;
		void Load(const nlohmann::json& basis, bool ecs_off = false);

		// get the bs'th Bspline (derivative dn)
		maths::complex bspline(double x, int bs, int dn = 0) const;                // x-before ecs rotation
		maths::complex bspline(maths::complex x, int bs, int dn = 0) const;               // x-after ecs rotation
		std::vector<maths::complex> getBSpline(const std::vector<double>& x, int bs, int dn = 0) const;// x-before ecs rotation

		// evaluate function expanded in Bspline basis at x - dn is the order of the derivative
        // x-before ecs rotation
		maths::complex FunctionEvaluate(double x, const std::vector<maths::complex>& fc, int dn = 0) const;
		std::vector<maths::complex> FunctionEvaluate(const std::vector<double>& x, const std::vector<maths::complex>& fc, int dn = 0) const;

		// integrate a function over two Bsplines over the whole grid
        // the parameter passed to 'f' is (complex) x
		maths::complex Integrate(int bs1, int bs2, int dn1 = 0, int dn2 = 0) const;
		maths::complex Integrate(int bs1, int bs2, std::function<maths::complex(maths::complex)> f, int dn1 = 0, int dn2 = 0) const;
		
        // integrate a function over restricted bounds
		maths::complex Integrate(double xmin, double xmax, int bs1, int bs2, std::function<maths::complex(maths::complex)> f, int dn1 = 0, int dn2 = 0) const;

		// helper functions
		const std::vector<double>& getGrid() const;
		int getNumBSplines() const;
		int getOrder() const;
		int getNumNodes() const;
		int whichInterval(double x) const;
		double getXmin() const;
		double getXmax() const;
		double getECS_R0() const;
		double getECS_Theta() const;
    };   




	
	struct LinearSequence : Sequence::Register<LinearSequence> {
		std::vector<double> GetGrid(double xmin, double xmax, int nodes);

		static Sequence::Ptr_t Create(const nlohmann::json& input);
		static bool Validate(const nlohmann::json& input);
		static std::string GetName();
	};
	struct ExponentialSequence : Sequence::Register<ExponentialSequence> {
		double g;
		
		static Sequence::Ptr_t Create(const nlohmann::json& input);
		static bool Validate(const nlohmann::json& input);
		static std::string GetName();

		std::vector<double> GetGrid(double xmin, double xmax, int nodes);
	};
	struct SinlikeSequence : Sequence::Register<SinlikeSequence> {
		double a;
		static Sequence::Ptr_t Create(const nlohmann::json& input);
		static bool Validate(const nlohmann::json& input);
		static std::string GetName();

		std::vector<double> GetGrid(double xmin, double xmax, int nodes);
	};
	struct ParabolicLinearSequence : Sequence::Register<ParabolicLinearSequence> {
		double x0;

		static Sequence::Ptr_t Create(const nlohmann::json& input);
		static bool Validate(const nlohmann::json& input);
		static std::string GetName();

		std::vector<double> GetGrid(double xmin, double xmax, int nodes);
	};
}

#endif