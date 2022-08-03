#include "objects/potentials/coulomb/coulomb_pot.h"
#include "common/utility/index_manip.h"
#include "common/utility/logger.h"
#include "common/utility/spherical_harmonics.h"
#include "common/bspline/bspline.h"

#include <iostream>
#include <functional>

using namespace tdse;
using namespace maths;


complex CoulombPotential::R_element(int i, int j, int l, double R0) {
		// std::pow(R0, -(complex)(l+1)) *
		// bspline::Basis::Integrate(0., R0, i+1, j+1, [l](complex r){
		// 	return std::pow(r, (double)l);
		// });// << " " <<
		// std::pow(R0, (complex)l) *
		// bspline::Basis::Integrate(R0, Infinity, i+1, j+1, [l](complex r){
		// 	return std::pow(r, -(double)(l+1));
		// });
		//<< std::endl;
// exit(-1);
	return 
		std::pow(R0, -(double)(l+1)) *
		bspline::Basis::Integrate(0., R0, i+1, j+1, [l](complex r){
			return std::pow(r, (double)l);
		}) 
		+
		std::pow(R0, (double)l) *
		bspline::Basis::Integrate(R0, Infinity, i+1, j+1, [l](complex r){
			return std::pow(r, -(double)(l+1));
		});
}

complex CoulombPotential::il1m1_jl2m2(
    int i, int l1, int m1, 
    int j, int l2, int m2, 
    const std::vector<BandedMatrix>& rmatrix) {

	complex total = 0;
    // this can be a different lmax
	// in the axial case is only one/two terms non-zero?
	// TripleIntergal has to satify the triangle inequality
	// |l1-l| <= l2 <= l1+l and
	// l1 + l + l2 must be even
	// so the only term that matters is l = abs(l1-l2)
    for (int l = 0; l <= _expansion_lmax; l++)
	//int l = abs(l1-l2);
	{
		complex term = 0;
		// for (int m = -l; m <= l; m++)
		int m = 0;
		{
			// if (std::abs(std::conj(Ylm(l, m, _theta0, _phi0)) *
			// 	TripleIntegral(l1, l, l2, m1, m, m2) ) != 0) {
			// 	std::cout <<  l1 << " " << l << " " << l2 << " ";
			// 	std::cout <<
			// 	std::conj(Ylm(l, m, _theta0, _phi0)) *
			// 		TripleIntegral(l1, l, l2, m1, m, m2) << std::endl; 
			// }
			
			term += std::conj(Ylm(l, m, _theta0, _phi0)) *
				TripleIntegral(l1, l, l2, m1, m, m2);
		}
		// std::cout << "l: " << l << " ";
		// std::cout << "term: " << term << " ";
		// std::cout << "std::conj(Ylm): " << std::conj(Ylm(l, m, _theta0, _phi0)) << " ";
		// std::cout << "TripleIntegral(l1, l, l2, m1, m, m2): " << TripleIntegral(l1, l, l2, m1, m, m2)<< std::endl;
		term *= rmatrix[l](i, j) * 4.*Pi/(2.*l + 1.);
		total += term;
	}
	
	return -_Z*total;
	
} 
BandedMatrix CoulombPotential::R_matrix(int N, int order, int l, double R0) {
	BandedMatrix rmatrix(N, 2*order-1);
	for (int i = 0; i < N; i++) {
		for (int j = i; j < std::min(N, i+order); j++) {
			rmatrix(i,j) = rmatrix(j, i) = R_element(i, j, l, R0);
		}
	}
	return rmatrix;
}

