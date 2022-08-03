#pragma once

#include <cstdint>
#include <vector>
#include <gsl/gsl_sf_coupling.h>

namespace Wigner3j {
    typedef double SCALAR;
    struct Fraction;

    Fraction ThreeJf(std::int8_t j1, std::int8_t j, std::int8_t j2, std::int8_t m1, std::int8_t m, std::int8_t m2);
    // SCALAR ThreeJd(std::int8_t j1, std::int8_t j, std::int8_t j2, std::int8_t m1, std::int8_t m, std::int8_t m2);
    inline SCALAR ThreeJd(std::int8_t j1, std::int8_t j, std::int8_t j2, std::int8_t m1, std::int8_t m, std::int8_t m2) {
        return gsl_sf_coupling_3j(2*j1, 2*j, 2*j2, 2*m1, 2*m, 2*m2);
    }

    struct FractionMinimal {
        int sign, num, denom;

        SCALAR Evaluate() const;
        operator Fraction() const;
    };

    struct Fraction {
        int sign;
        std::vector<int> num, denom;

        Fraction();
        Fraction(int i);
        Fraction& operator*= (int i);
        Fraction& operator/= (int i);
        Fraction& operator*= (const Fraction& o);
        Fraction& operator/= (const Fraction& o);
        Fraction operator* (int i) const;
        Fraction operator/ (int i) const;
        Fraction operator* (const Fraction& o) const;
        Fraction operator/ (const Fraction& o) const;
        Fraction& operator+= (const Fraction& b);
        Fraction& operator-= (const Fraction& b);
        Fraction operator+ (const Fraction& b) const;
        Fraction operator- (const Fraction& b) const;
        operator FractionMinimal() const;

        void Zero();
        void Simplify();
        SCALAR Evaluate();
        int Num() const;
        int Demom() const;
        void Contract();
        void Print() const;
    };

    Fraction Factorial(int N);
    int gcf(int a, int b);
    std::uint64_t hash(std::int8_t j1, std::int8_t j, std::int8_t j2, std::int8_t m1, std::int8_t m, std::int8_t m2);
}