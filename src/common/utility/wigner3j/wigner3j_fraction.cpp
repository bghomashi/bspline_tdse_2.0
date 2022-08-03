#include "common/utility/wigner3j/wigner3j.h"

#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>

using namespace Wigner3j;

SCALAR FractionMinimal::Evaluate() const {
    if (num == 0) return 0;
        
    return sign * (double)num / (double)denom;
}
FractionMinimal::operator Fraction() const {
    Fraction o;
    o.sign = sign;
    o.num.push_back(num);
    o.denom.push_back(denom);
    return o;
}

void Fraction::Zero() {
    sign = 1;
    num.clear();
    denom.clear();
}
void Fraction::Simplify() {
    // std::cout << "Simplify" << std::endl;
    std::sort(num.begin(), num.end(), std::less<int>());
    std::sort(denom.begin(), denom.end(), std::less<int>());

    // for each numerator is there a gcg
    // std::cout << "before: "; Print();
    for (auto& n : num) {
        if (n < 0) {
            sign *= -1;
            n *= -1;
        }
        if (n == 1) continue;
        
        for (auto& d : denom) {
            if (d < 0) {
                sign *= -1;
                d *= -1;
            }
            if (d == 1) continue;

            int fac;
            while (std::abs(fac = gcf(std::max(n, d), std::min(n, d))) != 1) {
                // std::cout << n << " " << d << " : " << fac << std::endl;

                n /= fac;
                d /= fac;
            }
            if (n == 1) break;
        }
    }
    while (num.size() > 1) {
        if (num[0] != 1) break;
        if (num[1] == 1) num.erase(num.begin());
        else break;
    }
    
    denom.erase(std::remove(denom.begin(), denom.end(), 1), denom.end());       // no 1's in the denominator

    // std::cout << "after: "; Print();
}
SCALAR Fraction::Evaluate() {
    if (num.size() == 0 && denom.size() == 0)
        return 0;
        
    Simplify();
    
    SCALAR total = 1;
    int minCommonFactors = std::min(num.size(), denom.size());
    int iFactor = 0;
    
    // first try to cancel large factors as much as possible
    for (iFactor = 0; iFactor < minCommonFactors; iFactor++)
        total *= (SCALAR)num[iFactor] / (SCALAR)denom[iFactor];

    // if there are remaining numerator factors
    for (int i = iFactor; i < num.size(); i++)
        total *= (SCALAR)num[i];
    // if there are remaining denom factors
    for (int i = iFactor; i < denom.size(); i++)
        total /= (SCALAR)denom[i];

    return total;
}
int Fraction::Num() const {
    SCALAR total = 1.;
    for (int i = 0; i < num.size(); i++)
        total *= (SCALAR)num[i];
    return total;
}
int Fraction::Demom() const {
    SCALAR total = 1.;
    for (int i = 0; i < denom.size(); i++)
        total /= (SCALAR)denom[i];
    return total;
}

Fraction::Fraction() : sign(1) {};
Fraction::Fraction(int i) : sign(i/std::abs(i)) {
    num.push_back(std::abs(i));
}
Fraction& Fraction::operator*= (int i) {
    sign *= (i < 0 ? -1 : 1);
    i *= (i < 0 ? -1 : 1);
    num.push_back(i);
    Simplify();
    return *this;
}
Fraction& Fraction::operator/= (int i) {
    sign *= (i < 0 ? -1 : 1);
    i *= (i < 0 ? -1 : 1);
    denom.push_back(i);
    Simplify();
    return *this;
}
Fraction& Fraction::operator*= (const Fraction& o) {
    sign *= o.sign;
    num.insert(num.end(), o.num.begin(), o.num.end());
    denom.insert(denom.end(), o.denom.begin(), o.denom.end());
    Simplify();
    
    return *this;
}
Fraction& Fraction::operator/= (const Fraction& o) {
    sign *= o.sign;
    denom.insert(denom.end(), o.num.begin(), o.num.end());
    num.insert(num.end(), o.denom.begin(), o.denom.end());
    Simplify();

    return *this;
}
// Fraction& operator=(const Fraction& o) {
//     sign = o.sign;
//     num = o.num;
//     denom = o.denom;
//     return *this;
// }
Fraction Fraction::operator* (int i) const {
    Fraction a = *this;
    a *= i;
    return a;
}
Fraction Fraction::operator/ (int i) const {
    Fraction a = *this;
    a /= i;
    return a;
}
Fraction Fraction::operator* (const Fraction& o) const {
    Fraction a = *this;
    a *= o;
    return a;
}
Fraction Fraction::operator/ (const Fraction& o) const {
    Fraction a = *this;
    a /= o;
    return a;
}
Fraction& Fraction::operator+= (const Fraction& b) {
    Fraction& a = *this;
    Fraction temp_a = *this;
    Fraction temp_b = b;

    if (b.num.empty()) {           // b = 0 so do nothing
        return *this;
    }
    if (a.num.empty()) {           // a = 0 so just copy b into a
        temp_a.sign = temp_b.sign;
        temp_a.denom.insert(temp_a.denom.end(), b.denom.begin(), b.denom.end());
        temp_a.num.insert(temp_a.num.end(), b.num.begin(), b.num.end());
    } else {
        // multiply 'a' (num and denom) by the denominator of 'b'. dont simplify
        temp_a.denom.insert(temp_a.denom.end(), b.denom.begin(), b.denom.end());
        temp_a.num.insert(temp_a.num.end(), b.denom.begin(), b.denom.end());
        // multiply 'b' (num and denom) by the denominator of 'a'. dont simplify
        temp_b.denom.insert(temp_b.denom.end(), a.denom.begin(), a.denom.end());
        temp_b.num.insert(temp_b.num.end(), a.denom.begin(), a.denom.end());
    
        // contract
        temp_a.Contract();
        temp_b.Contract();

        // both fractions should have the same denominator
        assert((temp_a.denom.empty() == temp_b.denom.empty()) && "both do not have denominators");
        assert(temp_a.denom.empty() || (temp_a.denom[0] == temp_b.denom[0]) && "Not the same denominator in addition");
    
        // combine the numerator
        // std::cout << temp_a.sign << "*" << temp_a.num[0] <<" + "<<temp_b.sign<<"*"<<temp_b.num[0] << std::endl;;
        temp_a.num[0] = temp_a.sign*temp_a.num[0] + temp_b.sign*temp_b.num[0];
        temp_a.sign = (temp_a.num[0] < 0 ? -1 : 1);
        temp_a.num[0] *= temp_a.sign;
    }

    // simplify if possible
    temp_a.Simplify();

    a = temp_a;

    return *this;
}
Fraction& Fraction::operator-= (const Fraction& b) {
    Fraction minusB = b*(-1);
    return ((*this) += minusB);
}
Fraction Fraction::operator+ (const Fraction& b) const {
    Fraction a = *this;
    a += b;
    return a;
}
Fraction Fraction::operator- (const Fraction& b) const {
    Fraction a = *this;
    a += b;
    return a;
}
void Fraction::Contract() {
    int tnum = 1, tdenom = 1;

    for (int i = 0; i < num.size(); i++)
        tnum *= num[i];
    num.clear();
    num.push_back(tnum);

    if (!denom.empty()) {
        for (int i = 0; i < denom.size(); i++)
            tdenom *= denom[i];
        denom.clear();
        denom.push_back(tdenom);
    }
}
void Fraction::Print() const {
    if (num.empty()) {
        std::cout << "0 ";
    } else {
        std::cout << (sign<0 ? "-" : "");
        for (int i = 0; i < num.size(); i++)
            std::cout << num[i] << " ";
    }
    if (!denom.empty()) {
        std::cout << "/ ";
        for (int i = 0; i < denom.size(); i++)
            std::cout << denom[i] << " ";
    }
    std::cout << std::endl;
}
Fraction::operator FractionMinimal() const {
    FractionMinimal o;
    o.sign = sign;
    if (num.empty())
        o.num = 0;
    else
        o.num = num[0];
    if (denom.empty())
        o.denom = 1;
    else
        o.denom = denom[0];
    return o;
}


Fraction Wigner3j::Factorial(int N) {
    Fraction a;
    if (N == 0 || N == 1)
        a.num.push_back(1);
    for (int i = 2; i <= N; i++)
        a.num.push_back(i);
        
    return a;
}