#include "common/utility/wigner3j/wigner3j.h"
#include <unordered_map>
#include <cstdint>
#include <cmath>

// #include <iostream>

using namespace Wigner3j;

static std::unordered_map<std::uint64_t, double> _3jDoubles;
static std::unordered_map<std::uint64_t, FractionMinimal> _3jFractions;

std::uint64_t Wigner3j::hash(std::int8_t j1, std::int8_t j, std::int8_t j2, std::int8_t m1, std::int8_t m, std::int8_t m2) {
    std::uint64_t hash = 0;
    // 
    std::int8_t* ptr = (std::int8_t*)&hash;
    *ptr = m2; ptr += 1;
    *ptr = m;  ptr += 1;
    *ptr = m1; ptr += 1;
    *ptr = j2; ptr += 1;
    *ptr = j;  ptr += 1;
    *ptr = j1;
    return hash;
}

Fraction Wigner3j::ThreeJf(std::int8_t j1, std::int8_t j, std::int8_t j2, 
                 std::int8_t m1, std::int8_t m, std::int8_t m2) {
    if (m1 + m + m2 != 0) return 0;
    if (std::abs(m1) > j1 || 
        std::abs(m) > j || 
        std::abs(m2) > j2 ) return 0;
    if (!(j2 >= std::abs(j1-j) && j1+j >= j)) return 0;

    auto it = _3jFractions.find(hash(j1, j, j2, m1, m, m2));
    if (it != _3jFractions.end())
        return it->second;

    int sign;
    int u = j - j1 + j2;
    int zMax = std::min(j2 - m2, u);
    Fraction w3j;
    Fraction D, d, term;

    D *= Factorial(j1 + j - j2);
    D *= Factorial(j1 - j + j2);
    D *= Factorial(-j1 + j + j2);
    D /= Factorial(j1 + j + j2 + 1);


    d *= Factorial(j2 - m2);
    d *= Factorial(j2 + m2);
    d /= Factorial(j - m);
    d /= Factorial(j + m);
    d /= Factorial(j1 - m - m2);
    d /= Factorial(j1 + m + m2);
    

    Fraction total;
    for (int z = 0; z <= zMax; z++) {
        if (j2+m2-u+z < 0) continue;

        sign = (2*j - j1 - m1 + z) % 2 == 0 ? 1 : -1;

        term.Zero();
        term *= sign;
        term *= Factorial(j + j2  - m - m2 - z);
        term *= Factorial(j1 + m + m2 + z);
        term /= Factorial(z);
        term /= Factorial(j2 - m2 - z);
        term /= Factorial(u - z);
        term /= Factorial(j2 + m2 - u + z);

        total += term;
    }

    w3j *= D;
    w3j *= d;
    w3j *= total;
    w3j *= total;
    // ------------
    w3j.sign = total.sign;
    // ------------
    w3j.Contract();

    // memoize
    _3jDoubles[hash(j1, j, j2, m1, m, m2)] = w3j.Evaluate();
    _3jFractions[hash(j1, j, j2, m1, m, m2)] = w3j;

    return w3j;
}


// double Wigner3j::ThreeJd( std::int8_t j1, std::int8_t j, std::int8_t j2, 
//                 std::int8_t m1, std::int8_t m, std::int8_t m2) {
//     auto it = _3jDoubles.find(hash(j1, j, j2, m1, m, m2));
//     if (it != _3jDoubles.end())
//         return it->second;

//     auto f = ThreeJf(j1, j, j2, m1, m, m2);
//     return f.Evaluate();
// }