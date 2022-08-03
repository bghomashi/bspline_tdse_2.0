#include "common/utility/wigner3j/wigner3j.h"
#include <cmath>

int Wigner3j::gcf(int a, int b) {
    int t;
    while (b != 0) {
        t = b;
        b = a % b;
        a = t;
    }
    return std::abs(a);
}