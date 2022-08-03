#include "common/bspline/bspline.h"

using namespace bspline;

int BSpline::getNumBSplines() const {
    int num = _numBSplines;
    if (_skipFirst) num--;
    if (_skipLast) num--;
    return num;
}
int BSpline::getOrder() const {
    return _order;
}
int BSpline::whichInterval(double x) const {
    for (int i = 0; i < (int)_grid.size(); i++) {
        if (_grid[i] <= x && x <= _grid[i + 1])
            return i;
    }
    return -1;
}

const std::vector<double>& BSpline::getGrid() const {
    return _grid;
}
void BSpline::setSkipFirst(bool flag) {
    _skipFirst = flag;
}
void BSpline::setSkipLast(bool flag) {
    _skipLast = flag;
}

double BSpline::getXmin() const {
    return _grid.front();
}
double BSpline::getXmax() const {
    return _grid.back();
}
double BSpline::getECS_R0() const {
    return _ecs.r0;
}
double BSpline::getECS_Theta() const {
    return _ecs.theta;
}
int BSpline::getNumNodes() const {
    return _nodes;
}