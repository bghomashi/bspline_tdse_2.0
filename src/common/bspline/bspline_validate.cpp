#include "common/bspline/bspline.h"
#include "common/utility/logger.h"

using namespace bspline;

bool BSpline::Validate(const nlohmann::json& basis) const {
    // has order?
    if (!(basis.contains("order") && basis["order"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: order");
        return false;
    }
    // has node sequence?
    if (!(basis.contains("node_sequence") && basis["node_sequence"].is_string())) {
        LOG_CRITICAL("basis must contain string entry: node_sequence");
        return false;
    } 

    if (!bspline::Sequence::Validate(basis["node_sequence"], basis)) {
        LOG_CRITICAL("Failed to contruct \"basis:node_sequence\".");
        return false;
    }
        
    if (!(basis.contains("num_nodes") && basis["num_nodes"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: num_nodes");
        return false;
    }
    if (!(basis.contains("x_min") && basis["x_min"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: x_min");
        return false;
    }
    if (!(basis.contains("x_max") && basis["x_max"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: x_max");
        return false;
    }
    if (basis.contains("ecs_r0") && (!basis["ecs_r0"].is_number() || 
            basis["ecs_r0"].get<double>() < 0 || basis["ecs_r0"].get<double>() > 1.0)) {
        LOG_CRITICAL("Optional entry \"ecs_r0\" in \"basis\" must be a number in the domain [0,1.0]");
        return false;
    }
    if (basis.contains("ecs_theta") && (!basis["ecs_theta"].is_number() || 
            basis["ecs_theta"].get<double>() < -M_PI || basis["ecs_theta"].get<double>() > M_PI)) {
        LOG_CRITICAL("Optional entry \"ecs_theta\" in \"basis\" must be a number in the domain [-Pi,Pi]");
        return false;
    }
    return true;
}
void BSpline::Load(const nlohmann::json& basis, bool ecs_off) {
    int order = basis["order"];
    int nodes = basis["num_nodes"];
    double xmin = basis["x_min"];
    double xmax = basis["x_max"];
    int lmax = basis["lmax"];
    int mmax = basis["mmax"];
    auto seq = bspline::Sequence::Create(basis["node_sequence"], basis);
    double ecs_r0 = 0.9;
    double ecs_theta = maths::Pi/4.;

    if (ecs_off) {
        ecs_r0 = 1.0;
        ecs_theta = 0.0;
    } else {
        // optional - ecs parameters
        if (basis.contains("ecs_r0")) ecs_r0 = basis["ecs_r0"];
        if (basis.contains("ecs_theta")) ecs_theta = basis["ecs_theta"];
    }
    Initialize(order, nodes, 
                      xmin, xmax, 
                      seq, 
                      bspline::ECS{ecs_r0, ecs_theta});
    setSkipFirst();     // enforce 0 boundary at r=0
    setSkipLast();      // enforce 0 boundary at r=xmax
}