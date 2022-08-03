#include "common/bspline/bspline.h"
#include "common/utility/json.hpp"
#include "common/utility/logger.h"

using namespace bspline;

inline static  std::string ToLower(std::string data) {
    std::transform(data.begin(), data.end(), data.begin(),
        [](unsigned char c){ return std::tolower(c); });
    return data;
}

// Sequence::Ptr_t Sequence::Create(const nlohmann::json& input) {
//     auto seq = input["node_sequence"].get<std::string>();
    
//     if (ToLower(input["node_sequence"]) == "linear") {
//         return Sequence::Ptr_t(LinearSequence::Create(input));
//     } else if (ToLower(input["node_sequence"]) == "exponential") {
//         return Sequence::Ptr_t(ExponentialSequence::Create(input));
//     } else if (ToLower(input["node_sequence"]) == "sinlike") {
//         return Sequence::Ptr_t(SinlikeSequence::Create(input));
//     } else if (ToLower(input["node_sequence"]) == "parabolic") {
//         return Sequence::Ptr_t(ParabolicLinearSequence::Create(input));
//     }

//     LOG_CRITICAL("node_sequence not supported!");
//     return nullptr;
// }
// bool Sequence::Validate(const nlohmann::json& input) {
//     if (!input.contains("node_sequence") || !input["node_sequence"].is_string()) {
//         LOG_CRITICAL("basis must contain string entry: node_sequence");
//         return false;
//     } 
//     if (ToLower(input["node_sequence"]) == "linear") {
//         return LinearSequence::Validate(input);
//     } else if (ToLower(input["node_sequence"]) == "exponential") {
//         return ExponentialSequence::Validate(input);
//     } else if (ToLower(input["node_sequence"]) == "sinlike") {
//         return SinlikeSequence::Validate(input);
//     } else if (ToLower(input["node_sequence"]) == "parabolic") {
//         return ParabolicLinearSequence::Validate(input);
//     } 
//     LOG_CRITICAL("node_sequence not supported!");
//     return false;
// }

Sequence::Ptr_t LinearSequence::Create(const nlohmann::json& input) {
    return Sequence::Ptr_t(new LinearSequence());
}
Sequence::Ptr_t ExponentialSequence::Create(const nlohmann::json& input) {
    auto seq = new ExponentialSequence();
    seq->g =input["parameter"].get<double>(); 
    return Sequence::Ptr_t(seq);
}
Sequence::Ptr_t ParabolicLinearSequence::Create(const nlohmann::json& input) {
    auto seq = new ParabolicLinearSequence();
    seq->x0 =input["parameter"].get<double>(); 
    return Sequence::Ptr_t(seq);
}
Sequence::Ptr_t SinlikeSequence::Create(const nlohmann::json& input) {
    auto seq = new SinlikeSequence();
    seq->a = input["parameter"].get<double>(); 
    return Sequence::Ptr_t(seq);
}

bool LinearSequence::Validate(const nlohmann::json& input) {
    return true;
}
bool ExponentialSequence::Validate(const nlohmann::json& input) {
    if (!(input.contains("parameter") && input["parameter"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: parameter");
        return false;
    }
    return true;
}
bool ParabolicLinearSequence::Validate(const nlohmann::json& input) {
    if (!(input.contains("parameter") && input["parameter"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: parameter");
        return false;
    }
    return true;
}
bool SinlikeSequence::Validate(const nlohmann::json& input) {
    if (!(input.contains("parameter") && input["parameter"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: parameter");
        return false;
    }
    return true;}

std::string LinearSequence::GetName() {
    return "linear";
}
std::string ExponentialSequence::GetName() {
    return "exponential";
}
std::string ParabolicLinearSequence::GetName() {
    return "parabolic";
}
std::string SinlikeSequence::GetName() {
    return "sinlike";
}

std::vector<double> LinearSequence::GetGrid(double xmin, double xmax, int nodes) {
    std::vector<double> grid(nodes);
    
    for (int i = 0; i < nodes; i++)
        grid[i] = xmin + ((xmax - xmin) * (double)i) / double(nodes - 1.);

    return grid;
}
std::vector<double> ExponentialSequence::GetGrid(double xmin, double xmax, int nodes) {
    std::vector<double> grid(nodes);
    for (int i = 0; i < nodes; i++)
        grid[i] = xmin + (xmax - xmin)*(exp(g*i/(nodes-1)) - 1.)/(exp(g) - 1.);
    return grid;
}
std::vector<double> ParabolicLinearSequence::GetGrid(double xmin, double xmax, int nodes) {
    std::vector<double> grid(nodes);

    double i0 = std::floor(2.*(nodes-1) / (1.+(xmax - xmin)/(x0 - xmin)));
    // double i0 = std::min(30, nodes);						// parameter
    double a0 = xmax / i0 /(2.*(nodes-1) - i0);
    double a1 = 2.*xmax /(2.*(nodes-1) - i0);
    double a2 = -xmax*i0 /(2.*(nodes-1) - i0);

    for (int i = 0; i < nodes; i++) {
        if (i < i0) 
            grid[i] = xmin + a0*i*i;
        else
            grid[i] = xmin + a2 + a1*i;
    }

    return grid;
}
std::vector<double> SinlikeSequence::GetGrid(double xmin, double xmax, int nodes) {
    std::vector<double> grid(nodes);
    for (int i = 0; i < nodes; i++)
        grid[i] = xmin + xmax*sin(M_PI/2. * pow(double(i)/double(nodes-1), a));
    return grid;
}