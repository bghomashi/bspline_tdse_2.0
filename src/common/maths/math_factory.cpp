#include "common/maths/math_factory.h"

using namespace maths;

static std::shared_ptr<IFactory> s_instance = nullptr;

bool Factory::Startup(int argc, char **args) {
    return s_instance->Startup(argc, args);
}
void Factory::Shutdown() {
    s_instance->Shutdown();
    s_instance = nullptr;
}

Vector Factory::CreateVector(int N) {
    return s_instance->CreateVector(N);
}
void Factory::DestroyVector(Vector& m) {
    s_instance->DestroyVector(m);
}

Matrix Factory::CreateMatrix(int rows, int cols, int numBands) {
    return s_instance->CreateMatrix(rows, cols, numBands);
}
void Factory::DestroyMatrix(Matrix& m) {
    s_instance->DestroyMatrix(m);
}

GMRESSolver Factory::CreateGMRESSolver(int restart_iter, int max_iter) {
    return s_instance->CreateGMRESSolver(restart_iter, max_iter);
}
void Factory::DestroyGMRESSolver(GMRESSolver& m) {
    s_instance->DestroyGMRESSolver(m);
}

EigenSolver Factory::CreateEigenSolver() {
    return s_instance->CreateEigenSolver();
}
void Factory::DestroyEigenSolver(EigenSolver& m) {
    s_instance->DestroyEigenSolver(m);
}

void Factory::SetInstance(IFactory* factory) {
    s_instance = std::shared_ptr<IFactory>(factory);
}