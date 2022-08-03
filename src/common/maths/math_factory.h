#ifndef __MATH_ENGINE_H__
#define __MATH_ENGINE_H__

#include "common/maths/math_common.h"
#include <memory>

namespace maths {
    class IFactory {
    public:
        virtual bool Startup(int argc, char **args) = 0;
        virtual void Shutdown() = 0;

        virtual Vector CreateVector(int N) = 0;
        virtual void DestroyVector(Vector& m) = 0;

        virtual Matrix CreateMatrix(int rows, int cols, int numBands) = 0;
        virtual void DestroyMatrix(Matrix& m) = 0;

        virtual GMRESSolver CreateGMRESSolver(int restart_iter = 500, int max_iter = 10000) = 0;
        virtual void DestroyGMRESSolver(GMRESSolver& m) = 0;

        virtual EigenSolver CreateEigenSolver() = 0;
        virtual void DestroyEigenSolver(EigenSolver& m) = 0;
    };

    class Factory {
    public:
        static bool Startup(int argc, char **args);
        static void Shutdown();

        static Vector CreateVector(int N);
        static void DestroyVector(Vector& m);

        static Matrix CreateMatrix(int rows, int cols, int numBands);
        static void DestroyMatrix(Matrix& m);

        static GMRESSolver CreateGMRESSolver(int restart_iter = 500, int max_iter = 10000);
        static void DestroyGMRESSolver(GMRESSolver& m);

        static EigenSolver CreateEigenSolver();
        static void DestroyEigenSolver(EigenSolver& m);

        static void SetInstance(IFactory* factory);
    };
}

#endif