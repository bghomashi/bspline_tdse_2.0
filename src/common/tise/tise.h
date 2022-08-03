#ifndef __TISE_H__
#define __TISE_H__

// class:       tise::TISE
// description: encapsulates eigenstate solver.
// usage:       Use the static API

#include "common/maths/vector.h"
#include "common/maths/math_factory.h"
#include "common/file_io/io_hdf5.h"
#include "common/objects/potential.h"

namespace tise {
    class TISE {
    protected:
        void AddPotential(tdse::Potential::Ptr_t p);

        void Solve();
        void SolveCentral();
        void SolveAxial();
        void _Execute();
        bool _Validate(const nlohmann::json& input) const;
        void _Load(const nlohmann::json& input);

        // list of potentials
        std::vector<tdse::Potential::Ptr_t> _potentials;

        // output
        io::HDF5 _tise_out;

        // numerical
        bool _ecs_on;
        double _tol;
        maths::EigenSolver _eigensolver;


        void ComputeAndOutputBoundStates(int nmax, maths::Matrix H0, maths::Matrix S, int l = -1);
        void ComputeAndOutputContinuumStates(int nmax, maths::Matrix H0, maths::Matrix S, int l = -1);
    public:
        typedef std::shared_ptr<TISE> Ptr_t;

        TISE();

        static void Execute();

        // load/validate the TDSE from the input file
        static bool Validate(const nlohmann::json& input);
        static void Load(const nlohmann::json& input);
    };
}


#endif