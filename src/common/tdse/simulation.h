#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "common/maths/vector.h"
#include "common/maths/math_factory.h"
#include "common/file_io/io_hdf5.h"

// #include "common/tdse/parameters.h"
#include "common/objects/pulse.h"
#include "common/objects/observable.h"
#include "common/objects/propagator.h"
#include "common/objects/initial_state.h"
#include "common/objects/potential.h"

namespace tdse {
    class Simulation {
        void WriteInitialState() const;
        void WriteFinalState() const;
        void Finish();

        void WriteParametersToTDSE() const;
        bool LoadLastCheckpoint();
        void LoadInitialState();
        void InitializePsi();
        void ComputeDuration();
        void ComputeFields();
        bool CompareTDSEH5wInput() const;


        
        bool ValidateBasis(const nlohmann::json& basis) const;
        void LoadBasis(const nlohmann::json& basis);


        void AddPulse(Pulse::Ptr_t p);
        void AddObservable(Observable::Ptr_t o);
        void AddPotential(Potential::Ptr_t p);



        void DoTimeSteps();
        void DoCheckpoint(int it);
        void DoObservables(int it);

    protected:
        //Parameters p;               // parameters from input file

        // list of pulses, potential, observables
        std::vector<Pulse::Ptr_t> _pulses;
        std::vector<Potential::Ptr_t> _potentials;
        std::vector<Observable::Ptr_t> _observables;
        bool _pol[maths::DimIndex::NUM];       // quick access if there is/is not polarization in x,y,z

        std::vector<double> _field[maths::DimIndex::NUM];

        // numerical setup information 
        int _dof;                       // total degrees of freedom
        std::vector<int> _Ms;           // all possible M values
        std::vector<int> _mRows;        // starting row for each M
        
        double _tmax;                   // time domain
        int _NT;                        // number of time steps
        int _checkpoints;               // how often to checkpoint
        bool _restarting, _do_propagate;// restarting from checkpoint?
        int _potentialSymmetry;         // most general symmetry of potential terms
        bool _laserAxial;               // laser symmetry
        
        InitialState _initialState;     // starting state
        Propagator::Ptr_t _propagator;

        maths::Vector _psi;
        int _startIteration;

        io::HDF5 _tdse_out;             // output file


        void _Execute();
        bool _Validate(const nlohmann::json& input) const;
        void _Load(const nlohmann::json& input);
        void CheckSymmetry();
    public:
        typedef std::shared_ptr<Simulation> Ptr_t;
        Simulation();


        static void Execute();

        // ----- getter functions to access the parameters -----
        static int GetDOF();
        static const std::vector<int>& GetMs();
        static const std::vector<int>& GetMRows(); 
        static double GetTmax();
        static double GetTimestep();
        static double GetTime();
        static int GetNumTimeSteps();
        static const bool* GetPolarization();
        static const std::vector<double>& GetField(int dim_index);
        static const std::vector<double>* GetFields();
        static const std::vector<Pulse::Ptr_t>& GetPulses();
        static const std::vector<Potential::Ptr_t>& GetPotentials();
        static maths::Vector& GetPsi();
        static std::string GetInitialStateFile();
        static int GetInitialStateNmax();
        static int GetPotentialSymmetry();
        static bool GetLaserSymmetry();


        // load/validate the TDSE from the input file
        static void Load(const nlohmann::json& input);
        static bool Validate(const nlohmann::json& input);
    };
}


#endif