// #ifndef __PETSC_PROFILER_H__
// #define __PETSC_PROFILER_H__


// // class:       PetscProfiler
// // description: Implements the profiler interface. The purpose of this class
// //              is to restrict access to the profiler resource to only the rank 0
// //              process.


// #include "common/utility/profiler.h"
// #include <string>

// class PetscProfiler : public Profiler {
// public:
//     void Push(const std::string& name);
//     void Pop(const std::string& name);
//     void Print();
//     bool PrintTo(const std::string& log_file);
// };

// #endif