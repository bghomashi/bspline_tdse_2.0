#pragma once

// class:       Profiler
// description: A simple class/interface to profile the application. Profiler
//              will essential keep a table of pushed values. When values are
//              "pushed" a counter is incremented and a timer will continue to
//              tick. "pop" the value to stop the timer.
//
//              This is to be used as a record for how long is spent in and how 
//              many times we visit a section of code.
//
//              ** Do not use the class directly. Instead use the static API 
//              defined at the end of the file **

#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include <memory>
#include "common/utility/timer.h"

// For profiling time spent in a particular function
#define ProfilerPush() Profile::Push(__FUNCTION__)
#define ProfilerPop() Profile::Pop(__FUNCTION__)


class Profiler {
public:
    typedef std::shared_ptr<Profiler> Ptr_t;
    
    // push/pop a value onto the table
    virtual void Push(const std::string& name);
    virtual void Pop(const std::string& name);

    // dump the entire table to a file.
    virtual bool PrintTo(const std::string& filename);
    // dump to standard output.
    virtual void Print();
};

// static API for the profiler. If a derived class is implemented set the
// profile to an instance of the derived class:
//
// Profile::SetProfiler(new DerivedProfiler());

namespace Profile {
    void Push(const std::string& name);
    void Pop(const std::string& name);

    bool PrintTo(const std::string& filename);
    void Print();
    
    void SetProfiler(Profiler* profiler);
};
