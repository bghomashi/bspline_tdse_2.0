#ifndef __OBSERVABLE_H__
#define __OBSERVABLE_H__

#include "common/utility/json.hpp"
#include "common/utility/factory.h"
#include <memory>
#include <string>
#include <unordered_map>

namespace tdse {
    class Observable : public utl::Factory<Observable> {
    protected:
        int _computePeriod;
    public:
        typedef std::shared_ptr<Observable> Ptr_t;

        Observable();
        void DoObservable(int it);

        void SetComputePeriod(int iterations);

        virtual void Flush() {};
        virtual int MemoryAllocated() const;
        virtual void Startup(int it) = 0;
        virtual void Shutdown() = 0;
        virtual void Compute(int it) = 0;
    };
}


#endif