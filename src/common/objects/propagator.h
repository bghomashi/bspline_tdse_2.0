#ifndef __PROPAGATOR_H__
#define __PROPAGATOR_H__

#include "common/utility/factory.h"
#include <memory>

namespace tdse {
    class Propagator : public utl::Factory<Propagator> {
    protected:
        double _dt;
        double _t;
    public:
        typedef std::shared_ptr<Propagator> Ptr_t;

        Propagator();

        double GetLastTime() const;
        double GetTimeStep() const;

        virtual void Initialize() = 0;
        virtual bool DoStep(int it) = 0;
        virtual void Finish() = 0;
    };
}


#endif