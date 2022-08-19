#ifndef __IONIZATION_TASK_H__
#define __IONIZATION_TASK_H__

#include "common/utility/json.hpp"
#include "common/utility/factory.h"
#include <memory>

class Ionization;
class IonizationTask : public utl::Factory<IonizationTask> {
protected:
public:
    typedef std::shared_ptr<IonizationTask> Ptr_t;

    Ionization* _parent;

    virtual void Execute() = 0;

};



#endif