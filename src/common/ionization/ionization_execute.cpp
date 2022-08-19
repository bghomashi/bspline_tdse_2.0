#include "common/ionization/ionization.h"
#include <iostream>
#include <fstream>

void Ionization::Execute() {
    ReadPopulations();
    CalculatePhaseshifts();
    
    for (auto& task : _ionization_tasks)
        task->Execute();
}