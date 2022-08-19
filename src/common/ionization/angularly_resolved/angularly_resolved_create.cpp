#include "common/ionization/angularly_resolved/angularly_resolved.h"

IonizationTask::Ptr_t AngularlyResolved::Create(const nlohmann::json& angularly_resolved) {
    AngularlyResolved* new_ar = new AngularlyResolved();

    auto& energy = angularly_resolved["energy"];
    auto& x = angularly_resolved["x"];
    auto& z = angularly_resolved["z"];
    new_ar->_output_filename = angularly_resolved["filename"];
    new_ar->_emin = energy["min"];
    new_ar->_emax = energy["max"];
    new_ar->_estep = energy["step"];
    new_ar->_xmin = x["min"];
    new_ar->_xmax = x["max"];
    new_ar->_xstep = x["step"];
    new_ar->_zmin = z["min"];
    new_ar->_zmax = z["max"];
    new_ar->_zstep = z["step"];

    return IonizationTask::Ptr_t(new_ar);
}
