#ifndef __DIPOLE_ACC_OBS_H__
#define __DIPOLE_ACC_OBS_H__

#include <vector>
#include "common/objects/observable.h"
#include "common/maths/math_common.h"
#include "common/file_io/io_ascii.h"
#include "common/utility/banded_matrix.h"

class DipoleObservable : public tdse::Observable::Register<DipoleObservable> {
    maths::Matrix _x[maths::DimIndex::NUM];
    maths::Vector _psi;            // shortcut to wavefunction
    maths::Vector _psiTemp;       // just from MatMult output

    std::string _outputFilename;
    io::ASCII _txtFile;

    void FillMatrixX(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows);
    void FillMatrixY(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows);
    void FillMatrixZ(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows);
    void BuildR(int N, int order, BandedMatrix& r);

    void FillBlock( int l1, int m1, int l2, int m2,
                    int lmax, int mmax, int N, int order,
                    const std::vector<int>& Ms,
                    const std::vector<int>& mRows,
                    const BandedMatrix& r,
                    std::function<maths::complex(int, int, int, int)> YlmYlm,
                    maths::Matrix m);
public:
    DipoleObservable();

    int MemoryAlloced() const;
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it);

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};

#endif