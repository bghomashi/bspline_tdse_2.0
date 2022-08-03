#ifndef __MPI_RANK_ZERO_ONLY_H__
#define __MPI_RANK_ZERO_ONLY_H__

#include <mpi.h>

template <typename T, typename Ret,  typename... Args1, typename... Args2>
inline typename std::enable_if<std::is_same<Ret, void>::value, Ret>::type
RankZeroOnly(T& ref, Ret(T::*func)(Args1...), Args2&&... args) {                
    // this function gets called if the return type of "func" is void
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        (ref.*func)(std::forward<Args2>(args)...);
}

template <typename T, typename Ret,  typename... Args1, typename... Args2>
inline typename std::enable_if<!std::is_same<Ret, void>::value, Ret>::type
RankZeroOnly(T& ref, Ret(T::*func)(Args1...), Args2&&... args) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        return (ref.*func)(std::forward<Args2>(args)...);
    return Ret();
}

template <typename T, typename Ret,  typename... Args1, typename... Args2>
inline typename std::enable_if<!std::is_same<Ret, void>::value, Ret>::type
RankZeroOnly(T& ref, Ret defaultRet, Ret(T::*func)(Args1...), Args2&&... args) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        return (ref.*func)(std::forward<Args2>(args)...);
    return defaultRet;
}






#endif