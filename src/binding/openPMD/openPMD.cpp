#include "openPMD.h"

#include <openPMD/openPMD.hpp>

namespace implementation
{
namespace
{
    auto access_c_to_cxx(openPMD_Access access) -> openPMD::Access
    {
        switch (access)
        {
        case openPMD_Access_create:
            return openPMD::Access::CREATE;
        case openPMD_Access_read_only:
            return openPMD::Access::READ_ONLY;
        }
        // unreachable
        return static_cast<openPMD::Access>(0);
    }

    template <typename... Args>
    auto Series_create(openPMD_Series *series_param, Args &&...constructor_args)
        -> int
    {
        auto series = reinterpret_cast<openPMD::Series **>(series_param);
        try
        {
            *series =
                new openPMD::Series(std::forward<Args>(constructor_args)...);
        }
        catch (...)
        {
            delete *series;
            return 1;
        }
        return 0;
    }
} // namespace
} // namespace implementation

extern "C"
{
    int openPMD_Series_create(
        char const *filename,
        openPMD_Access access,
        openPMD_Series *series_param)
    {
        return implementation::Series_create(
            series_param, filename, implementation::access_c_to_cxx(access));
    }

    int openPMD_Series_create_mpi(
        char const *filename,
        openPMD_Access access,
        MPI_Comm comm,
        openPMD_Series *series_param)
    {
        return implementation::Series_create(
            series_param,
            filename,
            implementation::access_c_to_cxx(access),
            comm);
    }

    int openPMD_Series_close(openPMD_Series series_param)
    {
        auto series = reinterpret_cast<openPMD::Series *>(series_param);
        series->close();
        delete series;
        return 0;
    }
}
