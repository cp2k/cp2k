#include "openPMD.h"

#include <openPMD/openPMD.hpp>

namespace
{
auto access_c_to_cxx(openPMD_Access_create access) -> openPMD::Access
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
} // namespace

extern "C"
{
    int openPMD_Series_create(
        char const *filename,
        openPMD_Access access,
        openPMD_Series *series_param)
    {
        auto series = reinterpret_cast<openPMD::Series **>(series_param);
        auto cxx_access = access_c_to_cxx(access);
        try
        {
            *series = new openPMD::Series(filename, cxx_access);
        }
        catch (...)
        {
            delete *series;
            return 1;
        }
        return 0;
    }

    int openPMD_Series_close(openPMD_Series series_param)
    {
        auto series = reinterpret_cast<openPMD::Series *>(series_param);
        series->close();
        delete series;
        return 0;
    }
}
