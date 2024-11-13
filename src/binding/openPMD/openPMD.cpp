#include "openPMD.h"

#include <openPMD/openPMD.hpp>

extern "C"
{
    int openPMD_Series_create(
        char const *filename, openPMD_Access access, openPMD_Series *series_param)
    {
        auto series = reinterpret_cast<openPMD::Series**>(series_param);
        auto cxx_access = [&]() {
            switch (access)
            {
            case openPMD_Access_create:
                return openPMD::Access::CREATE;
            case openPMD_Access_read_only:
                return openPMD::Access::READ_ONLY;
            }
            // unreachable
            return static_cast<openPMD::Access>(0);
        }();
        try
        {
            *series =
                new openPMD::Series(filename, cxx_access);
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
        auto series = reinterpret_cast<openPMD::Series*>(series_param);
        series->close();
        delete series;
        return 0;
    }
}
