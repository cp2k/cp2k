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
        catch (std::exception const &e)
        {
            std::cout << "[Series_create] Caught error: '" << e.what() << "'\n";
            delete *series;
            return 1;
        }
        catch (...)
        {
            std::cout << "[Series_create] Caught unknown error.\n";
            delete *series;
            return 1;
        }
        return 0;
    }

    template <
        typename From,
        typename To,
        typename FromOpaque,
        typename ToOpaque>
    int do_upcast(FromOpaque from_param, ToOpaque *to_param)
    {
        auto from = reinterpret_cast<From *>(from_param);
        auto to = reinterpret_cast<To **>(to_param);
        auto from_upcasted = static_cast<To *>(from);
        *to = from_upcasted;
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

    int openPMD_Series_write_Iteration(
        openPMD_Series series_param,
        openPMD_Iteration_Index_t const * index,
        openPMD_Iteration *iteration)
    {
        auto series = reinterpret_cast<openPMD::Series *>(series_param);
        auto &res_iteration =
            series->writeIterations()[openPMD::Iteration::IterationIndex_t(
                *index)];
        *reinterpret_cast<openPMD::Iteration **>(iteration) = &res_iteration;
        return 0;
    }

    int openPMD_Series_upcast_to_attributable(
        // in
        openPMD_Series series,
        // out
        openPMD_Attributable *attr)
    {
        return implementation::
            do_upcast<openPMD::Series, openPMD::Attributable>(series, attr);
    }

    int openPMD_Series_present(
        openPMD_Series series_param
    )
    {
        auto series = reinterpret_cast<openPMD::Series *>(series_param);
        bool res = series && *series;
        return res ? 1 : 0;
    }

    char const *openPMD_get_default_extension()
    {
        constexpr static char const *preferred_order[] = {
            "bp5", "bp4", "bp", "h5", "json"};
        auto const &available_extensions = openPMD::getFileExtensions();
        for (auto const &ext : preferred_order)
        {
            if (std::find(
                    available_extensions.begin(),
                    available_extensions.end(),
                    ext) != available_extensions.end())
            {
                return ext;
            }
        }
        return nullptr;
    }

    int openPMD_attributable_set_attribute_vec_int(
        openPMD_Attributable attr,
        char const *attr_name,
        const int *begin,
        int const *length)
    {
        auto attributable = reinterpret_cast<openPMD::Attributable *>(attr);
        attributable->setAttribute(
            attr_name, std::vector<int>{begin, begin + size_t(*length)});
        return 0;
    }

    int openPMD_Iteration_upcast_to_attributable(
        // in
        openPMD_Iteration iteration,
        // out
        openPMD_Attributable *attr)
    {
        return implementation::
            do_upcast<openPMD::Iteration, openPMD::Attributable>(
                iteration, attr);
    }
}
