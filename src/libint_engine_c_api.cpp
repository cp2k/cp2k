/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "libint_engine_c_api.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>

#if defined(__LIBINT)
#include <libint2/engine.h>
#include <libint2/initialize.h>
#endif

namespace {

// Keep exceptions and error storage on the C++ side of the ABI boundary.
template <typename Function>
int guarded_call(char *error, size_t capacity, Function &&function) noexcept {
  if (error && capacity)
    error[0] = '\0';
  try {
    function();
    return 0;
  } catch (const std::exception &exception) {
    if (error && capacity) {
      const auto length = std::min(capacity - 1, std::strlen(exception.what()));
      std::memcpy(error, exception.what(), length);
      error[length] = '\0';
    }
  } catch (...) {
    if (error && capacity) {
      const char message[] = "Unknown Libint Engine exception";
      const auto length = std::min(capacity - 1, sizeof(message) - 1);
      std::memcpy(error, message, length);
      error[length] = '\0';
    }
  }
  return 1;
}

#if defined(__LIBINT)
using libint2::BraKet;
using libint2::Engine;
using libint2::Operator;
using libint2::Shell;

static_assert(std::is_same<libint2::value_type, double>::value,
              "The CP2K Engine bridge requires double-precision Libint");

struct EngineContext {
  Engine engine;
  Engine::compute2_ptr_type compute = nullptr;
  std::array<Shell, 4> shells;
  std::array<bool, 4> ready{};
  int ncenters = 0;
  int ncomponents = 0;
};

EngineContext &get_context(void *handle) {
  if (!handle)
    throw std::invalid_argument("NULL Libint Engine handle");
  return *static_cast<EngineContext *>(handle);
}

template <Operator op, typename Parameters>
std::unique_ptr<EngineContext>
make_context(int ncenters, int derivative_order, int max_nprim, int max_l,
             double precision, const Parameters &parameters) {
  const std::array<BraKet, 3> brakets{BraKet::xs_xs, BraKet::xs_xx,
                                      BraKet::xx_xx};
  // Select once, avoiding generic compute dispatch and shell copies per
  // quartet.
  const Engine::compute2_ptr_type computes[3][2] = {
      {&Engine::compute2<op, BraKet::xs_xs, 0>,
       &Engine::compute2<op, BraKet::xs_xs, 1>},
      {&Engine::compute2<op, BraKet::xs_xx, 0>,
       &Engine::compute2<op, BraKet::xs_xx, 1>},
      {&Engine::compute2<op, BraKet::xx_xx, 0>,
       &Engine::compute2<op, BraKet::xx_xx, 1>}};
  auto context = std::make_unique<EngineContext>();
  context->engine = Engine(op, max_nprim, max_l, derivative_order, precision,
                           parameters, brakets[ncenters - 2]);
  context->compute = computes[ncenters - 2][derivative_order];
  context->ncenters = ncenters;
  context->ncomponents = derivative_order == 0 ? 1 : 3 * ncenters;
  for (int i = 0; i < ncenters; ++i) {
    auto &shell = context->shells[i];
    shell.alpha.reserve(max_nprim);
    shell.contr.resize(1);
    shell.contr[0].coeff.reserve(max_nprim);
    shell.max_ln_coeff.reserve(max_nprim);
  }
  return context;
}
#endif

} // namespace

extern "C" int cp2k_libint_engine_initialize(char *error,
                                             size_t error_capacity) {
  return guarded_call(error, error_capacity, [] {
#if defined(__LIBINT)
    libint2::initialize();
#else
    throw std::runtime_error("CP2K was built without Libint");
#endif
  });
}

extern "C" int cp2k_libint_engine_finalize(char *error, size_t error_capacity) {
  return guarded_call(error, error_capacity, [] {
#if defined(__LIBINT)
    libint2::finalize();
#else
    throw std::runtime_error("CP2K was built without Libint");
#endif
  });
}

extern "C" int cp2k_libint_engine_create(int operator_id, int ncenters,
                                         int derivative_order, int max_nprim,
                                         int max_l, const double parameters[3],
                                         double precision, void **handle,
                                         char *error, size_t error_capacity) {
  return guarded_call(error, error_capacity, [&] {
    if (!handle || *handle)
      throw std::invalid_argument("Engine creation requires a NULL handle");
#if defined(__LIBINT)
    if (!libint2::initialized())
      throw std::logic_error(
          "Libint must be initialized before creating Engines");
    if (ncenters < 2 || ncenters > 4 || derivative_order < 0 ||
        derivative_order > 1 || max_nprim < 1 || max_l < 0 ||
        !std::isfinite(precision) || precision < 0 || !parameters)
      throw std::invalid_argument("Invalid Libint Engine configuration");
    for (int i = 0; i < 3; ++i)
      if (!std::isfinite(parameters[i]))
        throw std::invalid_argument("Non-finite Libint operator parameter");
    if (parameters[0] < 0)
      throw std::invalid_argument("Negative Libint operator parameter");

    std::unique_ptr<EngineContext> context;
    switch (operator_id) {
    case CP2K_LIBINT_ENGINE_COULOMB:
      context = make_context<Operator::coulomb>(
          ncenters, derivative_order, max_nprim, max_l, precision,
          libint2::operator_traits<Operator::coulomb>::default_params());
      break;
    case CP2K_LIBINT_ENGINE_ERF_COULOMB:
      context = make_context<Operator::erf_coulomb>(ncenters, derivative_order,
                                                    max_nprim, max_l, precision,
                                                    parameters[0]);
      break;
    case CP2K_LIBINT_ENGINE_ERFC_COULOMB:
      context = make_context<Operator::erfc_coulomb>(ncenters, derivative_order,
                                                     max_nprim, max_l,
                                                     precision, parameters[0]);
      break;
    case CP2K_LIBINT_ENGINE_TRUNCATED_COULOMB:
      context = make_context<Operator::truncated_coulomb>(
          ncenters, derivative_order, max_nprim, max_l, precision,
          parameters[0]);
      break;
    case CP2K_LIBINT_ENGINE_GAUSSIAN:
      context = make_context<Operator::cgtg>(
          ncenters, derivative_order, max_nprim, max_l, precision,
          libint2::ContractedGaussianGeminal{{parameters[0], 1.0}});
      break;
    case CP2K_LIBINT_ENGINE_DELTA:
      context = make_context<Operator::delta>(
          ncenters, derivative_order, max_nprim, max_l, precision,
          libint2::operator_traits<Operator::delta>::default_params());
      break;
    case CP2K_LIBINT_ENGINE_ERFX_COULOMB:
      context = make_context<Operator::erfx_coulomb>(
          ncenters, derivative_order, max_nprim, max_l, precision,
          std::array<double, 3>{parameters[0], parameters[1], parameters[2]});
      break;
    default:
      throw std::invalid_argument("Unsupported Libint Engine operator");
    }
    *handle = context.release();
#else
    (void)operator_id;
    (void)ncenters;
    (void)derivative_order;
    (void)max_nprim;
    (void)max_l;
    (void)parameters;
    (void)precision;
    throw std::runtime_error("CP2K was built without Libint");
#endif
  });
}

extern "C" int cp2k_libint_engine_set_shell(void *handle, int slot, int l,
                                            int nprim, const double *exponents,
                                            const double *coefficients,
                                            const double center[3], char *error,
                                            size_t error_capacity) {
  return guarded_call(error, error_capacity, [&] {
#if defined(__LIBINT)
    auto &context = get_context(handle);
    if (slot < 1 || slot > context.ncenters || l < 0 ||
        static_cast<size_t>(l) > context.engine.max_l() || nprim < 1 ||
        static_cast<size_t>(nprim) > context.engine.max_nprim() || !exponents ||
        !coefficients || !center)
      throw std::invalid_argument("Invalid Libint shell dimensions or arrays");
    context.ready[slot - 1] = false;
    for (int i = 0; i < nprim; ++i)
      if (!std::isfinite(exponents[i]) || exponents[i] < 0 ||
          !std::isfinite(coefficients[i]))
        throw std::invalid_argument(
            "Invalid Libint shell exponent/coefficient");
    for (int i = 0; i < 3; ++i)
      if (!std::isfinite(center[i]))
        throw std::invalid_argument("Non-finite Libint shell center");

    auto &shell = context.shells[slot - 1];
    shell.alpha.assign(exponents, exponents + nprim);
    shell.O = {center[0], center[1], center[2]};
    shell.contr[0].l = l;
    shell.contr[0].pure = false;
    shell.contr[0].coeff.assign(coefficients, coefficients + nprim);
    // Shell's screening cache must match the raw, normalization-free
    // coefficients.
    shell.max_ln_coeff.resize(nprim);
    for (int i = 0; i < nprim; ++i)
      shell.max_ln_coeff[i] = coefficients[i] == 0
                                  ? -std::numeric_limits<double>::max()
                                  : std::log(std::abs(coefficients[i]));
    context.ready[slot - 1] = true;
#else
    (void)handle;
    (void)slot;
    (void)l;
    (void)nprim;
    (void)exponents;
    (void)coefficients;
    (void)center;
    throw std::runtime_error("CP2K was built without Libint");
#endif
  });
}

extern "C" int cp2k_libint_engine_compute(void *handle,
                                          const double *const **targets,
                                          size_t *ncart, int *ncomponents,
                                          char *error, size_t error_capacity) {
  if (targets)
    *targets = nullptr;
  if (ncart)
    *ncart = 0;
  if (ncomponents)
    *ncomponents = 0;
  return guarded_call(error, error_capacity, [&] {
    if (!targets || !ncart || !ncomponents)
      throw std::invalid_argument("NULL Libint Engine result descriptor");
#if defined(__LIBINT)
    auto &context = get_context(handle);
    size_t size = 1;
    for (int i = 0; i < context.ncenters; ++i) {
      if (!context.ready[i])
        throw std::logic_error(
            "Set every physical shell before Engine compute");
      const auto shell_size = context.shells[i].cartesian_size();
      if (size > std::numeric_limits<size_t>::max() / shell_size)
        throw std::overflow_error("Libint shell quartet size overflow");
      size *= shell_size;
    }
    const auto &unit = Shell::unit();
    const auto &a = context.shells[0];
    const auto &b = context.ncenters == 4 ? context.shells[1] : unit;
    const auto &c = context.shells[context.ncenters == 4 ? 2 : 1];
    const auto &d =
        context.ncenters == 2 ? unit : context.shells[context.ncenters - 1];
    const auto &result =
        (context.engine.*context.compute)(a, b, c, d, nullptr, nullptr);
    if (result.size() != static_cast<size_t>(context.ncomponents))
      throw std::runtime_error("Unexpected Libint Engine component count");
    *ncart = size;
    *ncomponents = context.ncomponents;
    *targets = result[0] ? result.data() : nullptr;
#else
    (void)handle;
    throw std::runtime_error("CP2K was built without Libint");
#endif
  });
}

extern "C" void cp2k_libint_engine_destroy(void **handle) {
  if (!handle)
    return;
#if defined(__LIBINT)
  delete static_cast<EngineContext *>(*handle);
#endif
  *handle = nullptr;
}
