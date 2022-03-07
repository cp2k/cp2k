/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef OFFLOAD_LIBRARY_H
#define OFFLOAD_LIBRARY_H

#if defined(__OFFLOAD_GRID) || defined(__OFFLOAD_PW) || defined(__OFFLOAD_DBM)

/* Check that __OFFLOAD_CUDA or __OFFLOAD_HIP are given. breaks the compilation
 * if not */

#if !defined(__OFFLOAD_CUDA) && !defined(__OFFLOAD_HIP)
#error                                                                         \
    "GPU support is activated for modules supporting it. To continue compilation please add -D__OFFLOAD_CUDA or -D__OFFLOAD_HIP to DFLAGS"
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif
/*******************************************************************************
 * \brief Returns the number of available devices.
 * \author Ole Schuett
 ******************************************************************************/
int offload_get_device_count(void);

/*******************************************************************************
 * \brief Selects the device to be used.
 * \author Ole Schuett
 ******************************************************************************/
void offload_set_device_id(int device_id);

/*******************************************************************************
 * \brief Returns the device to be used.
 * \author Ole Schuett
 ******************************************************************************/
int offload_get_device_id(void);

/*******************************************************************************
 * \brief Activates the device selected via offload_set_device_id()
 * \author Ole Schuett
 ******************************************************************************/
void offload_set_device(void);

/*******************************************************************************
 * \brief Starts a timing range.
 * \author Ole Schuett
 ******************************************************************************/
void offload_timeset(const char *message);

/*******************************************************************************
 * \brief Ends a timing range.
 * \author Ole Schuett
 ******************************************************************************/
void offload_timestop(void);

/*******************************************************************************
 * \brief Gets free and total device memory.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mem_info(size_t *free, size_t *total);

/*******************************************************************************
 * \brief Allocate pinned memory (or simple malloc when there is no gpu)
 ******************************************************************************/
int offload_host_malloc(void **ptr__, const size_t size__);

/*******************************************************************************
 * \brief free pinned memory (or simple free when there is no gpu)
 ******************************************************************************/
int offload_host_free(void *ptr__);

#ifdef __cplusplus
}
#endif
#endif
// EOF
