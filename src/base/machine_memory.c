/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#include <stdint.h>

#if defined(__APPLE__)
#include <mach/mach.h>
#include <sys/sysctl.h>
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/*******************************************************************************
 * \brief Physical memory where /proc/meminfo does not exist (macOS), in bytes:
 *        the total, the free pages, and the pages likely to be allocatable
 *        (free, inactive and purgeable, which the kernel reclaims on demand).
 *        All three stay zero where the system offers no such query, which
 *        m_memory_details reads as "not available".
 * \author Teodoro Laino
 ******************************************************************************/
void m_memory_details_c(int64_t *total, int64_t *free_bytes,
                        int64_t *likely_free) {
  *total = 0;
  *free_bytes = 0;
  *likely_free = 0;
#if defined(__APPLE__)
  uint64_t memsize = 0;
  size_t length = sizeof(memsize);
  if (sysctlbyname("hw.memsize", &memsize, &length, NULL, 0) == 0) {
    *total = (int64_t)memsize;
  }
  const mach_port_t host = mach_host_self();
  vm_size_t page_size = 0;
  vm_statistics64_data_t vm;
  mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
  if (host_page_size(host, &page_size) == KERN_SUCCESS &&
      host_statistics64(host, HOST_VM_INFO64, (host_info64_t)&vm, &count) ==
          KERN_SUCCESS) {
    *free_bytes = (int64_t)vm.free_count * (int64_t)page_size;
    *likely_free =
        (int64_t)(vm.free_count + vm.inactive_count + vm.purgeable_count) *
        (int64_t)page_size;
  }
  mach_port_deallocate(mach_task_self(), host);
#endif
}

#if defined(__cplusplus)
}
#endif
