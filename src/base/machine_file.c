/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/*******************************************************************************
 * \brief Create one directory, accepting an existing directory (or a symlink
 *        to one), but rejecting other existing filesystem objects.
 * \param path Directory name.
 * \return Zero on success, otherwise a POSIX error number.
 ******************************************************************************/
static int mkdir_one(const char *path) {
  if (mkdir(path, 0777) == 0)
    return 0;
  if (errno != EEXIST)
    return errno;
  struct stat info;
  if (stat(path, &info) != 0)
    return errno;
  return S_ISDIR(info.st_mode) ? 0 : ENOTDIR;
}

/*******************************************************************************
 * \brief Create a directory and missing parents without invoking a shell.
 *        Permissions follow the caller's umask. Paths are interpreted
 * literally.
 * \param path Directory name.
 * \return Zero on success, otherwise a POSIX error number.
 ******************************************************************************/
int cp_mkdir_p(const char *path) {
  if (path == NULL || path[0] == '\0')
    return EINVAL;
  char *buffer = malloc(strlen(path) + 1);
  if (buffer == NULL)
    return ENOMEM;
  strcpy(buffer, path);

  int error = 0;
  for (char *p = buffer + 1; *p != '\0'; p++) {
    if (*p != '/')
      continue;
    *p = '\0';
    error = mkdir_one(buffer);
    *p = '/';
    if (error != 0)
      break;
  }
  if (error == 0)
    error = mkdir_one(buffer);
  free(buffer);
  return error;
}

// EOF
