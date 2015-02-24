/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>

#include <stdlib.h>
#include <string.h>


int
cuda_error (cudaError_t cudaError)
{
  if (cudaError != cudaSuccess)
    {
      printf ("CUDA Error: %s\n", cudaGetErrorString (cudaError));
      return 1;
    }
  return 0;
};

extern "C" int
ma_get_ndevices_cu (int *n_devices)
{
  cudaError_t errDev;
  errDev = cudaGetDeviceCount (n_devices);

  if (cuda_error (errDev))
    return 1;

  return 0;
}

extern "C" int
ma_get_NUMAnode_cu (int cuDevice, int *NUMAnode)
{
  struct cudaDeviceProp devProp;
  cudaError_t errDev;
  FILE *gpuFile;
  char tmpFile[256], fileName[256];

  errDev = cudaGetDeviceProperties (&devProp, cuDevice);
  sprintf (fileName, "/sys/bus/pci/devices/%04x:%02x:%02x.0/",
	   devProp.pciDomainID, devProp.pciBusID, devProp.pciDeviceID);
  strcpy (tmpFile, fileName);
  strcat (tmpFile, "numa_node\0");
  gpuFile = fopen (tmpFile, "r");
  fscanf (gpuFile, "%d", NUMAnode);
  fclose (gpuFile);

  if (cuda_error (errDev))
    return 1;

  return 0;
}

extern "C" int
ma_get_cores_cu (int cuDevice, char *cores)
{
  struct cudaDeviceProp devProp;
  cudaError_t errDev;
  FILE *gpuFile;
  char tmpFile[256], fileName[256];

  errDev = cudaGetDeviceProperties (&devProp, cuDevice);
  sprintf (fileName, "/sys/bus/pci/devices/%04x:%02x:%02x.0/",
	   devProp.pciDomainID, devProp.pciBusID, devProp.pciDeviceID);
  strcpy (tmpFile, fileName);
  strcat (tmpFile, "local_cpulist\0");
  gpuFile = fopen (tmpFile, "r");
  fscanf (gpuFile, "%s", &cores);
  fclose (gpuFile);

  if (cuda_error (errDev))
    return 1;

  return 0;
}

extern "C" int
ma_get_closerDev_cu (int nodeId, int *devIds)
{
  struct cudaDeviceProp devProp;
  cudaError_t errDev;

  int i, j, devCount, gpuNUMAnode, *devDistances;
  FILE *gpuFile;
  char tmpFile[256], fileName[256];

  ma_get_ndevices_cu (&devCount);

  devDistances = (int *) malloc (devCount * sizeof (int));

  for (i = 0; i < devCount; ++i)
    {
      errDev = cudaGetDeviceProperties (&devProp, i);
      sprintf (fileName, "/sys/bus/pci/devices/%04x:%02x:%02x.0/",
	       devProp.pciDomainID, devProp.pciBusID, devProp.pciDeviceID);
      strcpy (tmpFile, fileName);
      strcat (tmpFile, "numa_node\0");
      gpuFile = fopen (tmpFile, "r");
      fscanf (gpuFile, "%d", &gpuNUMAnode);
      fclose (gpuFile);
      if (gpuNUMAnode == nodeId)
	devDistances[i] = 1;
      else
	devDistances[i] = 2;
    }

//closer devices 
  j = 0;
  for (i = 0; i < devCount; i++)
    if (devDistances[i] == 1)
      {
	devIds[j] = i;
	j++;
      }
//distant devices     
  for (i = 0; i < devCount; i++)
    if (devDistances[i] == 2)
      {
	devIds[j] = i;
	j++;
      }

  if (cuda_error (errDev))
    return 1;

  return 0;
}

extern "C" int
ma_get_uma_closerDev_cu (int *devIds)
{
  int errDev;
  int i, devCount;

  errDev = ma_get_ndevices_cu (&devCount);

  for (i = 0; i < devCount; i++)
    devIds[i] = i;

  return errDev;
}



extern "C" int
ma_get_core_cu (int coreId, int devId)
{
  struct cudaDeviceProp devProp;
  cudaError_t errDev;

  int i, devCount, core;
  FILE *gpuFile;
  char tmpFile[256], fileName[256];

  devId = -1;
  ma_get_ndevices_cu (&devCount);

  for (i = 0; i < devCount && devId == -1; ++i)
    {
      errDev = cudaGetDeviceProperties (&devProp, i);
      sprintf (fileName, "/sys/bus/pci/devices/%04x:%02x:%02x.0/",
	       devProp.pciDomainID, devProp.pciBusID, devProp.pciDeviceID);
      strcpy (tmpFile, fileName);
      strcat (tmpFile, "local_cpulist\0");
      gpuFile = fopen (tmpFile, "r");
      fscanf (gpuFile, "%d", &core);
      fclose (gpuFile);
      if (coreId == core)
	devId = i;
    }

  if (cuda_error (errDev))
    return 1;

  return 0;
}

extern "C" int
ma_get_nDevcu (int nodeId, int *count)
{
  struct cudaDeviceProp devProp;
  cudaError_t errDev;

  int i, j = 0, devCount, gpuNUMAnode;
  FILE *gpuFile;
  char tmpFile[256], fileName[256];

  ma_get_ndevices_cu (&devCount);

  for (i = 0; i < devCount; ++i)
    {
      errDev = cudaGetDeviceProperties (&devProp, i);
      sprintf (fileName, "/sys/bus/pci/devices/%04x:%02x:%02x.0/",
	       devProp.pciDomainID, devProp.pciBusID, devProp.pciDeviceID);
      strcpy (tmpFile, fileName);
      strcat (tmpFile, "numa_node\0");
      gpuFile = fopen (tmpFile, "r");
      fscanf (gpuFile, "%d", &gpuNUMAnode);
      fclose (gpuFile);
      if (gpuNUMAnode == nodeId)
	{
	  j++;
	}
    }

  *count = j;
  if (cuda_error (errDev))
    return 1;

  return 0;

}

extern "C" int
ma_get_cu (int nodeId, int devIds[])
{
  struct cudaDeviceProp devProp;
  cudaError_t errDev;

  int i, j = 0, devCount, gpuNUMAnode;
  FILE *gpuFile;
  char tmpFile[256], fileName[256];

  ma_get_ndevices_cu (&devCount);

  for (i = 0; i < devCount; ++i)
    {
      errDev = cudaGetDeviceProperties (&devProp, i);
      sprintf (fileName, "/sys/bus/pci/devices/%04x:%02x:%02x.0/",
	       devProp.pciDomainID, devProp.pciBusID, devProp.pciDeviceID);
      strcpy (tmpFile, fileName);
      strcat (tmpFile, "numa_node\0");
      gpuFile = fopen (tmpFile, "r");
      fscanf (gpuFile, "%d", &gpuNUMAnode);
      fclose (gpuFile);
      if (gpuNUMAnode == nodeId)
	{
	  devIds[j] = i;
	  j++;
	}
    }

  if (cuda_error (errDev))
    return 1;

  return 0;
}
