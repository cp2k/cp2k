#ifndef __MA_TOPOLOGY_H
#define __MA_TOPOLOGY_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#if defined __SEASTAR || __GEMINI

int extract_topology();
int remove_topology();
#endif
#endif
