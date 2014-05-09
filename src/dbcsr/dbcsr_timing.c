#include <stdio.h>
#include <stdlib.h>

#include <time.h>

static time_t first_sec = 0;
static time_t first_nano = 0;
double nano = 1.0 / 1.0E9;

clockid_t clk = CLOCK_MONOTONIC;

#ifdef __cplusplus
extern "C" {
#endif
double get_clk (){
    int error;
    struct timespec t;
    if (first_sec == 0)
      {
        error = clock_gettime (clk, &t);
        if (error)
  	return -1.0;
        first_sec = t.tv_sec;
        first_nano = t.tv_nsec;
        return ((double) (t.tv_sec - first_sec) +
  	      (double) ((double) (t.tv_nsec - first_nano) * nano));
      }
    else
      {
        error = clock_gettime (clk, &t);
        if (error)
  	return -1.0;
        return ((double) (t.tv_sec - first_sec) +
  	      (double) ((double) (t.tv_nsec - first_nano) * nano));
      }
  }
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
  double get_rate (){
    int error;
    struct timespec t;
    error = clock_getres (clk, &t);
    if (error)
      return -1.0;
    return (double) 1.0 / (((double) t.tv_sec) + ((double) t.tv_nsec) * nano);
  }
#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
extern "C" {
#endif
  double get_clk_ (){
    return get_clk ();
  }
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
  double get_rate_ (){
    return get_rate ();
  }
#ifdef __cplusplus
}
#endif
