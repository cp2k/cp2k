! Basic use statements and preprocessor macros
! should be included in the use statements

  USE base_hooks,                      ONLY: cp__a,&
                                             cp__b,&
                                             cp__w,&
                                             cp__h,&
                                             cp__l,&
                                             cp_abort,&
                                             cp_warn,&
                                             cp_hint,&
                                             timeset,&
                                             timestop

#if defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)
#define __OFFLOAD
#endif

! Check for OpenMP early on - ideally before the compiler fails with a cryptic message.
#if !defined(_OPENMP)
   "OpenMP is required. Please add the corresponding flag (eg. -fopenmp for GFortran) to your Fortran compiler flags."
#endif

! Dangerous: Full path can be arbitrarily long and might overflow Fortran line.
#if !defined(__SHORT_FILE__)
#define __SHORT_FILE__ __FILE__
#endif

#define __LOCATION__ cp__l(__SHORT_FILE__,__LINE__)
#define CPWARN(msg) CALL cp__w(__SHORT_FILE__,__LINE__,msg)
#define CPABORT(msg) CALL cp__b(__SHORT_FILE__,__LINE__,msg)
! In contrast to CPWARN, the warning counter is not increased
#define CPHINT(msg) CALL cp__h(__SHORT_FILE__,__LINE__,msg)

# define CPASSERT(cond) IF(.NOT.(cond))CALL cp__a(__SHORT_FILE__,__LINE__)

! The MARK_USED macro can be used to mark an argument/variable as used. It is intended to make
! it possible to switch on -Werror=unused-dummy-argument, but deal elegantly with, e.g.,
! library wrapper routines that take arguments only used if the library is linked in.
! This code should be valid for any Fortran variable, is always standard conforming,
! and will be optimized away completely by the compiler
#define MARK_USED(foo) IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(foo))==-1) EXIT ;  END DO ; ENDIF

! Calculate version number from 2 or 3 components. Can be used for comparison, e.g.,
! CPVERSION3(4, 9, 0) <= CPVERSION3(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__)
! CPVERSION(8, 0) <= CPVERSION(__GNUC__, __GNUC_MINOR__)
#define CPVERSION2(MAJOR, MINOR) ((MAJOR) * 10000 + (MINOR) * 100)
#define CPVERSION3(MAJOR, MINOR, UPDATE) (CPVERSION2(MAJOR, MINOR) + (UPDATE))
#define CPVERSION CPVERSION2

! Usage: CPVERSION_CHECK(9, 5, >, __GNUC__, __GNUC_MINOR__)
! On top of CPVERSION macro, test if MAJOR_TEST and MINOR_TEST are defined.
! Perform actual comparison according to COMP argument.
#define CPVERSION_CHECK(MAJOR_BASE, MINOR_BASE, COMP, MAJOR_TEST, MINOR_TEST) ( \
  (defined(MAJOR_TEST) && defined(MINOR_TEST)) && \
  (CPVERSION2(MAJOR_BASE, MINOR_BASE) COMP CPVERSION2(MAJOR_TEST, MINOR_TEST)))

! Avoid to default initialize type-components (default c'tor)
#if CPVERSION_CHECK(9, 5, >, __GNUC__, __GNUC_MINOR__) || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#define FTN_NO_DEFAULT_INIT
#endif

! gfortran before 8.3 complains about internal symbols not being specified in
! any data clause when using DEFAULT(NONE) and OOP procedures are called from
! within the parallel region.
#if CPVERSION_CHECK(8, 3, >, __GNUC__, __GNUC_MINOR__)
#define OMP_DEFAULT_NONE_WITH_OOP SHARED
#else
#define OMP_DEFAULT_NONE_WITH_OOP NONE
#endif
