c-----------------------------------------------------------------------------!
c   CP2K: A general program to perform molecular dynamics simulations         !
c   Copyright (C) 2000  CP2K developers group                                 !
c-----------------------------------------------------------------------------!
c
c     This is a fixed format file
c
c------------------------------------------------------------------------------
c
c     Holds External information for Message Passing Systems
c
c------------------------------------------------------------------------------

      MODULE parallel_include
#if defined(__parallel)
c
c     Include file for MPI
c
        INCLUDE 'mpif.h'
#endif
        
        PUBLIC
        
        INTEGER :: nothing
        
      END MODULE parallel_include
