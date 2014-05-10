!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE mltfftsg_tools
  USE fft_kinds,                       ONLY: dp

  PRIVATE
  PUBLIC :: matmov, zgetmo, get_cache_size, scaled

CONTAINS

! *****************************************************************************
  SUBROUTINE matmov ( n, m, a, lda, b, ldb )
    INTEGER                                  :: n, m, lda
    COMPLEX(dp)                              :: a( lda, * )
    INTEGER                                  :: ldb
    COMPLEX(dp)                              :: b( ldb, * )

    b ( 1:n , 1:m ) = a ( 1:n, 1:m )
  END SUBROUTINE matmov

! *****************************************************************************
  SUBROUTINE zgetmo ( a, lda, m, n, b, ldb )
    INTEGER                                  :: lda, m, n
    COMPLEX(dp)                              :: a( lda, n )
    INTEGER                                  :: ldb
    COMPLEX(dp)                              :: b( ldb, m )

    b ( 1:n, 1:m ) = TRANSPOSE ( a ( 1:m, 1:n ) )
  END SUBROUTINE zgetmo

! *****************************************************************************
  SUBROUTINE scaled ( n, sc, a )
    INTEGER                                  :: n
    REAL(dp)                                 :: sc
    COMPLEX(dp)                              :: a( n )

    CALL dscal(n,sc,a,1)

  END SUBROUTINE scaled

! *****************************************************************************
  FUNCTION get_cache_size ( )

    IMPLICIT NONE
    INTEGER :: get_cache_size

#if defined ( __T3E )
    get_cache_size = 1024*8
#elif defined ( __SX5 ) || defined ( __T90 )
    get_cache_size = 1024*128
#elif defined ( __ALPHA )
    get_cache_size = 1024*8
#elif defined ( __SGI )
    get_cache_size = 1024*4
#elif defined ( __POWER2 )
    get_cache_size = 1024*10
#elif defined ( __HP )
    get_cache_size = 1024*64
#else
    get_cache_size = 1024*2
#endif

  END FUNCTION get_cache_size

END MODULE mltfftsg_tools
