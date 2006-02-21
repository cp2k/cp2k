MODULE mltfftsg_tools

#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dbl = KIND(0.0)
#else
  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND ( 14, 200 )
#endif

  PRIVATE
  PUBLIC :: matmov, zgetmo, get_cache_size

CONTAINS
  
! *****************************************************

  SUBROUTINE matmov ( n, m, a, lda, b, ldb )
    IMPLICIT NONE
    INTEGER :: n, m, lda, ldb
    COMPLEX (dbl) :: a ( lda, * ), b ( ldb, * )
    b ( 1:n , 1:m ) = a ( 1:n, 1:m )
  END SUBROUTINE matmov
  
! *****************************************************

  SUBROUTINE zgetmo ( a, lda, m, n, b, ldb )
    IMPLICIT NONE
    INTEGER  :: lda, m, n, ldb
    COMPLEX(dbl) :: a ( lda, n ), b ( ldb, m )
    b ( 1:n, 1:m ) = TRANSPOSE ( a ( 1:m, 1:n ) )
  END SUBROUTINE zgetmo
  
! *****************************************************

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

! *****************************************************

END MODULE mltfftsg_tools
