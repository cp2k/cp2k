module mltfftsg_tools

  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND ( 14, 200 )

  private
  public matmov, zgetmo, get_cache_size

contains
  
  SUBROUTINE matmov ( n, m, a, lda, b, ldb )
    IMPLICIT NONE
    INTEGER :: n, m, lda, ldb
    COMPLEX (dbl) :: a ( lda, * ), b ( ldb, * )
    b ( 1:n , 1:m ) = a ( 1:n, 1:m )
  END SUBROUTINE matmov
  
  SUBROUTINE zgetmo ( a, lda, m, n, b, ldb )
    IMPLICIT NONE
    INTEGER  :: lda, m, n, ldb
    COMPLEX(dbl) :: a ( lda, n ), b ( ldb, m )
    b ( 1:n, 1:m ) = TRANSPOSE ( a ( 1:m, 1:n ) )
  END SUBROUTINE zgetmo
  
  FUNCTION get_cache_size ( ) RESULT ( ncache )
    IMPLICIT NONE
    INTEGER ncache
#if defined ( __T3E )
    ncache = 1024*8
#elif defined ( __SX5 ) || defined ( __T90 )
    ncache = 1024*128
#elif defined ( __ALPHA )
    ncache = 1024*8
#elif defined ( __SGI )
    ncache = 1024*4
#elif defined ( __POWER2 )
    ncache = 1024*10
#elif defined ( __HP )
    ncache = 1024*64
#else
    ncache = 1024*2
#endif
    
  END FUNCTION get_cache_size

END MODULE mltfftsg_tools
