#if !defined(__HAS_smm_dnn)
SUBROUTINE smm_dnn(M,N,K,A,B,C)
   INTEGER :: M,N,K
   REAL(KIND=KIND(0.0D0)) :: A(*), B(*), C(*)
   REAL(KIND=KIND(0.0D0)), PARAMETER :: one=1
   CALL DGEMM("n","n",M,N,K,one,A,M,B,K,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_dnt)
SUBROUTINE smm_dnt(M,N,K,A,B,C)
   INTEGER :: M,N,K
   REAL(KIND=KIND(0.0D0)) :: A(*), B(*), C(*)
   REAL(KIND=KIND(0.0D0)), PARAMETER :: one=1
   CALL DGEMM("n","t",M,N,K,one,A,M,B,N,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_dtn)
SUBROUTINE smm_dtn(M,N,K,A,B,C)
   INTEGER :: M,N,K
   REAL(KIND=KIND(0.0D0)) :: A(*), B(*), C(*)
   REAL(KIND=KIND(0.0D0)), PARAMETER :: one=1
   CALL DGEMM("t","n",M,N,K,one,A,K,B,K,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_dtt)
SUBROUTINE smm_dtt(M,N,K,A,B,C)
   INTEGER :: M,N,K
   REAL(KIND=KIND(0.0D0)) :: A(*), B(*), C(*)
   REAL(KIND=KIND(0.0D0)), PARAMETER :: one=1
   CALL DGEMM("t","t",M,N,K,one,A,K,B,N,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_snn)
SUBROUTINE smm_snn(M,N,K,A,B,C)
   INTEGER :: M,N,K
   REAL(KIND=KIND(0.0)) :: A(*), B(*), C(*)
   REAL(KIND=KIND(0.0)), PARAMETER :: one=1
   CALL SGEMM("n","n",M,N,K,one,A,M,B,K,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_snt)
SUBROUTINE smm_snt(M,N,K,A,B,C)
   INTEGER :: M,N,K
   REAL(KIND=KIND(0.0)) :: A(*), B(*), C(*)
   REAL(KIND=KIND(0.0)), PARAMETER :: one=1
   CALL SGEMM("n","t",M,N,K,one,A,M,B,N,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_stn)
SUBROUTINE smm_stn(M,N,K,A,B,C)
   INTEGER :: M,N,K
   REAL(KIND=KIND(0.0)) :: A(*), B(*), C(*)
   REAL(KIND=KIND(0.0)), PARAMETER :: one=1
   CALL SGEMM("t","n",M,N,K,one,A,K,B,K,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_stt)
SUBROUTINE smm_stt(M,N,K,A,B,C)
   INTEGER :: M,N,K
   REAL(KIND=KIND(0.0)) :: A(*), B(*), C(*)
   REAL(KIND=KIND(0.0)), PARAMETER :: one=1
   CALL SGEMM("t","t",M,N,K,one,A,K,B,N,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_znn)
SUBROUTINE smm_znn(M,N,K,A,B,C)
   INTEGER :: M,N,K
   COMPLEX(KIND=KIND(0.0D0)) :: A(*), B(*), C(*)
   COMPLEX(KIND=KIND(0.0D0)), PARAMETER :: one=1
   CALL ZGEMM("n","n",M,N,K,one,A,M,B,K,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_znt)
SUBROUTINE smm_znt(M,N,K,A,B,C)
   INTEGER :: M,N,K
   COMPLEX(KIND=KIND(0.0D0)) :: A(*), B(*), C(*)
   COMPLEX(KIND=KIND(0.0D0)), PARAMETER :: one=1
   CALL ZGEMM("n","t",M,N,K,one,A,M,B,N,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_ztn)
SUBROUTINE smm_ztn(M,N,K,A,B,C)
   INTEGER :: M,N,K
   COMPLEX(KIND=KIND(0.0D0)) :: A(*), B(*), C(*)
   COMPLEX(KIND=KIND(0.0D0)), PARAMETER :: one=1
   CALL ZGEMM("t","n",M,N,K,one,A,K,B,K,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_ztt)
SUBROUTINE smm_ztt(M,N,K,A,B,C)
   INTEGER :: M,N,K
   COMPLEX(KIND=KIND(0.0D0)) :: A(*), B(*), C(*)
   COMPLEX(KIND=KIND(0.0D0)), PARAMETER :: one=1
   CALL ZGEMM("t","t",M,N,K,one,A,K,B,N,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_cnn)
SUBROUTINE smm_cnn(M,N,K,A,B,C)
   INTEGER :: M,N,K
   COMPLEX(KIND=KIND(0.0)) :: A(*), B(*), C(*)
   COMPLEX(KIND=KIND(0.0)), PARAMETER :: one=1
   CALL CGEMM("n","n",M,N,K,one,A,M,B,K,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_cnt)
SUBROUTINE smm_cnt(M,N,K,A,B,C)
   INTEGER :: M,N,K
   COMPLEX(KIND=KIND(0.0)) :: A(*), B(*), C(*)
   COMPLEX(KIND=KIND(0.0)), PARAMETER :: one=1
   CALL CGEMM("n","t",M,N,K,one,A,M,B,N,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_ctn)
SUBROUTINE smm_ctn(M,N,K,A,B,C)
   INTEGER :: M,N,K
   COMPLEX(KIND=KIND(0.0)) :: A(*), B(*), C(*)
   COMPLEX(KIND=KIND(0.0)), PARAMETER :: one=1
   CALL CGEMM("t","n",M,N,K,one,A,K,B,K,one,C,M)
END SUBROUTINE
#endif
#if !defined(__HAS_smm_ctt)
SUBROUTINE smm_ctt(M,N,K,A,B,C)
   INTEGER :: M,N,K
   COMPLEX(KIND=KIND(0.0)) :: A(*), B(*), C(*)
   COMPLEX(KIND=KIND(0.0)), PARAMETER :: one=1
   CALL CGEMM("t","t",M,N,K,one,A,K,B,N,one,C,M)
END SUBROUTINE
#endif
