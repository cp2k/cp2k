# Builds the SMM library in 5 steps:
# 1. Builds programs to test tiny MM performance.
# 2. Runs programs to test tiny MM performance.
# 3. Builds programs to test small MM performance; needs data from
#    tiny MM performance
# 4. Runs programs to test small MM performance
# 5. Builds library; needs data from small MM performance

# Recommended targets:
# do.compile.tiny
# do.bench.tiny
# do.compile.small
# do.bench.small
# do.lib.test

.PHONY: all
all: do.lib.test

.PHONY: clean
clean:
	rm -fr run_tiny* run_small* *.x *.mod *.o

include config.mk

# Some dependency relaxation
ifeq ($(LOOSE),yes)
LOOSE=|
else
LOOSE=
endif
# Some dependency relaxation
ifeq ($(VERYLOOSE),yes)
VERYLOOSE=|
else
VERYLOOSE=
endif

# Builds the rules
TINY_SET = $(foreach m,$(dims_tiny),$(foreach n,$(dims_tiny),$(foreach k,$(dims_tiny),$m_$n_$k)))

SMALL_SET = $(foreach m,$(dims_small),$(foreach n,$(dims_small),$(foreach k,$(dims_small),$m_$n_$k)))

# More variable transformations
ifeq ($(data_type),1)
    filedatatype=d
    gemm=DGEMM
    strdat=REAL(KIND=KIND(0.0D0))
endif
ifeq ($(data_type),2)
    filedatatype=s
    gemm=SGEMM
    strdat=REAL(KIND=KIND(0.0))
endif
ifeq ($(data_type),3)
    filedatatype=z
    gemm=ZGEMM
    strdat=COMPLEX(KIND=KIND(0.0D0))
endif
ifeq ($(data_type),4)
    filedatatype=c
    gemm=CGEMM
    strdat=COMPLEX(KIND=KIND(0.0))
endif

ifeq ($(transpose_flavor),1)
    trtype=nn
    ta=N
    tb=N
    decl=A(M,K), B(K,N)
    lds=LDA=M ; LDB=K
endif
ifeq ($(transpose_flavor),2)
    trtype=tn
    ta=T
    tb=N
    decl=A(K,M), B(K,N)
    lds=LDA=K ; LDB=K
endif
ifeq ($(transpose_flavor),3)
    trtype=nt
    ta=N
    tb=T
    decl=A(M,K), B(N,K)
    lds=LDA=M ; LDB=N
endif
ifeq ($(transpose_flavor),4)
    trtype=tt
    ta=T
    tb=T
    decl=A(K,M), B(N,K)
    lds=LDA=K ; LDB=N
endif

filetrans=$(filedatatype)$(trtype)



# After completion all tiny benchmarks should have been compiled.
.PHONY: do.compile.tiny
do.compile.tiny: $(patsubst %,run_tiny_%/tiny_find.x,$(TINY_SET))

# After completion all small benchmarks should have been compiled.
.PHONY: do.compile.small
do.compile.small: $(patsubst %,run_small_%/small_find.x,$(SMALL_SET))

# After completion all tiny benchmarks should have been run and the
# result file generated.
.PHONY: do.bench.tiny
do.bench.tiny: tiny_gen_optimal.out

# After completion all tiny benchmarks should have been run and the
# result file generated.
.PHONY: do.bench.small
do.bench.small: small_gen_optimal.out

# Perform library test
.PHONY: do.lib.test
do.lib.test: test_smm_$(filetrans).out



# Generate tiny source files
.PRECIOUS : run_tiny_%/tiny_find.f90
run_tiny_%/tiny_find.f90: $(LOOSE) tiny_gen.x
	mkdir -p run_tiny_$*
	./tiny_gen.x $(subst _, ,$*) $(transpose_flavor) $(data_type) > $@

# Compile and link tiny source files
.PRECIOUS : run_tiny_%/tiny_find.x
run_tiny_%/tiny_find.x: run_tiny_%/tiny_find.f90
	cd run_tiny_$* ; $(target_compile) tiny_find.f90 -o tiny_find.x.tmp && mv tiny_find.x.tmp tiny_find.x

# Benchmark tiny
# Must be ordered to come after all compilations are done so that all
# processors are devoted to running jobs.
.PRECIOUS : run_tiny_%/tiny_find.out
run_tiny_%/tiny_find.out: run_tiny_%/tiny_find.x | do.compile.tiny
	cd run_tiny_$* ; ./tiny_find.x > tiny_find.out.tmp && mv tiny_find.out.tmp tiny_find.out


# Generate small source files
.PRECIOUS : run_small_%/small_find.f90
run_small_%/small_find.f90: $(VERYLOOSE) tiny_gen_optimal.out $(LOOSE) small_gen.x
	mkdir -p run_small_$*
	./small_gen.x $(subst _, ,$*) $(transpose_flavor) $(data_type) $(SIMD_size) "tiny_gen_optimal.out" > $@

# Compile and link small source files
.PRECIOUS : run_small_%/small_find.x
run_small_%/small_find.x: run_small_%/small_find.f90
	cd run_small_$* ; $(target_compile) small_find.f90 -o small_find.x.tmp $(blas_linking) && mv small_find.x.tmp small_find.x

# Benchmark small.
# Must be ordered to come after all compilations are done so that all
# processors are devoted to running jobs.
.PRECIOUS : run_small_%/small_find.out
run_small_%/small_find.out: run_small_%/small_find.x | do.compile.small
	cd run_small_$* ; ./small_find.x > small_find.tmp && mv small_find.tmp small_find.out


# Compile the helper/generation programs

small_gen.x: small_gen.f90 mults.o multrec_gen.o
	$(host_compile) -o $@ $^

tiny_gen.x: tiny_gen.f90 mults.o
	$(host_compile) -o $@ $^

lib_gen.x: lib_gen.f90 mults.o multrec_gen.o
	$(host_compile) -o $@ $^

multrec_gen.o: multrec_gen.f90 mults.mod
	$(host_compile) -o $@ -c $<

# Collection of tiny benchmark results into a single file.
ALL_TINY_BENCH = $(patsubst %,run_tiny_%/tiny_find.out,$(TINY_SET))
tiny_gen_optimal.out: $(VERYLOOSE) $(ALL_TINY_BENCH)
	@rm -f $@
	@echo $(foreach combo,$(TINY_SET),"$(subst _, ,$(combo))" "$(shell tail -n 1 run_tiny_$(combo)/tiny_find.out) @") | sed -e 's/@ /\n/g' -e 's/@//' > $@

# Collection of small benchmark results into a single file.
ALL_SMALL_BENCH = $(patsubst %,run_small_%/small_find.out,$(SMALL_SET))
small_gen_optimal.out: $(VERYLOOSE) $(ALL_SMALL_BENCH)
	@rm -f $@
	@echo $(foreach combo,$(SMALL_SET),"$(subst _, ,$(combo))" "$(shell tail -n 1 run_small_$(combo)/small_find.out) @") | sed -e 's/@ /\n/g' -e 's/@//' > $@


# Library building

libdir: | lib

lib:
	mkdir -p lib

SMALL_OBJECTS = $(patsubst %,lib/smm_m_$(filetrans)_%.o,$(SMALL_SET))
DRIVER = smm_$(filetrans)

# SMM Library
lib/libsmm_$(filetrans).a: lib/$(DRIVER).o $(SMALL_OBJECTS)
	cd lib ; ar -r $(notdir $@) $(foreach f,$^,$(notdir $f))

# Testing program
test_smm_$(filetrans).x: test_smm_$(filetrans).f90 lib/libsmm_$(filetrans).a
	$(target_compile) -o $@ $< -Llib -lsmm_$(filetrans) $(blas_linking)

# Main library entry point
lib/$(DRIVER).o: lib/$(DRIVER).f90
	cd lib ; $(target_compile) -c $(notdir $<)

# Small MM object
.PRECIOUS: lib/smm_%.o
lib/smm_m_$(filetrans)_%.o: lib/smm_m_$(filetrans)_%.f90
	cd lib ; $(target_compile) -c $(notdir $<)

# Small MM source
.PRECIOUS: lib/smm_m_$(filetrans)_%.f90
lib/smm_m_$(filetrans)_%.f90: lib_gen.x small_gen_optimal.out | libdir
	./lib_gen.x $(subst _, ,$*) $(transpose_flavor) $(data_type) $(SIMD_size) "small_gen_optimal.out" "tiny_gen_optimal.out" > $@.tmp
	@mv $@.tmp $@


# Library test
.PRECIOUS : test_smm_$(filetrans).out
test_smm_$(filetrans).out: test_smm_$(filetrans).x
	./$< | tee $@.tmp
	@grep 'BLAS and smm yield different results' ./test_smm_${filetrans}.out &> /dev/null && ( \
	echo "Library is miscompiled ... removing lib/libsmm_$(filetrans).a" ; \
	rm -f lib/libsmm_${filetrans}.a ; exit 1) || true
	@echo "Done... check performance looking at test_smm_$(filetrans).out"
	@echo "Final library can be linked as -L$(shell pwd)/lib -lsmm_$(filetrans)"
	@mv $@.tmp $@

# Generic compilation rule for utility programs
%.o %.mod: %.f90
	$(host_compile) -c $<


# Ugly code generation

numsize  = $(shell echo $(dims_small) | awk "{print NF;}")
maxsize = $(shell echo $(dims_small) | awk "{maxs=0; for(i=1; i<=NF; i++) {maxs=(\$$i>maxs)?\$$i:maxs;}; print maxs;}")
# eles is a mapping of block sizes to sequence number.
# E.g., 1 4 -> (/0,1,0,0,1/))
eles = $(shell echo $(dims_small) | awk "{\
maxs = 0;\
for (i=1;i<=NF; i++) {maxs=(\$$i>maxs)?\$$i:maxs;};\
order = 1;\
for(i=1; i<=NF; i++) map[\$$i]=order++;\
printf \"(/0\";\
for(i=1; i<=maxs; i++) printf \",%d\", map[i];\
printf \"/)\n\";}")

lib/$(DRIVER).f90: config.in libdir
	@printf "SUBROUTINE smm_$(filetrans)(M,N,K,A,B,C)\n INTEGER :: M,N,K,LDA,LDB\n $(strdat) :: C(M,N), $(decl)\n" > $@.tmp
	printf " INTEGER, PARAMETER :: indx(0:$(maxsize))=&\n $(eles)\n" >> $@.tmp
	@printf " INTEGER :: im,in,ik,itot\n" >> $@.tmp
	@printf " $(strdat), PARAMETER :: one=1\n" >> $@.tmp
	@printf " $(lds)\n" >> $@.tmp
	@printf " IF (M<=$(maxsize)) THEN\n   im=indx(M)\n ELSE\n   im=0\n ENDIF\n" >> $@.tmp
	@printf " IF (N<=$(maxsize)) THEN\n   in=indx(N)\n ELSE\n   in=0\n ENDIF\n" >> $@.tmp
	@printf " IF (K<=$(maxsize)) THEN\n   ik=indx(K)\n ELSE\n   ik=0\n ENDIF\n" >> $@.tmp
	@printf " itot=(ik*($(numsize)+1)+in)*($(numsize)+1)+im\n" >> $@.tmp
	@printf " SELECT CASE(itot)\n" >> $@.tmp
	@echo $(dims_small) | awk "{\
		cnt=0;\
		for (myk=0;myk<=NF;myk++)\
		for (myn=0;myn<=NF;myn++)\
		for (mym=0;mym<=NF;mym++) {\
			printf \" CASE(%d)\n\",cnt++;\
			prod=myk*myn*mym;\
			if (prod==0)\
				printf \"    GOTO 999\n\";\
			else\
				printf \"    CALL smm_$(filetrans)_%d_%d_%d(A,B,C)\n\",\$$mym,\$$myn,\$$myk;\
		};\
	}" >> $@.tmp
	@printf " END SELECT\n" >> $@.tmp
	@printf " RETURN\n" >> $@.tmp
	@printf "999 CONTINUE \n CALL $(gemm)(\'%s\',\'%s\',M,N,K,one,A,LDA,B,LDB,one,C,M)\n" $(ta) $(tb) >> $@.tmp
	printf "END SUBROUTINE smm_$(filetrans)\n" >> $@.tmp
	mv $@.tmp $@

#Testing program
test_smm_$(filetrans).f90: config.in
	printf "MODULE WTF\n" > $@.tmp
	@printf "  INTERFACE MYRAND\n" >> $@.tmp
	@printf "    MODULE PROCEDURE SMYRAND, DMYRAND, CMYRAND, ZMYRAND\n" >> $@.tmp
	@printf "  END INTERFACE\n" >> $@.tmp
	@printf "CONTAINS\n" >> $@.tmp
	@printf "  SUBROUTINE DMYRAND(A)\n" >> $@.tmp
	@printf "    REAL(KIND=KIND(0.0D0)), DIMENSION(:,:) :: A\n" >> $@.tmp
	@printf "    REAL(KIND=KIND(0.0)), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Aeq\n" >> $@.tmp
	@printf "    CALL RANDOM_NUMBER(Aeq)\n" >> $@.tmp
	@printf "    A=Aeq\n" >> $@.tmp
	@printf "  END SUBROUTINE\n" >> $@.tmp
	@printf "  SUBROUTINE SMYRAND(A)\n" >> $@.tmp
	@printf "    REAL(KIND=KIND(0.0)), DIMENSION(:,:) :: A\n" >> $@.tmp
	@printf "    REAL(KIND=KIND(0.0)), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Aeq\n" >> $@.tmp
	@printf "    CALL RANDOM_NUMBER(Aeq)\n" >> $@.tmp
	@printf "    A=Aeq\n" >> $@.tmp
	@printf "  END SUBROUTINE\n" >> $@.tmp
	@printf "  SUBROUTINE CMYRAND(A)\n" >> $@.tmp
	@printf "    COMPLEX(KIND=KIND(0.0)), DIMENSION(:,:) :: A\n" >> $@.tmp
	@printf "    REAL(KIND=KIND(0.0)), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Aeq,Beq\n" >> $@.tmp
	@printf "    CALL RANDOM_NUMBER(Aeq)\n" >> $@.tmp
	@printf "    CALL RANDOM_NUMBER(Beq)\n" >> $@.tmp
	@printf "    A=CMPLX(Aeq,Beq,KIND=KIND(0.0))\n" >> $@.tmp
	@printf "  END SUBROUTINE\n" >> $@.tmp
	@printf "  SUBROUTINE ZMYRAND(A)\n" >> $@.tmp
	@printf "    COMPLEX(KIND=KIND(0.0D0)), DIMENSION(:,:) :: A\n" >> $@.tmp
	@printf "    REAL(KIND=KIND(0.0)), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Aeq,Beq\n" >> $@.tmp
	@printf "    CALL RANDOM_NUMBER(Aeq)\n" >> $@.tmp
	@printf "    CALL RANDOM_NUMBER(Beq)\n" >> $@.tmp
	@printf "    A=CMPLX(Aeq,Beq,KIND=KIND(0.0D0))\n" >> $@.tmp
	@printf "  END SUBROUTINE\n" >> $@.tmp
	@printf "END MODULE\n" >> $@.tmp
	@printf "SUBROUTINE testit(M,N,K)\n" >> $@.tmp
	@printf "  USE WTF\n" >> $@.tmp
	@printf "  IMPLICIT NONE\n" >> $@.tmp
	@printf "  INTEGER :: M,N,K\n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "  $(strdat) :: C1(M,N), C2(M,N)\n" >> $@.tmp
	@printf "  $(strdat) :: $(decl)\n" >> $@.tmp
	@printf "  $(strdat), PARAMETER :: one=1\n" >> $@.tmp
	@printf "  INTEGER :: i,LDA,LDB\n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "  REAL(KIND=KIND(0.0D0)) :: flops,gflop\n" >> $@.tmp
	@printf "  REAL :: t1,t2,t3,t4\n" >> $@.tmp
	@printf "  INTEGER :: Niter\n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "  flops=2*REAL(M,KIND=KIND(0.0D0))*N*K\n" >> $@.tmp
	@printf "  gflop=1000.0D0*1000.0D0*1000.0D0\n" >> $@.tmp
	@printf "  ! assume we would like to do 5 Gflop for testing a subroutine\n" >> $@.tmp
	@printf "  Niter=MAX(1,CEILING(MIN(10000000.0D0,5*gflop/flops)))\n" >> $@.tmp
	@printf "  $(lds)\n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "  DO i=1,10\n" >> $@.tmp
	@printf "     CALL MYRAND(A)\n" >> $@.tmp
	@printf "     CALL MYRAND(B)\n" >> $@.tmp
	@printf "     CALL MYRAND(C1)\n" >> $@.tmp
	@printf "     C2=C1\n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "     CALL $(gemm)(\"$(ta)\",\"$(tb)\",M,N,K,one,A,LDA,B,LDB,one,C1,M) \n" >> $@.tmp
	@printf "     CALL smm_$(filetrans)(M,N,K,A,B,C2)\n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "     IF (MAXVAL(ABS(C2-C1))>100*EPSILON(REAL(1.0,KIND=KIND(A(1,1))))) THEN\n" >> $@.tmp
	@printf "        write(6,*) \"Matrix size\",M,N,K\n" >> $@.tmp
	@printf "        write(6,*) \"A=\",A\n" >> $@.tmp
	@printf "        write(6,*) \"B=\",B\n" >> $@.tmp
	@printf "        write(6,*) \"C1=\",C1\n" >> $@.tmp
	@printf "        write(6,*) \"C2=\",C2\n" >> $@.tmp
	@printf "        write(6,*) \"BLAS and smm yield different results : possible compiler bug... do not use the library ;-)\"\n" >> $@.tmp
	@printf "        STOP\n" >> $@.tmp
	@printf "     ENDIF\n" >> $@.tmp
	@printf "  ENDDO\n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "  A=0; B=0; C1=0 ; C2=0\n" >> $@.tmp
	@printf " \n" >> $@.tmp
	@printf "  CALL CPU_TIME(t1) \n" >> $@.tmp
	@printf "  DO i=1,Niter\n" >> $@.tmp
	@printf "     CALL $(gemm)(\"$(ta)\",\"$(tb)\",M,N,K,one,A,LDA,B,LDB,one,C1,M) \n" >> $@.tmp
	@printf "  ENDDO\n" >> $@.tmp
	@printf "  CALL CPU_TIME(t2) \n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "  CALL CPU_TIME(t3)\n" >> $@.tmp
	@printf "  DO i=1,Niter\n" >> $@.tmp
	@printf "     CALL smm_$(filetrans)(M,N,K,A,B,C2)\n" >> $@.tmp
	@printf "  ENDDO\n" >> $@.tmp
	@printf "  CALL CPU_TIME(t4)\n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "  WRITE(6,'(A,I5,I5,I5,A,F6.3,A,F6.3,A,F12.3,A)') \"Matrix size \",M,N,K, &\n" >> $@.tmp
	@printf "        \" smm: \",Niter*flops/(t4-t3)/gflop,\" Gflops. Linked blas: \",Niter*flops/(t2-t1)/gflop,&\n" >> $@.tmp
	@printf "        \" Gflops. Performance ratio: \",((t2-t1)/(t4-t3))*100,\"%%\"\n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "END SUBROUTINE \n" >> $@.tmp
	@printf "\n" >> $@.tmp
	@printf "PROGRAM tester\n" >> $@.tmp
	@printf "  IMPLICIT NONE\n" >> $@.tmp
	@echo $(dims_small) | awk " {\
		for (myk=1;myk<=NF;myk++)\
		for (myn=1;myn<=NF;myn++)\
		for (mym=1;mym<=NF;mym++)\
			printf \"   CALL testit(%d,%d,%d)\n\",\$$mym,\$$myn,\$$myk;\
	}" >> $@.tmp
	@printf "  ! checking 'peak' performance (and a size likely outside of the library)\n" >> $@.tmp
	@printf "  CALL testit(1000,1000,1000)\n" >> $@.tmp
	printf "END PROGRAM\n" >> $@.tmp
	mv $@.tmp $@
