  ! basic molecular_DOS analysis as done in the paper
  ! P. Hunt, M. Sprik, and R. Vuilleumier
  ! CPL 376 (2003) 68-74 

  ! major short comings : still assumes that just 1 molecular species is present
  INTEGER, PARAMETER :: dbl = 8 
  REAL(KIND=dbl), DIMENSION(:), POINTER     :: evals              ! eigenstate evals
  REAL(KIND=dbl), DIMENSION(:,:), POINTER   :: molevals           ! molecular evals alpha_gammaI
  REAL(KIND=dbl), DIMENSION(:,:), POINTER   :: molspread          ! molecular spread beta_gammaI
  REAL(KIND=dbl), DIMENSION(:,:,:), POINTER :: projection         ! D_kgammaI
  INTEGER, DIMENSION(:), POINTER            :: nmolstates_array

  REAL(KIND=dbl), DIMENSION(:), POINTER   :: molevals_m1        ! moment_1 molecular evals alpha_gammaI
  REAL(KIND=dbl), DIMENSION(:), POINTER   :: molevals_m2        ! moment_2 molecular evals alpha_gammaI
  REAL(KIND=dbl), DIMENSION(:), POINTER   :: molspread_m1       ! moment_1 molecular spread beta_gammaI
  REAL(KIND=dbl), DIMENSION(:,:), POINTER :: projection_m1    ! D_kgammaI

  REAL(KIND=dbl), PARAMETER :: au2ev = 2.72113834598679E+01_dbl

  INTEGER :: nstates,nmol,nmolstates_max,nmolstates
  INTEGER :: i,imol,imol_read,istate
  INTEGER :: unit_in,unit_out
  REAL(KIND=dbl) :: norm,ener
  CHARACTER(LEN=11) :: filename

  ! read in all data 
  unit_in=11
  OPEN(unit_in,FILE="Molecular_DOS")
  READ(unit_in,*) nstates
  WRITE(6,*) "nstates = ",nstates
  ALLOCATE(evals(nstates))
  DO i=1,nstates
     READ(unit_in,*) evals(i)
  ENDDO

  READ(unit_in,*) nmol,nmolstates_max
  WRITE(6,*) "nmol = ",nmol
  ALLOCATE(projection(nstates,nmolstates_max,nmol)) ! holds (psi_k phi_gammaI)**2
  ALLOCATE(molevals(nmolstates_max,nmol))
  ALLOCATE(molspread(nmolstates_max,nmol))
  ALLOCATE(nmolstates_array(nmol))
  molevals=0.0_dbl
 
  DO imol=1,nmol
     READ(unit_in,*) imol_read,nmolstates
     nmolstates_array(imol)=nmolstates
     IF (imol_read .NE. imol) THEN
        WRITE(6,*) "Error in input file format"
        STOP
     ENDIF
     DO i=1,nmolstates 
        READ(unit_in,*) molevals(i,imol) ! Eq. 5
     ENDDO 
     DO i=1,nmolstates
        DO k=1,nstates
           READ(unit_in,*) projection(k,i,imol)  ! Eq. 6
        ENDDO
     ENDDO
  ENDDO

  ! done reading, some sanity checking
  DO imol=1,nmol
     nmolstates=nmolstates_array(imol)
     DO i=1,nmolstates
        norm = 0.0_dbl
        ener = 0.0_dbl
        DO k=1,nstates
           norm=norm+projection(k,i,imol)   ! Eq. 7
           ener=ener+evals(k)*projection(k,i,imol) ! Eq. 8
        ENDDO
        IF (ABS(norm-1.0_dbl) .gt. 1.0E-6 ) THEN
           write(6,*) "Mol ",imol," State ",i," Warning : Eq. 7 off by ",norm-1.0_dbl
        ENDIF
        IF (ABS(ener-molevals(i,imol)) .gt. 1.0E-6 ) THEN
           write(6,*) "Mol ",imol," State ",i," Warning : Eq. 8 off by ",ener-molevals(i,imol)
        ENDIF
     ENDDO 
  ENDDO

  ! compute beta (Eq. 9-11)
  DO imol=1,nmol
     nmolstates=nmolstates_array(imol)
     write(6,*) "Molecule ",imol," alpha (eV)"," beta (eV)"
     DO i=1,nmolstates
        molspread(i,imol)=0.0_dbl
        DO k=1,nstates
           molspread(i,imol)=molspread(i,imol)+projection(k,i,imol)*(evals(k)-molevals(i,imol))**2
        ENDDO
        write(6,*) i,molevals(i,imol)*au2ev,sqrt(molspread(i,imol))*au2ev
     ENDDO
  ENDDO

  ! produce some averages / this is only sound if one molecule type is present !!!!!!!
  ALLOCATE(molevals_m1(nmolstates_max))
  ALLOCATE(molevals_m2(nmolstates_max))
  ALLOCATE(molspread_m1(nmolstates_max))
  ALLOCATE(projection_m1(nstates,nmolstates_max)) 
  molevals_m1=0.0_dbl
  molevals_m2=0.0_dbl
  molspread_m1=0.0_dbl
  projection_m1=0.0_dbl
  DO imol=1,nmol 
    DO k=1,nmolstates_array(imol)
       molevals_m1(k)=molevals_m1(k)+(1.0_dbl/nmol) * molevals(k,imol)                  ! Eq. 12
       molevals_m2(k)=molevals_m2(k)+(1.0_dbl/nmol) * molevals(k,imol) ** 2             ! Eq. 13
       projection_m1(:,k)=projection_m1(:,k)+(1.0_dbl/nmol) * projection(:,k,imol)      ! Eq. 14
       molspread_m1(k)=molspread_m1(k)+(1.0_dbl/nmol) * molspread(k,imol)               ! Eq. 19 (part 1)
    ENDDO
  ENDDO
  DO k=1,nmolstates_max
     molevals_m2(k)=molevals_m2(k) - molevals_m1(k)**2                                  ! Eq. 19 (part 2)
  ENDDO
  
  ! write some output
  write(6,'(A)') "         state         alpha         delta          beta            mu (eV, no squares)"
  DO k=1,nmolstates_max
     write(6,'(I14,4(F14.8))') k, molevals_m1(k)*au2ev, sqrt(molevals_m2(k))*au2ev,  &
                                  sqrt(molspread_m1(k))*au2ev,sqrt(molevals_m2(k)+molspread_m1(k))*au2ev
  ENDDO
  write(6,*) ""
  write(6,*) "Separated Molecular DOS written to files "
  DO i=1,nmolstates_max
     WRITE(filename,'(A7,I4.4)') "MolDOS_",i
     OPEN(11,FILE=filename) 
     DO k=1,nstates
        write(11,*) evals(k)*au2ev,projection_m1(k,i)
     ENDDO
     CLOSE(11)
  ENDDO

  END

