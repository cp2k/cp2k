  program dist_mat
  IMPLICIT NONE
  ! N atoms
  INTEGER, PARAMETER :: N=96*4 ! 128 water system
  ! P x Q cpus
  INTEGER, PARAMETER :: P=8
  INTEGER, PARAMETER :: Q=4
  INTEGER, PARAMETER :: NCPU=P*Q
  ! the number of atom groups, any number is allowed, but see below
  INTEGER, PARAMETER :: S=N/3

  ! locals
  INTEGER, DIMENSION(N,N) :: matrix
  INTEGER, DIMENSION(3,N)   :: particle_loc
  INTEGER, DIMENSION(0:P*Q-1) :: local_blocks
  INTEGER, DIMENSION(0:P*Q-1) :: vector_entries_read
  INTEGER, DIMENSION(0:P*Q-1) :: vector_entries_write
  INTEGER :: i,j,CI,CJ,gp,gq,k
  LOGICAL :: store, OK

  write(6,*) "# atoms",N
  write(6,*) "# cpus ",NCPU," P = ",P," Q = ",Q
  write(6,*) "# groups ",S

  !------------------------------------------------------------------------------------
  ! assign particles to groups, and groups to cpus
  ! this is where the real load balancing is, but in principle every atom can be its own group
  ! even a simple scheme will do for a first version
  !
  ! for later versions good heuristics for assigning atoms to groups and groups to cpus 
  ! might include:
  ! - the number of groups should possibly be larger than P*Q
  !       -> every cpu gets at least a few groups to work on
  ! - each group has about the same number of basis functions
  !       -> the surface of the full matrix is about equally partitioned
  ! - atoms that are spatially near could form a group
  !       -> communication useful sparsity will appear in the sparse matrix
  ! - a processor gets a bunch of groups that are near in space
  !       -> communication useful sparsity will appear in the sparse matrix
  !       -> some relation with the rs-grids
  !  one possible scheme :
  !    all atoms in a subcell (neighborliststyle) could be a group
  !------------------------------------------------------------------------------------
  ! each particle goes to a group 
  DO i=1,N
     particle_loc(1,i)=(S*i-1)/N 
  ENDDO
  ! and depending only on the group to a cpu  
  DO i=1,N
     particle_loc(2,i)= (particle_loc(1,i)*P)/S      ! the P-CPU that owns the group
     particle_loc(3,i)= MOD(particle_loc(1,i),Q)     ! the Q-CPU that owns the group
     ! write(6,*) i,particle_loc(:,i)
  ENDDO

  !------------------------------------------------------------------------------------
  ! distribute the matrix elements
  ! m_ij = -1 if not used due to the symmetry
  ! m_ij = icpu if used
  !------------------------------------------------------------------------------------
  matrix = -1
  DO i=1,N
   DO j=1,N
      CI=particle_loc(1,i)
      CJ=particle_loc(1,j)
      ! first label blocks
      store = .FALSE.
      IF (CI < CJ ) THEN
         IF (MOD(CI+CJ,2)==1) store=.TRUE. 
      ELSE
          IF (CI==CJ) THEN
             IF (i < j ) THEN
                 IF (MOD(i+j,2)==1) store=.TRUE. 
             ELSE
                 IF (MOD(i+j,2)==0) store=.TRUE. 
             ENDIF 
          ELSE
             IF (MOD(CI+CJ,2)==0) store=.TRUE. 
          ENDIF
      ENDIF 
      IF (store) THEN
          gp = particle_loc(2,j)
          gq = particle_loc(3,i)
          matrix(i,j)= gp+gq*P 
      ENDIF
   ENDDO
  ENDDO

  !------------------------------------------------------------------------------------
  ! compute properties of the distribution
  !------------------------------------------------------------------------------------

  ! check if this is indeed an acceptable matrix distribution. i.e. respecting symmetry
  OK = .TRUE.
  DO i=1,N
     if (matrix(i,i)==-1) OK= .FALSE.
     DO j=i+1,N 
        if (matrix(i,j)==-1 .AND. matrix(j,i)==-1) OK = .FALSE.
        if (matrix(i,j)/=-1 .AND. matrix(j,i)/=-1) OK = .FALSE.
     ENDDO
  ENDDO
  IF (OK) THEN
     write(6,*) "Matrix distribution respects symmetric nature"
  ELSE
     write(6,*) "Matrix distribution DOES NOT respect symmetric nature"
  ENDIF
  ! figure out the number of blocks on each cpu
  local_blocks=0
  DO i=1,N
     DO j=1,N 
        if (matrix(i,j).ge.0) local_blocks(matrix(i,j))=local_blocks(matrix(i,j))+1
     ENDDO
  ENDDO
  write(6,*) "number of local blocks  min ",MINVAL(local_blocks), &
                                    " max ",MAXVAL(local_blocks), " of ", (N*(N-1))/2
  DO i=0,P-1
   DO j=0,Q-1
     ! write(6,*) i*P+j,i,j,local_blocks(i*Q+j)
   ENDDO
  ENDDO

  ! figure out the numbers of columns occupied
  vector_entries_read=0
  DO k=0,P*Q-1
   DO j=1,N 
     store=.FALSE.
     DO i=1,N
        if (matrix(i,j).EQ.k) store = .TRUE.
     ENDDO
     if (store) vector_entries_read(k)=vector_entries_read(k)+1
   ENDDO
  ENDDO
  write(6,*) "number of entries read  min ",MINVAL(vector_entries_read), &
                                    " max ",MAXVAL(vector_entries_read), " of ", N
  DO i=0,P-1
   DO j=0,Q-1
     ! write(6,*) i+j*P,i,j,vector_entries_read(i+j*P)
   ENDDO
  ENDDO

  ! figure out the numbers of rows occupied
  vector_entries_write=0
  DO k=0,P*Q-1
   DO j=1,N 
     store=.FALSE.
     DO i=1,N
        if (matrix(j,i).EQ.k) store = .TRUE.
     ENDDO
     if (store) vector_entries_write(k)=vector_entries_write(k)+1
   ENDDO
  ENDDO
  write(6,*) "number of entries write min ",MINVAL(vector_entries_write), &
                                    " max ",MAXVAL(vector_entries_write), " of ", N
  DO i=0,P-1
   DO j=0,Q-1
     ! write(6,*) i+j*P,i,j,vector_entries_write(i+j*P)
   ENDDO
  ENDDO

  !------------------------------------------------------------------------------------
  ! write the local matrix entries to fort.(10+icpu), nice for gnuplot
  !------------------------------------------------------------------------------------
  DO k=0,P*Q-1
     DO i=1,N
       DO j=1,N 
          if (matrix(i,j)==k) write(10+k,*) i,j
       ENDDO
     ENDDO
  ENDDO
  end program dist_mat
