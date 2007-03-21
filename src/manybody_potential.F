!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/manybody_potential [1.0] *
!!
!!   NAME
!!     manybody_potential
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE manybody_potential

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,pbc
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE pair_potential_types,            ONLY: ea_type,&
                                             tersoff_type,&
                                             pair_potential_pp_type,&
                                             pair_potential_single_type
  USE particle_types,                  ONLY: particle_type
  USE splines,                         ONLY: init_spline,&
                                             init_splinexy,&
                                             potential_s,&
                                             spline_data_create,&
                                             spline_data_release,&
                                             spline_data_type,&
                                             spline_environment_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_force_eam
  PUBLIC :: density_nonbond, energy_manybody
  PUBLIC :: force_nonbond_manybody
  PUBLIC :: destroy_newlist,get_newlist
  PUBLIC :: rho_f, rho_s
  PUBLIC :: spline_manybody_control
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_potential'

CONTAINS

  !******************************************************************************

  FUNCTION manybody_pot ( pot, r )

    TYPE(pair_potential_single_type), &
      INTENT(IN)                             :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: manybody_pot

    IF (pot%type==ea_type) THEN
       ! No splining necessary for EAM.  Zero splines
       manybody_pot = 0._dp
    ELSEIF (pot%type==tersoff_type) THEN
       WRITE(*,*) "ERROR: lsdjkfjklsdflkdsfklf"
       STOP
    ELSE
       CALL stop_program("manybody_potential%ener_pot","Potential type unknown")
    END IF

  END FUNCTION manybody_pot

  !******************************************************************************

  ! creates the splines for the manybody contribution

  SUBROUTINE spline_manybody_control ( spline_env, potparm, eps_spline, max_energy, rlow_nb,&
       emax_spline, iw, iw2, error)

    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), INTENT(IN)                :: eps_spline, max_energy, &
                                                rlow_nb, emax_spline
    INTEGER, INTENT(IN)                      :: iw, iw2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_manybody_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, jx, kx, n, nkx, &
                                                npoints, nsize, nspline, ntype
    LOGICAL                                  :: check_bounds, failure
    REAL(KIND=dp)                            :: diffmax, dx2, e, ediff, &
                                                hicut, locut, x, x2, xdum, &
                                                xdum1, xsav
    TYPE(spline_data_type), POINTER          :: spline_data

    failure = .FALSE.
    n = 0
    ntype = SIZE ( potparm%pot, 1 )
    n = 0
    DO i = 1, ntype ! i:  first atom type
       DO j = 1, i  ! j:  second atom type
          n = n + 1
          nspline = n
          SELECT CASE ( potparm%pot ( i, j ) % pot%type )
          CASE ( ea_type )
             spline_data => spline_env % spline_data_p ( n ) % spline_data
             locut = rlow_nb
             hicut = SQRT ( potparm%pot ( i, j )%pot % rcutsq )
             IF (ABS(hicut) <= 1.0E-15_dp) hicut = 1.0E3_dp
             ! Find the real locut according emax_spline
             dx2 = (hicut-locut)/10000.0_dp
             x   = hicut
             DO jx = 1, 10000
                e = manybody_pot ( potparm%pot ( i, j )%pot, x )
                IF (ABS(e) > emax_spline) THEN
                   locut = x
                   EXIT
                END IF
                x = x - dx2
             END DO
             npoints = 20
             IF (jx == 10001) npoints = 2  ! Zero spline
             DO WHILE (.TRUE.)

                CALL init_splinexy( spline_data, npoints+1 )

                dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
                x2  = 1.0_dp/hicut**2
                DO jx = 1, npoints + 1 ! jx: loop over 1/distance**2
                   x = SQRT ( 1.0_dp/x2 )
                   e = manybody_pot ( potparm%pot ( i, j )%pot, x )
                   spline_data % x(jx) = x2
                   spline_data % y(jx) = e
                   x2 = x2 + dx2
                END DO

                CALL init_spline(spline_data,dx=dx2)

                ! This is the check for required accuracy on spline setup
                dx2 = (hicut - locut)/REAL(5*npoints+1,KIND=dp)
                x2  = locut + dx2
                diffmax = 0.0_dp
                xsav = hicut
                DO jx = 1, 5*npoints
                   x = x2
                   e = manybody_pot ( potparm%pot ( i, j )%pot, x )
                   IF (e<max_energy) THEN
                      xdum1   = ABS(e-potential_s(spline_env%spline_data_p(n)%spline_data,x*x,xdum))
                      diffmax = MAX(diffmax,xdum1)
                      xsav = MIN(x,xsav)
                   END IF
                   x2 = x2 + dx2
                END DO
                IF (npoints > 5E4) THEN
                   IF (iw>0) THEN
                      WRITE(iw,'(A,I5,A,F15.9)') "SPLINE_INFO| Number of points: ",&
                           npoints," obtained accuracy ",diffmax
                   END IF
                   WRITE(*,'(A)')"MM SPLINE: no convergence on required accuracy (adjust EPS_SPLINE and rerun)"
                   CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                END IF
                IF (diffmax > eps_spline ) THEN
                   npoints   = CEILING(1.2_dp*REAL(npoints,KIND=dp))
                ELSE
                   EXIT
                END IF
             END DO

             !check spline existence
             nsize = SIZE(spline_env%spline_data_p(n)%spline_data%y)
             DO kx = 1, n-1
                nkx = SIZE(spline_env%spline_data_p(kx)%spline_data%y)
                check_bounds = (spline_env%spline_data_p(kx)%spline_data%x1 == &
                     spline_env%spline_data_p(n )%spline_data%x1).AND.&
                     (spline_env%spline_data_p(kx)%spline_data%xn == &
                     spline_env%spline_data_p(n )%spline_data%xn)
                IF (nsize==nkx) THEN
                   ediff = MAXVAL(ABS(spline_env%spline_data_p(n)%spline_data%y-&
                        spline_env%spline_data_p(kx)%spline_data%y))
                   IF (ediff <= eps_spline) THEN
                      CALL spline_data_release(spline_env%spline_data_p(n)%spline_data,&
                           error)
                      CALL spline_data_create ( spline_env%spline_data_p(n)%spline_data,error=error)
                      n = n - 1
                      nspline = kx
                      EXIT
                   END IF
                END IF
             END DO
             IF (iw>0) THEN
                WRITE(iw,*)
                WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Spline number:",i*(i-1)/2+j
                WRITE(iw,'("SPLINE_INFO|",5X,A,I7,A,2I7)')"Unique Spline number",nspline," Kinds involved:",i,j
                WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of points used in the Splines:",npoints
                WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy requested:",eps_spline
                WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy achieved: ",diffmax
                WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX  used for spline (bohr):               ",&
                     locut,hicut
                WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX used to achieve spline accuracy (bohr):",&
                     xsav,hicut
                WRITE(iw,*)
             END IF
             ! print spline data on file
             IF (nspline==n) THEN
                IF (iw2>0) THEN
                   dx2 = (hicut - locut)/REAL(5000+1,KIND=dp)
                   x2  = locut + dx2
                   DO jx = 1, 5000
                      x = x2
                      e = manybody_pot ( potparm%pot ( i, j )%pot, x )
                      WRITE(2000+n,*) x, e, potential_s(spline_env%spline_data_p(n)%spline_data,x*x,xdum)
                      x2 = x2 + dx2
                   END DO
                   WRITE(2000+n,*)
                END IF
             END IF

          CASE default
             ! Nullify spline_data not involved in interaction
             CALL spline_data_release ( spline_env % spline_data_p ( n ) % spline_data ,error=error)
             n = n - 1
             nspline = n
          END SELECT
          spline_env % spltab ( i, j ) = nspline
          spline_env % spltab ( j, i ) = nspline

       END DO ! j:  second atom type
    END DO ! i:  first atom type

    IF (iw>0) THEN
       WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of manybody potential Splines allocated:",&
            n
    END IF

  END SUBROUTINE spline_manybody_control

  !******************************************************************************

  SUBROUTINE rho_s ( spline_table, rijsq, energy, fscalar )


    TYPE(spline_data_type), POINTER          :: spline_table
    REAL(KIND=dp), INTENT(IN)                :: rijsq
    REAL(KIND=dp), INTENT(OUT)               :: energy, fscalar

!------------------------------------------------------------------------------

    energy = potential_s ( spline_table, rijsq, fscalar )

  END SUBROUTINE rho_s

  !******************************************************************************

  SUBROUTINE rho_f ( rijsq, potparm, iatom, jatom, energy )


    REAL(KIND=dp), INTENT(IN)                :: rijsq
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), INTENT(OUT)               :: energy

    REAL(KIND=dp)                            :: rij

!------------------------------------------------------------------------------

    rij = SQRT ( rijsq )

    ! non-bond energy
    energy = manybody_pot ( potparm %pot( iatom, jatom )%pot, rij )

  END SUBROUTINE rho_f

  !******************************************************************************







  SUBROUTINE density_nonbond ( fist_nonbond_env, particle_set, cell, para_env, error )
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'density_nonbond'

    INTEGER                                  :: atom_a, atom_b, handle, &
                                                ikind, ilist, ipair, &
                                                iparticle, istat, jkind, &
                                                nkinds, nparticle
    LOGICAL                                  :: do_eam
    REAL(KIND=dp)                            :: flops, rab2, rab2_max
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rho
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc

    CALL timeset ( routineN, 'I', 'Mflops', handle )
    flops = 0.0_dp
    do_eam = .FALSE.
    CALL fist_nonbond_env_get ( fist_nonbond_env, nonbonded = nonbonded, &
         potparm = potparm, r_last_update=r_last_update,&
         r_last_update_pbc=r_last_update_pbc,error=error)
    nkinds = SIZE ( potparm%pot, 1 )
    DO ikind = 1, nkinds
       DO jkind = ikind, nkinds
          IF (potparm%pot( ikind, jkind )%pot % TYPE== ea_type) do_eam = .TRUE.
       END DO
    END DO
    ! Only if EAM potential are present
    IF (do_eam) THEN
       nparticle = SIZE ( particle_set )
       ALLOCATE ( rho ( nparticle ), STAT = istat )
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "rho",nparticle*dp_size)
       rho ( : ) = 0._dp
       ! Starting the force loop
       DO ilist=1,nonbonded%nlists
          neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
          IF (neighbor_kind_pair%npairs ==0) CYCLE
          ikind = neighbor_kind_pair%ikind
          jkind = neighbor_kind_pair%jkind
          SELECT CASE ( potparm%pot( ikind, jkind )%pot % TYPE )
          CASE ( ea_type)
             rab2_max = potparm%pot(ikind,jkind)%pot%rcutsq
             CALL matvec_3x3(cell_v, cell%hmat,REAL(neighbor_kind_pair%cell_vector,KIND=dp))
             DO ipair = 1, neighbor_kind_pair%npairs
                atom_a = neighbor_kind_pair%list(1,ipair)
                atom_b = neighbor_kind_pair%list(2,ipair)
                rab = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
                rab = rab+cell_v
                rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                IF (rab2 <= rab2_max) THEN
                   pot => potparm %pot( ikind, jkind )%pot
                   CALL get_rho_eam ( pot, rab2, atom_a, atom_b, rho  )
                END IF
             END DO
          END SELECT
       END DO
       CALL mp_sum( rho, para_env%group )
       DO iparticle=1,nparticle
          particle_set(iparticle)%rho = rho(iparticle)
       END DO

       DEALLOCATE ( rho, STAT = istat )
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "deallocate rho")
    END IF
    CALL timestop ( flops, handle )

  END SUBROUTINE density_nonbond
  !******************************************************************************
  SUBROUTINE energy_manybody ( fist_nonbond_env, atomic_kind_set, &
       local_particles, particle_set, cell, pot_manybody, para_env, error )

    ! computes the embedding contribution to the energy

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), INTENT(INOUT)                  :: pot_manybody
    TYPE(cp_para_env_type), POINTER,OPTIONAL :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'energy_manybody', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, index, index1,index2, &
                                                iparticle, iparticle_local, &
                                                jparticle, kparticle, &
                                                istat, nparticle, &
                                                nparticle_local
    REAL(dp)                                 :: embed, pot_loc, qr
    REAL(dp), POINTER                        :: fembed( : )
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot

    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    INTEGER, DIMENSION(:, :), POINTER        :: list
    INTEGER                                  :: iatom,natom
    INTEGER                                  :: ipair,npairs
    INTEGER                                  :: ilist,nnlists
    INTEGER, POINTER                         :: newlist(:,:),nl(:)


    REAL(KIND=dp)                            :: rij(3),drij
    REAL(dp)                                 :: a_ij,b_ij,b_ji
    REAL(dp)                                 :: f_C, f_R, f_A
    REAL(dp)                                 :: RpD
    INTEGER                                  :: handle

    CALL timeset ( routineN, handle )

    CALL fist_nonbond_env_get ( fist_nonbond_env, potparm = potparm ,error=error)
    DO ikind = 1, SIZE ( atomic_kind_set )
       pot => potparm %pot ( ikind, ikind ) % pot
       SELECT CASE ( pot % type )
       CASE DEFAULT
          CYCLE
       CASE ( ea_type )

          nparticle = SIZE ( particle_set )
          ALLOCATE ( fembed ( nparticle ), STAT = istat )
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
               "fembed",nparticle*dp_size)
          fembed ( : ) = 0._dp
          ! computation of embedding function and energy
          nparticle_local = local_particles%n_el(ikind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(ikind)%array(iparticle_local)
             index = INT ( particle_set ( iparticle ) % rho / pot % eam % drhoar) + 1
             IF ( index > 4999 ) index = 5000
             qr = particle_set ( iparticle ) % rho  -  pot % eam % rhoval ( index )

             embed = pot % eam % frho ( index ) + &
                  qr * pot % eam % frhop ( index )
             fembed ( iparticle ) = pot % eam % frhop ( index ) + &
                  qr * ( pot % eam % frhop ( index + 1 ) -    &
                  pot % eam % frhop ( index ) ) / pot % eam %  drhoar

             pot_manybody = pot_manybody + embed
          END DO
          ! communicate data
          CALL mp_sum( fembed, para_env%group )
          DO iparticle=1,nparticle
             particle_set(iparticle)%f_embed = fembed(iparticle)
          END DO

          DEALLOCATE ( fembed, STAT = istat )
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
               "deallocate fembed")

       CASE ( tersoff_type )

         RpD = pot%tersoff%bigR + pot%tersoff%bigD

         CALL fist_nonbond_env_get(fist_nonbond_env,nonbonded=nonbonded,error=error)

         natom = SIZE(particle_set)
    
         CALL get_newlist(fist_nonbond_env,newlist,nl,natom,error)

         DO ilist=1,nonbonded%nlists
            neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
            npairs=neighbor_kind_pair%npairs
            IF (npairs ==0) CYCLE
            list  =>neighbor_kind_pair%list
            DO ipair = 1, npairs
               iparticle = list(1,ipair)
               jparticle = list(2,ipair)

               pot_loc = 0.0_dp
               rij(:) = pbc(particle_set(jparticle)%r(:),particle_set(iparticle)%r(:),cell)
               drij = SQRT(DOT_PRODUCT(rij,rij))
               IF(drij>RpD) CYCLE
               !a_ij   = ter_a_ij(pot,drij)
               !a_ij   = 1.0_dp
               b_ij   = ter_b_ij(pot,particle_set,cell,iparticle,jparticle,newlist,nl)
               b_ji   = ter_b_ij(pot,particle_set,cell,jparticle,iparticle,newlist,nl)
               f_C    = ter_f_C(pot,drij)
               f_A    = ter_f_A(pot,drij)
               f_R    = ter_f_R(pot,drij)
               pot_loc = f_C*(f_R + 0.5_dp*(b_ij+b_ji)*f_A)
               pot_manybody = pot_manybody + pot_loc

            END DO
        END DO

        CALL destroy_newlist(newlist,nl,error)

       END SELECT
    END DO
    CALL timestop ( handle )
  END SUBROUTINE energy_manybody


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_newlist(fist_nonbond_env,newlist,nl,natom,error)
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    INTEGER, POINTER                         :: newlist(:,:),nl(:)
    INTEGER                                  :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    INTEGER, DIMENSION(:, :), POINTER        :: list
    INTEGER                                  :: ilist,nlists
    INTEGER                                  :: ipair,npairs
    INTEGER                                  :: iparticle,jparticle

    CALL fist_nonbond_env_get(fist_nonbond_env,nonbonded=nonbonded,error=error)

    NULLIFY(newlist,nl)
    ALLOCATE(newlist(20,natom))
    ALLOCATE(nl(natom))
    newlist(:,:)=0
    nl(:)=0
    
    DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs=neighbor_kind_pair%npairs
       IF (npairs ==0) CYCLE
       list  =>neighbor_kind_pair%list
       DO ipair = 1, npairs
          iparticle = list(1,ipair)
          jparticle = list(2,ipair)

          nl(iparticle) = nl(iparticle) + 1
          newlist(nl(iparticle),iparticle) = jparticle
          IF(nl(iparticle)>20) STOP
          nl(jparticle) = nl(jparticle) + 1
          newlist(nl(jparticle),jparticle) = iparticle
          IF(nl(jparticle)>20) STOP

       END DO
    END DO
  END SUBROUTINE get_newlist

  SUBROUTINE destroy_newlist(newlist,nl,error)
    INTEGER, POINTER                         :: newlist(:,:),nl(:)
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: istat

    IF(ASSOCIATED(nl)) THEN
      DEALLOCATE(nl,STAT=istat)
      !CP ERROR CHECK HERE
    END IF
    IF(ASSOCIATED(newlist)) THEN
      DEALLOCATE(newlist,STAT=istat)
      !CP ERROR CHECK HERE
    END IF

  END SUBROUTINE destroy_newlist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_f_C(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_C', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_f_C
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: bigR,bigD
    REAL(KIND=dp)                                       :: RpD,RmD
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )
    bigR = pot%tersoff%bigR
    bigD = pot%tersoff%bigD
    RmD = pot%tersoff%bigR-pot%tersoff%bigD
    RpD = pot%tersoff%bigR+pot%tersoff%bigD
    ter_f_C = 0.0_dp
    IF(r < RmD) ter_f_C = 1.0_dp
    IF(r > RpD) ter_f_C = 0.0_dp
    IF( (r < RpD) .AND. (r > RmD) ) THEN
      ter_f_C = 0.5_dp *(1.0_dp - SIN(0.5_dp*PI*(r-bigR)/(bigD)))
    END IF
    !CALL timestop ( handle )
  END FUNCTION ter_f_C

  FUNCTION ter_f_C_d(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_C_d', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_f_C_d
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: bigR,bigD
    REAL(KIND=dp)                                       :: RpD,RmD
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    bigR = pot%tersoff%bigR
    bigD = pot%tersoff%bigD
    RmD = pot%tersoff%bigR-pot%tersoff%bigD
    RpD = pot%tersoff%bigR+pot%tersoff%bigD
    ter_f_C_d = 0.0_dp
    IF(r < RmD) ter_f_C_d = 0.0_dp
    IF(r > RpD) ter_f_C_d = 0.0_dp
    IF( (r < RpD) .AND. (r > RmD) ) THEN
      ter_f_C_d = (0.25_dp*PI/bigD) * COS(0.5_dp*PI*(r-bigR)/(bigD)) / r
    END IF
    !CALL timestop ( handle )
  END FUNCTION ter_f_C_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_f_R(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_R', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_f_R
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    A = pot%tersoff%A
    lambda1 = pot%tersoff%lambda1
    ter_f_R = 0.0_dp
    ter_f_R = A * EXP(-lambda1*r)
    !CALL timestop ( handle )
  END FUNCTION ter_f_R

  FUNCTION ter_f_R_d(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_R_d', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_f_R_d
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: f_R
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    A = pot%tersoff%A
    lambda1 = pot%tersoff%lambda1
    f_R = A * EXP(-lambda1*r)
    ter_f_R_d = 0.0_dp
    ter_f_R_d = lambda1 * f_R / r
    !CALL timestop ( handle )
  END FUNCTION ter_f_R_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_f_A(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_A', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_f_A
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    B = pot%tersoff%B
    lambda2 = pot%tersoff%lambda2
    ter_f_A = 0.0_dp
    ter_f_A = - B * EXP(-lambda2*r)
    !CALL timestop ( handle )
  END FUNCTION ter_f_A

  FUNCTION ter_f_A_d(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_A_d', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_f_A_d
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    B = pot%tersoff%B
    lambda2 = pot%tersoff%lambda2
    ter_f_A_d = 0.0_dp
    ter_f_A_d = - B *lambda2 * EXP(-lambda2*r) / r
    !CALL timestop ( handle )
  END FUNCTION ter_f_A_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_a_ij(pot,r)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: r
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_a_ij', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_a_ij
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    n = pot%tersoff%n
    alpha = pot%tersoff%alpha
    ter_a_ij = 0.0_dp
    !Note alpha = 0.0_dp for the parameters in the paper so using simplified term
    !ter_a_ij = (1.0_dp+(alpha*ter_n_ij(pot,iparticle,jparticle,r))**n)**(-0.5_dp/n)
    ter_a_ij = 1.0_dp
    !CALL timestop ( handle )
  END FUNCTION ter_a_ij
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_b_ij(pot,particle_set,cell,iparticle,jparticle,newlist,nl)
    TYPE(pair_potential_single_type), POINTER           :: pot
    TYPE(particle_type), INTENT(IN)                     :: particle_set( : )
    TYPE(cell_type), POINTER                            :: cell
    INTEGER, INTENT(IN)                                 :: iparticle, jparticle
    INTEGER, POINTER                                    :: newlist(:,:),nl(:)
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_b_ij', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_b_ij
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: zeta_ij
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    n = pot%tersoff%n
    beta = pot%tersoff%beta
    ter_b_ij = 0.0_dp
    zeta_ij = ter_zeta_ij(pot,particle_set,cell,iparticle,jparticle,newlist,nl)
    ter_b_ij = (1.0_dp+(beta*zeta_ij)**n)**(-0.5_dp/n)
    !CALL timestop ( handle )
  END FUNCTION ter_b_ij

  FUNCTION ter_b_ij_d(pot,particle_set,cell,iparticle,jparticle,newlist,nl)
    TYPE(pair_potential_single_type), POINTER           :: pot
    TYPE(particle_type), INTENT(IN)                     :: particle_set( : )
    TYPE(cell_type), POINTER                            :: cell
    INTEGER, INTENT(IN)                                 :: iparticle, jparticle
    INTEGER, POINTER                                    :: newlist(:,:),nl(:)
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_b_ij_d', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_b_ij_d
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: beta_n,zeta_ij
    REAL(KIND=dp)                                       :: zeta_ij_nm1,zeta_ij_n
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    n           = pot%tersoff%n
    beta        = pot%tersoff%beta
    beta_n      = beta**n
    zeta_ij     = ter_zeta_ij(pot,particle_set,cell,iparticle,jparticle,newlist,nl)
    zeta_ij_nm1 = zeta_ij**(n-1.0_dp)
    zeta_ij_n   = zeta_ij**(n)

    ter_b_ij_d = 0.0_dp
    ter_b_ij_d = -0.5_dp * beta_n * zeta_ij_nm1 * &
                   ((1.0_dp + beta_n * zeta_ij_n)**((-0.5_dp/n) - 1.0_dp))
    !CALL timestop ( handle )
  END FUNCTION ter_b_ij_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_zeta_ij(pot,particle_set,cell,iparticle,jparticle,newlist,nl)
    TYPE(pair_potential_single_type), POINTER           :: pot
    TYPE(particle_type), INTENT(IN)                     :: particle_set( : )
    TYPE(cell_type), POINTER                            :: cell
    INTEGER, INTENT(IN)                                 :: iparticle, jparticle
    INTEGER, POINTER                                    :: newlist(:,:),nl(:)
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_zeta_ij', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_zeta_ij
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: RpD
    INTEGER                                             :: kparticle, nparticle
    REAL(KIND=dp)                                       :: rij(3),rik(3),drij,drik
    REAL(KIND=dp)                                       :: costheta,sintheta
    REAL(KIND=dp)                                       :: f_C,gterm,expterm
    INTEGER                                             :: handle
    INTEGER :: ipair,npair

    !CALL timeset ( routineN, handle )

    ter_zeta_ij = 0.0_dp

    n = pot%tersoff%n
    lambda3 = pot%tersoff%lambda3
    RpD = pot%tersoff%bigR + pot%tersoff%bigD
    rij(:) = pbc(particle_set(jparticle)%r(:),particle_set(iparticle)%r(:),cell)
    drij = SQRT(DOT_PRODUCT(rij,rij))
    ter_zeta_ij = 0.0_dp
    !npair = nl(iparticle)
    !DO ipair=1,npair
    !  kparticle = newlist(ipair,iparticle)
    nparticle = SIZE ( particle_set )
    DO kparticle=1,nparticle
      IF((iparticle==kparticle).OR.(jparticle==kparticle)) CYCLE
      rik(:) = pbc(particle_set(kparticle)%r(:),particle_set(iparticle)%r(:),cell)
      drik = SQRT(DOT_PRODUCT(rik,rik))
      IF(drik>RpD) CYCLE
      costheta = DOT_PRODUCT(rij,rik) / (drij*drik)
      IF(costheta<-1.0_dp) costheta=-1.0_dp
      IF(costheta>+1.0_dp) costheta=+1.0_dp
      f_C = ter_f_C(pot,drik)
      gterm = ter_g(pot,costheta)
      expterm = exp((lambda3*(drij-drik))**3)
      ter_zeta_ij = ter_zeta_ij + f_C*gterm*expterm
    END DO
    !CALL timestop ( handle )
  END FUNCTION ter_zeta_ij

  SUBROUTINE ter_zeta_ij_d(pot,particle_set,cell,iparticle,jparticle,f_nonbond,ptens_nonbond,prefactor,newlist,nl)
    TYPE(pair_potential_single_type), POINTER           :: pot
    TYPE(particle_type), INTENT(IN)                     :: particle_set( : )
    TYPE(cell_type), POINTER                            :: cell
    INTEGER, INTENT(IN)                                 :: iparticle, jparticle
    REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)       :: f_nonbond,ptens_nonbond
    REAL(KIND=dp), INTENT(IN)                           :: prefactor
    INTEGER, POINTER                                    :: newlist(:,:),nl(:)
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_zeta_ij_d', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: RpD
    INTEGER                                             :: kparticle, nparticle
    REAL(KIND=dp), DIMENSION(3)                         :: rij,rik,rij_hat,rik_hat
    REAL(KIND=dp)                                       :: costheta
    REAL(KIND=dp)                                       :: drij,drik
    REAL(KIND=dp), DIMENSION(3)                         :: dri,drj,drk
    REAL(KIND=dp), DIMENSION(3)                         :: dcosdri,dcosdrj,dcosdrk
    REAL(KIND=dp)                                       :: f_C,f_C_d
    REAL(KIND=dp)                                       :: gterm,gterm_d
    REAL(KIND=dp)                                       :: expterm,expterm_d
    REAL(KIND=dp)                                       :: crapi,crapj,crapk
    INTEGER                                             :: handle
    INTEGER :: ipair,npair

    !CALL timeset ( routineN, handle )

    n = pot%tersoff%n
    lambda3 = pot%tersoff%lambda3
    RpD = pot%tersoff%bigR + pot%tersoff%bigD

    rij(:) = pbc(particle_set(jparticle)%r(:),particle_set(iparticle)%r(:),cell)
    drij = SQRT(DOT_PRODUCT(rij,rij))
    rij_hat(:) = rij(:) / drij

    !npair = nl(iparticle)
    !DO ipair=1,npair
    !  kparticle = newlist(ipair,iparticle)
    nparticle = SIZE ( particle_set )
    DO kparticle=1,nparticle
      IF((iparticle==kparticle).OR.(jparticle==kparticle)) CYCLE
      rik(:) = pbc(particle_set(kparticle)%r(:),particle_set(iparticle)%r(:),cell)
      drik = SQRT(DOT_PRODUCT(rik,rik))
      rik_hat(:) = rik(:) / drik

      IF(drik>RpD) CYCLE

      costheta = DOT_PRODUCT(rij,rik) / (drij*drik)
      IF(costheta<-1.0_dp) costheta=-1.0_dp
      IF(costheta>+1.0_dp) costheta=+1.0_dp

      CALL ter_costheta_d(particle_set,cell,iparticle,jparticle,kparticle,&
                          dcosdri,dcosdrj,dcosdrk)

      f_C       = ter_f_C(pot,drik)
      f_C_d     = ter_f_C_d(pot,drik)
      gterm     = ter_g(pot,costheta)
      gterm_d   = ter_g_d(pot,costheta) !still need d(costheta)/dR term
      expterm   = exp((lambda3*(drij-drik))**3)
      expterm_d = (3.0_dp)*(lambda3**3)*((drij-drik)**2)*expterm

      dri =  f_C_d * gterm   * expterm   * (rik) &
           + f_C   * gterm_d * expterm   * (dcosdri) &
           + f_C   * gterm   * expterm_d * (-rij_hat+rik_hat)

      !No f_C_d component for Rj
      drj =  f_C   * gterm_d * expterm   * (dcosdrj) &
           + f_C   * gterm   * expterm_d * (rij_hat)

      drk =  f_C_d * gterm   * expterm   * (-rik) &
           + f_C   * gterm_d * expterm   * (dcosdrk) &
           + f_C   * gterm   * expterm_d * (-rik_hat)

      f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + prefactor*dri(1)
      f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + prefactor*dri(2)
      f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + prefactor*dri(3)
  
      f_nonbond(1,jparticle) = f_nonbond(1,jparticle) + prefactor*drj(1)
      f_nonbond(2,jparticle) = f_nonbond(2,jparticle) + prefactor*drj(2)
      f_nonbond(3,jparticle) = f_nonbond(3,jparticle) + prefactor*drj(3)

      f_nonbond(1,kparticle) = f_nonbond(1,kparticle) + prefactor*drk(1)
      f_nonbond(2,kparticle) = f_nonbond(2,kparticle) + prefactor*drk(2)
      f_nonbond(3,kparticle) = f_nonbond(3,kparticle) + prefactor*drk(3)

      ptens_nonbond(1,1) = ptens_nonbond(1,1) - prefactor*(rij(1)*drj(1) + rik(1)*drk(1))
      ptens_nonbond(1,2) = ptens_nonbond(1,2) - prefactor*(rij(1)*drj(2) + rik(1)*drk(2))
      ptens_nonbond(1,3) = ptens_nonbond(1,3) - prefactor*(rij(1)*drj(3) + rik(1)*drk(3))

      ptens_nonbond(2,1) = ptens_nonbond(2,1) - prefactor*(rij(2)*drj(1) + rik(2)*drk(1))
      ptens_nonbond(2,2) = ptens_nonbond(2,2) - prefactor*(rij(2)*drj(2) + rik(2)*drk(2))
      ptens_nonbond(2,3) = ptens_nonbond(2,3) - prefactor*(rij(2)*drj(3) + rik(2)*drk(3))

      ptens_nonbond(3,1) = ptens_nonbond(3,1) - prefactor*(rij(3)*drj(1) + rik(3)*drk(1))
      ptens_nonbond(3,2) = ptens_nonbond(3,2) - prefactor*(rij(3)*drj(2) + rik(3)*drk(2))
      ptens_nonbond(3,3) = ptens_nonbond(3,3) - prefactor*(rij(3)*drj(3) + rik(3)*drk(3))

    END DO
    !CALL timestop ( handle )
  END SUBROUTINE ter_zeta_ij_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ter_g(pot,costheta)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: costheta
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_g', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_g
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: c2,d2
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    c = pot%tersoff%c
    d = pot%tersoff%d
    h = pot%tersoff%h
    c2 = c*c
    d2 = d*d
    ter_g = 0.0_dp
    ter_g = 1.0_dp + (c2/d2) - (c2)/(d2+(h-costheta)**2)
    !CALL timestop ( handle )
  END FUNCTION ter_g

  FUNCTION ter_g_d(pot,costheta)
    TYPE(pair_potential_single_type), POINTER           :: pot
    REAL(KIND=dp), INTENT(IN)                           :: costheta
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_g_d', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: ter_g_d
    REAL(KIND=dp)                                       :: A,B
    REAL(KIND=dp)                                       :: lambda1,lambda2,lambda3
    REAL(KIND=dp)                                       :: alpha,beta
    REAL(KIND=dp)                                       :: n,c,d,h
    REAL(KIND=dp)                                       :: c2,d2,hc
    REAL(KIND=dp)                                       :: sintheta
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    c = pot%tersoff%c
    d = pot%tersoff%d
    h = pot%tersoff%h
    c2 = c*c
    d2 = d*d
    hc = h-costheta

    sintheta = SQRT(1.0 - costheta**2)

    ter_g_d = 0.0_dp

    ! Still need d(theta)/dR
    !ter_g_d = 2.0_dp * c**2 * sintheta * (h - costheta) / &
    !          ( d**2 + (h-costheta)**2 )**2

    ! Still need d(costheta)/dR
    ter_g_d = (-2.0_dp * c2 * hc)/(d2 + hc**2)**2
    !CALL timestop ( handle )
  END FUNCTION ter_g_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ter_costheta_d(particle_set,cell,iparticle,jparticle,kparticle,dri,drj,drk)
    TYPE(pair_potential_single_type), POINTER           :: pot
    TYPE(particle_type), INTENT(IN)                     :: particle_set( : )
    TYPE(cell_type), POINTER                            :: cell
    INTEGER, INTENT(IN)                                 :: iparticle,jparticle,kparticle
    REAL(KIND=dp), INTENT(OUT)                          :: dri(3),drj(3),drk(3)
    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_costheta_d', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                                       :: rij(3),drij,rij_hat(3)
    REAL(KIND=dp)                                       :: rik(3),drik,rik_hat(3)
    REAL(KIND=dp)                                       :: costheta
    INTEGER                                             :: handle

    !CALL timeset ( routineN, handle )

    !Note that the center atom is at atom i not the usual atom j

    rij(:) = pbc(particle_set(jparticle)%r(:),particle_set(iparticle)%r(:),cell)
    drij = SQRT(DOT_PRODUCT(rij,rij))
    rij_hat(:) = rij(:)/drij

    rik(:) = pbc(particle_set(kparticle)%r(:),particle_set(iparticle)%r(:),cell)
    drik = SQRT(DOT_PRODUCT(rik,rik))
    rik_hat(:) = rik(:)/drik

    costheta = DOT_PRODUCT(rij,rik) / (drij*drik)
    IF(costheta<-1.0_dp) costheta=-1.0_dp
    IF(costheta>+1.0_dp) costheta=+1.0_dp

    drj(:) = (1.0_dp/(drij)) * ( rik_hat(:) - costheta*rij_hat(:) )
    drk(:) = (1.0_dp/(drik)) * ( rij_hat(:) - costheta*rik_hat(:) )
    dri(:) = - ( drj(:) + drk(:) )
    !CALL timestop ( handle )
  END SUBROUTINE ter_costheta_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !******************************************************************************
  SUBROUTINE get_rho_eam ( potparm, rab2, atom_a, atom_b, rho )
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    REAL(dp), INTENT(IN)                     :: rab2
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    REAL(dp), INTENT(INOUT)                  :: rho( : )

    INTEGER                                  :: index
    REAL(dp)                                 :: qq, rab, rhoi, rhoj

! Computation of rhotot

    rab = SQRT ( rab2 )
    index = INT ( rab / potparm % eam % drar ) + 1

    IF ( index > 5000 ) THEN
       index = 5000
    ELSEIF ( index < 1 ) THEN
       index = 1
    ENDIF

    qq = rab - potparm % eam % rval ( index )

    rhoi = potparm % eam % rho ( index ) + qq * potparm % eam % rhop ( index )
    rhoj = potparm % eam % rho ( index ) + qq * potparm % eam % rhop ( index )

    rho ( atom_a ) = rho ( atom_a ) + rhoi
    rho ( atom_b ) = rho ( atom_b ) + rhoj
  END SUBROUTINE get_rho_eam
  !******************************************************************************
  SUBROUTINE get_force_eam ( rab2, potparm, particle_set, atom_a, atom_b, f_eam )
    REAL(dp), INTENT(IN)                     :: rab2
    TYPE(pair_potential_single_type), &
      POINTER                                :: potparm
    TYPE(particle_type), INTENT(IN)          :: particle_set( : )
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    REAL(dp), INTENT(OUT)                    :: f_eam

    INTEGER                                  :: index
    REAL(dp)                                 :: denspi, denspj, fcp, qq, rab

! compute scaler part of the force

    f_eam = 0._dp
    rab = SQRT ( rab2 )
    index = INT ( rab / potparm % eam % drar ) + 1

    IF ( index > 5000 ) THEN
       index = 5000
    ELSEIF ( index < 1 ) THEN
       index = 1
    ENDIF
    qq = rab - potparm % eam % rval ( index )

    denspi = potparm % eam % rhop ( index ) +  &
         qq * ( potparm % eam % rhop ( index + 1 ) - &
         potparm % eam % rhop ( index ) ) / potparm % eam % drar
    denspj = potparm % eam % rhop ( index ) + &
         qq *( potparm % eam % rhop ( index + 1 ) - &
         potparm % eam % rhop ( index ) ) / potparm % eam % drar

    fcp = denspj * particle_set ( atom_a ) % f_embed + &
         denspi * particle_set ( atom_b ) % f_embed
    f_eam =  fcp / rab
  END SUBROUTINE get_force_eam

  !******************************************************************************

  SUBROUTINE get_force_tersoff(pot,particle_set,cell,iparticle,jparticle,f_nonbond,&
                               ptens_nonbond,newlist,nl)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(particle_type), INTENT(IN)          :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond,ptens_nonbond
    INTEGER, POINTER                         :: newlist(:,:),nl(:)
    REAL(KIND=dp)                            :: rij(3),drij,rij_hat(3)
    REAL(KIND=dp)                            :: RpD,RmD
    REAL(KIND=dp)                            :: f_R1, f_R2
    REAL(KIND=dp)                            :: f_A1, f_A2
    REAL(KIND=dp)                            :: b_ij,b_ij_d
    REAL(KIND=dp)                            :: b_ji,b_ji_d
    REAL(KIND=dp)                            :: f_A,f_A_d
    REAL(KIND=dp)                            :: f_C,f_C_d
    REAL(KIND=dp)                            :: f_R,f_R_d
    REAL(KIND=dp)                            :: prefactor



    rij(:) = pbc(particle_set(jparticle)%r(:),particle_set(iparticle)%r(:),cell)
    drij = SQRT(DOT_PRODUCT(rij,rij))
    rij_hat(:) = rij(:) / drij 

    RpD = pot%tersoff%bigR + pot%tersoff%bigD
    IF(drij>RpD) RETURN

    b_ij  = 0.5_dp * ter_b_ij(pot,particle_set,cell,iparticle,jparticle,newlist,nl)
    b_ij_d= 0.5_dp * ter_b_ij_d(pot,particle_set,cell,iparticle,jparticle,newlist,nl)
    b_ji  = 0.5_dp * ter_b_ij(pot,particle_set,cell,jparticle,iparticle,newlist,nl)
    b_ji_d= 0.5_dp * ter_b_ij_d(pot,particle_set,cell,jparticle,iparticle,newlist,nl)
    f_A   = ter_f_A(pot,drij)
    f_A_d = ter_f_A_d(pot,drij)
    f_C   = ter_f_C(pot,drij)
    f_C_d = ter_f_C_d(pot,drij)
    f_R   = ter_f_R(pot,drij)
    f_R_d = ter_f_R_d(pot,drij)

    ! Lets do the easy one first, the repulsive term
    ! Note a_ij = 1.0_dp so just going to ignore it...
    f_R1 = f_C_d * f_R
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_R1*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_R1*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_R1*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_R1*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_R1*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_R1*rij(3)

    ptens_nonbond(1,1) = ptens_nonbond(1,1) + f_R1*rij(1)*rij(1)
    ptens_nonbond(1,2) = ptens_nonbond(1,2) + f_R1*rij(1)*rij(2)
    ptens_nonbond(1,3) = ptens_nonbond(1,3) + f_R1*rij(1)*rij(3)
    ptens_nonbond(2,1) = ptens_nonbond(2,1) + f_R1*rij(2)*rij(1)
    ptens_nonbond(2,2) = ptens_nonbond(2,2) + f_R1*rij(2)*rij(2)
    ptens_nonbond(2,3) = ptens_nonbond(2,3) + f_R1*rij(2)*rij(3)
    ptens_nonbond(3,1) = ptens_nonbond(3,1) + f_R1*rij(3)*rij(1)
    ptens_nonbond(3,2) = ptens_nonbond(3,2) + f_R1*rij(3)*rij(2)
    ptens_nonbond(3,3) = ptens_nonbond(3,3) + f_R1*rij(3)*rij(3)

    f_R2 = f_C * f_R_d
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_R2*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_R2*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_R2*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_R2*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_R2*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_R2*rij(3)

    ptens_nonbond(1,1) = ptens_nonbond(1,1) + f_R2*rij(1)*rij(1)
    ptens_nonbond(1,2) = ptens_nonbond(1,2) + f_R2*rij(1)*rij(2)
    ptens_nonbond(1,3) = ptens_nonbond(1,3) + f_R2*rij(1)*rij(3)
    ptens_nonbond(2,1) = ptens_nonbond(2,1) + f_R2*rij(2)*rij(1)
    ptens_nonbond(2,2) = ptens_nonbond(2,2) + f_R2*rij(2)*rij(2)
    ptens_nonbond(2,3) = ptens_nonbond(2,3) + f_R2*rij(2)*rij(3)
    ptens_nonbond(3,1) = ptens_nonbond(3,1) + f_R2*rij(3)*rij(1)
    ptens_nonbond(3,2) = ptens_nonbond(3,2) + f_R2*rij(3)*rij(2)
    ptens_nonbond(3,3) = ptens_nonbond(3,3) + f_R2*rij(3)*rij(3)

    ! Lets do the f_A1 piece derivative of F_C
    f_A1 = f_C_d * (b_ij+b_ji) * f_A
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_A1*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_A1*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_A1*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_A1*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_A1*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_A1*rij(3)

    ptens_nonbond(1,1) = ptens_nonbond(1,1) + f_A1*rij(1)*rij(1)
    ptens_nonbond(1,2) = ptens_nonbond(1,2) + f_A1*rij(1)*rij(2)
    ptens_nonbond(1,3) = ptens_nonbond(1,3) + f_A1*rij(1)*rij(3)
    ptens_nonbond(2,1) = ptens_nonbond(2,1) + f_A1*rij(2)*rij(1)
    ptens_nonbond(2,2) = ptens_nonbond(2,2) + f_A1*rij(2)*rij(2)
    ptens_nonbond(2,3) = ptens_nonbond(2,3) + f_A1*rij(2)*rij(3)
    ptens_nonbond(3,1) = ptens_nonbond(3,1) + f_A1*rij(3)*rij(1)
    ptens_nonbond(3,2) = ptens_nonbond(3,2) + f_A1*rij(3)*rij(2)
    ptens_nonbond(3,3) = ptens_nonbond(3,3) + f_A1*rij(3)*rij(3)

    ! Lets do the f_A2 piece derivative of F_A
    f_A2 = f_C * (b_ij+b_ji) * f_A_d
    f_nonbond(1,iparticle) = f_nonbond(1,iparticle) + f_A2*rij(1)
    f_nonbond(2,iparticle) = f_nonbond(2,iparticle) + f_A2*rij(2)
    f_nonbond(3,iparticle) = f_nonbond(3,iparticle) + f_A2*rij(3)
    f_nonbond(1,jparticle) = f_nonbond(1,jparticle) - f_A2*rij(1)
    f_nonbond(2,jparticle) = f_nonbond(2,jparticle) - f_A2*rij(2)
    f_nonbond(3,jparticle) = f_nonbond(3,jparticle) - f_A2*rij(3)

    ptens_nonbond(1,1) = ptens_nonbond(1,1) + f_A2*rij(1)*rij(1)
    ptens_nonbond(1,2) = ptens_nonbond(1,2) + f_A2*rij(1)*rij(2)
    ptens_nonbond(1,3) = ptens_nonbond(1,3) + f_A2*rij(1)*rij(3)
    ptens_nonbond(2,1) = ptens_nonbond(2,1) + f_A2*rij(2)*rij(1)
    ptens_nonbond(2,2) = ptens_nonbond(2,2) + f_A2*rij(2)*rij(2)
    ptens_nonbond(2,3) = ptens_nonbond(2,3) + f_A2*rij(2)*rij(3)
    ptens_nonbond(3,1) = ptens_nonbond(3,1) + f_A2*rij(3)*rij(1)
    ptens_nonbond(3,2) = ptens_nonbond(3,2) + f_A2*rij(3)*rij(2)
    ptens_nonbond(3,3) = ptens_nonbond(3,3) + f_A2*rij(3)*rij(3)

    ! Lets do the f_A3 piece derivative of b_ij
    prefactor = f_C * b_ij_d * f_A    ! Note need to do d(Zeta_ij)/dR
    CALL ter_zeta_ij_d(pot,particle_set,cell,iparticle,jparticle,f_nonbond,ptens_nonbond,prefactor,newlist,nl)
    prefactor = f_C * b_ji_d * f_A    ! Note need to do d(Zeta_ij)/dR
    CALL ter_zeta_ij_d(pot,particle_set,cell,jparticle,iparticle,f_nonbond,ptens_nonbond,prefactor,newlist,nl)

  END SUBROUTINE get_force_tersoff

  !******************************************************************************

  SUBROUTINE force_nonbond_manybody ( fist_nonbond_env, particle_set, cell,  &
       pot_nonbond, f_nonbond, ptens_nonbond, error )

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: pot_nonbond
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond, ptens_nonbond
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_nonbond_manybody', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, handle, &
                                                ikind, ilist, inode, ipair, &
                                                jkind, nkinds, nlist, nnode, &
                                                npairs
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: is_eam, is_tersoff, lnl
    REAL(KIND=dp) :: energy, f_eam, fr(3), fr_com(3), fscalar, ptens11, &
      ptens12, ptens13, ptens21, ptens22, ptens23, ptens31, ptens32, ptens33, &
      rab(3), rab2, rab2_max
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot

    INTEGER, POINTER                         :: newlist(:,:),nl(:)
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc
    TYPE(spline_data_type), POINTER          :: spline_data
    REAL(KIND=dp), DIMENSION(3)              :: cvi,cell_v
    

    CALL timeset ( routineN, handle )

    lnl = .FALSE.

    CALL fist_nonbond_env_get(fist_nonbond_env,nonbonded=nonbonded,potparm=potparm,&
              natom_types=nkinds,r_last_update_pbc=r_last_update_pbc,error=error)

    ! initializing the potential energy, pressure tensor and force
    ptens11 = 0.0_dp ; ptens12 = 0.0_dp ; ptens13 = 0.0_dp
    ptens21 = 0.0_dp ; ptens22 = 0.0_dp ; ptens23 = 0.0_dp
    ptens31 = 0.0_dp ; ptens32 = 0.0_dp ; ptens33 = 0.0_dp
    ! starting the force loop
! incorrect omp (i.e. threads can write to the variable an other thread might be reading (e.g. pot_nonbonded))
! !$OMP parallel do default(none) &
! !$OMP shared(nkinds) &
! !$OMP private(ikind,jkind) &
! !$OMP shared(potparm) &
! !$OMP private(is_eam) &
! !$OMP shared(nonbonded) &
! !$OMP private(rab,rab2_max) &
! !$OMP private(ilist,nlist,atom_a,atom_b,npairs,list,inode,nnode) &
! !$OMP shared(cell) &
! !$OMP shared(pot_nonbond) &
! !$OMP private(f_eam,pot,fr,fscalar,energy) &
! !$OMP shared(f_nonbond,ptens_nonbond,r_last_update) &
! !$OMP private(rab2) &
! !$OMP private(neighbor_kind_pair) &
! !$OMP shared(particle_set,r_last_update_pbc) &
! !$OMP reduction(+:ptens11,ptens12,ptens13,ptens21,ptens22,ptens23, &
! !$OMP             ptens31,ptens32,ptens33) &
! !$OMP private(ipair) 

    DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs=neighbor_kind_pair%npairs
       IF (npairs ==0) CYCLE
       ikind = neighbor_kind_pair%ikind
       jkind = neighbor_kind_pair%jkind
       list  =>neighbor_kind_pair%list
       cvi   = neighbor_kind_pair%cell_vector

       is_eam=(potparm %pot( ikind, jkind ) %pot % TYPE == ea_type)
       is_tersoff=(potparm %pot( ikind, jkind ) %pot % TYPE == tersoff_type)

       IF((.NOT.is_eam).AND.(.NOT.is_tersoff)) CYCLE

       IF((.NOT.lnl).AND.(is_tersoff)) THEN
         CALL get_newlist(fist_nonbond_env,newlist,nl,SIZE(particle_set),error)
         lnl = .TRUE.
       END IF

       spline_data => potparm % pot( ikind, jkind ) %pot % pair_spline_data
       rab2_max = potparm%pot(ikind,jkind)%pot%rcutsq
       CALL matvec_3x3(cell_v, cell%hmat,cvi)

       DO ipair = 1, npairs
          atom_a = list(1,ipair)
          atom_b = list(2,ipair)

          !set this outside the potential type in case need multiple potentials
          fr(:) = 0.0_dp

          !Do everything necessary for EAM here
          IF (is_eam) THEN
            rab = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
            rab = rab+cell_v
            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

            IF (rab2 <= rab2_max) THEN
              energy =  potential_s (spline_data,rab2,fscalar)
              pot_nonbond = pot_nonbond + energy
              fr(1) = fscalar*rab(1)
              fr(2) = fscalar*rab(2)
              fr(3) = fscalar*rab(3)

              pot => potparm %pot( ikind, jkind )%pot
              CALL get_force_eam ( rab2, pot, particle_set, atom_a, atom_b, f_eam )
              fr(1) =  fr(1) - f_eam * rab ( 1 )
              fr(2) =  fr(2) - f_eam * rab ( 2 )
              fr(3) =  fr(3) - f_eam * rab ( 3 )

              f_nonbond(1,atom_a) = f_nonbond(1,atom_a) - fr(1)
              f_nonbond(2,atom_a) = f_nonbond(2,atom_a) - fr(2)
              f_nonbond(3,atom_a) = f_nonbond(3,atom_a) - fr(3)

              f_nonbond(1,atom_b) = f_nonbond(1,atom_b) + fr(1)
              f_nonbond(2,atom_b) = f_nonbond(2,atom_b) + fr(2)
              f_nonbond(3,atom_b) = f_nonbond(3,atom_b) + fr(3)
              ptens11 = ptens11 + rab(1)*fr(1)
              ptens21 = ptens21 + rab(2)*fr(1)
              ptens31 = ptens31 + rab(3)*fr(1)
              ptens12 = ptens12 + rab(1)*fr(2)
              ptens22 = ptens22 + rab(2)*fr(2)
              ptens32 = ptens32 + rab(3)*fr(2)
              ptens13 = ptens13 + rab(1)*fr(3)
              ptens23 = ptens23 + rab(2)*fr(3)
              ptens33 = ptens33 + rab(3)*fr(3)
            ENDIF
          END IF
          !Do everything necessary for Tersoff here
          IF (is_tersoff) THEN
            pot => potparm %pot( ikind, jkind )%pot
            CALL get_force_tersoff ( pot, particle_set, cell, &
                                     atom_a, atom_b, f_nonbond, ptens_nonbond, &
                                     newlist,nl )
          ENDIF
       END DO
    END DO
    ptens_nonbond(1,1) = ptens_nonbond(1,1) + ptens11
    ptens_nonbond(1,2) = ptens_nonbond(1,2) + ptens12
    ptens_nonbond(1,3) = ptens_nonbond(1,3) + ptens13
    ptens_nonbond(2,1) = ptens_nonbond(2,1) + ptens21
    ptens_nonbond(2,2) = ptens_nonbond(2,2) + ptens22
    ptens_nonbond(2,3) = ptens_nonbond(2,3) + ptens23
    ptens_nonbond(3,1) = ptens_nonbond(3,1) + ptens31
    ptens_nonbond(3,2) = ptens_nonbond(3,2) + ptens32
    ptens_nonbond(3,3) = ptens_nonbond(3,3) + ptens33

    IF(lnl) THEN
      CALL destroy_newlist(newlist,nl,error)
      lnl = .TRUE.
    END IF

    CALL timestop ( handle )

  END SUBROUTINE force_nonbond_manybody

  !******************************************************************************


END MODULE manybody_potential

