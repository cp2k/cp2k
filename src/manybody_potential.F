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
!!     CJM, I-Feng W. Kuo, Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     Efficient tersoff implementation and general "lifting" of manybody_potential module
!!
!!   SOURCE
!******************************************************************************

MODULE manybody_potential

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
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
                                             pair_potential_pp_type,&
                                             pair_potential_single_type,&
                                             tersoff_type
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
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_force_eam
  PUBLIC :: density_nonbond, energy_manybody
  PUBLIC :: force_nonbond_manybody
  PUBLIC :: rho_f, rho_s
  PUBLIC :: spline_manybody_control
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_potential'

CONTAINS

!!****f* manybody_potential/manybody_pot [1.0] *
!!
!!   NAME
!!     manybody_pot
!!
!!   FUNCTION
!!     Implements possible potential expression for manybody potential
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo, Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION manybody_pot ( pot, r )

    TYPE(pair_potential_single_type), &
      INTENT(IN)                             :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: manybody_pot

    IF (pot%type==ea_type) THEN
       ! No splining necessary for EAM.  Zero splines
       manybody_pot = 0._dp
    ELSEIF (pot%type==tersoff_type) THEN
       CALL stop_program("manybody_potential%ener_pot",&
            "We should never reach this point when doing Tersoff!")
    ELSE
       CALL stop_program("manybody_potential%ener_pot","Potential type unknown")
    END IF

  END FUNCTION manybody_pot

!!****f* manybody_potential/spline_manybody_control [1.0] *
!!
!!   NAME
!!     spline_manybody_control
!!
!!   FUNCTION
!!     creates the splines for the manybody contribution
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo, Teodoro Laino [2006-2007]
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
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

!!****f* manybody_potential/rho_s [1.0] *
!!
!!   NAME
!!     rho_s
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE rho_s ( spline_table, rijsq, energy, fscalar )


    TYPE(spline_data_type), POINTER          :: spline_table
    REAL(KIND=dp), INTENT(IN)                :: rijsq
    REAL(KIND=dp), INTENT(OUT)               :: energy, fscalar

    energy = potential_s ( spline_table, rijsq, fscalar )

  END SUBROUTINE rho_s

!!****f* manybody_potential/rho_f [1.0] *
!!
!!   NAME
!!     rho_f
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE rho_f ( rijsq, potparm, iatom, jatom, energy )

    REAL(KIND=dp), INTENT(IN)                :: rijsq
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), INTENT(OUT)               :: energy

    REAL(KIND=dp)                            :: rij

    rij = SQRT ( rijsq )

    ! non-bond energy
    energy = manybody_pot ( potparm %pot( iatom, jatom )%pot, rij )

  END SUBROUTINE rho_f

!!****f* manybody_potential/density_nonbond [1.0] *
!!
!!   NAME
!!     density_nonbond
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
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

!!****f* manybody_potential/density_nonbond [1.0] *
!!
!!   NAME
!!     density_nonbond
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     CJM, I-Feng W. Kuo, Teodoro Laino 
!!
!!   MODIFICATION HISTORY
!!     tlaino [2007] - New algorithm for tersoff potential
!!
!!*** **********************************************************************
  SUBROUTINE energy_manybody ( fist_nonbond_env, atomic_kind_set, &
       local_particles, particle_set, cell, pot_manybody, para_env, error )

    ! computes the embedding contribution to the energy

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), INTENT(INOUT)                  :: pot_manybody
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'energy_manybody', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, ikind, ilist, index, ipair, iparticle, &
      iparticle_local, istat, iunique, jkind, jlist, junique, m, nloc_size, &
      npairs, npairs_2, nparticle, nparticle_local, nunique, stat
    INTEGER, DIMENSION(:), POINTER           :: work_list
    INTEGER, DIMENSION(:, :), POINTER        :: full_loc_list, list, sort_list
    LOGICAL                                  :: failure, is_tersoff
    REAL(dp)                                 :: b_ij, embed, f_A, f_C, f_R, &
                                                pot_loc, qr, rab2_max
    REAL(dp), POINTER                        :: fembed( : )
    REAL(KIND=dp)                            :: drij, rij(3)
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, cell_v_2, cvi
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: loc_cell_v
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair, &
                                                neighbor_kind_pair_2
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc

    CALL timeset ( routineN, handle )
    CALL fist_nonbond_env_get ( fist_nonbond_env, r_last_update_pbc=r_last_update_pbc,&
         potparm = potparm ,error=error)
    ! EAM requires a single loop 
    DO ikind = 1, SIZE ( atomic_kind_set )
       pot => potparm %pot ( ikind, ikind ) % pot
       SELECT CASE ( pot % TYPE )
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
       END SELECT
    END DO
    ! Other manybody potential
    CALL fist_nonbond_env_get(fist_nonbond_env,nonbonded=nonbonded,error=error)
    DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs=neighbor_kind_pair%npairs
       IF (npairs ==0) CYCLE
       ikind = neighbor_kind_pair%ikind
       jkind = neighbor_kind_pair%jkind
       list  =>neighbor_kind_pair%list
       cvi   = neighbor_kind_pair%cell_vector
       pot => potparm %pot ( ikind, jkind ) % pot
       is_tersoff=(pot % TYPE == tersoff_type)
       IF(.NOT.is_tersoff) CYCLE

       rab2_max = pot%rcutsq
       CALL matvec_3x3(cell_v, cell%hmat,cvi)
       pot => potparm %pot ( ikind, jkind ) % pot

       IF (npairs /=0) THEN
          ALLOCATE(sort_list(2,npairs),&
               work_list (npairs),stat=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
               "allocate sort_list")       
          sort_list = list(:,1:npairs)
          ! Sort the list of neighbors, this increases the efficiency for single 
          ! potential contributions
          CALL sort(sort_list(1,:),npairs,work_list)
          DO ipair = 1, npairs
             work_list(ipair)=sort_list(2,work_list(ipair))
          END DO
          sort_list(2,:) = work_list
          ! find number of unique elements of array index 1
          nunique = 1
          DO ipair = 1, npairs-1
             IF (sort_list(1,ipair+1)/=sort_list(1,ipair)) nunique = nunique + 1
          END DO
          ipair   = 1
          junique = sort_list(1,ipair)
          DO iunique = 1, nunique
             atom_a = junique
             ! Create a local list (atom_a local only) of all possible neighbors
             ! i.e. kind independent
             ! (1) First get array size
             nloc_size = 0
             DO jlist=1,nonbonded%nlists
                neighbor_kind_pair_2 => nonbonded%neighbor_kind_pairs(jlist)
                npairs_2=neighbor_kind_pair_2%npairs
                IF ((npairs_2 /=0).AND.(neighbor_kind_pair_2%ikind==ikind)) THEN
                   nloc_size = nloc_size + COUNT(neighbor_kind_pair_2%list(1,1:npairs_2)==atom_a)
                END IF
             END DO
             ! (2) Allocate Array and fill with values
             ALLOCATE(full_loc_list(2,nloc_size),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(loc_cell_v(3,nloc_size),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             nloc_size = 0
             DO jlist=1,nonbonded%nlists
                neighbor_kind_pair_2 => nonbonded%neighbor_kind_pairs(jlist)
                npairs_2=neighbor_kind_pair_2%npairs
                IF ((npairs_2 /=0).AND.(neighbor_kind_pair_2%ikind==ikind)) THEN
                   CALL matvec_3x3(cell_v_2, cell%hmat, REAL(neighbor_kind_pair_2%cell_vector,KIND=dp))
                   DO m = 1, npairs_2
                      IF (neighbor_kind_pair_2%list(1,m)==atom_a) THEN
                         nloc_size = nloc_size + 1
                         full_loc_list(1:2,nloc_size) = neighbor_kind_pair_2%list(1:2,m)
                         loc_cell_v   (1:3,nloc_size) = cell_v_2
                      END IF
                   END DO
                END IF
             END DO
             DO WHILE (ipair<=npairs)
                IF (sort_list(1,ipair)/=junique) EXIT
                atom_b = sort_list(2,ipair)
                ! Energy terms
                pot_loc = 0.0_dp
                rij(:) = r_last_update_pbc(atom_b)%r(:)-r_last_update_pbc(atom_a)%r(:)+cell_v
                drij = DOT_PRODUCT(rij,rij)
                ipair = ipair + 1
                IF(drij>rab2_max) CYCLE
                drij = SQRT(drij)
                b_ij   = ter_b_ij(pot,r_last_update_pbc,atom_a,atom_b,full_loc_list,loc_cell_v,cell_v)
                f_C    = ter_f_C (pot,drij)
                f_A    = ter_f_A (pot,drij)
                f_R    = ter_f_R (pot,drij)
                pot_loc = f_C*(f_R + b_ij*f_A)
                pot_manybody = pot_manybody + 0.5_dp*pot_loc
             END DO
             IF (ipair<=npairs) junique = sort_list(1,ipair)
             ! Destroy local list
             DEALLOCATE(full_loc_list,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(loc_cell_v,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END DO
          DEALLOCATE(sort_list,work_list,stat=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
               "deallocate sort_list")
       END IF
    END DO
    CALL timestop ( handle )
  END SUBROUTINE energy_manybody

!!****f* manybody_potential/ter_f_C [1.0] *
!!
!!   NAME
!!     ter_f_C
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_f_C(pot,r)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_C

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_C', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: bigD, bigR, RmD, RpD

    bigR = pot%tersoff%bigR
    bigD = pot%tersoff%bigD
    RmD  = pot%tersoff%bigR-pot%tersoff%bigD
    RpD  = pot%tersoff%bigR+pot%tersoff%bigD
    ter_f_C = 0.0_dp
    IF(r < RmD) ter_f_C = 1.0_dp
    IF(r > RpD) ter_f_C = 0.0_dp
    IF( (r < RpD) .AND. (r > RmD) ) THEN
      ter_f_C = 0.5_dp *(1.0_dp - SIN(0.5_dp*PI*(r-bigR)/(bigD)))
    END IF
  END FUNCTION ter_f_C

!!****f* manybody_potential/ter_f_C_d [1.0] *
!!
!!   NAME
!!     ter_f_C_d
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_f_C_d(pot,r)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_C_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_C_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: bigD, bigR, RmD, RpD

    bigR = pot%tersoff%bigR
    bigD = pot%tersoff%bigD
    RmD  = pot%tersoff%bigR-pot%tersoff%bigD
    RpD  = pot%tersoff%bigR+pot%tersoff%bigD
    ter_f_C_d = 0.0_dp
    IF(r < RmD) ter_f_C_d = 0.0_dp
    IF(r > RpD) ter_f_C_d = 0.0_dp
    IF( (r < RpD) .AND. (r > RmD) ) THEN
      ter_f_C_d = (0.25_dp*PI/bigD) * COS(0.5_dp*PI*(r-bigR)/(bigD)) / r
    END IF

  END FUNCTION ter_f_C_d

!!****f* manybody_potential/ter_f_R [1.0] *
!!
!!   NAME
!!     ter_f_R
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_f_R(pot,r)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_R

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_R', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: A, lambda1

    A       = pot%tersoff%A
    lambda1 = pot%tersoff%lambda1
    ter_f_R = 0.0_dp
    ter_f_R = A * EXP(-lambda1*r)

  END FUNCTION ter_f_R

!!****f* manybody_potential/ter_f_R_d [1.0] *
!!
!!   NAME
!!     ter_f_R_d
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_f_R_d(pot,r)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_R_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_R_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: A, f_R, lambda1

    A = pot%tersoff%A
    lambda1 = pot%tersoff%lambda1
    f_R = A * EXP(-lambda1*r)
    ter_f_R_d = 0.0_dp
    ter_f_R_d = lambda1 * f_R / r

  END FUNCTION ter_f_R_d

!!****f* manybody_potential/ter_f_A [1.0] *
!!
!!   NAME
!!     ter_f_A
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_f_A(pot,r)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_A

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_A', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: B, lambda2

    B = pot%tersoff%B
    lambda2 = pot%tersoff%lambda2
    ter_f_A = 0.0_dp
    ter_f_A = - B * EXP(-lambda2*r)

  END FUNCTION ter_f_A

!!****f* manybody_potential/ter_f_A_d [1.0] *
!!
!!   NAME
!!     ter_f_A_d
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_f_A_d(pot,r)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_f_A_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_f_A_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: B, lambda2

    B = pot%tersoff%B
    lambda2 = pot%tersoff%lambda2
    ter_f_A_d = 0.0_dp
    ter_f_A_d = - B *lambda2 * EXP(-lambda2*r) / r

  END FUNCTION ter_f_A_d

!!****f* manybody_potential/ter_a_ij [1.0] *
!!
!!   NAME
!!     ter_a_ij
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_a_ij(pot,r)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: r
    REAL(KIND=dp)                            :: ter_a_ij

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_a_ij', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: alpha, n

    n = pot%tersoff%n
    alpha = pot%tersoff%alpha
    ter_a_ij = 0.0_dp
    !Note alpha = 0.0_dp for the parameters in the paper so using simplified term
    !ter_a_ij = (1.0_dp+(alpha*ter_n_ij(pot,iparticle,jparticle,r))**n)**(-0.5_dp/n)
    ter_a_ij = 1.0_dp

  END FUNCTION ter_a_ij

!!****f* manybody_potential/ter_b_ij [1.0] *
!!
!!   NAME
!!     ter_b_ij
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo, Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_b_ij(pot,r_last_update_pbc,iparticle,jparticle,full_loc_list, loc_cell_v, cell_v)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    INTEGER, DIMENSION(:, :), POINTER        :: full_loc_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    REAL(KIND=dp)                            :: ter_b_ij

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_b_ij', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: beta, n, zeta_ij

    n = pot%tersoff%n
    beta = pot%tersoff%beta
    ter_b_ij = 0.0_dp
    zeta_ij = ter_zeta_ij(pot,r_last_update_pbc,iparticle,jparticle,full_loc_list, loc_cell_v, cell_v)
    ter_b_ij = (1.0_dp+(beta*zeta_ij)**n)**(-0.5_dp/n)

  END FUNCTION ter_b_ij

!!****f* manybody_potential/ter_b_ij_d [1.0] *
!!
!!   NAME
!!     ter_b_ij_d
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo, Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_b_ij_d(pot,r_last_update_pbc,iparticle,jparticle,full_loc_list, loc_cell_v, cell_v)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    INTEGER, DIMENSION(:, :), POINTER        :: full_loc_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    REAL(KIND=dp)                            :: ter_b_ij_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_b_ij_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: beta, beta_n, n, zeta_ij, &
                                                zeta_ij_n, zeta_ij_nm1

    n           = pot%tersoff%n
    beta        = pot%tersoff%beta
    beta_n      = beta**n
    zeta_ij     = ter_zeta_ij(pot,r_last_update_pbc,iparticle,jparticle,full_loc_list, loc_cell_v, cell_v)
    zeta_ij_nm1 = zeta_ij**(n-1.0_dp)
    zeta_ij_n   = zeta_ij**(n)

    ter_b_ij_d = 0.0_dp
    ter_b_ij_d = -0.5_dp * beta_n * zeta_ij_nm1 * &
                   ((1.0_dp + beta_n * zeta_ij_n)**((-0.5_dp/n) - 1.0_dp))

  END FUNCTION ter_b_ij_d

!!****f* manybody_potential/ter_zeta_ij [1.0] *
!!
!!   NAME
!!     ter_zeta_ij
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo, Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     Using a local list of neighbors - [tlaino] 2007
!!
!!*** **********************************************************************
  FUNCTION ter_zeta_ij(pot,r_last_update_pbc,iparticle,jparticle,full_loc_list, loc_cell_v, cell_v)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    INTEGER, DIMENSION(:, :), POINTER        :: full_loc_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    REAL(KIND=dp)                            :: ter_zeta_ij

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_zeta_ij', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilist, kparticle
    REAL(KIND=dp)                            :: cell_v_2(3), costheta, drij, &
                                                drik, expterm, f_C, gterm, &
                                                lambda3, n, rab2_max, rij(3), &
                                                rik(3)

    ter_zeta_ij = 0.0_dp
    n = pot%tersoff%n
    lambda3 = pot%tersoff%lambda3
    rab2_max = pot%rcutsq
    rij(:) = r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v
    drij = SQRT(DOT_PRODUCT(rij,rij))
    ter_zeta_ij = 0.0_dp
    DO ilist = 1, SIZE(full_loc_list,2)
       kparticle = full_loc_list(2,ilist)
       IF (kparticle==jparticle) CYCLE
       cell_v_2 = loc_cell_v(:,ilist)
       rik(:) = r_last_update_pbc(kparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v_2
       drik = DOT_PRODUCT(rik,rik)
       IF(drik>rab2_max) CYCLE
       drik = SQRT(drik)
       costheta = DOT_PRODUCT(rij,rik) / (drij*drik)
       IF(costheta<-1.0_dp) costheta=-1.0_dp
       IF(costheta>+1.0_dp) costheta=+1.0_dp
       f_C = ter_f_C(pot,drik)
       gterm = ter_g(pot,costheta)
       expterm = EXP((lambda3*(drij-drik))**3)
       ter_zeta_ij = ter_zeta_ij + f_C*gterm*expterm
    END DO

  END FUNCTION ter_zeta_ij

!!****f* manybody_potential/ter_zeta_ij_d [1.0] *
!!
!!   NAME
!!     ter_zeta_ij_d
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo, Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!      Using a local list of neighbors - [tlaino] 2007
!!
!!*** **********************************************************************
  SUBROUTINE ter_zeta_ij_d(pot,r_last_update_pbc,iparticle,jparticle,f_nonbond,ptens_nonbond,prefactor,&
       full_loc_list, loc_cell_v, cell_v)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond, ptens_nonbond
    REAL(KIND=dp), INTENT(IN)                :: prefactor
    INTEGER, DIMENSION(:, :), POINTER        :: full_loc_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: loc_cell_v
    REAL(KIND=dp), DIMENSION(3)              :: cell_v

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_zeta_ij_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilist, kparticle, nparticle
    REAL(KIND=dp)                            :: costheta, drij, drik, &
                                                expterm, expterm_d, f_C, &
                                                f_C_d, gterm, gterm_d, &
                                                lambda3, n, rab2_max
    REAL(KIND=dp), DIMENSION(3)              :: cell_v_2, dcosdri, dcosdrj, &
                                                dcosdrk, dri, drj, drk, rij, &
                                                rij_hat, rik, rik_hat

    n = pot%tersoff%n
    lambda3 = pot%tersoff%lambda3
    rab2_max = pot%rcutsq

    rij(:) = r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v
    drij = SQRT(DOT_PRODUCT(rij,rij))
    rij_hat(:) = rij(:) / drij

    nparticle = SIZE ( r_last_update_pbc )
    DO ilist = 1, SIZE (full_loc_list,2)
       kparticle = full_loc_list(2,ilist)
       IF (kparticle==jparticle) CYCLE
       cell_v_2 = loc_cell_v(:,ilist)
       rik(:) = r_last_update_pbc(kparticle)%r(:)-r_last_update_pbc(iparticle)%r(:)+cell_v_2
       drik = DOT_PRODUCT(rik,rik)

       IF(drik>rab2_max) CYCLE
       drik = SQRT(drik)
       rik_hat(:) = rik(:) / drik
       costheta = DOT_PRODUCT(rij,rik) / (drij*drik)
       IF(costheta<-1.0_dp) costheta=-1.0_dp
       IF(costheta>+1.0_dp) costheta=+1.0_dp

       dcosdrj(:) = (1.0_dp/(drij)) * ( rik_hat(:) - costheta*rij_hat(:) )
       dcosdrk(:) = (1.0_dp/(drik)) * ( rij_hat(:) - costheta*rik_hat(:) )
       dcosdri(:) = - ( dcosdrj(:) + dcosdrk(:) )

       f_C       = ter_f_C(pot,drik)
       f_C_d     = ter_f_C_d(pot,drik)
       gterm     = ter_g(pot,costheta)
       gterm_d   = ter_g_d(pot,costheta) !still need d(costheta)/dR term
       expterm   = EXP((lambda3*(drij-drik))**3)
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
  END SUBROUTINE ter_zeta_ij_d

!!****f* manybody_potential/ter_g [1.0] *
!!
!!   NAME
!!     ter_g
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_g(pot,costheta)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: costheta
    REAL(KIND=dp)                            :: ter_g

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_g', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: c, c2, d, d2, h

    c = pot%tersoff%c
    d = pot%tersoff%d
    h = pot%tersoff%h
    c2 = c*c
    d2 = d*d
    ter_g = 0.0_dp
    ter_g = 1.0_dp + (c2/d2) - (c2)/(d2+(h-costheta)**2)

  END FUNCTION ter_g

!!****f* manybody_potential/ter_g_d [1.0] *
!!
!!   NAME
!!     ter_g_d
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  FUNCTION ter_g_d(pot,costheta)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: costheta
    REAL(KIND=dp)                            :: ter_g_d

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ter_g_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: c, c2, d, d2, h, hc, sintheta

    c = pot%tersoff%c
    d = pot%tersoff%d
    h = pot%tersoff%h
    c2 = c*c
    d2 = d*d
    hc = h-costheta

    sintheta = SQRT(1.0 - costheta**2)

    ter_g_d = 0.0_dp
    ! Still need d(costheta)/dR
    ter_g_d = (-2.0_dp * c2 * hc)/(d2 + hc**2)**2
  END FUNCTION ter_g_d

!!****f* manybody_potential/get_rho_eam [1.0] *
!!
!!   NAME
!!     get_rho_eam
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
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

!!****f* manybody_potential/get_force_eam [1.0] *
!!
!!   NAME
!!     get_force_eam
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
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

!!****f* manybody_potential/get_force_tersoff [1.0] *
!!
!!   NAME
!!     get_force_tersoff
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo, Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!      Using a local list of neighbors - [tlaino] 2007
!!
!!*** **********************************************************************
  SUBROUTINE get_force_tersoff(pot,r_last_update_pbc, cell_v, full_loc_list, loc_cell_v,&
       iparticle,jparticle,f_nonbond,ptens_nonbond)
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    REAL(KIND=dp), DIMENSION(3)              :: cell_v
    INTEGER, DIMENSION(:, :), POINTER        :: full_loc_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: loc_cell_v
    INTEGER, INTENT(IN)                      :: iparticle, jparticle
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond, ptens_nonbond

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_force_tersoff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp) :: b_ij, b_ij_d, drij, f_A, f_A1, f_A2, f_A_d, f_C, f_C_d, &
      f_R, f_R1, f_R2, f_R_d, fac, prefactor, rij(3), rij_hat(3)

    CALL timeset ( routineN, handle )
    rij(:) = r_last_update_pbc(jparticle)%r(:)-r_last_update_pbc(iparticle)%r(:) + cell_v
    drij = SQRT(DOT_PRODUCT(rij,rij))
    rij_hat(:) = rij(:) / drij 

    fac = - 0.5_dp
    b_ij  = ter_b_ij  (pot,r_last_update_pbc,iparticle,jparticle,full_loc_list, loc_cell_v,cell_v)
    b_ij_d= ter_b_ij_d(pot,r_last_update_pbc,iparticle,jparticle,full_loc_list, loc_cell_v,cell_v)
    f_A   = ter_f_A   (pot,drij)
    f_A_d = ter_f_A_d (pot,drij)
    f_C   = ter_f_C   (pot,drij)
    f_C_d = ter_f_C_d (pot,drij)
    f_R   = ter_f_R   (pot,drij)
    f_R_d = ter_f_R_d (pot,drij)

    ! Lets do the easy one first, the repulsive term
    ! Note a_ij = 1.0_dp so just going to ignore it...
    f_R1 = f_C_d * f_R * fac
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

    f_R2 = f_C * f_R_d * fac
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
    f_A1 = f_C_d * b_ij * f_A * fac
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
    f_A2 = f_C * b_ij * f_A_d * fac
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
    prefactor = f_C * b_ij_d * f_A * fac  ! Note need to do d(Zeta_ij)/dR
    CALL ter_zeta_ij_d(pot,r_last_update_pbc,iparticle,jparticle,f_nonbond,ptens_nonbond,prefactor,&
         full_loc_list, loc_cell_v,cell_v)

    CALL timestop ( handle )
  END SUBROUTINE get_force_tersoff

!!****f* manybody_potential/force_nonbond_manybody [1.0] *
!!
!!   NAME
!!     force_nonbond_manybody
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     I-Feng W. Kuo, Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     Fast implementation of the tersoff potential - [tlaino] 2007
!!
!!*** **********************************************************************
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

    INTEGER :: atom_a, atom_b, handle, ikind, ilist, ipair, istat, iunique, &
      jkind, jlist, junique, m, nkinds, nloc_size, npairs, npairs_2, nunique, &
      stat
    INTEGER, DIMENSION(:), POINTER           :: work_list
    INTEGER, DIMENSION(:, :), POINTER        :: full_loc_list, list, sort_list
    LOGICAL                                  :: failure, is_eam, is_tersoff
    REAL(KIND=dp) :: energy, f_eam, fr(3), fscalar, ptens11, ptens12, &
      ptens13, ptens21, ptens22, ptens23, ptens31, ptens32, ptens33, rab(3), &
      rab2, rab2_max, rtmp(3)
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, cell_v_2, cvi
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: loc_cell_v
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair, &
                                                neighbor_kind_pair_2
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update_pbc
    TYPE(spline_data_type), POINTER          :: spline_data

    failure = .FALSE.
    CALL timeset ( routineN, handle )

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

       spline_data => potparm % pot( ikind, jkind ) %pot % pair_spline_data
       rab2_max = potparm%pot(ikind,jkind)%pot%rcutsq
       CALL matvec_3x3(cell_v, cell%hmat,cvi)

       IF (is_eam) THEN
          DO ipair = 1, npairs
             atom_a = list(1,ipair)
             atom_b = list(2,ipair)
             !set this outside the potential type in case need multiple potentials
             fr(:) = 0.0_dp
             !Do everything necessary for EAM here
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
          END DO
       END IF
       !Do everything necessary for Tersoff here
       IF (is_tersoff) THEN
          IF (npairs /=0) THEN
             ALLOCATE(sort_list(2,npairs),&
                  work_list (npairs),stat=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                  "allocate sort_list")       
             sort_list = list(:,1:npairs)
             ! Sort the list of neighbors, this increases the efficiency for single 
             ! potential contributions
             CALL sort(sort_list(1,:),npairs,work_list)
             DO ipair = 1, npairs
                work_list(ipair)=sort_list(2,work_list(ipair))
             END DO
             sort_list(2,:) = work_list
             ! find number of unique elements of array index 1
             nunique = 1
             DO ipair = 1, npairs-1
                IF (sort_list(1,ipair+1)/=sort_list(1,ipair)) nunique = nunique + 1
             END DO
             ipair   = 1
             junique = sort_list(1,ipair)
             DO iunique = 1, nunique
                atom_a = junique
                ! Create a local list (atom_a local only) of all possible neighbors
                ! i.e. kind independent
                ! (1) First get array size
                nloc_size = 0
                DO jlist=1,nonbonded%nlists
                   neighbor_kind_pair_2 => nonbonded%neighbor_kind_pairs(jlist)
                   npairs_2=neighbor_kind_pair_2%npairs
                   IF ((npairs_2 /=0).AND.(neighbor_kind_pair_2%ikind==ikind)) THEN
                      nloc_size = nloc_size + COUNT(neighbor_kind_pair_2%list(1,1:npairs_2)==atom_a)
                   END IF
                END DO
                ! (2) Allocate Array and fill with values
                ALLOCATE(full_loc_list(2,nloc_size),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(loc_cell_v(3,nloc_size),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                nloc_size = 0
                DO jlist=1,nonbonded%nlists
                   neighbor_kind_pair_2 => nonbonded%neighbor_kind_pairs(jlist)
                   npairs_2=neighbor_kind_pair_2%npairs
                   IF ((npairs_2 /=0).AND.(neighbor_kind_pair_2%ikind==ikind)) THEN
                      CALL matvec_3x3(cell_v_2, cell%hmat, REAL(neighbor_kind_pair_2%cell_vector,KIND=dp))
                      DO m = 1, npairs_2
                         IF (neighbor_kind_pair_2%list(1,m)==atom_a) THEN
                            nloc_size = nloc_size + 1
                            full_loc_list(1:2,nloc_size) = neighbor_kind_pair_2%list(1:2,m)
                            loc_cell_v   (1:3,nloc_size) = cell_v_2
                         END IF
                      END DO
                   END IF
                END DO
                DO WHILE (ipair<=npairs)
                   IF (sort_list(1,ipair)/=junique) EXIT
                   atom_b = sort_list(2,ipair)
                   ! Derivative terms
                   rtmp = r_last_update_pbc(atom_b)%r(:)-r_last_update_pbc(atom_a)%r(:) + cell_v
                   pot => potparm %pot( ikind, jkind )%pot
                   ipair = ipair + 1
                   IF (DOT_PRODUCT(rtmp,rtmp)<=pot%rcutsq) THEN
                      CALL get_force_tersoff ( pot, r_last_update_pbc, cell_v, full_loc_list, loc_cell_v,&
                           atom_a, atom_b, f_nonbond, ptens_nonbond)
                   END IF
                END DO
                IF (ipair<=npairs) junique = sort_list(1,ipair)
                ! Destroy local list
                DEALLOCATE(full_loc_list,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(loc_cell_v,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END DO
             DEALLOCATE(sort_list,work_list,stat=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                  "deallocate sort_list")
          ENDIF
       END IF
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

    CALL timestop ( handle )
  END SUBROUTINE force_nonbond_manybody

END MODULE manybody_potential

