!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE nose
  
  USE constraint, ONLY : roll
  USE eigenvalueproblems, ONLY : diagonalise
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
!*apsi USE lapack, ONLY : rs
  USE molecule_types, ONLY : particle_node_type, molecule_structure_type
  USE mp, ONLY : mp_bcast, mp_sum
  USE stop_program, ONLY : stop_prg
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : gasdev
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: lnhc, lnhcp, lnhcpf, initial_npt, initial_nhc
  PUBLIC :: extended_parameters_type, nhc_info_type, npt_info_type
  PUBLIC :: yoshida_coef
  PUBLIC :: particle_kinetic_type, particle_scale_type
  
  TYPE npt_info_type
     REAL ( dbl ) :: eps
     REAL ( dbl ) :: v
     REAL ( dbl ) :: f
     REAL ( dbl ) :: mass
  END TYPE npt_info_type
  
  TYPE nhc_info_type
     REAL ( dbl ) :: eta
     REAL ( dbl ) :: v
     REAL ( dbl ) :: f
     REAL ( dbl ) :: nkt
     REAL ( dbl ) :: mass
  END TYPE nhc_info_type
  
  TYPE particle_kinetic_type
     REAL ( dbl ), POINTER :: point
  END TYPE particle_kinetic_type
  
  TYPE particle_scale_type
     REAL ( dbl ), POINTER :: point
  END TYPE particle_scale_type
  
  TYPE extended_parameters_type
     CHARACTER ( LEN = 10 ) :: dis_type
     INTEGER :: nyosh, nc, nhc_len, num_nhc
     REAL ( dbl ), POINTER :: dt_yosh ( : )
     TYPE (nhc_info_type ), POINTER :: nvt ( :, : )
     TYPE (npt_info_type ), POINTER :: npt ( :, : )
     TYPE (particle_kinetic_type ), POINTER :: p_kin ( :, : )
     TYPE (particle_scale_type ), POINTER :: p_scale ( :, : )
     REAL ( dbl ), POINTER :: v_scale ( : )
     REAL ( dbl ), POINTER :: s_kin ( : )
  END TYPE extended_parameters_type
  
CONTAINS

!******************************************************************************

SUBROUTINE lnhc ( nhcp, pnode, group )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
  INTEGER, INTENT ( IN ) :: group
  
! Locals
  INTEGER :: i, nnodes, inc, iyosh, n, nhc, handle, nx1, nx2
  REAL ( dbl ) :: scale, flops
  
!------------------------------------------------------------------------------
  
  CALL timeset('LNHC','I','Mflops',handle)
  flops = 0.0_dbl
  
  nx1 = SIZE ( nhcp % nvt, 1 )
  nx2 = SIZE ( nhcp % nvt, 2 )
  
! get force on first thermostat for all the chains in the system.
  nhcp % v_scale = 1.0_dbl
  nhcp % s_kin = 0.0_dbl
  nnodes = SIZE ( pnode )
  DO i = 1, nnodes
     nhcp % p_kin(1,i) % point = nhcp % p_kin(1,i) % point &
          + pnode(i) % p % prop % mass*pnode(i) % p % v(1)*pnode(i) % p % v(1)
     nhcp % p_kin(2,i) % point = nhcp % p_kin(2,i) % point &
          + pnode(i) % p % prop % mass*pnode(i) % p % v(2)*pnode(i) % p % v(2)
     nhcp % p_kin(3,i) % point = nhcp % p_kin(3,i) % point &
          + pnode(i) % p % prop % mass*pnode(i) % p % v(3)*pnode(i) % p % v(3)
  END DO
#if defined(__parallel)
! for replicated thermostats and distributed molecules, we have to
! sum s_kin (pointed to by p_kin%point) over all processors
  IF ( nhcp % dis_type == 'INTER_REP' ) THEN
     CALL mp_sum ( nhcp % s_kin, group )
  END IF
#endif
  flops = flops + REAL ( nnodes * 9, dbl )

! force on the first bead in every thermostat chain
  DO n = 1, nhcp % num_nhc
     nhcp % nvt(1,n) % f = (nhcp % s_kin(n)-nhcp % nvt(1,n) % nkt)/ &
          nhcp % nvt(1,n) % mass
  END DO
  flops = flops + REAL ( nhcp % num_nhc * 2, dbl )

! perform multiple time stepping using Yoshida
  NCLOOP: DO inc = 1, nhcp % nc
     YOSH: DO iyosh = 1, nhcp % nyosh
        
! update velocity on the last thermostat in the chain    ! O1
        nhcp % nvt(nhcp % nhc_len,:) % v = nhcp % nvt(nhcp % nhc_len,:) % v + &
             nhcp % nvt(nhcp % nhc_len,:) % f*0.25_dbl*nhcp % dt_yosh(iyosh)
        flops = flops + REAL ( 3 * nx2, dbl )
        
! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
        DO n = 1, nhcp % num_nhc
           DO nhc = nhcp % nhc_len - 1, 1, -1
              scale = exp(-0.125_dbl*nhcp % nvt(nhc+1,n) % v*nhcp % dt_yosh(iyosh) &
                   )
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v + &
                   nhcp % nvt(nhc,n) % f*0.25_dbl*nhcp % dt_yosh(iyosh) ! shift
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
           END DO
        END DO
        flops = flops &
             + REAL ( nhcp % num_nhc * (nhcp % nhc_len - 1 ) * 8, dbl )

! the core of the operator ----- START------
! update nhc positions
        nhcp % nvt ( :, : ) % eta = nhcp % nvt ( :, : ) % eta + &
             0.5_dbl*nhcp % nvt ( :, : ) % v*nhcp % dt_yosh(iyosh)
        flops = flops + REAL(2*nx1*nx2)
! now accumulate the scale factor for particle velocities
        DO n = 1, nhcp % num_nhc
           nhcp % v_scale(n) = nhcp % v_scale(n)*exp(-0.5_dbl*nhcp % dt_yosh( &
                iyosh)*nhcp % nvt(1,n) % v)
        END DO
        flops = flops + REAL(nhcp % num_nhc*4)
! the core of the operator ------ END ------

! update the force on first thermostat again (since particle velocities 
! have changed)
        DO n = 1, nhcp % num_nhc
           nhcp % nvt(1,n) % f = (nhcp % s_kin(n)*nhcp % v_scale(n)*nhcp % v_scale(n &
                )-nhcp % nvt(1,n) % nkt)/nhcp % nvt(1,n) % mass
        END DO
        flops = flops + REAL(nhcp % num_nhc*4)

! update velocity of other thermostats on chain (from 1 to nhc_len-1)  ! O2
        DO nhc = 1, nhcp % nhc_len - 1
           DO n = 1, nhcp % num_nhc
              scale = exp(-0.125_dbl*nhcp % nvt(nhc+1,n) % v*nhcp % dt_yosh(iyosh) &
                   )
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v + &
                   nhcp % nvt(nhc,n) % f*0.25_dbl*nhcp % dt_yosh(iyosh) ! shift
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
           END DO
! updating the forces on all the thermostats
           DO n = 1, nhcp % num_nhc
              nhcp % nvt(nhc+1,n) % f = (nhcp % nvt(nhc,n) % mass*nhcp % nvt(nhc,n) % v &
                   *nhcp % nvt(nhc,n) % v-nhcp % nvt(nhc+1,n) % nkt)/ &
                   nhcp % nvt(nhc+1,n) % mass
           END DO
        END DO
        flops = flops + REAL(nhcp % nhc_len-1)*REAL(nhcp % num_nhc*12)
! update velocity on last thermostat                                    ! O1
        nhcp % nvt(nhcp % nhc_len,:) % v = nhcp % nvt(nhcp % nhc_len,:) % v + &
             nhcp % nvt(nhcp % nhc_len,:) % f*0.25_dbl*nhcp % dt_yosh(iyosh)
        flops = flops + REAL(3*nx2)
     END DO YOSH
  END DO NCLOOP

! now scale the particle velocities
  DO i = 1, nnodes
     pnode(i) % p % v(1) = pnode(i) % p % v(1)*nhcp % p_scale(1,i) % point
     pnode(i) % p % v(2) = pnode(i) % p % v(2)*nhcp % p_scale(2,i) % point
     pnode(i) % p % v(3) = pnode(i) % p % v(3)*nhcp % p_scale(3,i) % point
  END DO
  flops = flops + REAL ( nnodes*3, dbl )
  
  flops = flops * 1.E-6_dbl
  CALL timestop(flops,handle)
  
END SUBROUTINE lnhc

!******************************************************************************

SUBROUTINE lnhcp ( nhcp, molecule, pnode, pv_virial, pv_kin, pv_const, pext, &
     deth, nfree, flag, group )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  TYPE ( molecule_structure_type ), INTENT ( INOUT ), DIMENSION ( : ) :: &
       molecule
  TYPE ( particle_node_type ), INTENT ( INOUT ), DIMENSION ( : ) :: pnode
  REAL ( dbl ), DIMENSION (3,3), INTENT ( IN ) :: pv_virial, pv_kin, pv_const
  REAL ( dbl ), INTENT ( IN ) :: pext, deth
  INTEGER, INTENT ( IN ) :: nfree
  LOGICAL, INTENT ( IN ) :: flag
  INTEGER, INTENT ( IN ) :: group
  
! Locals
  INTEGER :: i, ii, nnodes, inc, iyosh, n, nhc, handle, nx1, nx2, isos
  REAL ( dbl ) :: scale, flops, baro_kin, aa, infree, fdotr, kin
! roll stuff
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( :, : ), SAVE :: veps, veta, mass
  
!------------------------------------------------------------------------------
  
  CALL timeset('LNHC','I','Mflops', handle )
  flops = 0.0_dbl
  
  infree = 1.0_dbl / REAL ( nfree, dbl )
  
  nx1 = size(nhcp % nvt(:,1))
  nx2 = size(nhcp % nvt(1,:))
  
! get force on barostat
  kin = 0.0_dbl
  fdotr = 0.0_dbl
  DO ii = 1, 3
     kin = kin + pv_kin(ii,ii)
     fdotr = fdotr + pv_virial(ii,ii) + pv_const(ii,ii)
  END DO
  nhcp % npt ( :, : ) % f = (1._dbl+(3._dbl*infree))*kin + fdotr - &
       3._dbl*pext*deth

! get force on first thermostat for all the chains in the system.
  nhcp % v_scale = 1.0_dbl
  nhcp % s_kin = 0.0_dbl
  nnodes = size(pnode)
  DO i = 1, nnodes
     nhcp % p_kin(1,i) % point = nhcp % p_kin(1,i) % point + &
          pnode(i) % p % prop % mass*pnode(i) % p % v(1)*pnode(i) % p % v(1)
     nhcp % p_kin(2,i) % point = nhcp % p_kin(2,i) % point + &
          pnode(i) % p % prop % mass*pnode(i) % p % v(2)*pnode(i) % p % v(2)
     nhcp % p_kin(3,i) % point = nhcp % p_kin(3,i) % point + &
          pnode(i) % p % prop % mass*pnode(i) % p % v(3)*pnode(i) % p % v(3)
  END DO
  
#if defined(__parallel)
! for replicated thermostats and distributed molecules, we have to
! sum s_kin (pointed to by p_kin%point) over all processors
  IF (nhcp % dis_type=='INTER_REP') THEN
     CALL mp_sum(nhcp % s_kin,group)
  END IF
#endif
  
  flops = flops + REAL ( nnodes * 9, dbl )
  
! force on the first bead in every thermostat chain attached to system
  DO n = 1, nhcp % num_nhc - 1
     nhcp % nvt(1,n) % f = (nhcp % s_kin(n)-nhcp % nvt(1,n) % nkt)/ &
          nhcp % nvt(1,n) % mass
  END DO
  
! force on the first bead attached to barostat
  baro_kin = nhcp % npt(1,1) % mass*nhcp % npt(1,1) % v*nhcp % npt(1,1) % v
  nhcp % nvt(1,nhcp % num_nhc) % f = (baro_kin-nhcp % nvt(1,nhcp % num_nhc) % nkt)/ &
       nhcp % nvt(1,nhcp % num_nhc) % mass

  flops = flops + REAL(nhcp % num_nhc*2)

! perform multiple time stepping using Yoshida
  NCLOOP: DO inc = 1, nhcp % nc
     YOSH:     DO iyosh = 1, nhcp % nyosh

! update velocity on the last thermostat in the chain    ! O1
        nhcp % nvt(nhcp % nhc_len,:) % v = nhcp % nvt(nhcp % nhc_len,:) % v + &
             nhcp % nvt(nhcp % nhc_len,:) % f*0.25_dbl*nhcp % dt_yosh(iyosh)
        flops = flops + REAL(3*nx2)

! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
        DO n = 1, nhcp % num_nhc
           DO nhc = nhcp % nhc_len - 1, 1, -1
              scale = exp(-0.125_dbl*nhcp % nvt(nhc+1,n) % v*nhcp % dt_yosh(iyosh) &
                   )
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v + &
                   nhcp % nvt(nhc,n) % f*0.25_dbl*nhcp % dt_yosh(iyosh) ! shift
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
           END DO
        END DO
! updating veps:
        aa = exp(-.125_dbl*nhcp % dt_yosh(iyosh)*nhcp % nvt(1,nhcp % num_nhc) % v)
        nhcp % npt ( :, : ) % v = nhcp % npt ( :, : ) % v*aa*aa + &
             .25_dbl*nhcp % dt_yosh(iyosh)*nhcp % npt ( :, : ) % f*aa/ &
             nhcp % npt ( :, : ) % mass
!
        flops = flops + REAL(nhcp % num_nhc*(nhcp % nhc_len-1)*8)

! the core of the operator ----- START------
! update nhc positions
        nhcp % nvt ( :, : ) % eta = nhcp % nvt ( :, : ) % eta + &
             0.5_dbl*nhcp % nvt ( :, : ) % v*nhcp % dt_yosh(iyosh)
        flops = flops + REAL(2*nx1*nx2)
! now accumulate the scale factor for particle velocities
        DO n = 1, nhcp % num_nhc - 1
           nhcp % v_scale(n) = nhcp % v_scale(n)*exp(-0.5_dbl*nhcp % dt_yosh( &
                iyosh)*(nhcp % nvt(1,n) % v+(1._dbl+3._dbl*infree)*nhcp % npt(1, &
                1) % v))
        END DO
        flops = flops + REAL(nhcp % num_nhc*4)
! the core of the operator ------ END ------

! when constraints, call roll
!
        IF (flag) THEN
           isos = 0
           IF ( .NOT. ALLOCATED ( veta)) ALLOCATE (veta(size(nhcp % nvt, &
                1),size(nhcp % nvt,2)),STAT=isos)
           IF ( .NOT. ALLOCATED ( veps)) ALLOCATE (veps(1,1),STAT=isos)
           IF ( .NOT. ALLOCATED ( mass)) ALLOCATE (mass(1,1),STAT=isos)
           IF ( isos /= 0 ) CALL stop_prg ( 'lnhcp','failed to allocate')
           veta ( :, : ) = nhcp % nvt ( :, : ) % v
           veps(1,1) = nhcp % npt(1,1) % v
           mass(1,1) = nhcp % npt(1,1) % mass
           CALL roll(molecule,veps,veta,mass,nhcp % dt_yosh(iyosh), &
                nhcp % num_nhc,nfree)
        END IF
!
! update the force on first thermostat again (since particle velocities 
! have changed) connected to the system
        kin = 0._dbl
        DO n = 1, nhcp % num_nhc - 1
           nhcp % nvt(1,n) % f = (nhcp % s_kin(n)*nhcp % v_scale(n)*nhcp % v_scale(n &
                )-nhcp % nvt(1,n) % nkt)/nhcp % nvt(1,n) % mass
           kin = kin + nhcp % s_kin(n)*nhcp % v_scale(n)*nhcp % v_scale(n)
        END DO
        flops = flops + REAL(nhcp % num_nhc*4)

        nhcp % npt ( :, : ) % f = (1._dbl+(3._dbl*infree))*kin + fdotr - &
             3._dbl*pext*deth
! updating veps:
        aa = exp(-.125_dbl*nhcp % dt_yosh(iyosh)*nhcp % nvt(1,nhcp % num_nhc) % v)
        nhcp % npt ( :, : ) % v = nhcp % npt ( :, : ) % v*aa*aa + &
             .25_dbl*nhcp % dt_yosh(iyosh)*nhcp % npt ( :, : ) % f*aa/ &
             nhcp % npt ( :, : ) % mass
! update force on first thermostat connected to barostat
        baro_kin = nhcp % npt(1,1) % mass*nhcp % npt(1,1) % v*nhcp % npt(1,1) % v
        nhcp % nvt(1,nhcp % num_nhc) % f = (baro_kin-nhcp % nvt(1,nhcp % num_nhc) % &
             nkt)/nhcp % nvt(1,nhcp % num_nhc) % mass

! update velocity of other thermostats on chain (from 1 to nhc_len-1)  ! O2
        DO nhc = 1, nhcp % nhc_len - 1
           DO n = 1, nhcp % num_nhc
              scale = exp(-0.125_dbl*nhcp % nvt(nhc+1,n) % v*nhcp % dt_yosh(iyosh) &
                   )
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v + &
                   nhcp % nvt(nhc,n) % f*0.25_dbl*nhcp % dt_yosh(iyosh) ! shift
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
           END DO
! updating the forces on all the thermostats
           DO n = 1, nhcp % num_nhc
              nhcp % nvt(nhc+1,n) % f = (nhcp % nvt(nhc,n) % mass*nhcp % nvt(nhc,n) % v &
                   *nhcp % nvt(nhc,n) % v-nhcp % nvt(nhc+1,n) % nkt)/ &
                   nhcp % nvt(nhc+1,n) % mass
           END DO
        END DO
        flops = flops + REAL(nhcp % nhc_len-1)*REAL(nhcp % num_nhc*12)
! update velocity on last thermostat                                    ! O1
        nhcp % nvt(nhcp % nhc_len,:) % v = nhcp % nvt(nhcp % nhc_len,:) % v + &
             nhcp % nvt(nhcp % nhc_len,:) % f*0.25_dbl*nhcp % dt_yosh(iyosh)
        flops = flops + REAL(3*nx2)
     END DO YOSH
  END DO NCLOOP

! now scale the particle velocities
  DO i = 1, nnodes
     pnode(i) % p % v(1) = pnode(i) % p % v(1)*nhcp % p_scale(1,i) % point
     pnode(i) % p % v(2) = pnode(i) % p % v(2)*nhcp % p_scale(2,i) % point
     pnode(i) % p % v(3) = pnode(i) % p % v(3)*nhcp % p_scale(3,i) % point
  END DO
  flops = flops + REAL ( nnodes * 3, dbl )
  
  flops = flops * 1.0E-6_dbl
  CALL timestop ( flops, handle )
  
END SUBROUTINE lnhcp

!******************************************************************************

SUBROUTINE lnhcpf ( nhcp, molecule, pnode, pv_virial, pv_kin, pv_const, pext, &
     deth, nfree, flag, group )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( IN ) :: pv_virial, pv_kin
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( IN ) :: pv_const
  REAL ( dbl ), INTENT ( IN ) :: pext, deth
  INTEGER, INTENT ( IN ) :: nfree
  LOGICAL, INTENT ( IN ) :: flag
  INTEGER, INTENT ( IN ) :: group
  
! Locals
  INTEGER :: i, ipart, ii, jj, nnodes, inc, iyosh, n, nhc, handle, nx1, nx2
  INTEGER :: isos
  REAL ( dbl ) :: scale, flops, baro_kin, aa, infree, trvg, kin
  REAL ( dbl ), DIMENSION (3,3) :: u, unit, pv_kin_loc
  REAL ( dbl ), DIMENSION (3) :: uv, e, e_val
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( :, : ), SAVE :: veps, veta, mass
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'LNHC', 'I', 'Mflops', handle )
  flops = 0.0_dbl
  
! initializing locals
  unit = 0.0_dbl
  unit(1,1) = 1.0_dbl
  unit(2,2) = 1.0_dbl
  unit(3,3) = 1.0_dbl
  infree = 1.0_dbl / REAL ( nfree, dbl )
  
  nx1 = size(nhcp % nvt(:,1))
  nx2 = size(nhcp % nvt(1,:))
  
! get force on barostat
  kin =  0._dbl
  DO ii = 1, 3
     kin = kin + pv_kin(ii,ii)
  END DO
  nhcp % npt ( :, : ) % f = pv_kin ( :, : ) + infree*kin*unit ( :, : ) + &
       pv_virial ( :, : ) + pv_const ( :, : ) - pext*deth*unit ( :, : )

! get force on first thermostat for all the chains in the system.
  nhcp % v_scale = 1.0_dbl
  nhcp % s_kin = 0.0_dbl
  nnodes = size(pnode)
  DO i = 1, nnodes
     nhcp % p_kin(1,i) % point = nhcp % p_kin(1,i) % point + &
          pnode(i) % p % prop % mass*pnode(i) % p % v(1)*pnode(i) % p % v(1)
     nhcp % p_kin(2,i) % point = nhcp % p_kin(2,i) % point + &
          pnode(i) % p % prop % mass*pnode(i) % p % v(2)*pnode(i) % p % v(2)
     nhcp % p_kin(3,i) % point = nhcp % p_kin(3,i) % point + &
          pnode(i) % p % prop % mass*pnode(i) % p % v(3)*pnode(i) % p % v(3)
  END DO
#if defined(__parallel)
! for replicated thermostats and distributed molecules, we have to
! sum s_kin (pointed to by p_kin % point) over all processors
  IF (nhcp % dis_type=='INTER_REP') THEN
     CALL mp_sum(nhcp % s_kin,group)
  END IF
#endif
  flops = flops + REAL(nnodes*9)

! force on the first bead in every thermostat chain attached to system
  DO n = 1, nhcp % num_nhc - 1
     nhcp % nvt(1,n) % f = (nhcp % s_kin(n)-nhcp % nvt(1,n) % nkt)/ &
          nhcp % nvt(1,n) % mass
  END DO
! force on the first bead attached to barostat
  baro_kin = 0._dbl
  DO ii = 1, 3
     DO jj = 1, 3
        baro_kin = baro_kin*nhcp % npt(ii,jj) % v*nhcp % npt(ii,jj) % v* &
             nhcp % npt(ii,jj) % mass
     END DO
  END DO
  nhcp % nvt(1,nhcp % num_nhc) % f = (baro_kin-nhcp % nvt(1,nhcp % num_nhc) % nkt)/ &
       nhcp % nvt(1,nhcp % num_nhc) % mass

  flops = flops + REAL(nhcp % num_nhc*2)

! perform multiple time stepping using Yoshida
  NCLOOP: DO inc = 1, nhcp % nc
     YOSH:     DO iyosh = 1, nhcp % nyosh

! update velocity on the last thermostat in the chain    ! O1
        nhcp % nvt(nhcp % nhc_len,:) % v = nhcp % nvt(nhcp % nhc_len,:) % v + &
             nhcp % nvt(nhcp % nhc_len,:) % f*0.25_dbl*nhcp % dt_yosh(iyosh)
        flops = flops + REAL(3*nx2)

! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
        DO n = 1, nhcp % num_nhc
           DO nhc = nhcp % nhc_len - 1, 1, -1
              scale = exp(-0.125_dbl*nhcp % nvt(nhc+1,n) % v*nhcp % dt_yosh(iyosh) &
                   )
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v + &
                   nhcp % nvt(nhc,n) % f*0.25_dbl*nhcp % dt_yosh(iyosh) ! shift
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
           END DO
        END DO
! updating veps:
        aa = exp(-.125_dbl*nhcp % dt_yosh(iyosh)*nhcp % nvt(1,nhcp % num_nhc) % v)
        nhcp % npt ( :, : ) % v = nhcp % npt ( :, : ) % v*aa*aa + &
             .25_dbl*nhcp % dt_yosh(iyosh)*nhcp % npt ( :, : ) % f*aa/ &
             nhcp % npt ( :, : ) % mass
! computing the trace
        trvg = infree*(nhcp % npt(1,1) % v+nhcp % npt(2,2) % v+nhcp % npt(3,3) % v)
!
        flops = flops + REAL(nhcp % num_nhc*(nhcp % nhc_len-1)*8)
! the core of the operator ----- START------
! update nhc positions
        nhcp % nvt ( :, : ) % eta = nhcp % nvt ( :, : ) % eta + &
             0.5_dbl*nhcp % nvt ( :, : ) % v*nhcp % dt_yosh(iyosh)
        flops = flops + REAL(2*nx1*nx2)
!
! accumulate the matrix part of the scale factor on the velocities
! Get eigenvectors and eigenvalues of nhcp % npt % v
!
!*apsi            CALL rs(nhcp % npt % v,e_val,u)
        CALL diagonalise ( matrix = nhcp % npt % v, mysize = 3, &
             storageform = "UPPER", eigenvalues = e_val, eigenvectors = u )

        e = exp(-.5_dbl*nhcp % dt_yosh(iyosh)*e_val)
! now get the scalar part of the scale factor 
        DO n = 1, nhcp % num_nhc - 1
           nhcp % v_scale(n) = exp(-0.5_dbl*nhcp % dt_yosh(iyosh)*(nhcp % nvt(1, &
                n) % v+trvg))
        END DO
        nhcp % s_kin = 0._dbl
        pv_kin_loc = 0._dbl
        kin = 0._dbl
        DO ipart = 1, nnodes
! update the particles with the matrix contribution
           uv = matmul(transpose(u),pnode(ipart) % p % v)
           uv = uv*e
           pnode(ipart) % p % v = matmul(u,uv)
! update velocites to compute kinetic energy dependent pieces
           pnode(ipart) % p % v(1) = pnode(ipart) % p % v(1)* &
                nhcp % p_scale(1,ipart) % point
           pnode(ipart) % p % v(2) = pnode(ipart) % p % v(2)* &
                nhcp % p_scale(2,ipart) % point
           pnode(ipart) % p % v(3) = pnode(ipart) % p % v(3)* &
                nhcp % p_scale(3,ipart) % point
           nhcp % p_kin(1,ipart) % point = nhcp % p_kin(1,ipart) % point + &
                pnode(ipart) % p % prop % mass*pnode(ipart) % p % v(1)* &
                pnode(ipart) % p % v(1)
           nhcp % p_kin(2,ipart) % point = nhcp % p_kin(2,ipart) % point + &
                pnode(ipart) % p % prop % mass*pnode(ipart) % p % v(2)* &
                pnode(ipart) % p % v(2)
           nhcp % p_kin(3,ipart) % point = nhcp % p_kin(3,ipart) % point + &
                pnode(ipart) % p % prop % mass*pnode(ipart) % p % v(3)* &
                pnode(ipart) % p % v(3)
           DO ii = 1, 3
              kin = kin + pnode(ipart) % p % prop % mass*pnode(ipart) % p % v(ii)* &
                   pnode(ipart) % p % v(ii)
              DO jj = 1, 3
                 pv_kin_loc(ii,jj) = pv_kin_loc(ii,jj) + &
                      pnode(ipart) % p % prop % mass*pnode(ipart) % p % v(ii)* &
                      pnode(ipart) % p % v(jj)
              END DO
           END DO
        END DO
#if defined(__parallel)
! for replicated thermostats and distributed molecules, we have to
! sum s_kin (pointed to by p_kin % point) over all processors
        IF (nhcp % dis_type=='INTER_REP') THEN
           CALL mp_sum(nhcp % s_kin,group)
        END IF
#endif
        flops = flops + REAL(nhcp % num_nhc*4)
! the core of the operator ------ END ------

! when constraints, call roll
!
        IF (flag) THEN
           IF ( .NOT. ALLOCATED ( veta)) ALLOCATE (veta(size(nhcp % nvt, &
                1),size(nhcp % nvt,2)),STAT=isos)
           IF ( .NOT. ALLOCATED ( veps)) ALLOCATE (veps(3,3),STAT=isos)
           IF ( .NOT. ALLOCATED ( mass)) ALLOCATE (mass(3,3),STAT=isos)
           IF ( isos /= 0 ) CALL stop_prg ( 'lnhcpf','failed to allocate')
           veta ( :, : ) = nhcp % nvt ( :, : ) % v
           mass ( :, : ) = nhcp % npt ( :, : ) % mass
           CALL roll(molecule,veps,veta,mass,nhcp % dt_yosh(iyosh), &
                nhcp % num_nhc,nfree,e_val,u)
        END IF
!
! update the force on first thermostat again (since particle velocities 
! have changed) connected to the system
        DO n = 1, nhcp % num_nhc - 1
           nhcp % nvt(1,n) % f = (nhcp % s_kin(n)-nhcp % nvt(1,n) % nkt)/ &
                nhcp % nvt(1,n) % mass
        END DO
        flops = flops + REAL(nhcp % num_nhc*4)

        nhcp % npt ( :, : ) % f = pv_kin_loc ( :, : ) + infree*kin*unit ( :, : ) + &
             pv_virial ( :, : ) + pv_const ( :, : ) - pext*deth*unit ( :, : )
! updating veps:
        aa = exp(-.125_dbl*nhcp % dt_yosh(iyosh)*nhcp % nvt(1,nhcp % num_nhc) % v)
        nhcp % npt ( :, : ) % v = nhcp % npt ( :, : ) % v*aa*aa + &
             .25_dbl*nhcp % dt_yosh(iyosh)*nhcp % npt ( :, : ) % f*aa/ &
             nhcp % npt ( :, : ) % mass
! update force on first thermostat connected to barostat
        baro_kin = 0._dbl
        DO ii = 1, 3
           DO jj = 1, 3
              baro_kin = baro_kin*nhcp % npt(ii,jj) % v*nhcp % npt(ii,jj) % v* &
                   nhcp % npt(ii,jj) % mass
           END DO
        END DO
        nhcp % nvt(1,nhcp % num_nhc) % f = (baro_kin-nhcp % nvt(1,nhcp % num_nhc) % &
             nkt)/nhcp % nvt(1,nhcp % num_nhc) % mass

! update velocity of other thermostats on chain (from 1 to nhc_len-1)  ! O2
        DO nhc = 1, nhcp % nhc_len - 1
           DO n = 1, nhcp % num_nhc
              scale = exp(-0.125_dbl*nhcp % nvt(nhc+1,n) % v*nhcp % dt_yosh(iyosh) &
                   )
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v + &
                   nhcp % nvt(nhc,n) % f*0.25_dbl*nhcp % dt_yosh(iyosh) ! shift
              nhcp % nvt(nhc,n) % v = nhcp % nvt(nhc,n) % v*scale ! scale
           END DO
! updating the forces on all the thermostats
           DO n = 1, nhcp % num_nhc
              nhcp % nvt(nhc+1,n) % f = (nhcp % nvt(nhc,n) % mass*nhcp % nvt(nhc,n) % v &
                   *nhcp % nvt(nhc,n) % v-nhcp % nvt(nhc+1,n) % nkt)/ &
                   nhcp % nvt(nhc+1,n) % mass
           END DO
        END DO
        flops = flops + REAL(nhcp % nhc_len-1)*REAL(nhcp % num_nhc*12)
! update velocity on last thermostat                                    ! O1
        nhcp % nvt(nhcp % nhc_len,:) % v = nhcp % nvt(nhcp % nhc_len,:) % v + &
             nhcp % nvt(nhcp % nhc_len,:) % f*0.25_dbl*nhcp % dt_yosh(iyosh)
        flops = flops + REAL(3*nx2)
     END DO YOSH
  END DO NCLOOP

  flops = flops + REAL(nnodes*3)
!
  flops = flops*1.E-6_dbl
  CALL timestop(flops,handle)
END SUBROUTINE lnhcpf
!******************************************************************************

SUBROUTINE initial_nhc(nhcp,tau_nhc,temp_ext, globenv )
! Initializes the NHC velocities to the Maxwellian distribution
  IMPLICIT NONE
  TYPE (extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  REAL ( dbl ), INTENT ( IN ) :: tau_nhc, temp_ext
  INTEGER :: i, j, number
  REAL ( dbl ) :: akin, temp, v

! first initializing the mass of the nhc variables
  nhcp % nvt ( :, : ) % mass = nhcp % nvt ( :, : ) % nkt*tau_nhc*tau_nhc
! initializing velocities
  DO i = 1, nhcp % num_nhc
     DO j = 1, nhcp % nhc_len
        v = gasdev(globenv % idum)
#if defined(__parallel)
        IF (nhcp % dis_type=='INTER_REP') CALL mp_bcast(v,globenv % source, &
             globenv % group)
#endif
        nhcp % nvt(j,i) % v = v
     END DO
  END DO
  akin = 0.0_dbl
  DO i = 1, nhcp % num_nhc
     DO j = 1, nhcp % nhc_len
        akin = akin + 0.5_dbl*(nhcp % nvt(j,i) % mass*nhcp % nvt(j,i) % v*nhcp % nvt &
             (j,i) % v)
     END DO
  END DO
  number = nhcp % num_nhc
#if defined(__parallel)
  IF (nhcp % dis_type=='DIS_REP') CALL mp_sum(akin,globenv % group)
  IF (nhcp % dis_type=='DIS_REP') CALL mp_sum(number,globenv % group)
#endif

  temp = 2.0_dbl*akin/REAL(number)

! scale velocities to get the correct initial temperature
  DO i = 1, nhcp % num_nhc
     DO j = 1, nhcp % nhc_len
        nhcp % nvt(j,i) % v = sqrt(temp_ext/temp)*nhcp % nvt(j,i) % v
! dbg
!            nhcp % nvt(j,i) % v = 0._dbl
! dbg
        nhcp % nvt(j,i) % eta = 0.0_dbl
     END DO
  END DO

! initializing all of the forces on the thermostats
  DO i = 1, nhcp % num_nhc
     DO j = 2, nhcp % nhc_len
        nhcp % nvt(j,i) % f = nhcp % nvt(j-1,i) % mass*nhcp % nvt(j-1,i) % v* &
             nhcp % nvt(j-1,i) % v - nhcp % nvt(j,i) % nkt
        nhcp % nvt(j,i) % f = nhcp % nvt(j,i) % f/nhcp % nvt(j,i) % mass
     END DO
  END DO
  
END SUBROUTINE initial_nhc

!******************************************************************************

! Initializes the NHC velocities to the Maxwellian distribution

SUBROUTINE initial_npt ( nhcp, tau_cell, temp_ext, deth, nfree, ensemble, &
     globenv )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  REAL ( dbl ), INTENT ( IN ) :: tau_cell, temp_ext, deth
  CHARACTER ( LEN = 20 ), INTENT ( IN ) :: ensemble
  INTEGER, INTENT ( IN ) :: nfree
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  
! Locals
  INTEGER :: i, j, number, isos
  REAL ( dbl ) :: akin, temp, v
  
!------------------------------------------------------------------------------
  
! first initializing the mass of the nhc variables
  SELECT CASE ( ensemble )
  CASE ( 'NPT_I')
     ALLOCATE ( nhcp % npt ( 1, 1 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_prg ( 'initial_npt', &
          'failed to allocate nhcp % npt')
     nhcp % npt ( :, : ) % mass = REAL ( nfree + 3, dbl ) &
          * temp_ext * tau_cell ** 2
     
  CASE ( 'NPT_F')
     ALLOCATE ( nhcp % npt ( 3, 3 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_prg ( 'initial_npt', &
          'failed to allocate nhcp % npt')
     nhcp % npt ( :, : ) % mass = REAL ( nfree + 3, dbl ) &
          * temp_ext * tau_cell ** 2 / 3.0_dbl
  END SELECT
  
! initializing velocities
  DO i = 1, size ( nhcp % npt,1)
     DO j = 1, size ( nhcp % npt,2)
        v = gasdev ( globenv % idum)
#if defined(__parallel)
        CALL mp_bcast ( v, globenv % source, globenv % group )
#endif
        nhcp % npt(j,i) % v = v
     END DO
  END DO
  
  akin = 0.0_dbl
  DO i = 1, size(nhcp % npt,1)
     DO j = 1, size(nhcp % npt,2)
        akin = akin + 0.5_dbl*(nhcp % npt(j,i) % mass*nhcp % npt(j,i) % v*nhcp % npt &
             (j,i) % v)
     END DO
  END DO
  number = SIZE ( nhcp % npt, 1 ) * SIZE ( nhcp % npt, 2 )
  
  temp = 2.0_dbl * akin / REAL ( number, dbl )
  
! scale velocities to get the correct initial temperature
  DO i = 1, size(nhcp % npt,1)
     DO j = 1, size(nhcp % npt,2)
        nhcp % npt(j,i) % v = sqrt(temp_ext/temp)*nhcp % npt(j,i) % v
! dbg
!            nhcp % npt(j,i) % v = 0._dbl
! dbg
     END DO
  END DO
  
! initialize the position
  IF ( ensemble == 'NPT_I' ) nhcp % npt ( :, : ) % eps = LOG ( deth ) / 3.0_dbl
  
END SUBROUTINE initial_npt

!******************************************************************************

! this lists the coefficients for the Yoshida method (higher
!    order integrator used in NVT)

SUBROUTINE yoshida_coef ( nhcp, dt )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  REAL ( dbl ), INTENT ( IN ) :: dt
  
! Locals
  REAL ( dbl ), DIMENSION (nhcp % nyosh) :: yosh_wt
  
!------------------------------------------------------------------------------
  
  SELECT CASE (nhcp % nyosh)
  CASE (1)
     yosh_wt(1) = 1.0_dbl
  CASE (3)
     yosh_wt(1) = 1.0_dbl/(2.0_dbl-(2.0_dbl)**(1.0_dbl/3.0_dbl))
     yosh_wt(2) = 1.0_dbl - 2.0_dbl*yosh_wt(1)
     yosh_wt(3) = yosh_wt(1)
  CASE (5)
     yosh_wt(1) = 1.0_dbl/(4.0_dbl-(4.0_dbl)**(1.0_dbl/3.0_dbl))
     yosh_wt(2) = yosh_wt(1)
     yosh_wt(4) = yosh_wt(1)
     yosh_wt(5) = yosh_wt(1)
     yosh_wt(3) = 1.0_dbl - 4.0_dbl*yosh_wt(1)
  CASE (7)
     yosh_wt(1) = .78451361047756_dbl
     yosh_wt(2) = .235573213359357_dbl
     yosh_wt(3) = -1.17767998417887_dbl
     yosh_wt(4) = 1.0_dbl - 2.0_dbl*(yosh_wt(1)+yosh_wt(2)+yosh_wt(3))
     yosh_wt(5) = yosh_wt(3)
     yosh_wt(6) = yosh_wt(2)
     yosh_wt(7) = yosh_wt(1)
  END SELECT
  
  nhcp % dt_yosh = dt * yosh_wt / REAL ( nhcp % nc, dbl )
  
END SUBROUTINE yoshida_coef

!******************************************************************************

END MODULE nose

