!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/scp_dispersion
!!
!!   NAME
!!     scp_dispersion
!!
!!   FUNCTION
!!     Perform the calculation of the dispersion interaction
!!
!!   AUTHOR
!!      CJM
!!
!!   SOURCE
!******************************************************************************

MODULE  scp_dispersion
!******************************************************************************
!!****f* scp_dispersion/scp_calculate_dispersion [1.0] *
!!
!!   NAME
!!     scp_calculate_dispersion
!!
!!   FUNCTION
!!     Analytic calculation of the dispersion energy and forces 
!!     for SCP.   
!!
!!   NOTE
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!     CJM
!!
!!
!!*** **********************************************************************

  SUBROUTINE scp_calculate_dispersion ( qs_env, calculate_forces, error, debug )

    TYPE(qs_environment_type), POINTER                 :: qs_env
    LOGICAL, INTENT ( IN ), OPTIONAL                   :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)                 :: error
    LOGICAL, INTENT ( IN ), OPTIONAL                   :: debug

    CHARACTER(LEN=*), PARAMETER :: routineN = "scp_calculate_dispersion"

    TYPE(atomic_kind_type), DIMENSION(:), POINTER       :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER                 :: local_particles
    TYPE(atomic_kind_type),       POINTER      :: atomic_kind
    TYPE(gto_basis_set_type),     POINTER      :: aux_basis
    TYPE(cp_para_env_type), POINTER            :: para_env
    TYPE ( aux_coeff_set_type ), POINTER       :: aux_coeff_set
    TYPE ( aux_coeff_type ), POINTER           :: local_coeffs
    TYPE ( scp_energy_type ), POINTER          :: energy
    REAL(kind=dp)                              :: e_scp

    REAL(dp), DIMENSION(:,:),       ALLOCATABLE  :: work, gcca, gccb,fcoeff
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE  :: vv
    REAL(dp), DIMENSION ( :, : ), POINTER      :: vab

    REAL(dp), DIMENSION ( :, : ), POINTER      :: sphi, zet, rpgf
    REAL(dp), DIMENSION ( : ), POINTER         :: zeta, zetb, rpgfa, rpgfb
    REAL(dp), DIMENSION(:),   ALLOCATABLE      :: ff

    INTEGER :: handle, iatom,  i, istat
    INTEGER :: icoa, icob, ncoa, ncob
    INTEGER :: iparticle_local, nparticle_local, ikind, nkind
    INTEGER :: iseta, isetb, sgfa, sgfb
    INTEGER :: ncotot, maxco, nset
    INTEGER :: offseta, offsetb, nza, nzb, ldvab1, ldvab2
    INTEGER :: maxco_global, maxlgto, maxsgf_set, output_unit
    INTEGER, DIMENSION ( : ), POINTER :: l_max, l_min, npgf, nsgf_set
    INTEGER, DIMENSION ( :, : ), POINTER :: first_sgf

    LOGICAL :: ionode, my_debug, forces
    


!   ---------------------------------------------------------------------------
    my_debug = .FALSE.
    IF ( PRESENT ( debug ) ) my_debug = debug
    forces = .FALSE.
    IF ( PRESENT ( calculate_forces ) ) forces = calculate_forces

    NULLIFY( atomic_kind, aux_basis, aux_coeff_set)
    NULLIFY( energy, local_coeffs, para_env )

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                    cell=cell, scp_env = scp_env, &
                    para_env=para_env,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,error=error)

    CALL get_scp_env(scp_env=scp_env, para_env=para_env, & 
                     energy=energy, aux_coeff_set=aux_coeff_set, &
                     error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco_global,maxlgto=maxlgto)

    nkind = SIZE(atomic_kind_set,1)
    natom = SIZE(particle_set,1)

    maxder = 0
    IF ( forces ) THEN

      CALL timeset(routineN//" (forces)","I"," ",handle)

      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
           "atom_of_kind",natom*int_size)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

      CALL get_scp_env(scp_env=scp_env, force=force,error=error)

      maxder = 1
      maxlgto1 = maxlgto + maxder
    ELSE

      CALL timeset(routineN,"I"," ",handle)
      maxlgto1 = maxlgto

    END IF



    CALL timeset(routineN,"I"," ",handle)

    ALLOCATE (work (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "work",maxco_global*1*dp_size)
    ALLOCATE (gcca (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "gcca",maxco_global*1*dp_size)
    ALLOCATE (gccb (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "gccb",maxco_global*1*dp_size)

    ALLOCATE (work (maxco_global,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "work",maxco_global*1*dp_size)
    ALLOCATE(vab(1,1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vab",1*1*dp_size)
    ldvab1 = 1
    ldvab2 = 1
    vab = 0.0_dp

    ALLOCATE(vv(ncoset(maxlgto1),ncoset(maxlgto),maxlgto1+maxlgto+1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vv",maxco_global*maxco_global*(maxlgto1+maxlgto+1)*dp_size)
    ALLOCATE(ff(0:maxlgto1+maxlgto),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "ff", (maxlgto1+maxlgto+1)*dp_size)

    e_disp = 0._dp
    DO ikind = 1,nkind

      NULLIFY(atomic_kind)
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,aux_basis_set=aux_basis )
                           

! Check to see if the atom is an SCP atom
      IF (.NOT.ASSOCIATED(aux_basis)) CYCLE
      NULLIFY(first_sgfa,la_max,la_min,npgfa,nsgf_seta)
      NULLIFY(rpgfa,sphia,zeta)
      CALL get_gto_basis_set(gto_basis_set=aux_basis,&
                            first_sgf=first_sgfa,&
                            lmax=la_max,&
                            lmin=la_min,&
                            maxco=maxcoa,&
                            npgf=npgfa,&
                            nset=nseta,&
                            nsgf_set=nsgf_seta,&
                            pgf_radius=rpgfa,&
                            sphi=sphia,&
                            zet=zeta)
      ncotota = maxcoa * nseta
      DO jkind = 1, nkind 
        NULLIFY(atomic_kind)
        atomic_kind => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,aux_basis_set=aux_basis )
                           

! Check to see if the atom is an SCP atom
        IF (.NOT.ASSOCIATED(aux_basis)) CYCLE
        NULLIFY(first_sgfb,lb_max,lb_min,npgfb,nsgf_setb)
        NULLIFY(rpgfb,sphib,zetb)
        CALL get_gto_basis_set(gto_basis_set=aux_basis,&
                              first_sgf=first_sgfb,&
                              lmax=lb_max,&
                              lmin=lb_min,&
                              maxco=maxcob,&
                              npgf=npgfb,&
                              nset=nsetb,&
                              nsgf_set=nsgf_setb,&
                              pgf_radius=rpgfb,&
                              sphi=sphib,&
                              zet=zetb)
        ncototb = maxcob * nsetb

        ab = ikind + nkind*(jkind-1)

        NULLIFY(sab_orb_neighbor_list)
        IF(ASSOCIATED(sab_orb(ab)%neighbor_list_set)) THEN
          CALL get_neighbor_list_set(neighbor_list_set=&
               sab_orb(ab)%neighbor_list_set,nlist=nlist)

          sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)
        ELSE
          nlist = 0
        END IF

        DO ilist = 1,nlist

          sab_orb_neighbor_list_local => sab_orb_neighbor_list
          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom,nnode=nnode)

          ra(:) = pbc(particle_set(iatom)%r,cell)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab)
            rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

            rb(:) = pbc(particle_set(jatom)%r,cell)

            IF(iatom == jatom) THEN
              CYCLE
            ELSE
              ieqj = 2.0_dp
            END IF

            offseta = 0
            DO iseta = 1, nseta
              ncoa = npgf ( iseta ) * ncoset ( l_max ( iseta ) )
              sgfa = first_sgf ( 1, iseta )
              nza = npgf ( iseta )
              rpgfa => rpgf ( 1:nza, iseta )
              zeta => zet ( 1 : nza, iseta )
! Assigning the contraction coefficients
              DO i=1, nsgf_set ( iseta )
                 work ( i, 1 ) = 1.0_dp
              END DO
              CALL dgemm ( "N", "N", ncoa, 1, nsgf_seta ( iseta ), &
                   1.0_dp, sphia ( 1, sgfa ), SIZE ( sphia, 1 ), &
                   work ( 1, 1 ), SIZE ( work, 1 ), &
                   0.0_dp, gcca ( 1, 1 ), SIZE ( gcca, 1 ) )
              offsetb = 0
              DO isetb = 1, nset
                 ncob = npgf ( isetb ) * ncoset ( l_max ( isetb ) )
                 sgfb = first_sgf ( 1, isetb )
                 nzb = npgf ( isetb )
                 rpgfb => rpgf ( 1 : nzb, isetb )
                 zetb => zet ( 1 : nzb, isetb )
! Assigning the contraction coefficients
                 DO i=1, nsgf_set ( isetb )
                    work ( i, 1 ) = 1.0_dp
                 END DO
                 CALL dgemm ( "N", "N", ncob, 1, nsgf_setb ( isetb ), &
                      1.0_dp, sphib ( 1, sgfb ), SIZE ( sphib, 1 ), &
                      work ( 1, 1 ), SIZE ( work, 1 ), &
                      0.0_dp, gccb ( 1, 1 ), SIZE ( gccb, 1 ) )
! Allcating the arrays for the integrals
                 IF ( ncoa > ldvab1 .OR.  ncob > ldvab2 ) THEN
                   CALL reallocate ( vab , 1, ncoa, 1, ncob )
                   ldvab1 = ncoa
                   ldvab2 = ncob
                 ELSE
                   vab = 0._dp 
                 ENDIF 

                 IF ( forces ) THEN

                  ncoap = npgf ( iseta ) * ncoset ( l_max ( iseta ) + 1 )
                  ALLOCATE(vab_plus(ncoap,ncob), I_tmp(ncoap), STAT=istat)
                  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                  "vab_plus,I_tmp",ncoap*ncob*dp_size+ncoap*dp_size)
            
                   CALL coulomb2( l_max ( iseta ) + 1, npgfa ( iseta ), zeta, rpgfa, l_min ( iseta ), &
                                  l_max ( isetb ), npgfb ( isetb ), zetb, rpgfb, l_min ( isetb ), &
                                  rab, rab2, vab, vv, ff ( 0: ), maxder, vac_plus = vab_plus )
                   I_tmp = 0.0_dp
                   DO icoa = 1, ncoa ! may need a stride
                     DO icob = 1, ncob ! may need a stride
             ! **** need to omit (s|s) ****
                       I_tmp ( icoa ) = I_tmp ( icoa ) + vab ( icoa, icob ) * gccb ( icob )
                       e_disp = e_disp + gcca ( icoa )**2 * vab ( icoa, icob )**2 * gccb ( icob )**2
                     END DO
                   END DO
                   DO idir = 1,3
                     DO icoa = 1, ncoa
                       la ( 1:3 ) = indco ( 1:3, icoa )
                       lap ( 1:3 ) = la ( 1:3 )
                       lam ( 1:3 ) = la ( 1:3 )
                       lap ( idir ) = la ( idir ) + 1
                       lam ( idir ) = lam ( idir ) - 1
                       icoap = coset ( lap ( 1 ), lap ( 2 ), lap ( 3 ) )
             ! **** need to omit (s|s) ****
                       IF ( lam ( idir ) < 0 ) CYCLE
                       icoam = coset ( lam ( 1 ), lam ( 2 ), lam ( 3 ) )
                       force_a ( idir ) = force_a ( idir ) -  &
                                          2._dp * I_tmp ( icoa ) * gcca ( icoa ) *  &
                                          ( zeta ( ) * gcca ( icoa ) * I_tmp ( icoap ) + &
                                           la ( idir ) * gcca ( icoa ) * I_tmp ( icoam ) )
                     END DO
                   END DO ! idir
!             ***  sum up in the force arrey  ***
                   i_atom = atom_of_kind(iatom)
                   j_atom = atom_of_kind(jatom)
                   DO idir = 1,3
                     force(ikind)%scp_disp(idir,i_atom) = &
                       force(ikind)%scp_disp(idir,i_atom) - force_a(idir)
                     force(jkind)%scp_disp(idir,j_atom) = &
                        force(jkind)%scp_disp(idir,j_atom) + force_a(idir)
                   END DO
                   DEALLOCATE (vab_plus, I_tmp, STAT=istat )
                   IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                   "vab_plus.I_tmp")
                 ELSE
                   CALL coulomb2( l_max ( iseta ), npgfa ( iseta ), zeta, rpgfa, l_min ( iseta ), &
                                  l_max ( isetb ), npgfb ( isetb ), zetb, rpgfb, l_min ( isetb ), &
                                  rab, rab2, vab, vv, ff ( 0: ) )

                   DO icoa = 1, ncoa ! may need a stride
                     DO icob = 1, ncob ! may need a stride
             ! **** need to omit (s|s) ****
                       I_tmp ( icoa ) = I_tmp ( icoa ) + vab ( icoa, icob ) 
                       e_disp = e_disp + vab ( icoa, icob ) * vab ( icoa, icob )
                     END DO
                   END DO

                ENDIF

              END DO ! isetb
            END DO ! iseta
            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          END DO ! inode
          sab_orb_neighbor_list => next(sab_orb_neighbor_list)
        END DO ! ilist 
      END DO ! jkind
    END DO  ! ikind
    e_disp = 0.5_dp * e_disp

    CALL mp_sum(e_disp,para_env%group)
    energy % e_scp_disp = e_disp
   
    DEALLOCATE(vab,vv,ff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
          "vab.vv.ff")

    DEALLOCATE(work,gcca,gccb,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
          "work.gcca.gccb")

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE scp_calculate_dispersion

!******************************************************************************

