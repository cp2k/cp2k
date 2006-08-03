!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****h* cp2k/restraint_methods *
!!
!!   NAME
!!     restraint_methods
!!
!!   FUNCTION
!!     methods used for restraints
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     ikuo
!!     
!!
!!   MODIFICATION HISTORY
!!     07.2006 [ikuo]
!!
!!   SOURCE
!****************************************************************************
MODULE restraint_methods

  USE beta_gamma_psi,                  ONLY: erf,&
                                             gamln,&
                                             psi
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_output_handling,              ONLY: cp_iter_string,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE restraint_types,                 ONLY: restraint_env_type
  USE kinds,                           ONLY: dp,default_string_length
  USE mathconstants,                   ONLY: fac,&
                                             maxfac,&
                                             oorootpi
  USE mol_force,                       ONLY: force_bends,&
                                             force_bonds,&
                                             force_torsions
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE restraint_types,                 ONLY: fp_type,&
                                             harmonic_type,&
                                             dist_type,&
                                             angle_type,&
                                             torsion_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: restraint_eval

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'restraint_methods'

CONTAINS

!!****f* restraint_methods/restraint_eval *
!!
!!   NAME
!!    restraint_eval
!!
!!   FUNCTION
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [ikuo]
!!
!!*** **********************************************************************
  SUBROUTINE restraint_eval(restraint_env, subsys, cell, error)
    TYPE(restraint_env_type), POINTER        :: restraint_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'restraint_eval', &
      routineP = moduleN//':'//routineN

    restraint_env%energy = 0.0_dp
    restraint_env%harmonic_env%energy = 0.0_dp
    restraint_env%dist_env%energy = 0.0_dp
    restraint_env%angle_env%energy = 0.0_dp
    restraint_env%torsion_env%energy = 0.0_dp

    IF(restraint_env%fp_env%in_use) THEN
      CALL fp_eval(restraint_env%fp_env, subsys, cell, error)
      restraint_env%energy = restraint_env%energy + restraint_env%fp_env%energy
    END IF

    IF(restraint_env%harmonic_env%in_use) THEN
      CALL harmonic_eval(restraint_env%harmonic_env, subsys, cell, error)
      restraint_env%energy = restraint_env%energy + restraint_env%harmonic_env%energy
    END IF

    IF(restraint_env%dist_env%in_use) THEN
      CALL dist_eval(restraint_env%dist_env, subsys, cell, error)
      restraint_env%energy = restraint_env%energy + restraint_env%dist_env%energy
    END IF

    IF(restraint_env%angle_env%in_use) THEN
      CALL angle_eval(restraint_env%angle_env, subsys, cell, error)
      restraint_env%energy = restraint_env%energy + restraint_env%angle_env%energy
    END IF

    IF(restraint_env%torsion_env%in_use) THEN
      CALL torsion_eval(restraint_env%torsion_env, subsys, cell, error)
      restraint_env%energy = restraint_env%energy + restraint_env%torsion_env%energy
    END IF

  END SUBROUTINE restraint_eval



!!****f* restraint_methods/fp_eval *
!!
!!   NAME
!!    fp_eval
!!
!!   FUNCTION
!!    computest the forces and the energy due to the flexible potential & bias, 
!!    and writes the weights file
!!
!!   NOTES
!!    currently limited to 1 subsystem
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     04.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE fp_eval(fp_env, subsys, cell, error)
    TYPE(fp_type), POINTER                   :: fp_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fp_eval', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=15)                        :: tmpstr
    INTEGER                                  :: handle, i, icenter, &
                                                iparticle, output_unit
    LOGICAL                                  :: failure, zero_weight
    REAL(KIND=dp)                            :: c, dcdr, kT, r, rab(3), sf, &
                                                strength
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_list_type), POINTER        :: particles_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles

    CALL timeset(routineN,handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(fp_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(fp_env%in_use,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    ! limitation can be easily removed
    CPPrecondition(SIZE(subsys,1)==1,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys(1)%subsys,particles=particles_list,error=error)
    particles=>particles_list%els

    ! compute the force due to the reflecting walls
    ! and count the distribution in discrete and contiguous ways
    zero_weight=.FALSE.
    fp_env%restraint_energy=0.0_dp
    icenter=fp_env%central_atom
    strength=fp_env%strength
    fp_env%i1=0 ; fp_env%i2=0 ; fp_env%o1=0 ; fp_env%o2=0
    fp_env%ri1=0.0_dp ; fp_env%ri2=0.0_dp ; fp_env%ro1=0.0_dp ; fp_env%ro2=0.0_dp
    fp_env%energy=0.0_dp

    ! inner particles 
    DO i=1,SIZE(fp_env%inner_atoms)
       iparticle=fp_env%inner_atoms(i)
       rab=particles(iparticle)%r-particles(icenter)%r
       rab=pbc(rab,cell)
       r=SQRT(SUM(rab**2))
       ! constraint wall  (they feel to outer wall)
       IF (r>fp_env%outer_radius) THEN
          zero_weight=.TRUE.
          fp_env%restraint_energy=fp_env%restraint_energy+0.5_dp*strength*(r-fp_env%outer_radius)**2
          sf=strength*(r-fp_env%outer_radius)/r
          particles(iparticle)%f=particles(iparticle)%f-sf*rab
          particles(icenter)%f  =particles(icenter)%f  +sf*rab
       ENDIF
       ! count the distribution
       IF (r>fp_env%inner_radius) THEN
          fp_env%i2=fp_env%i2+1 
       ELSE
          fp_env%i1=fp_env%i1+1 
       ENDIF
       ! smooth count the distribution
       CALL smooth_count(r,fp_env%inner_radius,fp_env%smooth_width,c,dcdr)
       fp_env%ri1=fp_env%ri1+c
       fp_env%ri2=fp_env%ri2+(1.0_dp-c)
    ENDDO

    ! outer particles 
    DO i=1,SIZE(fp_env%outer_atoms)
       iparticle=fp_env%outer_atoms(i)
       rab=particles(iparticle)%r-particles(icenter)%r
       rab=pbc(rab,cell)
       r=SQRT(SUM(rab**2))
       ! constraint wall (they feel the inner wall)
       IF (r<fp_env%inner_radius) THEN
          zero_weight=.TRUE.
          fp_env%restraint_energy=fp_env%restraint_energy+ &
                              0.5_dp*strength*(r-fp_env%inner_radius)**2
          sf=strength*(r-fp_env%inner_radius)/r
          particles(iparticle)%f=particles(iparticle)%f-sf*rab
          particles(icenter)%f  =particles(icenter)%f  +sf*rab
       ENDIF
       ! count the distribution
       IF (r>fp_env%outer_radius) THEN
          fp_env%o2=fp_env%o2+1 
       ELSE
          fp_env%o1=fp_env%o1+1 
       ENDIF
       ! smooth count the distribution
       CALL smooth_count(r,fp_env%outer_radius,fp_env%smooth_width,c,dcdr)
       fp_env%ro1=fp_env%ro1+c
       fp_env%ro2=fp_env%ro2+(1.0_dp-c)
    ENDDO
    fp_env%energy=fp_env%energy+fp_env%restraint_energy

    ! the combinatorial weight
    i=fp_env%i2+fp_env%o1
    CPPrecondition(i<=maxfac,cp_failure_level,routineP,error,failure)
    fp_env%comb_weight=(fac(fp_env%i2)*fac(fp_env%o1))/fac(i)

    ! we can add the bias potential now.
    ! this bias has the form
    ! kT * { ln[(o1+i2)!] - ln[o1!] - ln[i2!] }
    ! where the smooth counts are used for o1 and i2
    fp_env%bias_energy=0.0_dp
    IF (fp_env%bias) THEN
       kT=fp_env%temperature
       fp_env%bias_energy=kT*(gamln(fp_env%ro1+fp_env%ri2+1) - &
                              gamln(fp_env%ro1+1)-gamln(fp_env%ri2+1))

       ! and add the corresponding forces
       ! inner particles
       DO i=1,SIZE(fp_env%inner_atoms)
          iparticle=fp_env%inner_atoms(i)
          rab=particles(iparticle)%r-particles(icenter)%r
          rab=pbc(rab,cell)
          r=SQRT(SUM(rab**2))
          CALL smooth_count(r,fp_env%inner_radius,fp_env%smooth_width,c,dcdr)
          sf=kT*(psi(fp_env%ro1+fp_env%ri2+1)-psi(fp_env%ri2+1))*(-dcdr)/r
          particles(iparticle)%f=particles(iparticle)%f-sf*rab
          particles(icenter)%f  =particles(icenter)%f  +sf*rab
       ENDDO
       ! outer particles
       DO i=1,SIZE(fp_env%outer_atoms)
          iparticle=fp_env%outer_atoms(i)
          rab=particles(iparticle)%r-particles(icenter)%r
          rab=pbc(rab,cell)
          r=SQRT(SUM(rab**2))
          CALL smooth_count(r,fp_env%outer_radius,fp_env%smooth_width,c,dcdr)
          sf=kT*(psi(fp_env%ro1+fp_env%ri2+1)-psi(fp_env%ro1+1))*(dcdr)/r
          particles(iparticle)%f=particles(iparticle)%f-sf*rab
          particles(icenter)%f  =particles(icenter)%f  +sf*rab
       ENDDO
    ENDIF
    fp_env%energy=fp_env%energy+fp_env%bias_energy
    fp_env%bias_weight=EXP(fp_env%bias_energy/kT)

    ! if this configuration is a valid one, compute its weight
    IF (zero_weight) THEN
        fp_env%weight=0.0_dp 
    ELSE
        fp_env%weight=fp_env%comb_weight*fp_env%bias_weight
    ENDIF

    ! put weights and other info on file
    logger => cp_error_get_logger(error)
    output_unit=cp_print_key_unit_nr(logger,fp_env%print_section,"", &
                                     extension=".weights",error=error)
    IF (output_unit>0) THEN
       tmpstr=cp_iter_string(logger%iter_info, fp_env%print_section)
       WRITE(output_unit,'(T2,A15,6(1X,F16.10),4(1X,I4),4(1X,F16.10))') &
         tmpstr, &
         fp_env%weight, fp_env%comb_weight, fp_env%bias_weight, &
         fp_env%energy, fp_env%restraint_energy, fp_env%bias_energy, &
         fp_env%i1,fp_env%i2,fp_env%o1,fp_env%o2, &
         fp_env%ri1,fp_env%ri2,fp_env%ro1,fp_env%ro2
    ENDIF

    CALL cp_print_key_finished_output(output_unit,logger,fp_env%print_section,&
                                      "",error=error)

    CALL timestop(handle)

  END SUBROUTINE fp_eval

!!****f* restraint_methods/smooth_count *
!!
!!   NAME
!!     smooth_count
!!
!!   FUNCTION
!!     counts in a smooth way (error function with width=width)
!!     if r is closer than r1. Returns 1.0 for the count=c if r<<r1
!!     and the derivative wrt r dcdr
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     04.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE smooth_count(r,r1,width,c,dcdr)
    REAL(KIND=dp), INTENT(IN)                :: r, r1, width
    REAL(KIND=dp), INTENT(OUT)               :: c, dcdr

    REAL(KIND=dp)                            :: arg

    arg=(r1-r)/width

    c   = (1.0_dp + erf(arg))/2.0_dp
    dcdr= (-oorootpi/width)*EXP(-arg**2)

  END SUBROUTINE


!!****f* restraint_methods/harmonic_eval *
!!
!!   NAME
!!    harmonic_eval
!!
!!   FUNCTION
!!    computest the forces and the energy due to the harmonic potential
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [ikuo]
!!    
!!
!!*** **********************************************************************
  SUBROUTINE harmonic_eval(harmonic_env,subsys,cell,error)
    TYPE(harmonic_type), POINTER             :: harmonic_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harmonic_eval', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: handle
    CHARACTER(len=default_string_length)     :: name
    LOGICAL                                  :: failure
    TYPE(particle_list_type), POINTER        :: particles_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    INTEGER :: ires,nres,iatom
    REAL(KIND=dp) ::rab(3),k(3),r0,energy,fscalar


    CALL timeset(routineN,handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(harmonic_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(harmonic_env%in_use,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    ! limitation can be easily removed
    CPPrecondition(SIZE(subsys,1)==1,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys(1)%subsys,particles=particles_list,error=error)
    particles=>particles_list%els

    nres = harmonic_env % nrestraint

    IF(harmonic_env%first_time) THEN
      DO ires=1,nres
        iatom=harmonic_env%atom_a(ires)
        harmonic_env%r_orig(:,ires) = particles(iatom)%r(:)
      END DO
      harmonic_env%first_time = .FALSE.
    ELSE
      name = "HARMONIC"
      DO ires=1,nres
        iatom=harmonic_env%atom_a(ires)
        k(:)=0.0_dp
        k(1)=harmonic_env%k
        r0=0.0_dp
        rab=particles(iatom)%r-harmonic_env%r_orig(:,iatom)
        CALL force_bonds(name,rab,r0,k,energy,fscalar)
        particles(iatom)%f(1)=particles(iatom)%f(1)-rab(1)*fscalar
        particles(iatom)%f(2)=particles(iatom)%f(2)-rab(2)*fscalar
        particles(iatom)%f(3)=particles(iatom)%f(3)-rab(3)*fscalar
        harmonic_env%energy = harmonic_env%energy + energy
      END DO
    END IF

    CALL timestop(handle)



  END SUBROUTINE harmonic_eval

!!****f* restraint_methods/dist_eval *
!!
!!   NAME
!!    dist_eval
!!
!!   FUNCTION
!!    computest the forces and the energy due to the dist potential
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [ikuo]
!!    
!!
!!*** **********************************************************************
  SUBROUTINE dist_eval(dist_env,subsys,cell,error)
    TYPE(dist_type), POINTER             :: dist_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dist_eval', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: name
    LOGICAL                                  :: failure
    INTEGER                                  :: handle
    TYPE(particle_list_type), POINTER        :: particles_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    INTEGER :: ires,nres,iatom,jatom
    REAL(KIND=dp) ::rab(3),k(3),r0,energy,fscalar

    CALL timeset(routineN,handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(dist_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(dist_env%in_use,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    ! limitation can be easily removed
    CPPrecondition(SIZE(subsys,1)==1,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys(1)%subsys,particles=particles_list,error=error)
    particles=>particles_list%els

    name = "HARMONIC"

    nres = dist_env % nrestraint
    DO ires=1,nres
      iatom=dist_env%atom_a(ires)
      jatom=dist_env%atom_b(ires)
      k=dist_env%k(ires)
      r0=dist_env%r0(ires)
      rab=particles(iatom)%r-particles(jatom)%r
      CALL force_bonds(name,rab,r0,k,energy,fscalar)
      particles(iatom)%f(1)=particles(iatom)%f(1)-rab(1)*fscalar
      particles(iatom)%f(2)=particles(iatom)%f(2)-rab(2)*fscalar
      particles(iatom)%f(3)=particles(iatom)%f(3)-rab(3)*fscalar
      particles(jatom)%f(1)=particles(jatom)%f(1)+rab(1)*fscalar
      particles(jatom)%f(2)=particles(jatom)%f(2)+rab(2)*fscalar
      particles(jatom)%f(3)=particles(jatom)%f(3)+rab(3)*fscalar
      dist_env%energy = dist_env%energy + energy
    END DO

    CALL timestop(handle)

  END SUBROUTINE dist_eval


!!****f* restraint_methods/angle_eval *
!!
!!   NAME
!!    angle_eval
!!
!!   FUNCTION
!!    computest the forces and the energy due to the angle potential
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [ikuo]
!!    
!!
!!*** **********************************************************************
  SUBROUTINE angle_eval(angle_env,subsys,cell,error)
    TYPE(angle_type), POINTER             :: angle_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'angle_eval', &
      routineP = moduleN//':'//routineN
    CHARACTER(len=default_string_length)     :: name
    LOGICAL                                  :: failure

    INTEGER                                  :: handle
    TYPE(particle_list_type), POINTER        :: particles_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles

    INTEGER :: index_a,index_b,index_c
    INTEGER :: ires,nres
    REAL(KIND=dp) :: b12(3),b32(3),d12,id12,d32,id32,dist,theta,theta0,k
    REAL(KIND=dp) :: energy,fscalar
    REAL(KIND=dp) :: g1(3),g2(3),g3(3)




    CALL timeset(routineN,handle)


    failure=.FALSE.

    CPPrecondition(ASSOCIATED(angle_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(angle_env%in_use,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    ! limitation can be easily removed
    CPPrecondition(SIZE(subsys,1)==1,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys(1)%subsys,particles=particles_list,error=error)
    particles=>particles_list%els

    name = "HARMONIC"

    nres = angle_env % nrestraint
    DO ires=1,nres
        index_a = angle_env%atom_a(ires)
        index_b = angle_env%atom_b(ires)
        index_c = angle_env%atom_c(ires)
        b12 = particles ( index_a ) % r - particles ( index_b ) % r
        b32 = particles ( index_c ) % r - particles ( index_b ) % r
        d12 = SQRT(DOT_PRODUCT(b12,b12))
        id12 = 1.0_dp/d12
        d32 = SQRT(DOT_PRODUCT(b32,b32))
        id32 = 1.0_dp/d32
        dist = DOT_PRODUCT(b12,b32)
        theta = (dist*id12*id32)
        IF(theta<-1.0_dp) theta=-1.0_dp
        IF(theta>+1.0_dp) theta=+1.0_dp
        theta = ACOS(theta)
        k=angle_env%k(ires)
        theta0=angle_env%theta0(ires)
        CALL force_bends( name, &
                         b12,b32,d12,d32,id12,id32,dist,theta, &
                         theta0,k,&
                         g1,g2,g3,energy,fscalar)
        angle_env%energy = angle_env%energy + energy
        particles ( index_a )% f(1) = particles ( index_a ) % f(1) + g1(1) * fscalar
        particles ( index_a )% f(2) = particles ( index_a ) % f(2) + g1(2) * fscalar
        particles ( index_a )% f(3) = particles ( index_a ) % f(3) + g1(3) * fscalar
        particles ( index_b )% f(1) = particles ( index_b ) % f(1) + g2(1) * fscalar
        particles ( index_b )% f(2) = particles ( index_b ) % f(2) + g2(2) * fscalar
        particles ( index_b )% f(3) = particles ( index_b ) % f(3) + g2(3) * fscalar
        particles ( index_c )% f(1) = particles ( index_c ) % f(1) + g3(1) * fscalar
        particles ( index_c )% f(2) = particles ( index_c ) % f(2) + g3(2) * fscalar
        particles ( index_c )% f(3) = particles ( index_c ) % f(3) + g3(3) * fscalar
    END DO



    CALL timestop(handle)

  END SUBROUTINE angle_eval


!!****f* restraint_methods/torsion_eval *
!!
!!   NAME
!!    torsion_eval
!!
!!   FUNCTION
!!    computest the forces and the energy due to the torsion potential
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [ikuo]
!!    
!!
!!*** **********************************************************************
  SUBROUTINE torsion_eval(torsion_env,subsys,cell,error)
    TYPE(torsion_type), POINTER             :: torsion_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'torsion_eval', &
      routineP = moduleN//':'//routineN
    CHARACTER(len=default_string_length)     :: name
    LOGICAL                                  :: failure

    INTEGER                                  :: handle
    TYPE(particle_list_type), POINTER        :: particles_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(0:3)            :: a
        

    INTEGER :: index_a,index_b,index_c,index_d
    INTEGER :: ires,nres,imul
    REAL(KIND=dp) :: t12(3),t32(3),t43(3),s32,is32,dist1,dist2,u1(3),u2(3),w1,w2,iw1,iw2,k,phi0
    REAL(KIND=dp) :: gt1(3),gt2(3),gt3(3),gt4(3)
    REAL(KIND=dp) :: energy,fscalar



    CALL timeset(routineN,handle)


    failure=.FALSE.

    CPPrecondition(ASSOCIATED(torsion_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(torsion_env%in_use,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    ! limitation can be easily removed
    CPPrecondition(SIZE(subsys,1)==1,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys(1)%subsys,particles=particles_list,error=error)
    particles=>particles_list%els

    name = "CHARMM"
    a(:)=0.0_dp

    nres = torsion_env % nrestraint
    DO ires=1,nres


        index_a = torsion_env%atom_a(ires)
        index_b = torsion_env%atom_b(ires)
        index_c = torsion_env%atom_c(ires)
        index_d = torsion_env%atom_d(ires)
        t12 = particles ( index_a ) % r - particles ( index_b ) % r
        t32 = particles ( index_c ) % r - particles ( index_b ) % r
        t43 = particles ( index_d ) % r - particles ( index_c ) % r
        s32 = SQRT(DOT_PRODUCT(t32,t32))
        is32 = 1.0_dp/s32
        dist1 = DOT_PRODUCT(t12,t32)
        dist2 = DOT_PRODUCT(t43,t32)
        u1 = t12 - dist1*t32*is32**2
        u2 = t43 - dist2*t32*is32**2
        w1 = SQRT(DOT_PRODUCT(u1,u1))
        w2 = SQRT(DOT_PRODUCT(u2,u2))
        iw1 = 1.0_dp/w1
        iw2 = 1.0_dp/w2
        k=torsion_env%k(ires)
        phi0=COS(torsion_env%phi0(ires))
        imul=torsion_env%imul(ires)
            CALL force_torsions(name, &
                              is32,dist1,dist2,u1,u2,iw1,iw2, &
                              a, &
                              k, &
                              phi0, &
                              imul, &
                              gt1,gt2,gt3,gt4,energy,fscalar)
            torsion_env%energy = torsion_env%energy + energy
            particles ( index_a )% f(1) = particles ( index_a ) % f(1) + gt1(1) * fscalar
            particles ( index_a )% f(2) = particles ( index_a ) % f(2) + gt1(2) * fscalar
            particles ( index_a )% f(3) = particles ( index_a ) % f(3) + gt1(3) * fscalar
            particles ( index_b )% f(1) = particles ( index_b ) % f(1) + gt2(1) * fscalar
            particles ( index_b )% f(2) = particles ( index_b ) % f(2) + gt2(2) * fscalar
            particles ( index_b )% f(3) = particles ( index_b ) % f(3) + gt2(3) * fscalar
            particles ( index_c )% f(1) = particles ( index_c ) % f(1) + gt3(1) * fscalar
            particles ( index_c )% f(2) = particles ( index_c ) % f(2) + gt3(2) * fscalar
            particles ( index_c )% f(3) = particles ( index_c ) % f(3) + gt3(3) * fscalar
            particles ( index_d )% f(1) = particles ( index_d ) % f(1) + gt4(1) * fscalar
            particles ( index_d )% f(2) = particles ( index_d ) % f(2) + gt4(2) * fscalar
            particles ( index_d )% f(3) = particles ( index_d ) % f(3) + gt4(3) * fscalar


    END DO

    CALL timestop(handle)

  END SUBROUTINE torsion_eval


END MODULE restraint_methods
