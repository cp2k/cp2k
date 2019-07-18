!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

PROGRAM parallel_rng_types_TEST
   USE message_passing, ONLY: mp_world_finalize, &
                              mp_world_init, &
                              mp_environ
   USE kinds, ONLY: dp
   USE machine, ONLY: m_walltime, &
                      default_output_unit
   USE parallel_rng_types, ONLY: GAUSSIAN, &
                                 UNIFORM, &
                                 check_rng, &
                                 rng_stream_type, &
                                 rng_stream_type_from_record, &
                                 rng_name_length, &
                                 rng_record_length

   IMPLICIT NONE

   INTEGER                          :: mpi_comm, i, nsamples, nargs, stat, mepos, num_pe
   LOGICAL                          :: ionode
   REAL(KIND=dp)                    :: t, tend, tmax, tmin, tstart, tsum, tsum2
   TYPE(rng_stream_type)            :: rng_stream
   CHARACTER(len=32)                :: arg

   nsamples = 1000
   nargs = command_argument_count()

   IF (nargs .GT. 1) &
      ERROR STOP "Usage: parallel_rng_types_TEST [<int:nsamples>]"

   IF (nargs == 1) THEN
      CALL get_command_argument(1, arg)
      READ (arg, *, iostat=stat) nsamples
      IF (stat /= 0) &
         ERROR STOP "Usage: parallel_rng_types_TEST [<int:nsamples>]"
   END IF

   CALL mp_world_init(mpi_comm)
   CALL mp_environ(numtask=num_pe, taskid=mepos, groupid=mpi_comm)
   ionode = mepos == 0

   CALL check_rng(default_output_unit, ionode)

   ! Check performance

   IF (ionode) THEN
      WRITE (UNIT=default_output_unit, FMT="(/,/,T2,A,I10,A)") &
         "Check distributions using", nsamples, " random numbers:"
   END IF

   ! Test uniform distribution [0,1]

   rng_stream = rng_stream_type(name="Test uniform distribution [0,1]", &
                                distribution_type=UNIFORM, &
                                extended_precision=.TRUE.)

   IF (ionode) &
      CALL rng_stream%write(default_output_unit)

   tmax = -HUGE(0.0_dp)
   tmin = +HUGE(0.0_dp)
   tsum = 0.0_dp
   tsum2 = 0.0_dp

   tstart = m_walltime()
   DO i = 1, nsamples
      t = rng_stream%next()
      tsum = tsum + t
      tsum2 = tsum2 + t*t
      IF (t > tmax) tmax = t
      IF (t < tmin) tmin = t
   END DO
   tend = m_walltime()

   IF (ionode) THEN
      CALL rng_stream%write(default_output_unit, write_all=.TRUE.)
      WRITE (UNIT=default_output_unit, FMT="(/,(T4,A,F12.6))") &
         "Minimum: ", tmin, &
         "Maximum: ", tmax, &
         "Average: ", tsum/REAL(nsamples, KIND=dp), &
         "Variance:", tsum2/REAL(nsamples, KIND=dp), &
         "Time [s]:", tend - tstart
   END IF

   ! Test normal Gaussian distribution

   rng_stream = rng_stream_type(name="Test normal Gaussian distribution", &
                                distribution_type=GAUSSIAN, &
                                extended_precision=.TRUE.)

   IF (ionode) &
      CALL rng_stream%write(default_output_unit)

   tmax = -HUGE(0.0_dp)
   tmin = +HUGE(0.0_dp)
   tsum = 0.0_dp
   tsum2 = 0.0_dp

   tstart = m_walltime()
   DO i = 1, nsamples
      t = rng_stream%next()
      tsum = tsum + t
      tsum2 = tsum2 + t*t
      IF (t > tmax) tmax = t
      IF (t < tmin) tmin = t
   END DO
   tend = m_walltime()

   IF (ionode) THEN
      CALL rng_stream%write(default_output_unit)
      WRITE (UNIT=default_output_unit, FMT="(/,(T4,A,F12.6))") &
         "Minimum: ", tmin, &
         "Maximum: ", tmax, &
         "Average: ", tsum/REAL(nsamples, KIND=dp), &
         "Variance:", tsum2/REAL(nsamples, KIND=dp), &
         "Time [s]:", tend - tstart
   END IF

   IF (ionode) &
      CALL dump_reload_check()

   CALL mp_world_finalize()

CONTAINS
! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dump_reload_check()
      TYPE(rng_stream_type)            :: rng_stream
      CHARACTER(len=rng_record_length) :: rng_record
      REAL(KIND=dp), DIMENSION(3, 2)   :: ig, ig_orig, cg, cg_orig, bg, bg_orig
      CHARACTER(len=rng_name_length)   :: name, name_orig
      CHARACTER(len=*), PARAMETER      :: serialized_string = &
         "qtb_rng_gaussian                         1 F T F   0.0000000000000000E+00&
         &                12.0                12.0                12.0&
         &                12.0                12.0                12.0&
         &                12.0                12.0                12.0&
         &                12.0                12.0                12.0&
         &                12.0                12.0                12.0&
         &                12.0                12.0                12.0"

      WRITE (UNIT=default_output_unit, FMT="(/,/,T2,A)") &
         "Checking dump and load round trip:"

      rng_stream = rng_stream_type(name="Roundtrip for normal Gaussian distrib", &
                                   distribution_type=GAUSSIAN, &
                                   extended_precision=.TRUE.)

      CALL rng_stream%advance(7, 42)
      CALL rng_stream%get(ig=ig_orig, cg=cg_orig, bg=bg_orig, name=name_orig)
      CALL rng_stream%dump(rng_record)

      rng_stream = rng_stream_type_from_record(rng_record)
      CALL rng_stream%get(ig=ig, cg=cg, bg=bg, name=name)

      IF (ANY(ig /= ig_orig) .OR. ANY(cg /= cg_orig) .OR. ANY(bg /= bg_orig) &
          .OR. (name /= name_orig)) &
         ERROR STOP "Stream dump and load roundtrip failed"

      WRITE (UNIT=default_output_unit, FMT="(T4,A)") &
         "Roundtrip successful"

      WRITE (UNIT=default_output_unit, FMT="(/,/,T2,A)") &
         "Checking dumped format:"

      ig(:, :) = 12.0_dp
      rng_stream = rng_stream_type(name="qtb_rng_gaussian", &
                                   distribution_type=GAUSSIAN, &
                                   extended_precision=.TRUE., &
                                   seed=ig)

      CALL rng_stream%dump(rng_record)

      WRITE (UNIT=default_output_unit, FMT="(T4,A10,A433)") &
         "EXPECTED:", serialized_string

      WRITE (UNIT=default_output_unit, FMT="(T4,A10,A433)") &
         "GENERATED:", rng_record

      IF (rng_record /= serialized_string) &
         ERROR STOP "Serialized record does not match the expected output"

      WRITE (UNIT=default_output_unit, FMT="(T4,A)") &
         "Serialized record matches the expected output"

   END SUBROUTINE

END PROGRAM parallel_rng_types_TEST
! vim: set ts=3 sw=3 tw=132 :
