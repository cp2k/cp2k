! Test driver for the same numerical module compiled into CP2K.
program wilson_driver
  use kinds, only: dp
  use topology_wilson, only: wilson_step, wilson_spectrum, z2_from_wcc, surface_resolved
  implicit none
  integer :: n, points, loops, i, j, l, p, status, invariant
  complex(dp), allocatable :: product(:, :), overlap(:, :)
  real(dp), allocatable :: wcc(:, :)
  real(dp) :: re, im, singular, berry
  read (*, *) n, points, loops
  allocate(product(n,n), overlap(n,n), wcc(n,loops))
  do l = 1, loops
    product = (0.0_dp, 0.0_dp)
    do i = 1, n
      product(i,i) = (1.0_dp, 0.0_dp)
    end do
    do p = 1, points
      do j = 1, n
        do i = 1, n
          read (*, *) re, im
          overlap(i,j) = cmplx(re, im, dp)
        end do
      end do
      call wilson_step(product, overlap, singular, status, 1.e-10_dp)
      if (status /= 0) then
        write (*,*) 'ERROR', status
        stop 1
      end if
    end do
    call wilson_spectrum(product, wcc(:,l), berry, status)
    if (status /= 0) stop 2
    write (*,'(A,*(1X,ES24.16))') 'WCC', wcc(:,l)
  end do
  if (loops > 1) then
    call z2_from_wcc(wcc, invariant, status, 1.e-7_dp)
    write (*,*) 'Z2', invariant, status
    write (*,*) 'RESOLVED', surface_resolved(wcc)
  end if
end program wilson_driver
