! npy_maxdiff A.npy B.npy [TOL] -- max |A - B| of two (512 x 1024 by default) .npy grids, read with the
! shipped reader (fortran/src/npy_io_mod.f90).  Prints the value; exit status 0 if it is < TOL
! (default 1e-9), 1 otherwise.  Options --nlat N --nlon M before the file names for other grids.
program npy_maxdiff
   use iso_fortran_env, only: real64
   use npy_io_mod, only: npy_read_r8
   implicit none
   real(real64), allocatable :: a(:), b(:)
   real(real64) :: tol, d
   integer :: nlat, nlon, i, n, ios
   character(len=4096) :: arg, fa, fb
   character(len=256) :: msg
   nlat = 512; nlon = 1024; tol = 1.0e-9_real64; fa = ''; fb = ''
   i = 1; n = command_argument_count()
   do while (i <= n)
      call get_command_argument(i, arg)
      if (trim(arg) == '--nlat') then
         i = i + 1; call get_command_argument(i, arg); read(arg, *) nlat
      else if (trim(arg) == '--nlon') then
         i = i + 1; call get_command_argument(i, arg); read(arg, *) nlon
      else if (len_trim(fa) == 0) then
         fa = arg
      else if (len_trim(fb) == 0) then
         fb = arg
      else
         read(arg, *) tol
      end if
      i = i + 1
   end do
   allocate(a(nlat * nlon), b(nlat * nlon))
   call npy_read_r8(trim(fa), a, size(a), ios, errmsg=msg)
   if (ios /= 0) then; print '(A)', trim(msg); stop 2; end if
   call npy_read_r8(trim(fb), b, size(b), ios, errmsg=msg)
   if (ios /= 0) then; print '(A)', trim(msg); stop 2; end if
   d = maxval(abs(a - b))
   print '(A,ES10.3,A,ES10.3,A,F7.1,A)', 'max |A - B| = ', d, ' m (tolerance ', tol, ' m; max |B| ', maxval(abs(b)), ' m)'
   if (d >= tol) stop 1
end program npy_maxdiff
