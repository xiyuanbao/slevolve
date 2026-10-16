!=======================================================================
! make_synthetic_case -- Fortran twin of tools/make_synthetic_case.py: writes the same synthetic Earth
! (Earth-like random topography, a polar and a mid-latitude Vialov ice sheet melting over 8 steps) as
! an INPUT_FILES tree, with the shipped .npy reader/writer (fortran/src/npy_io_mod.f90).
!
!   make_synthetic_case --outdir DIR [--nglv 512] [--love FILE] [--seed FILE] [--topography FILE]
!
!   --outdir      INPUT_FILES directory to fill (DIR/grid, icemodel, earthmodel, others are created)
!   --nglv        number of Gauss-Legendre latitudes N (grid N x 2N; default 512)
!   --love        Love-number file to copy into DIR/earthmodel/ (default lovenum/prem_512.l96C.ump3.lm1
!                 relative to the current directory)
!   --seed        topography seed, 512 x 1024 .npy (default reference/synthetic/topography_seed.npy
!                 relative to the current directory)
!   --topography  (N, 2N) .npy bedrock grid to use instead of the synthetic topography
!
! The topography is not synthesised here: the random field, its hypsometric transfer and the fractal
! detail are built once in Python (build_topography_seed in tools/make_synthetic_case.py, with the
! shipped spherical-harmonic transform and numpy's random generator), rounded to whole millimetres and
! shipped as reference/synthetic/topography_seed.npy.  Both generators read that file (and sample it
! at the nearest grid point when N /= 512), so the topographies are identical by construction.  The
! ice (two Vialov profiles H (1 - (r/R)^(4/3))^(3/8), only where the topography is above -500 m) is
! computed here as in the Python script and rounded to whole millimetres, which makes the files
! identical.
! Build (tests/test_fortran_synthetic.sh does this):
!   gfortran -O2 -ffree-line-length-none fortran/src/npy_io_mod.f90 tools/make_synthetic_case.f90 -o make_synthetic_case
!=======================================================================
program make_synthetic_case
   use iso_fortran_env, only: int64, real64
   use npy_io_mod, only: read_grid, write_grid
   implicit none
   integer, parameter :: dp = real64, NSTEP = 8, NS = 512
   integer(int64), parameter :: SEED = 20261013_int64
   real(dp), parameter :: PI = 3.141592653589793_dp, D2R = PI / 180.0_dp, ICE_MIN_BED = -500.0_dp
   real(dp), parameter :: POLAR(0:NSTEP) = [3000._dp, 2625._dp, 2250._dp, 1875._dp, 1500._dp, 1125._dp, 750._dp, 375._dp, 0._dp]
   real(dp), parameter :: MIDLAT(0:NSTEP) = [2500._dp, 2400._dp, 2150._dp, 1800._dp, 1350._dp, 900._dp, 500._dp, 200._dp, 0._dp]
   ! ice sheets: centre latitude, longitude, radius (deg); polar sheet centred on the south pole
   real(dp), parameter :: R1 = 17.0_dp, LAT2 = 56.0_dp, LON2 = 132.0_dp, R2 = 12.0_dp
   integer :: nglv, nlon, nargs, i, j, k, ios, u, jj, ii
   character(len=4096) :: arg, outdir, love, topofile, seedfile
   character(len=32) :: kstr
   real(dp), allocatable :: lat(:), lon(:), topo(:,:), ice(:,:), polar_s(:,:), mid_s(:,:), x(:), w(:)
   real(dp), allocatable :: f0(:,:), seedbuf(:,:), xs(:), ws(:), lats(:)
   integer, allocatable :: rows(:)
   real(dp) :: landfrac, p, p0, dl, cc, d, best
   logical :: ex

   nglv = 512; outdir = ''; love = 'lovenum/prem_512.l96C.ump3.lm1'; topofile = ''
   seedfile = 'reference/synthetic/topography_seed.npy'
   nargs = command_argument_count(); i = 1
   do while (i <= nargs)
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('--outdir');     i = i + 1; call get_command_argument(i, outdir)
      case ('--nglv');       i = i + 1; call get_command_argument(i, arg); read(arg, *) nglv
      case ('--love');       i = i + 1; call get_command_argument(i, love)
      case ('--seed');       i = i + 1; call get_command_argument(i, seedfile)
      case ('--topography'); i = i + 1; call get_command_argument(i, topofile)
      case default; call usage(); stop 2
      end select
      i = i + 1
   end do
   if (len_trim(outdir) == 0) then
      call usage(); stop 2
   end if
   nlon = 2 * nglv

   ! ---- grid ------------------------------------------------------------
   allocate(x(nglv), w(nglv), lat(nglv), lon(nlon))
   call gauleg(nglv, x, w)                   ! ascending nodes, as numpy leggauss
   do j = 1, nglv
      lat(j) = asin(x(nglv + 1 - j)) * (180.0_dp / PI)   ! north to south
   end do
   do i = 1, nlon
      lon(i) = real(i - 1, dp) * 360.0_dp / real(nlon, dp)
   end do

   ! ---- topography ------------------------------------------------------
   allocate(topo(nglv, nlon))
   if (len_trim(topofile) > 0) then
      allocate(f0(nlon, nglv))
      call read_grid(f0, trim(strip_npy(topofile)), ios)
      if (ios /= 0) then; print '(A)', 'cannot read '//trim(topofile); stop 1; end if
      topo = transpose(f0)
   else
      ! seed(lat j, lon i) = seedbuf(i, j): the .npy payload is row-major, longitude fastest
      allocate(seedbuf(2 * NS, NS), xs(NS), ws(NS), lats(NS), rows(nglv))
      inquire(file=trim(seedfile), exist=ex)
      if (.not. ex) then; print '(A)', 'topography seed not found: '//trim(seedfile)//' (--seed)'; stop 1; end if
      call read_grid(seedbuf, trim(strip_npy(seedfile)), ios)
      if (ios /= 0) then; print '(A)', 'cannot read '//trim(seedfile); stop 1; end if
      call gauleg(NS, xs, ws)
      do j = 1, NS
         lats(j) = asin(xs(NS + 1 - j)) * (180.0_dp / PI)
      end do
      do j = 1, nglv                       ! nearest seed row in latitude (first one on a tie, as argmin)
         rows(j) = 1; best = abs(lats(1) - lat(j))
         do jj = 2, NS
            if (abs(lats(jj) - lat(j)) < best) then; best = abs(lats(jj) - lat(j)); rows(j) = jj; end if
         end do
      end do
      do i = 1, nlon                       ! nearest seed column: (2 i NS + N) div (2 N), i from 0
         ii = mod((2 * (i - 1) * NS + nglv) / (2 * nglv), 2 * NS) + 1
         do j = 1, nglv
            topo(j, i) = seedbuf(ii, rows(j))
         end do
      end do
   end if

   ! ---- ice: Vialov profiles --------------------------------------------
   allocate(ice(nglv, nlon), polar_s(nglv, nlon), mid_s(nglv, nlon))
   p0 = LAT2 * D2R
   do i = 1, nlon
      dl = (lon(i) - LON2) * D2R
      do j = 1, nglv
         polar_s(j, i) = vialov(90.0_dp + lat(j), R1)
         p = lat(j) * D2R
         cc = sin(p) * sin(p0) + cos(p) * cos(p0) * cos(dl)
         d = acos(min(max(cc, -1.0_dp), 1.0_dp)) / D2R
         mid_s(j, i) = vialov(d, R2)
      end do
   end do

   ! ---- write -----------------------------------------------------------
   ! Fields are held as f(lat, lon); the .npy payload (and the solver's grids) are in row-major
   ! order, longitude fastest, so they are written as reshape(transpose(f), [nglv, nlon]).
   call execute_command_line('mkdir -p "'//trim(outdir)//'/grid" "'//trim(outdir)//'/icemodel" "'// &
        trim(outdir)//'/earthmodel" "'//trim(outdir)//'/others"')
   write(kstr, '(I0)') nglv
   open(newunit=u, file=trim(outdir)//'/grid/GLlat_'//trim(kstr)//'.txt', status='replace')
   do j = 1, nglv; write(u, '(ES16.7)') lat(j); end do
   close(u)
   open(newunit=u, file=trim(outdir)//'/grid/GLlon_'//trim(kstr)//'.txt', status='replace')
   do i = 1, nlon; write(u, '(ES16.7)') lon(i); end do
   close(u)
   call write_grid(reshape(transpose(topo), [nglv, nlon]), trim(outdir)//'/others/synth_topo', ios); call chk(ios, 'synth_topo')
   do k = 0, NSTEP
      ice = POLAR(k) * polar_s + MIDLAT(k) * mid_s
      where (topo <= ICE_MIN_BED) ice = 0.0_dp      ! ice on land or shelf only
      call round_mm(ice)
      write(kstr, '(I0)') k
      call write_grid(reshape(transpose(ice), [nglv, nlon]), trim(outdir)//'/icemodel/synth_ice_'//trim(kstr), ios); call chk(ios, 'synth_ice')
   end do
   inquire(file=trim(love), exist=ex)
   if (.not. ex) then; print '(A)', 'Love-number file not found: '//trim(love)//' (--love)'; stop 1; end if
   call copy_file(trim(love), trim(outdir)//'/earthmodel/'//trim(basename(love)))

   ! ---- report ----------------------------------------------------------
   landfrac = 0.0_dp
   do i = 1, nlon
      do j = 1, nglv
         if (topo(j, i) >= 0.0_dp) landfrac = landfrac + w(nglv + 1 - j)
      end do
   end do
   landfrac = landfrac / (sum(w) * nlon)
   print '(A,I0,A)', 'seed ', SEED, ': random topography (reference/synthetic/topography_seed.npy), two Vialov ice sheets'
   print '(A,F0.1,A,I0,A)', '  polar ice sheet: centre -90.0 N 0.0 E, radius ', R1, ' deg, H0 ', nint(POLAR(0)), ' m'
   print '(A,F0.1,A,F0.1,A,F0.1,A,I0,A)', '  mid-latitude ice sheet: centre ', LAT2, ' N ', LON2, ' E, radius ', R2, &
        ' deg, H0 ', nint(MIDLAT(0)), ' m'
   print '(A,9I7)', '  step:           ', (k, k = 0, NSTEP)
   print '(A,9I7)', '  polar ice (m):  ', nint(POLAR)
   print '(A,9I7)', '  mid-lat ice (m):', nint(MIDLAT)
   print '(A,I0,A,I0,A,F5.3,A)', 'wrote '//trim(outdir)//': grid ', nglv, 'x', nlon, ', land fraction ', landfrac, &
        ', 9 ice files synth_ice_0..8, topography others/synth_topo.npy, earthmodel/'//trim(basename(love))

contains

   subroutine usage()
      print '(A)', 'usage: make_synthetic_case --outdir DIR [--nglv N] [--love FILE] [--seed FILE] [--topography FILE]'
   end subroutine usage

   ! Vialov profile shape (1 - (d/R)^(4/3))^(3/8) inside d < R, 0 outside (as vialov() in the Python script)
   real(dp) function vialov(dd, R)
      real(dp), intent(in) :: dd, R
      vialov = max(1.0_dp - (min(dd, R) / R)**(4.0_dp / 3.0_dp), 0.0_dp)**0.375_dp
   end function vialov

   ! round to whole millimetres (np.round(1000 a) / 1000 in the Python script)
   subroutine round_mm(a)
      real(dp), intent(inout) :: a(:,:)
      a = anint(a * 1000.0_dp) / 1000.0_dp
   end subroutine round_mm

   ! Gauss-Legendre nodes (ascending) and weights on [-1, 1] by Newton iteration
   subroutine gauleg(n, xx, ww)
      integer, intent(in) :: n
      real(dp), intent(out) :: xx(n), ww(n)
      integer :: m, it, l
      real(dp) :: z, z1, p1, p2, p3, pp
      do m = 1, (n + 1) / 2
         z = cos(PI * (real(m, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))
         do it = 1, 100
            p1 = 1.0_dp; p2 = 0.0_dp
            do l = 1, n
               p3 = p2; p2 = p1
               p1 = ((2.0_dp * l - 1.0_dp) * z * p2 - (l - 1.0_dp) * p3) / l
            end do
            pp = n * (z * p1 - p2) / (z * z - 1.0_dp)
            z1 = z; z = z1 - p1 / pp
            if (abs(z - z1) < 1.0e-15_dp) exit
         end do
         xx(m) = -z; xx(n + 1 - m) = z
         ww(m) = 2.0_dp / ((1.0_dp - z * z) * pp * pp); ww(n + 1 - m) = ww(m)
      end do
   end subroutine gauleg

   subroutine chk(st, what)
      integer, intent(in) :: st
      character(len=*), intent(in) :: what
      if (st /= 0) then; print '(A)', 'cannot write '//what; stop 1; end if
   end subroutine chk

   function strip_npy(s) result(r)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: r
      integer :: n
      r = s; n = len_trim(s)
      if (n > 4) then
         if (s(n-3:n) == '.npy') r = s(1:n-4)
      end if
   end function strip_npy

   function basename(s) result(r)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: r
      r = s(index(trim(s), '/', back=.true.) + 1:)
   end function basename

   subroutine copy_file(src, dst)
      character(len=*), intent(in) :: src, dst
      character(len=:), allocatable :: buf
      integer(int64) :: sz
      integer :: ui, uo
      open(newunit=ui, file=src, access='stream', form='unformatted', status='old', action='read')
      inquire(unit=ui, size=sz)
      allocate(character(len=sz) :: buf)
      read(ui) buf; close(ui)
      open(newunit=uo, file=dst, access='stream', form='unformatted', status='replace', action='write')
      write(uo) buf; close(uo)
   end subroutine copy_file

end program make_synthetic_case
