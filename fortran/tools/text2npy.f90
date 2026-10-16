!=======================================================================
! text2npy -- convert grids in the legacy list-directed grid text format
! (any number of values per line, NLAT*NLON values in stream order) into .npy files.
!
!   text2npy [--f4] [--nlat 512] [--nlon 1024] [--outdir DIR] FILE [FILE ...]
!   text2npy [--f4] ... --dir DIR --prefix PFX --range K0 K1
!           (converts DIR/PFX<K0> ... DIR/PFX<K1>, the <prefix><step> ice-file naming)
!
! Output: <outdir>/<basename>.npy  (default outdir = directory of the input),
! shape (NLAT, NLON), fortran_order False, payload in text stream order
! (numpy a[row, col] = value col of text row for 1024-per-line files).
! Parsing is the solver's own list-directed read into real(8), so the
! values are identical to what read_sl produces.
!   --f4  store '<f4' when every value of the file is exactly representable
!         in float32 (lossless); otherwise fall back to '<f8' and say so.
!=======================================================================
program text2npy
   use iso_fortran_env, only: int64, real64
   use npy_io_mod
   implicit none
   integer :: nlat, nlon, nargs, i, u, ios, p, k0, k1, k, nconv, nf4
   logical :: want_f4, use_f4
   character(len=4096) :: arg, outdir, dir, pfx, fin, fout, base, kstr
   character(len=4096), allocatable :: files(:)
   real(real64), allocatable :: a(:)
   integer(int64) :: c0, c1, cr

   nlat = 512; nlon = 1024; want_f4 = .false.; outdir = ''; dir = ''; pfx = ''
   k0 = -1; k1 = -2
   nargs = command_argument_count()
   allocate(files(nargs)); nconv = 0
   i = 1
   do while (i <= nargs)
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('--f4');     want_f4 = .true.
      case ('--f8');     want_f4 = .false.
      case ('--nlat');   i = i + 1; call get_command_argument(i, arg); read(arg,*) nlat
      case ('--nlon');   i = i + 1; call get_command_argument(i, arg); read(arg,*) nlon
      case ('--outdir'); i = i + 1; call get_command_argument(i, outdir)
      case ('--dir');    i = i + 1; call get_command_argument(i, dir)
      case ('--prefix'); i = i + 1; call get_command_argument(i, pfx)
      case ('--range')
         i = i + 1; call get_command_argument(i, arg); read(arg,*) k0
         i = i + 1; call get_command_argument(i, arg); read(arg,*) k1
      case ('-h', '--help')
         call usage(); stop
      case default
         nconv = nconv + 1; files(nconv) = arg
      end select
      i = i + 1
   enddo
   if (k1 >= k0) then
      deallocate(files); allocate(files(nconv + k1 - k0 + 1))
      nconv = 0
      do k = k0, k1
         write(kstr,'(I0)') k
         nconv = nconv + 1
         files(nconv) = trim(dir)//'/'//trim(pfx)//trim(kstr)
      enddo
   endif
   if (nconv == 0) then
      call usage(); stop 1
   endif
   allocate(a(int(nlat,int64)*nlon))
   call system_clock(c0, cr)
   nf4 = 0
   do i = 1, nconv
      fin = files(i)
      open(newunit=u, file=trim(fin), form='formatted', access='sequential', &
           status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(0,*) 'text2npy: cannot open ', trim(fin); stop 1
      endif
      read(u, *, iostat=ios) a
      close(u)
      if (ios /= 0) then
         write(0,*) 'text2npy: short or bad file ', trim(fin); stop 1
      endif
      p = index(fin, '/', back=.true.)
      base = fin(p+1:)
      if (len_trim(outdir) > 0) then
         fout = trim(outdir)//'/'//trim(base)//'.npy'
      else
         fout = trim(fin)//'.npy'
      endif
      use_f4 = .false.
      if (want_f4) then
         use_f4 = npy_f4_exact(a, size(a))
         if (.not. use_f4) write(0,'(A)') 'text2npy: '//trim(base)// &
              ' is not exactly representable in float32 -> stored as <f8'
      endif
      if (use_f4) then
         call npy_write_r8(fout, a, size(a), [int(nlat,int64), int(nlon,int64)], ios, dtype='f4')
         nf4 = nf4 + 1
      else
         call npy_write_r8(fout, a, size(a), [int(nlat,int64), int(nlon,int64)], ios)
      endif
      if (ios /= 0) then
         write(0,*) 'text2npy: write failed ', trim(fout); stop 1
      endif
   enddo
   call system_clock(c1)
   write(*,'(A,I0,A,I0,A,F8.2,A)') 'text2npy: converted ', nconv, ' files (', nf4, ' as <f4) in ', &
        real(c1-c0)/real(cr), ' s'
contains
   subroutine usage()
      write(*,'(A)') 'usage: text2npy [--f4] [--nlat N] [--nlon M] [--outdir DIR] FILE...'
      write(*,'(A)') '       text2npy [--f4] [--outdir DIR] --dir DIR --prefix PFX --range K0 K1'
   end subroutine usage
end program text2npy
