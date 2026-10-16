program make_cache
  ! Convert the text grid inputs into flat little-endian fp64 stream files
  ! (<name>.bin, same element order as the text file).  The parsing uses the
  ! very same list-directed read as the solver, so the cached values are
  ! identical to what the text path would have produced.
  !
  !   make_cache <nglv> <outdir> <file> [<file> ...]
  !
  ! writes <outdir>/<basename>.bin for each input file.
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: nglv, nargs, i, u, ios, p
  character(len=1024) :: arg, outdir, fin, fout, base
  real(dp), allocatable :: a(:,:)
  integer :: c0, c1, cr

  nargs = command_argument_count()
  if (nargs < 3) stop 'usage: make_cache NGLV OUTDIR FILE [FILE...]'
  call get_command_argument(1, arg); read(arg,*) nglv
  call get_command_argument(2, outdir)
  allocate(a(nglv, 2*nglv))
  call system_clock(c0, cr)
  do i = 3, nargs
     call get_command_argument(i, fin)
     open(newunit=u, file=trim(fin), form='formatted', access='sequential', &
          status='old', action='read', iostat=ios)
     if (ios /= 0) then
        write(*,*) 'cannot open ', trim(fin); stop 1
     endif
     read(u,*) a
     close(u)
     base = fin
     p = index(base, '/', back=.true.)
     if (p > 0) base = base(p+1:)
     fout = trim(outdir)//'/'//trim(base)//'.bin'
     open(newunit=u, file=trim(fout), form='unformatted', access='stream', &
          status='replace', action='write')
     write(u) a
     close(u)
  enddo
  call system_clock(c1)
  write(*,'(A,I6,A,F8.2,A)') 'cached ', nargs-2, ' files in ', real(c1-c0)/real(cr), ' s'
end program make_cache
