!=======================================================================
! npy2text -- write .npy grids back in the legacy grid text formats.
!
!   npy2text [--format es16.9|row-e16.7|row-e.18|row-E18.10] [--ncol 1024]
!           [--outdir DIR] FILE.npy [FILE.npy ...]
!
! Output name: <outdir>/<basename without .npy> (default: next to input).
! Values are written in C-flat ("stream") order, i.e. the order the solver's
! read_sl / write_sl use.  Formats (C printf notation):
!   es16.9      %16.9E, one value per line   (the ES16.9E2 one-value-per-line output format of write_sl)
!   row-e16.7   %16.7e, NCOL per line, no separator  (OSF archive topo/tgrid files)
!   row-e.18    %.18e,  NCOL per line, ' ' separator (np.savetxt default; OSF ice files)
!   row-E18.10  %18.10E, NCOL per line, ' ' separator (nadi-case ice files)
! Any supported dtype is read (f4 promoted exactly to f8 before printing).
!=======================================================================
program npy2text
   use iso_fortran_env, only: int64, real64
   use npy_io_mod
   implicit none
   integer :: nargs, i, j, u, ios, p, ncol, nfile, w, pos
   character(len=16) :: fmtname
   character(len=4096) :: arg, outdir, fin, fout, base
   character(len=4096), allocatable :: files(:)
   character(len=32) :: s
   character(len=:), allocatable :: line
   real(real64), allocatable :: a(:)
   type(npy_header) :: h
   character(len=256) :: msg

   fmtname = 'es16.9'; ncol = 1024; outdir = ''
   nargs = command_argument_count()
   allocate(files(nargs)); nfile = 0
   i = 1
   do while (i <= nargs)
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('--format'); i = i + 1; call get_command_argument(i, arg); fmtname = arg
      case ('--ncol');   i = i + 1; call get_command_argument(i, arg); read(arg,*) ncol
      case ('--outdir'); i = i + 1; call get_command_argument(i, outdir)
      case ('-h', '--help')
         write(*,'(A)') 'usage: npy2text [--format es16.9|row-e16.7|row-e.18|row-E18.10] [--ncol N] [--outdir DIR] FILE.npy...'
         stop
      case default
         nfile = nfile + 1; files(nfile) = arg
      end select
      i = i + 1
   enddo
   select case (trim(fmtname))
   case ('es16.9', 'row-e16.7', 'row-e.18', 'row-E18.10')
   case default
      write(0,*) 'npy2text: unknown format ', trim(fmtname); stop 2
   end select
   allocate(character(len=ncol*27) :: line)
   do i = 1, nfile
      fin = files(i)
      call npy_read_header(fin, h, ios, msg)
      if (ios /= 0) then
         write(0,*) 'npy2text: ', trim(msg); stop 1
      endif
      if (allocated(a)) deallocate(a)
      allocate(a(h%nelem))
      call npy_read_r8(fin, a, size(a), ios, errmsg=msg)
      if (ios /= 0) then
         write(0,*) 'npy2text: ', trim(msg); stop 1
      endif
      p = index(fin, '/', back=.true.)
      base = fin(p+1:)
      if (len_trim(base) > 4) then
         if (base(len_trim(base)-3:len_trim(base)) == '.npy') base = base(1:len_trim(base)-4)
      endif
      if (len_trim(outdir) > 0) then
         fout = trim(outdir)//'/'//trim(base)
      else
         fout = fin(1:p)//trim(base)
      endif
      open(newunit=u, file=trim(fout), form='formatted', access='sequential', &
           status='replace', action='write', recl=ncol*27+8)
      if (trim(fmtname) == 'es16.9') then
         write(u,'(ES16.9E2)') a
      else
         do j = 1, size(a)
            if (mod(j-1, ncol) == 0) pos = 0
            select case (trim(fmtname))
            case ('row-e16.7')
               write(s,'(ES16.7E2)') a(j); call lower_e(s); w = 16
               line(pos+1:pos+w) = s(1:w); pos = pos + w
            case ('row-e.18')
               write(s,'(ES26.18E2)') a(j); s = adjustl(s); call lower_e(s); w = len_trim(s)
               if (mod(j-1, ncol) /= 0) then
                  line(pos+1:pos+1) = ' '; pos = pos + 1
               endif
               line(pos+1:pos+w) = s(1:w); pos = pos + w
            case ('row-E18.10')
               write(s,'(ES18.10E2)') a(j); w = 18
               if (mod(j-1, ncol) /= 0) then
                  line(pos+1:pos+1) = ' '; pos = pos + 1
               endif
               line(pos+1:pos+w) = s(1:w); pos = pos + w
            end select
            if (mod(j, ncol) == 0 .or. j == size(a)) write(u,'(A)') line(1:pos)
         enddo
      endif
      close(u)
   enddo
   write(*,'(A,I0,A)') 'npy2text: wrote ', nfile, ' files ('//trim(fmtname)//')'
contains
   subroutine lower_e(t)
      character(len=*), intent(inout) :: t
      integer :: q
      q = index(t, 'E')
      if (q > 0) t(q:q) = 'e'
   end subroutine lower_e
end program npy2text
