!=======================================================================
! npy_io_mod -- dependency-free reader/writer for NumPy .npy files
! (format versions 1.0, 2.0 and 3.0; little-endian payloads).
!
! Supported dtypes
!   read : '<f8' '<f4' '|i1' '<i1' '|u1' '|b1' '<i4' '<i8' (-> real(8)),
!          '<c16' (-> complex(8))
!   write: '<f8', '<f4' (real(8) data demoted; caller should check
!          npy_f4_exact first), '|i1' (int8 masks), '<c16'
!
! Element-order convention ("stream order")
! -----------------------------------------
! Every read routine returns the elements in the array's *C-order
! flattening*, i.e. the order numpy.ravel(a, order='C') gives.  For a
! file written with fortran_order=False this is the payload as stored;
! for fortran_order=True the payload is transposed on the fly.  The
! values land in the caller's buffer in its storage order.
!
! For the 512x1024 Gauss-Legendre grids this is exactly the order of the text files
! and of the legacy .bin caches, so
!     read(u,*) data_slm            (text)
!     read(u)   data_slm            (.bin stream)
!     call npy_read_r8(file, data_slm, size(data_slm), ...)   (.npy)
! fill data_slm identically.  The writers store the caller's buffer in
! storage order under the header the caller supplies:
!   * grids:  shape=(nglv, 2*nglv), fortran_order=.false.  -> numpy sees
!             a[lat_row, lon], the text-file line/column layout;
!   * spectral fields dS(0:N,0:N): shape=(N+1,N+1), fortran_order=.true.
!             -> numpy a[l, m] == dS(l, m).
! Both choices keep the payload identical to the legacy .bin files.
!
! The header is padded with spaces and terminated by '\n' so that the
! payload starts on a 64-byte boundary (numpy >= 1.26 convention; any
! 16-byte-aligned header is accepted on read).
!=======================================================================
module npy_io_mod
   use iso_fortran_env, only: int8, int16, int32, int64, real32, real64
   implicit none
   private
   ! minimal solver-facing interface (all a solver needs; no other deps)
   public :: read_grid, write_grid, read_ice_cache, write_snapshot
   ! low-level layer
   public :: npy_header, npy_read_header, npy_read_r8, npy_read_c16, &
             npy_write_r8, npy_write_c16, npy_write_i1, npy_f4_exact, &
             npy_header_string

   interface write_snapshot
      module procedure write_snapshot_r8, write_snapshot_c16
   end interface write_snapshot

   integer, parameter :: maxdim = 8

   type :: npy_header
      character(len=8) :: descr = ''
      logical          :: fortran_order = .false.
      integer          :: ndim = 0
      integer(int64)   :: shape(maxdim) = 0
      integer(int64)   :: nelem = 0
      integer(int64)   :: data_offset = 0      ! 0-based byte offset of the payload
      integer          :: major = 0, minor = 0
   end type npy_header

contains

   !====================================================================
   ! Minimal interface.  PATH is the file name without extension, as the
   ! legacy text files are named (e.g. 'OUTPUT_SLM/tgrid5').
   !====================================================================

   ! Read a 2-D grid: PATH.npy (any real dtype, f4 promoted exactly), else
   ! PATH.bin (legacy headerless little-endian f8 stream), else the text
   ! file PATH read list-directed (the legacy grid text format).  All three
   ! fill A in the same (stream) order.  SOURCE returns 'npy'|'bin'|'text'.
   subroutine read_grid(a, path, ios, source, errmsg)
      real(real64), intent(inout) :: a(:,:)
      character(len=*), intent(in) :: path
      integer, intent(out) :: ios
      character(len=*), intent(out), optional :: source, errmsg
      character(len=256) :: msg
      character(len=4) :: src
      logical :: ex
      integer :: u
      msg = ''
      inquire(file=trim(path)//'.npy', exist=ex)
      if (ex) then
         src = 'npy'
         call npy_read_r8(trim(path)//'.npy', a, size(a), ios, errmsg=msg)
      else
         inquire(file=trim(path)//'.bin', exist=ex)
         if (ex) then
            src = 'bin'
            open(newunit=u, file=trim(path)//'.bin', form='unformatted', access='stream', &
                 status='old', action='read', convert='little_endian', iostat=ios)
            if (ios == 0) then
               read(u, iostat=ios) a
               close(u)
            endif
         else
            src = 'text'
            open(newunit=u, file=trim(path), form='formatted', access='sequential', &
                 status='old', action='read', iostat=ios)
            if (ios == 0) then
               read(u, *, iostat=ios) a
               close(u)
            endif
         endif
         if (ios /= 0) msg = 'cannot read '//trim(path)//' ('//trim(src)//')'
      endif
      if (present(source)) source = src
      if (present(errmsg)) errmsg = msg
   end subroutine read_grid

   ! Ice inputs: identical lookup to read_grid; caches are normally '<f4'
   ! (written by text2npy --f4 when lossless).
   subroutine read_ice_cache(a, path, ios, source, errmsg)
      real(real64), intent(inout) :: a(:,:)
      character(len=*), intent(in) :: path
      integer, intent(out) :: ios
      character(len=*), intent(out), optional :: source, errmsg
      character(len=256) :: msg
      character(len=4) :: src
      call read_grid(a, path, ios, src, msg)
      if (present(source)) source = src
      if (present(errmsg)) errmsg = msg
   end subroutine read_ice_cache

   ! Write PATH.npy: '<f8' (or '<f4' with dtype='f4', only if lossless,
   ! else f8), shape (n1,n2), fortran_order=False, payload = storage order.
   subroutine write_grid(a, path, ios, dtype)
      real(real64), intent(in) :: a(:,:)
      character(len=*), intent(in) :: path
      integer, intent(out) :: ios
      character(len=*), intent(in), optional :: dtype
      character(len=3) :: dt
      dt = 'f8'
      if (present(dtype)) then
         if (trim(dtype) == 'f4') then
            if (npy_f4_exact(a, size(a))) dt = 'f4'
         endif
      endif
      call npy_write_r8(trim(path)//'.npy', a, size(a), &
           [int(size(a,1),int64), int(size(a,2),int64)], ios, dtype=dt)
   end subroutine write_grid

   subroutine write_snapshot_r8(a, path, ios)
      real(real64), intent(in) :: a(:,:)
      character(len=*), intent(in) :: path
      integer, intent(out) :: ios
      call write_grid(a, path, ios)
   end subroutine write_snapshot_r8

   ! complex spectral field z(0:N,0:N): '<c16', fortran_order=True, so
   ! numpy a[l,m] == z(l,m); payload = storage order.
   subroutine write_snapshot_c16(z, path, ios)
      complex(real64), intent(in) :: z(:,:)
      character(len=*), intent(in) :: path
      integer, intent(out) :: ios
      call npy_write_c16(trim(path)//'.npy', z, size(z), &
           [int(size(z,1),int64), int(size(z,2),int64)], ios, fortran_order=.true.)
   end subroutine write_snapshot_c16

   !====================================================================
   ! low-level layer
   !====================================================================

   !--------------------------------------------------------------------
   ! header
   !--------------------------------------------------------------------
   subroutine npy_read_header(filename, h, ios, errmsg)
      character(len=*), intent(in) :: filename
      type(npy_header), intent(out) :: h
      integer, intent(out) :: ios
      character(len=*), intent(out), optional :: errmsg
      character(len=6) :: magic
      integer(int8) :: ver(2)
      integer(int16) :: hl2
      integer(int32) :: hl4
      integer(int64) :: hlen
      character(len=:), allocatable :: hs
      character(len=256) :: msg
      integer :: u
      msg = ''
      open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
           status='old', action='read', convert='little_endian', iostat=ios)
      if (ios /= 0) then
         msg = 'cannot open '//trim(filename); goto 900
      endif
      read(u, iostat=ios) magic, ver
      if (ios /= 0 .or. magic /= char(147)//'NUMPY') then
         msg = 'not a .npy file: '//trim(filename); ios = 1; close(u); goto 900
      endif
      h%major = ver(1); h%minor = ver(2)
      if (h%major == 1) then
         read(u, iostat=ios) hl2
         hlen = iand(int(hl2, int64), 65535_int64)
         h%data_offset = 10 + hlen
      else if (h%major == 2 .or. h%major == 3) then
         read(u, iostat=ios) hl4
         hlen = iand(int(hl4, int64), 4294967295_int64)
         h%data_offset = 12 + hlen
      else
         msg = 'unsupported .npy version'; ios = 2; close(u); goto 900
      endif
      allocate(character(len=hlen) :: hs)
      read(u, iostat=ios) hs
      close(u)
      if (ios /= 0) then
         msg = 'truncated .npy header'; goto 900
      endif
      call parse_dict(hs, h, ios, msg)
900   continue
      if (present(errmsg)) errmsg = msg
   end subroutine npy_read_header

   subroutine parse_dict(hs, h, ios, msg)
      character(len=*), intent(in) :: hs
      type(npy_header), intent(inout) :: h
      integer, intent(out) :: ios
      character(len=*), intent(inout) :: msg
      integer :: p, q, r, k
      character(len=64) :: tok
      ios = 0
      ! descr
      p = index(hs, "'descr'")
      if (p == 0) then
         ios = 3; msg = 'header has no descr'; return
      endif
      q = p + 7 + index(hs(p+7:), "'")          ! first char after the opening quote
      r = q - 1 + index(hs(q:), "'")            ! closing quote
      h%descr = hs(q:r-1)
      ! fortran_order
      p = index(hs, "'fortran_order'")
      if (p == 0) then
         ios = 3; msg = 'header has no fortran_order'; return
      endif
      q = p + 15 + index(hs(p+15:), ':')
      tok = adjustl(hs(q:min(len(hs), q+10)))
      h%fortran_order = (tok(1:4) == 'True')
      ! shape
      p = index(hs, "'shape'")
      if (p == 0) then
         ios = 3; msg = 'header has no shape'; return
      endif
      q = p + index(hs(p:), '(')
      r = q - 1 + index(hs(q:), ')')
      h%ndim = 0
      k = q
      do while (k < r)
         p = index(hs(k:r-1), ',')
         if (p == 0) then
            tok = hs(k:r-1); k = r
         else
            tok = hs(k:k+p-2); k = k + p
         endif
         if (len_trim(tok) == 0) cycle
         h%ndim = h%ndim + 1
         if (h%ndim > maxdim) then
            ios = 4; msg = 'too many dimensions'; return
         endif
         read(tok, *, iostat=ios) h%shape(h%ndim)
         if (ios /= 0) then
            msg = 'bad shape entry'; return
         endif
      enddo
      h%nelem = 1
      do k = 1, h%ndim
         h%nelem = h%nelem * h%shape(k)
      enddo
   end subroutine parse_dict

   ! Build the padded header block (magic .. '\n') for a given dict.
   function npy_header_string(descr, fortran_order, shp) result(block)
      character(len=*), intent(in) :: descr
      logical, intent(in) :: fortran_order
      integer(int64), intent(in) :: shp(:)
      character(len=:), allocatable :: block
      character(len=:), allocatable :: d
      character(len=32) :: num
      integer :: k, pre, tot, pad
      integer(int32) :: hl
      d = "{'descr': '"//trim(descr)//"', 'fortran_order': "
      if (fortran_order) then
         d = d//'True'
      else
         d = d//'False'
      endif
      d = d//", 'shape': ("
      do k = 1, size(shp)
         write(num, '(I0)') shp(k)
         if (k > 1) d = d//' '
         d = d//trim(num)//','
      enddo
      if (size(shp) > 1) d = d(1:len(d)-1)       ! (a, b) ; 1-d keeps (a,)
      d = d//'), }'
      pre = 10
      if (len(d) + 1 + 10 > 65535) pre = 12
      tot = ((pre + len(d) + 1 + 63) / 64) * 64
      pad = tot - pre - len(d) - 1
      d = d//repeat(' ', pad)//char(10)
      hl = len(d)
      if (pre == 10) then
         block = char(147)//'NUMPY'//char(1)//char(0)//char(iand(hl,255))//char(iand(ishft(hl,-8),255))//d
      else
         block = char(147)//'NUMPY'//char(2)//char(0)//char(iand(hl,255))//char(iand(ishft(hl,-8),255)) &
                 //char(iand(ishft(hl,-16),255))//char(iand(ishft(hl,-24),255))//d
      endif
   end function npy_header_string

   !--------------------------------------------------------------------
   ! C-order flat index of the element stored at F-order flat index k
   ! (both 0-based) for an array of the given shape.
   !--------------------------------------------------------------------
   pure integer(int64) function f2c_index(k, ndim, shp) result(c)
      integer(int64), intent(in) :: k
      integer, intent(in) :: ndim
      integer(int64), intent(in) :: shp(:)
      integer(int64) :: rem, idx(maxdim)
      integer :: d
      rem = k
      do d = 1, ndim
         idx(d) = mod(rem, shp(d)); rem = rem / shp(d)
      enddo
      c = 0
      do d = 1, ndim
         c = c * shp(d) + idx(d)
      enddo
   end function f2c_index

   !--------------------------------------------------------------------
   ! Read any supported real/int dtype into real(8) buf(n), C-flat order.
   ! If expect_shape is given the file shape must match it exactly.
   !--------------------------------------------------------------------
   subroutine npy_read_r8(filename, buf, n, ios, expect_shape, hdr, errmsg)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: n
      real(real64), intent(out) :: buf(n)
      integer, intent(out) :: ios
      integer(int64), intent(in), optional :: expect_shape(:)
      type(npy_header), intent(out), optional :: hdr
      character(len=*), intent(out), optional :: errmsg
      type(npy_header) :: h
      character(len=256) :: msg
      real(real32), allocatable :: f4(:)
      real(real64), allocatable :: f8(:)
      integer(int8), allocatable :: i1(:)
      integer(int32), allocatable :: i4(:)
      integer(int64), allocatable :: i8(:)
      integer :: u
      integer(int64) :: k
      call npy_read_header(filename, h, ios, msg)
      if (ios /= 0) goto 900
      if (present(hdr)) hdr = h
      if (h%nelem /= int(n, int64)) then
         ios = 5; write(msg,'(A,I0,A,I0)') 'element count ', h%nelem, ' /= buffer ', n; goto 900
      endif
      if (present(expect_shape)) then
         if (size(expect_shape) /= h%ndim) then
            ios = 5; msg = 'rank mismatch'; goto 900
         endif
         if (any(expect_shape /= h%shape(1:h%ndim))) then
            ios = 5; msg = 'shape mismatch'; goto 900
         endif
      endif
      open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
           status='old', action='read', convert='little_endian', iostat=ios)
      if (ios /= 0) goto 900
      select case (trim(h%descr))
      case ('<f8')
         if (h%fortran_order .and. h%ndim > 1) then
            allocate(f8(n)); read(u, pos=h%data_offset+1, iostat=ios) f8
            do k = 0, n-1
               buf(f2c_index(k, h%ndim, h%shape) + 1) = f8(k+1)
            enddo
         else
            read(u, pos=h%data_offset+1, iostat=ios) buf
         endif
      case ('<f4')
         allocate(f4(n)); read(u, pos=h%data_offset+1, iostat=ios) f4
         if (h%fortran_order .and. h%ndim > 1) then
            do k = 0, n-1
               buf(f2c_index(k, h%ndim, h%shape) + 1) = real(f4(k+1), real64)
            enddo
         else
            buf = real(f4, real64)
         endif
      case ('|i1', '<i1', '|u1', '|b1')
         allocate(i1(n)); read(u, pos=h%data_offset+1, iostat=ios) i1
         if (trim(h%descr) == '|u1') then
            f8 = real(iand(int(i1, int32), 255), real64)
         else
            f8 = real(i1, real64)
         endif
         call place(f8)
      case ('<i4')
         allocate(i4(n)); read(u, pos=h%data_offset+1, iostat=ios) i4
         f8 = real(i4, real64); call place(f8)
      case ('<i8')
         allocate(i8(n)); read(u, pos=h%data_offset+1, iostat=ios) i8
         f8 = real(i8, real64); call place(f8)
      case default
         ios = 6; msg = 'unsupported dtype '//trim(h%descr)
      end select
      close(u)
900   continue
      if (ios /= 0 .and. len_trim(msg) == 0) msg = 'read error in '//trim(filename)
      if (present(errmsg)) errmsg = msg
   contains
      subroutine place(src)
         real(real64), intent(in) :: src(:)
         integer(int64) :: kk
         if (h%fortran_order .and. h%ndim > 1) then
            do kk = 0, n-1
               buf(f2c_index(kk, h%ndim, h%shape) + 1) = src(kk+1)
            enddo
         else
            buf = src
         endif
      end subroutine place
   end subroutine npy_read_r8

   subroutine npy_read_c16(filename, buf, n, ios, hdr, errmsg)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: n
      complex(real64), intent(out) :: buf(n)
      integer, intent(out) :: ios
      type(npy_header), intent(out), optional :: hdr
      character(len=*), intent(out), optional :: errmsg
      type(npy_header) :: h
      character(len=256) :: msg
      complex(real64), allocatable :: t(:)
      integer :: u
      integer(int64) :: k
      call npy_read_header(filename, h, ios, msg)
      if (ios == 0 .and. trim(h%descr) /= '<c16') then
         ios = 6; msg = 'expected <c16, got '//trim(h%descr)
      endif
      if (ios == 0 .and. h%nelem /= int(n, int64)) then
         ios = 5; msg = 'element count mismatch'
      endif
      if (ios == 0) then
         if (present(hdr)) hdr = h
         open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
              status='old', action='read', convert='little_endian', iostat=ios)
         if (ios == 0) then
            if (h%fortran_order .and. h%ndim > 1) then
               allocate(t(n)); read(u, pos=h%data_offset+1, iostat=ios) t
               do k = 0, n-1
                  buf(f2c_index(k, h%ndim, h%shape) + 1) = t(k+1)
               enddo
            else
               read(u, pos=h%data_offset+1, iostat=ios) buf
            endif
            close(u)
         endif
      endif
      if (present(errmsg)) errmsg = msg
   end subroutine npy_read_c16

   !--------------------------------------------------------------------
   ! writers: payload = buf in storage order, header as given
   !--------------------------------------------------------------------
   subroutine npy_write_r8(filename, buf, n, shp, ios, fortran_order, dtype)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: n
      real(real64), intent(in) :: buf(n)
      integer(int64), intent(in) :: shp(:)
      integer, intent(out) :: ios
      logical, intent(in), optional :: fortran_order
      character(len=*), intent(in), optional :: dtype     ! 'f8' (default) | 'f4'
      logical :: fo
      character(len=3) :: dt
      integer :: u
      fo = .false.; if (present(fortran_order)) fo = fortran_order
      dt = 'f8'; if (present(dtype)) dt = dtype
      if (product(shp) /= int(n, int64)) then
         ios = 5; return
      endif
      open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
           status='replace', action='write', convert='little_endian', iostat=ios)
      if (ios /= 0) return
      if (trim(dt) == 'f4') then
         write(u, iostat=ios) npy_header_string('<f4', fo, shp)
         write(u, iostat=ios) real(buf, real32)
      else
         write(u, iostat=ios) npy_header_string('<f8', fo, shp)
         write(u, iostat=ios) buf
      endif
      close(u)
   end subroutine npy_write_r8

   subroutine npy_write_c16(filename, buf, n, shp, ios, fortran_order)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: n
      complex(real64), intent(in) :: buf(n)
      integer(int64), intent(in) :: shp(:)
      integer, intent(out) :: ios
      logical, intent(in), optional :: fortran_order
      logical :: fo
      integer :: u
      fo = .false.; if (present(fortran_order)) fo = fortran_order
      open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
           status='replace', action='write', convert='little_endian', iostat=ios)
      if (ios /= 0) return
      write(u, iostat=ios) npy_header_string('<c16', fo, shp)
      write(u, iostat=ios) buf
      close(u)
   end subroutine npy_write_c16

   subroutine npy_write_i1(filename, buf, n, shp, ios, fortran_order)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: n
      integer(int8), intent(in) :: buf(n)
      integer(int64), intent(in) :: shp(:)
      integer, intent(out) :: ios
      logical, intent(in), optional :: fortran_order
      logical :: fo
      integer :: u
      fo = .false.; if (present(fortran_order)) fo = fortran_order
      open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
           status='replace', action='write', convert='little_endian', iostat=ios)
      if (ios /= 0) return
      write(u, iostat=ios) npy_header_string('|i1', fo, shp)
      write(u, iostat=ios) buf
      close(u)
   end subroutine npy_write_i1

   ! .true. when every value survives real(8) -> real(4) -> real(8) exactly
   ! (NaN never does; -0.0 does).
   logical function npy_f4_exact(buf, n)
      integer, intent(in) :: n
      real(real64), intent(in) :: buf(n)
      integer :: i
      npy_f4_exact = .true.
      do i = 1, n
         if (real(real(buf(i), real32), real64) /= buf(i)) then
            npy_f4_exact = .false.; return
         endif
      enddo
   end function npy_f4_exact

end module npy_io_mod
