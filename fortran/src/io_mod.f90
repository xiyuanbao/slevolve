!=======================================================================
! I/O layer.
!
! The legacy solver spends ~40 % of its wall time in formatted text
! I/O: 13 MB ice files parsed with list-directed reads, five 8.9 MB
! ES16.9E2 grids written per step, a 7.9 MB complex dS field, and a 53 MB
! hereditary-state checkpoint written and read back every step.
!
! This module provides drop-in replacements that
!   * read grids from a flat little-endian fp64 stream cache when one
!     exists (same element order as the text file: Fortran array order),
!     falling back to the original text reader otherwise;
!   * write grids as flat fp64 streams, and optionally also as the
!     original text files (SLEVOLVE_TEXT=1).
!
! Runtime switches come from the environment so the same binary can be
! benchmarked in every mode:
!   SLEVOLVE_TEXT   = 0|1   also emit the legacy text outputs   (default 0)
!   SLEVOLVE_BIN    = 0|1   emit binary outputs                 (default 1)
!   SLEVOLVE_BININ  = 0|1   prefer binary inputs when present    (default 1)
!   SLEVOLVE_CKPT   = K     write the exact-memory checkpoint every K
!                           steps (0 = never; the state stays resident)
!
! Snapshot mode (added 2026-09-23 for archive comparisons):
!   SLEVOLVE_SNAP   = all | none | list   restrict the per-step GRID outputs
!                     (tgrid, ocean, beta, dS_converged, R, G) to the listed
!                     steps.  list = comma-separated steps and ranges, e.g.
!                     "1,2,5,10,50,100,150,200,300,400,550" or "1-20,50".
!                     Step 0 is always added to a list.  "none" writes no
!                     grid at all.  Unset = legacy behaviour (every field,
!                     every step).  The final-step pred_pres_topo_/tgrid0_
!                     copies are not written in snapshot mode.
!   SLEVOLVE_SNAP_FIELDS = tgrid[,ocean,beta,dS,R,G,horiz,Psi] | all   (default tgrid)
! Snapshot files use the normal names (tgrid<N>.bin, ...): flat
! little-endian stream in Fortran array order, identical to the .bin ice
! caches; OUTPUT_SLM/snap_index.txt lists step, field, file, dtype, shape.
! SLEVOLVE_TEXT=1 still writes the text form, restricted to the same
! steps/fields when SLEVOLVE_SNAP is set.  Scalar per-step text files
! (numiter, elapsed_wall_time, TPW<N>, areas, gmsl, ice_volume) are
! unaffected.
!
! NumPy .npy is the internal binary format.
!   SLEVOLVE_FMT = npy (default) | bin   binary OUTPUT container.  npy writes
!                  <name>.npy: grids '<f8' shape (n1,n2) fortran_order=False,
!                  complex dS '<c16' shape (n1,n2) fortran_order=True; the
!                  payload is identical to the legacy headerless .bin.
!   Inputs: <path>.npy (any of f8/f4/i1..., f4 promoted exactly), then
!   <path>.bin (legacy LE f8 stream), then the text file <path>.
!   See fortran/IO_FORMATS.md and npy_io_mod.f90.
!=======================================================================
module io_mod
   use npy_io_mod, only: read_grid, write_grid, write_snapshot
   implicit none
   public

   logical, save :: pp_text_out  = .false.
   logical, save :: pp_bin_out   = .true.
   logical, save :: pp_bin_in    = .true.
   logical, save :: pp_npy_out   = .true.     ! SLEVOLVE_FMT=npy (default) | bin
   integer, save :: pp_ckpt_every = 0
   ! SLEVOLVE_COMPAT bit flags (verification only: re-enable legacy per-step work to
   ! check that keeping state resident does not change results; default 0):
   !   1  re-read tinit_0 / cxy0 / beta0 / cxy from disk every step
   !   2  re-parse the Earth model every step
   !   4  save + reload the hereditary-state checkpoint every step
   !   8  recompute t0lm every step
   integer, save :: pp_compat = 0
   logical, save :: pp_opts_read = .false.
   ! transform batching and per-step re-zeroing.
   !   SLEVOLVE_BATCH = 0|1  batch mutually independent SH transforms (default 1)
   !   SLEVOLVE_ZERO  = full|min  how much of the working state to re-zero at
   !                     the top of every timestep (default min)
   !   SLEVOLVE_SHTPROF = 0|1 print the internal transform profile at the end
   logical, save :: pp_batch_opt = .true.
   logical, save :: pp_zero_full = .false.
   logical, save :: pp_shtprof   = .false.

   ! transform-count reductions.
   !   SLEVOLVE_OPT = all | none | comma list of
   !     diag      the four post-loop GMSL/area analyses -> degree-0 quadrature
   !     dsl       deltaslxy = dslxy + conserv*K00 instead of a synthesis
   !     lazyrg    skip the R synthesis when neither R nor G is written this step
   !     icebatch  grounded-ice analysis rides in the first-guess sweep
   !     rt        one analysis of (rO - tO) instead of rO and tO separately
   !     merge     rO analysis and the next ocean-function re-guess in one sweep
   !               (conserv from the degree-0 quadrature of rO; implies dsl)
   !   default all.  SLEVOLVE_OPT=none disables all six reductions.
   logical, save :: opt_diag = .true., opt_dsl = .true., opt_lazyrg = .true.
   logical, save :: opt_icebatch = .true., opt_rt = .true., opt_merge = .true.

   ! optional iteration accelerations (default OFF).
   !   SLEVOLVE_INNER = off (default) | extrap | anderson
   !     extrap    initial guess of the step-n ocean-load increment
   !               dS_n = 2 dS_{n-1} - dS_{n-2} (converged increments of the two
   !               previous steps of this process); first two steps unchanged
   !     anderson  extrap + Anderson mixing of the spectral dS iterate from the
   !               3rd inner iteration on, ocean function frozen from iteration 2
   !   SLEVOLVE_AA_M  = Anderson depth m (default 3, 1..5)
   integer, save :: inner_mode = 0
   integer, save :: inner_aa_m  = 3
   logical, save :: inner_refresh = .false.   ! SLEVOLVE_AA_MASK=refresh: do not exit on a stale frozen mask
   logical, save :: inner_log   = .false.
   ! SLEVOLVE_OUTER = picard | secant   (unset = no outer-iteration files in snapshot mode)
   !   When set: the outer-iteration files pred_pres_topo_<k>, tgrid0_<k> are
   !   written (binary, plus text with SLEVOLVE_TEXT=1) even in snapshot mode,
   !   and read back with pp_read_grid (binary preferred).  secant: from outer
   !   iteration 3 on, per-cell secant update of the initial topography.
   logical, save :: outer_set = .false.
   ! SLEVOLVE_REUSE = off (default) | warm   cross-outer-iteration warm start:
   !   pass ITERSL writes OUTPUT_SLM/reuse_<ITERSL>.bin (converged dS_n and ocean
   !   mask of every step n) and, if OUTPUT_SLM/reuse_<ITERSL-1>.bin exists,
   !   starts step n's inner loop from the stored dS_n and mask.
   integer, save :: warm_mode = 0
   integer, save :: outer_mode = 0          ! 0 picard, 1 secant     ! SLEVOLVE_INNER_LOG=1: per-iteration xi, |F|, mask changes

   ! SLEVOLVE_HORIZ = 0 (default) | 1 | 2   optional horizontal (tangential)
   !   displacement read-out from the same modal states as R.
   !     1  write the spectral scalar Psi_lm (u_h = grad_h Psi) as horiz<N>.npy
   !        ('<c16', (L+1,L+1), fortran_order=True, like dS_converged)
   !     2  also synthesize Psi on the grid: Psi<N>.npy
   !   The gradient (u_theta, u_phi) is formed in post-processing
   !   (python/tools/horiz_grad.py).  0 leaves every other output unchanged.
   integer, save :: pp_horiz = 0

   ! snapshot mode
   integer, parameter :: pp_snap_maxstep = 100000
   logical, save :: pp_snap_on   = .false.     ! SLEVOLVE_SNAP set
   logical, save :: pp_snap_all  = .false.
   logical, save :: pp_snap_step(0:pp_snap_maxstep) = .false.
   integer, save :: pp_snap_fields = 1         ! bit mask, see pp_field_bit
   logical, save :: pp_index_started = .false.

   ! wall-clock accounting for the I/O layer (seconds)
   real(8), save :: pp_io_read = 0.0d0, pp_io_write = 0.0d0

contains

   subroutine pp_read_options()
      character(len=64) :: v
      integer :: st, k
      if (pp_opts_read) return
      pp_opts_read = .true.
      call get_environment_variable('SLEVOLVE_TEXT', v, status=st)
      if (st == 0) pp_text_out = (trim(v) == '1' .or. trim(v) == 'T' .or. trim(v) == 'true')
      call get_environment_variable('SLEVOLVE_BIN', v, status=st)
      if (st == 0) pp_bin_out = .not. (trim(v) == '0' .or. trim(v) == 'F' .or. trim(v) == 'false')
      call get_environment_variable('SLEVOLVE_FMT', v, status=st)
      if (st == 0) pp_npy_out = (trim(adjustl(v)) /= 'bin')
      call get_environment_variable('SLEVOLVE_BININ', v, status=st)
      if (st == 0) pp_bin_in = .not. (trim(v) == '0' .or. trim(v) == 'F' .or. trim(v) == 'false')
      call get_environment_variable('SLEVOLVE_COMPAT', v, status=st)
      if (st == 0) then
         read(v,*,iostat=k) pp_compat
         if (k /= 0) pp_compat = 0
      endif
      call get_environment_variable('SLEVOLVE_CKPT', v, status=st)
      if (st == 0) then
         read(v,*,iostat=k) pp_ckpt_every
         if (k /= 0) pp_ckpt_every = 0
      endif
      call get_environment_variable('SLEVOLVE_BATCH', v, status=st)
      if (st == 0) pp_batch_opt = .not. (trim(v) == '0' .or. trim(v) == 'F' .or. trim(v) == 'false')
      call get_environment_variable('SLEVOLVE_ZERO', v, status=st)
      if (st == 0) then
         if (trim(v) == 'full') then
            pp_zero_full = .true.
         else if (trim(v) == 'min') then
            pp_zero_full = .false.
         else
            print *,'SLEVOLVE_ZERO must be full or min'; stop
         endif
      endif
      call get_environment_variable('SLEVOLVE_SHTPROF', v, status=st)
      if (st == 0) pp_shtprof = (trim(v) == '1')
      call get_environment_variable('SLEVOLVE_HORIZ', v, status=st)
      if (st == 0) then
         read(v,*,iostat=k) pp_horiz
         if (k /= 0 .or. pp_horiz < 0 .or. pp_horiz > 2) then
            print *,'SLEVOLVE_HORIZ must be 0, 1 or 2'; stop
         endif
      endif
      call pp_read_snap_options()
      call opt_read_options()
      call inner_read_options()
   end subroutine pp_read_options

   subroutine inner_read_options()
      character(len=64) :: v
      integer :: st, k
      call get_environment_variable('SLEVOLVE_INNER', v, status=st)
      if (st == 0) then
         select case (trim(adjustl(v)))
         case ('', 'off', 'none', 'picard'); inner_mode = 0
         case ('extrap');   inner_mode = 1
         case ('anderson'); inner_mode = 2
         case default
            print *,'SLEVOLVE_INNER must be off, extrap or anderson'; stop
         end select
      endif
      call get_environment_variable('SLEVOLVE_OUTER', v, status=st)
      if (st == 0) then
         outer_set = .true.
         select case (trim(adjustl(v)))
         case ('picard'); outer_mode = 0
         case ('secant'); outer_mode = 1
         case default
            print *,'SLEVOLVE_OUTER must be picard or secant'; stop
         end select
      endif
      call get_environment_variable('SLEVOLVE_REUSE', v, status=st)
      if (st == 0) then
         select case (trim(adjustl(v)))
         case ('', 'off', 'none'); warm_mode = 0
         case ('warm'); warm_mode = 1
         case default
            print *,'SLEVOLVE_REUSE must be off or warm'; stop
         end select
      endif
      call get_environment_variable('SLEVOLVE_AA_MASK', v, status=st)
      if (st == 0) then
         select case (trim(adjustl(v)))
         case ('frozen');  inner_refresh = .false.
         case ('refresh'); inner_refresh = .true.
         case default
            print *,'SLEVOLVE_AA_MASK must be frozen or refresh'; stop
         end select
      endif
      call get_environment_variable('SLEVOLVE_INNER_LOG', v, status=st)
      if (st == 0) inner_log = (trim(v) == '1')
      call get_environment_variable('SLEVOLVE_AA_M', v, status=st)
      if (st == 0) then
         read(v,*,iostat=k) inner_aa_m
         if (k /= 0 .or. inner_aa_m < 1 .or. inner_aa_m > 5) then
            print *,'SLEVOLVE_AA_M must be 1..5'; stop
         endif
      endif
      if (inner_mode == 2 .and. .not. (opt_merge .and. opt_rt)) then
         print *,'SLEVOLVE_INNER=anderson requires SLEVOLVE_OPT with merge and rt'; stop
      endif
   end subroutine inner_read_options

   subroutine opt_read_options()
      character(len=256) :: v, tok
      integer :: st, n, p, q
      call get_environment_variable('SLEVOLVE_OPT', v, status=st)
      if (st /= 0) return                        ! unset: all on
      v = adjustl(v)
      if (trim(v) == 'all') return
      opt_diag = .false.; opt_dsl = .false.; opt_lazyrg = .false.
      opt_icebatch = .false.; opt_rt = .false.; opt_merge = .false.
      if (trim(v) == 'none' .or. len_trim(v) == 0) return
      p = 1
      n = len_trim(v)
      do while (p <= n)
         q = index(v(p:n), ',')
         if (q == 0) then
            tok = v(p:n); p = n + 1
         else
            tok = v(p:p+q-2); p = p + q
         endif
         tok = adjustl(tok)
         select case (trim(tok))
         case ('diag');     opt_diag = .true.
         case ('dsl');      opt_dsl = .true.
         case ('lazyrg');   opt_lazyrg = .true.
         case ('icebatch'); opt_icebatch = .true.
         case ('rt');       opt_rt = .true.
         case ('merge');    opt_merge = .true.
         case ('');         continue
         case default
            print *,'SLEVOLVE_OPT: unknown option "'//trim(tok)//'"'; stop
         end select
      enddo
      if (opt_merge) opt_dsl = .true.            ! merge needs deltasl on the grid
   end subroutine opt_read_options

   subroutine pp_read_snap_options()
      character(len=4096) :: v, tok
      integer :: st, n, p, q, a, b, i, k1, k2
      call get_environment_variable('SLEVOLVE_SNAP', v, length=n, status=st)
      if (st == 1) return                       ! unset: legacy output
      if (st /= 0) then
         print *,'SLEVOLVE_SNAP too long (max 4096 chars)'; stop
      endif
      pp_snap_on = .true.
      v = adjustl(v)
      if (trim(v) == 'all') then
         pp_snap_all = .true.
      else if (trim(v) == 'none' .or. len_trim(v) == 0) then
         continue
      else
         pp_snap_step(0) = .true.
         p = 1
         n = len_trim(v)
         do while (p <= n)
            q = index(v(p:n), ',')
            if (q == 0) then
               tok = v(p:n); p = n + 1
            else
               tok = v(p:p+q-2); p = p + q
            endif
            tok = adjustl(tok)
            if (len_trim(tok) == 0) cycle
            i = index(tok, '-')
            if (i > 1) then
               read(tok(1:i-1),*,iostat=k1) a
               read(tok(i+1:),*,iostat=k2) b
            else
               read(tok,*,iostat=k1) a
               b = a; k2 = 0
            endif
            if (k1 /= 0 .or. k2 /= 0 .or. a < 0 .or. b < a .or. b > pp_snap_maxstep) then
               print *,'SLEVOLVE_SNAP: bad entry "'//trim(tok)//'"'; stop
            endif
            pp_snap_step(a:b) = .true.
         enddo
      endif
      call get_environment_variable('SLEVOLVE_SNAP_FIELDS', v, status=st)
      if (st == 0) then
         v = adjustl(v)
         if (trim(v) == 'all') then
            pp_snap_fields = 255
         else
            pp_snap_fields = 0
            p = 1
            n = len_trim(v)
            do while (p <= n)
               q = index(v(p:n), ',')
               if (q == 0) then
                  tok = v(p:n); p = n + 1
               else
                  tok = v(p:p+q-2); p = p + q
               endif
               tok = adjustl(tok)
               if (len_trim(tok) == 0) cycle
               i = pp_field_bit(trim(tok))
               if (trim(tok) == 'dS') i = 8
               if (i == 0) then
                  print *,'SLEVOLVE_SNAP_FIELDS: unknown field "'//trim(tok)//'"'; stop
               endif
               pp_snap_fields = ior(pp_snap_fields, i)
            enddo
         endif
      endif
   end subroutine pp_read_snap_options

   pure integer function pp_field_bit(filename)
      character(len=*), intent(in) :: filename
      select case (trim(filename))
      case ('tgrid');        pp_field_bit = 1
      case ('ocean');        pp_field_bit = 2
      case ('beta');         pp_field_bit = 4
      case ('dS_converged'); pp_field_bit = 8
      case ('R');            pp_field_bit = 16
      case ('G');            pp_field_bit = 32
      case ('horiz');        pp_field_bit = 64
      case ('Psi');          pp_field_bit = 128
      case default;          pp_field_bit = 0
      end select
   end function pp_field_bit

   !--------------------------------------------------------------------
   ! In snapshot mode: is <filename><suffix> a requested snapshot?  The
   ! step number is the numeric suffix every per-step grid carries.
   ! Always .true. outside snapshot mode.
   !--------------------------------------------------------------------
   logical function pp_snap_want(filename, suffix, step)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: suffix
      integer, intent(out) :: step
      integer :: bit, ios
      step = -1
      pp_snap_want = .true.
      if (.not. pp_snap_on) return
      pp_snap_want = .false.
      bit = pp_field_bit(filename)
      if (bit == 0 .or. iand(bit, pp_snap_fields) == 0) return
      if (.not. present(suffix)) return
      read(suffix,*,iostat=ios) step
      if (ios /= 0 .or. step < 0 .or. step > pp_snap_maxstep) return
      pp_snap_want = pp_snap_all .or. pp_snap_step(step)
   end function pp_snap_want

   subroutine pp_snap_index(filepath, filename, suffix, step, dtype, n1, n2)
      character(len=*), intent(in) :: filepath, filename, suffix, dtype
      integer, intent(in) :: step, n1, n2
      integer :: u
      if (.not. pp_snap_on) return
      if (.not. pp_index_started) then
         open(newunit=u, file=trim(filepath)//'snap_index.txt', form='formatted', &
              status='replace', action='write')
         if (pp_npy_out) then
            write(u,'(A)') '# SLEvolve snapshot index.  Files are NumPy .npy (numpy.load).  Grids: <f8,'
            write(u,'(A)') '# shape (n1,n2), fortran_order=False, a[row,col] in the text-file stream order'
            write(u,'(A)') '# (row = Gauss-Legendre latitude N->S, col = longitude).  dS: <c16, shape'
            write(u,'(A)') '# (n1,n2), fortran_order=True, a[l,m].  See docs/src/formats.md.'
         else
            write(u,'(A)') '# SLEvolve snapshot index.  Files are raw little-endian streams in Fortran'
            write(u,'(A)') '# array order (first index fastest), no header, same layout as the .bin ice'
            write(u,'(A)') '# caches: numpy.fromfile(f, "<f8") (or "<c16") gives the flat array whose order'
            write(u,'(A)') '# equals the text output read with numpy.fromfile(sep=" ").'
         endif
         write(u,'(A)') '# step field file dtype n1 n2'
         pp_index_started = .true.
      else
         open(newunit=u, file=trim(filepath)//'snap_index.txt', form='formatted', &
              status='old', position='append', action='write')
      endif
      write(u,'(I6,1X,A,1X,A,1X,A,2(1X,I6))') step, trim(filename), &
           trim(filename)//trim(suffix)//merge('.npy', '.bin', pp_npy_out), trim(dtype), n1, n2
      close(u)
   end subroutine pp_snap_index

   subroutine pp_report(unit_num)
      integer, intent(in) :: unit_num
      call pp_read_options()
      write(unit_num,'(A,L2,A,L2,A,L2,A,I5,A,A)') ' SLEvolve io: text_out=', pp_text_out, &
           '  bin_out=', pp_bin_out, '  bin_in=', pp_bin_in, '  ckpt_every=', pp_ckpt_every, &
           '  fmt=', merge('npy', 'bin', pp_npy_out)
      write(unit_num,'(A,I3)') ' SLEvolve compat flags: ', pp_compat
      write(unit_num,'(A,6(A,L2))') ' SLEvolve opts:', ' diag=', opt_diag, ' dsl=', opt_dsl, &
           ' lazyrg=', opt_lazyrg, ' icebatch=', opt_icebatch, ' rt=', opt_rt, ' merge=', opt_merge
      if (inner_mode /= 0) write(unit_num,'(A,I2,A,I2)') ' SLEvolve inner acceleration: mode=', inner_mode, &
           '  (1=extrap 2=anderson)  anderson depth=', inner_aa_m
      if (pp_snap_on) then
         write(unit_num,'(A,L2,A,I6,A,I3)') ' SLEvolve snapshot mode: all=', pp_snap_all, &
              '  nsteps listed=', count(pp_snap_step), '  field mask=', pp_snap_fields
      endif
   end subroutine pp_report

   pure function pp_path(filename, filepath, suffix) result(p)
      character(len=*), intent(in) :: filename, filepath
      character(len=*), intent(in), optional :: suffix
      character(len=1024) :: p
      if (present(suffix)) then
         p = trim(filepath)//trim(filename)//trim(suffix)
      else
         p = trim(filepath)//trim(filename)
      endif
   end function pp_path

   !--------------------------------------------------------------------
   ! Read a real grid.  Prefers <path>.bin (flat fp64 stream) over <path>.
   !--------------------------------------------------------------------
   subroutine pp_read_grid(data_slm, filename, filepath, suffix)
      real, dimension(:,:), intent(inout) :: data_slm
      character(len=*), intent(in) :: filename, filepath
      character(len=*), intent(in), optional :: suffix
      character(len=1024) :: path
      integer :: u, ios
      character(len=256) :: msg
      real(8) :: t0, t1
      call pp_read_options()
      call pp_clock(t0)
      path = pp_path(filename, filepath, suffix)
      if (pp_bin_in) then
         call read_grid(data_slm, path, ios, errmsg=msg)     ! .npy, then .bin, then text
         if (ios /= 0) then
            print *, 'pp_read_grid: ', trim(msg); stop 1
         endif
      else
         open(newunit=u, file=trim(path), form='formatted', access='sequential', &
              status='old', action='read')
         read(u,*) data_slm
         close(u)
      endif
      call pp_clock(t1); pp_io_read = pp_io_read + (t1-t0)
   end subroutine pp_read_grid

   !--------------------------------------------------------------------
   ! Write a real grid: binary stream always (unless disabled), plus the
   ! original ES16.9E2 text file when SLEVOLVE_TEXT=1.
   !--------------------------------------------------------------------
   subroutine pp_write_grid(data_slm, filename, filepath, suffix, force)
      real, dimension(:,:), intent(in) :: data_slm
      character(len=*), intent(in) :: filename, filepath
      character(len=*), intent(in), optional :: suffix
      logical, intent(in), optional :: force
      character(len=1024) :: path
      integer :: u, step
      real(8) :: t0, t1
      logical :: frc
      call pp_read_options()
      frc = .false.
      if (present(force)) frc = force
      if (.not. frc .and. .not. pp_snap_want(filename, suffix, step)) return
      if (frc) step = -1
      call pp_clock(t0)
      path = pp_path(filename, filepath, suffix)
      if (pp_bin_out .and. pp_npy_out) then
         call write_snapshot(data_slm, path, u)
         if (u /= 0) then
            print *, 'pp_write_grid: cannot write ', trim(path)//'.npy'; stop 1
         endif
         if (pp_snap_on .and. .not. frc) call pp_snap_index(filepath, filename, suffix, step, 'f8', size(data_slm,1), size(data_slm,2))
      else if (pp_bin_out) then
         open(newunit=u, file=trim(path)//'.bin', form='unformatted', access='stream', &
              status='replace', action='write')
         write(u) data_slm
         close(u)
         if (pp_snap_on .and. .not. frc) call pp_snap_index(filepath, filename, suffix, step, 'f8', size(data_slm,1), size(data_slm,2))
      endif
      if (pp_text_out) then
         open(newunit=u, file=trim(path), form='formatted', access='sequential', &
              status='replace', action='write')
         write(u,'(ES16.9E2)') data_slm
         close(u)
      endif
      call pp_clock(t1); pp_io_write = pp_io_write + (t1-t0)
   end subroutine pp_write_grid

   !--------------------------------------------------------------------
   ! Complex spectral field (dS_converged).  Text form is ES14.7E2 with
   ! real and imaginary parts as consecutive records, exactly as the
   ! baseline writes it.
   !--------------------------------------------------------------------
   subroutine pp_write_cplx(z, filename, filepath, suffix, fmt)
      complex, dimension(:,:), intent(in) :: z
      character(len=*), intent(in) :: filename, filepath, fmt
      character(len=*), intent(in), optional :: suffix
      character(len=1024) :: path
      integer :: u, step
      real(8) :: t0, t1
      call pp_read_options()
      if (.not. pp_snap_want(filename, suffix, step)) return
      call pp_clock(t0)
      path = pp_path(filename, filepath, suffix)
      if (pp_bin_out .and. pp_npy_out) then
         call write_snapshot(z, path, u)
         if (u /= 0) then
            print *, 'pp_write_cplx: cannot write ', trim(path)//'.npy'; stop 1
         endif
         if (pp_snap_on) call pp_snap_index(filepath, filename, suffix, step, 'c16', size(z,1), size(z,2))
      else if (pp_bin_out) then
         open(newunit=u, file=trim(path)//'.bin', form='unformatted', access='stream', &
              status='replace', action='write')
         write(u) z
         close(u)
         if (pp_snap_on) call pp_snap_index(filepath, filename, suffix, step, 'c16', size(z,1), size(z,2))
      endif
      if (pp_text_out) then
         open(newunit=u, file=trim(path), form='formatted', access='sequential', &
              status='replace', action='write')
         write(u,fmt) z
         close(u)
      endif
      call pp_clock(t1); pp_io_write = pp_io_write + (t1-t0)
   end subroutine pp_write_cplx

   subroutine pp_clock(t)
      real(8), intent(out) :: t
      integer(8) :: c, r
      call system_clock(c, r)
      t = real(c,8)/real(r,8)
   end subroutine pp_clock

end module io_mod
