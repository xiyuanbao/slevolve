! run_outer_mem.f90 -- driver running the legacy outer (initial-topography)
! iteration INSIDE ONE PROCESS:
!   ./exe LAST_ITER K0 K1
! runs outer passes itersl = K0..K1 (Picard: tinit_0 <- tinit_0 + (T_obs - T_pred)).
! Pass K0 reads OUTPUT_SLM/pred_pres_topo_<K0-1>, tgrid0_<K0-1> if K0 > 1 (the
! seed, as in a run of run_outer_plus); later passes take them from memory.  The
! Earth model, Legendre table, FFTW plans, hereditary-memory allocation, ice
! grids (SLEVOLVE_ICECACHE, default 1) and the SLEVOLVE_REUSE=warm store stay
! resident.  pred_pres_topo_<k>.bin / tgrid0_<k>.bin are still written
! (diagnostics).  Per pass: OUTPUT_SLM/numiter is copied to numiter_<k> and
! one line '== pass k wall ... s' is printed.
program run_outer_mem
  use sl_model_mod
  use user_specs_mod, only: dt1
  use spharmt, only: sht_profile
  use io_mod, only: pp_shtprof, pp_read_options
  implicit none
  integer :: iter,last,k,k0,k1,ios,st
  character(len=32)::arg, v
  integer(8) :: c0, c1, cr, cs
  real(8) :: ti0, tq0, tt0, tc0
  if (command_argument_count() /= 3) call usage()
  call get_command_argument(1,arg); read(arg,*,iostat=ios) last; if(ios/=0 .or. last<1) call usage()
  call get_command_argument(2,arg); read(arg,*,iostat=ios) k0; if(ios/=0 .or. k0<1) call usage()
  call get_command_argument(3,arg); read(arg,*,iostat=ios) k1; if(ios/=0 .or. k1<k0) call usage()
  call pp_read_options()
  if (pp_shtprof) sht_profile = .true.
  call get_environment_variable('SLEVOLVE_ICECACHE', v, status=st)
  if (st == 0) pp_icecache = .not. (trim(v) == '0')
  pp_mem_outer = .true.
  call system_clock(cs, cr)
  call sl_solver_checkpoint(k0,dt1)          ! planet, SHT env, Legendre/FFTW: once
  call system_clock(c0)
  write(6,'(A,F9.3,A)') '== one-time setup (checkpoint: planet, SHT tables, FFTW plans): ', real(c0-cs,8)/real(cr,8), ' s'
  do k = k0, k1
    call system_clock(c0)
    ti0 = pp_t_init; tq0 = pp_t_quant; tt0 = pp_t_topo; tc0 = pp_t_ice
    if (k > k0) call sl_outer_reset()
    call sl_timewindow(0)
    call sl_allocate_and_initialize_array
    call sl_solver_init(k,0.0)
    call sl_deallocate_array
    do iter=1,last
      call sl_timewindow(iter)
      call sl_allocate_and_initialize_array
      call sl_solver(k,iter,dt1,0.0)
      call sl_deallocate_array
    enddo
    call system_clock(c1)
    write(arg,'(I0)') k
    call execute_command_line('cp OUTPUT_SLM/numiter numiter_'//trim(arg))
    write(6,'(A,I3,A,F9.3,A,4F8.3)') '== pass', k, ' wall', real(c1-c0,8)/real(cr,8), &
         ' s ; init, quantize, truetopo, ice-read s:', pp_t_init-ti0, pp_t_quant-tq0, pp_t_topo-tt0, pp_t_ice-tc0
    flush(6)
  enddo
  call system_clock(c1)
  write(6,'(A,F9.3,A)') '== process total ', real(c1-cs,8)/real(cr,8), ' s'
  call sl_print_counts(6)
contains
  subroutine usage()
    write(0,'(A)') 'usage: ./exe LAST K0 K1   (driver run_outer_mem: outer passes K0..K1 of steps 1..LAST in one process)'
    call exit(2)
  end subroutine usage
end program run_outer_mem
