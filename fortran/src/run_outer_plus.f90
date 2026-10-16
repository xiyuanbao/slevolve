! run_outer_plus.f90 -- driver for one OUTER iteration (the legacy topography
! iteration of Han et al. 2022), one pass per process:
!   ./exe LAST_ITER ITERSL
! ITERSL = 1 starts from the present topography; ITERSL = k > 1 reads
! OUTPUT_SLM/pred_pres_topo_<k-1> and tgrid0_<k-1> (and, with SLEVOLVE_OUTER=secant
! and k >= 3, also <k-2>).  dtime = user_specs_mod::dt1.
program run_outer_plus
  use sl_model_mod
  use user_specs_mod, only: dt1
  use spharmt, only: sht_profile
  use io_mod, only: pp_shtprof, pp_read_options
  implicit none
  integer :: iter,last,itersl,ios,nargs
  character(len=32)::arg
  nargs = command_argument_count()
  if(nargs<1 .or. nargs>2) call usage()
  call get_command_argument(1,arg); read(arg,*,iostat=ios) last; if(ios/=0 .or. last<1) call usage()
  itersl = 1
  if(nargs==2) then
     call get_command_argument(2,arg); read(arg,*,iostat=ios) itersl; if(ios/=0 .or. itersl<1) call usage()
  endif
  call pp_read_options()
  if (pp_shtprof) sht_profile = .true.
  call sl_solver_checkpoint(itersl,dt1)
  call sl_timewindow(0)
  call sl_allocate_and_initialize_array
  call sl_solver_init(itersl,0.0)
  call sl_deallocate_array
  do iter=1,last
    call sl_timewindow(iter)
    call sl_allocate_and_initialize_array
    call sl_solver(itersl,iter,dt1,0.0)
    call sl_deallocate_array
  enddo
  call sl_print_counts(6)
contains
  subroutine usage()
    write(0,'(A)') 'usage: ./exe LAST [ITERSL]   (driver run_outer_plus: one outer pass ITERSL (default 1) of steps 1..LAST)'
    call exit(2)
  end subroutine usage
end program run_outer_plus
