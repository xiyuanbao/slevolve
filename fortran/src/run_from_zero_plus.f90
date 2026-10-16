program run_from_zero
  use sl_model_mod
  use spharmt, only: sht_profile, sht_times, sht_counts
  use io_mod, only: pp_shtprof, pp_read_options
  use user_specs_mod, only: dt1
  implicit none
  integer :: iter,last,ios
  character(len=32)::arg
  if(command_argument_count()/=1) call usage()
  call get_command_argument(1,arg); read(arg,*,iostat=ios) last; if(ios/=0 .or. last<1) call usage()
  call pp_read_options()
  if (pp_shtprof) sht_profile = .true.
  call sl_solver_checkpoint(1,dt1)
  call sl_timewindow(0)
  call sl_allocate_and_initialize_array
  call sl_solver_init(1,0.0)
  call sl_deallocate_array
  do iter=1,last
    call sl_timewindow(iter)
    call sl_allocate_and_initialize_array
    call sl_solver(1,iter,dt1,0.0)
    call sl_deallocate_array
  enddo
  call sl_print_counts(6)
  if (pp_shtprof) then
    write(6,'(A)') '== SHT internal profile (s / calls) =='
    write(6,'(A,F9.3,I8)') '  grid pack        ', sht_times(1), sht_counts(1)
    write(6,'(A,F9.3,I8)') '  fft forward      ', sht_times(2), sht_counts(2)
    write(6,'(A,F9.3,I8)') '  interleave+wts   ', sht_times(3), sht_counts(3)
    write(6,'(A,F9.3,I8)') '  legendre analysis', sht_times(4), sht_counts(4)
    write(6,'(A,F9.3,I8)') '  legendre synth   ', sht_times(6), sht_counts(6)
    write(6,'(A,F9.3,I8)') '  de-interleave    ', sht_times(7), sht_counts(7)
    write(6,'(A,F9.3,I8)') '  fft inverse      ', sht_times(8), sht_counts(8)
    write(6,'(A,F9.3,I8)') '  grid unpack      ', sht_times(5), sht_counts(5)
  endif
contains
  subroutine usage()
    write(0,'(A)') 'usage: ./exe LAST   (driver run_from_zero_plus: forward run of steps 1..LAST; the step length is the build-time DT)'
    call exit(2)
  end subroutine usage
end program run_from_zero
