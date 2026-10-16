!=======================================================================================================MAIN MODULE=====!
!_______________________________________________________________________________________________________________________!
module sl_model_mod

   use spharmt
   use user_specs_mod
   use planets_mod
   use sl_io_mod
   use modal_state_mod, only: dp, response_spectrum, response_memory_state, lm_index, extrap_coeffs
   use io_mod

   implicit none
   private

   public :: sl_timewindow, set_planet, sl_solver_checkpoint, sl_solver_init, &
             sl_solver, sl_allocate_and_initialize_array, sl_deallocate_array, &
             sl_set_unit_num, sl_call_readnl, sl_print_counts, sl_outer_reset ! public module
   public :: iterstr ! public variables

   !============================================  Variables for the time window============================================|
   integer :: L_TW, TW_nfiles                     !    L_TW; length of a full time window                                  |
                                                  !    TW_nfiles: total number of ice files in a full time window          |
   integer :: ncalls, TW_nmelt, difference, ice   !    TW_nmelt: total number of melting episodes in the TW                |
   integer :: iter_TW                             !    nmelt_TW:  index variable for melting episodes in the TW            |
   integer :: Travel, Travel_total                !    number of travels that the timewindow has made                      |
                                                  !    Travel_total: total number of travelling                            |
   integer :: masked_iceload                      !    elements in mask and iceload arrays                                 |
   integer, dimension(:), allocatable :: mask, iceload, icefiles!                                                          |
   integer, dimension(:), allocatable :: TIMEWINDOW        !   TW array                                                    |
   integer :: err                                 !   I/O error                                                            |
   integer :: dummy                               ! dummy variable                                                         |
   integer, dimension(4) :: int_dt, Rdt, Ndt, Ldt         !   Values of internal timesteps of the time window              |
                                             !    Rdt: ratio between the smallest internal timestep (dt1)                  |
                                             !      to other internal time steps                                           |
                                             !    Ndt: total number of each internal timesteps                             |
                                             !    Ldt: an array for the time covered by each internal timestep             |
                                             !                                                                             |
   !=======================================================================================================================|

   !======================================== Variables for dynamic arrays==================================================|
   real, dimension(:), allocatable :: times        ! Timesteps of ice model (years)                                        |
   real, dimension(:,:,:), allocatable :: icexy    ! Spatial inputs of ice model                                           |
   real, dimension(:,:,:), allocatable :: sl       ! Big arrays of sea level change                                        |
   complex, dimension(:,:,:), allocatable:: dicestar,dS,deltaicestar,deltaS ! Big arrays of changes in loads               |
                                                                                          !  used in Love number viscous   |
                                                                                          !  response                      |
   real, dimension(:,:,:), allocatable :: rr,gg                         !  R and G (radial displacement and geoid change)  |
   complex, dimension(:,:,:), allocatable :: dlambda, deltalambda       ! Big arrays of changes in rotational driving      |
   real, dimension(:,:,:), allocatable :: dil                           ! Big array of changes in IL                       |
   real, dimension(:,:), allocatable :: dm                              ! Big array of changes in m                        |
   real, dimension(:), allocatable :: lovebetatt, lovebetattrr          ! Used in Love number calculations                 |
   real, dimension(:,:), allocatable :: lovebeta                        !                                                  |
   real, dimension(:,:), allocatable ::   lovebetarr                    ! Used in Love number calculations                 |
   !=======================================================================================================================|

   !========================= Variables for topography correction (i.e. outer iteration)===================================|
   real :: current_time                            ! time (in years) since the start of the simulation                     |
   real, dimension(:,:), allocatable :: tinit_0         ! topography at the very beginning of the simulation (i.e. time0)  |
   real, dimension(:,:), allocatable :: tinit_0_last    ! tinit_0 from the previous outer-iteration                        |
   real, dimension(:,:), allocatable :: pred_pres_topo  ! predicted present topography at current outer-iteration loop     |
   real, dimension(:,:), allocatable :: init_topo_corr  ! correction applied to compute tinit_0 at current outer-iter loop |
   character(6) :: iterstr                         ! String for timestep number for reading/writing files                  |
   !=======================================================================================================================|
   ! Inputs
   real, dimension(:,:), allocatable :: truetopo        ! Present-day topography
   real, dimension(:,:), allocatable :: rprime,r,s      ! Love numbers
   real, dimension(:), allocatable :: ke,he                ! Love numbers
   real, dimension(:,:), allocatable :: rprimeT,rT      ! Love numbers (tidal)
   real, dimension(:), allocatable :: kTE, hTE             ! Love numbers (tidal)

   ! grid lat-lon
   real, dimension(:), allocatable :: latgrid
   real, dimension(:), allocatable :: longrid

   ! Model calculations
   real :: gmslc_ocn, gmslc_ocnBeta                         ! Global mean sea-level change over time-evolving ocean areas
   real, dimension(:,:), allocatable :: glw_matrix               ! weigtht in the Gaussian-Legendre grid
   real, dimension(:,:), allocatable :: deltaslxy, dslxy         ! Total sea level change
   real, dimension(:,:), allocatable :: slcxy_ocn, slcxy_ocnBeta ! Total sea level change with topography masked out
                                                            !  total spatially heterogeneous sea level change
   real, dimension(:,:), allocatable :: icestarxy                ! Grounded ice thickness
   real, dimension(:,:), allocatable :: beta, cstarxy, cstar0    ! Grounded ice mask, ice-free ocean function
   real, dimension(:,:), allocatable :: tOxy, rOxy, tTxy         ! Projections used to calculate loads and shoreline migration
   complex, dimension(:,:), allocatable :: slclm_ocn,slclm_ocnBeta,clm,cstarlm,oldcstarlm,&
                                            tOlm,rOlm,load_inc_lm,old_load_inc_lm,icestarlm,dicestarlm,&
                                            deltaicestarlm,oldicestarlm,icestar0, &
                                            t0lm,oldt0lm,tTlm,oldtTlm,sl_inc_lm,sl_tot_lm  ! Above, in spectral domain

   real :: conserv                                          ! Uniform geoid shift (ΔΦ/g)
   real :: ttl, ekhl                                        ! Used in Love number calculations
   complex, dimension(:,:), allocatable :: viscous          ! Used in Love number calculations
   real :: xi, zeta                                         ! Convergence checks
   real :: ice_volume                                       ! ice volume
   ! For calculating R and G separately
   real, dimension(:,:), allocatable :: rrxy, drrxy_computed
   complex, dimension(:,:), allocatable :: rrlm, dgglm, drrlm_computed

   complex :: viscousrr

   ! Rotation calculations
   ! real, parameter :: --------------------------------->     ! Fluid Love number is defined in the planetary modules
   complex, dimension(0:2,0:2) :: lambda, oldlambda, lambda0   ! Rotational driving
   real, dimension(3) :: mm, oldm, sum_m                       ! Rotational perturbation vector
   real, dimension(3,3) :: il, oldil, sum_il                     ! Load perturbations to MOI
   complex, dimension(0:2,0:2) :: dsl_rot, rr_rot              ! Sea level change from rotational feedbacks
   real :: ekhTE                                               ! Used in Love number calculations
   real :: betatt, betattprime                                 ! Used in Love number calculations
   complex :: viscoustt,viscousttrr                            ! Used in Love number calculations

   ! Miscellaneous variables
   integer :: ninner                                           ! Iteration of inner loop
   integer :: i,j,k,l,m,n,nn                                   ! Do-loop indices
   character(6) :: numstr, numstr2                             ! String for timestep number for reading/writing files
   integer :: counti, countf,countrate                         ! Computation timing
   real :: counti_cpu, countf_cpu
   type(sphere) :: spheredat                                   ! SH transform data to be passed to subroutines

   ! For reading the Love numbers
   integer, dimension(:), allocatable :: legord,nmod,nmodes
   integer :: ll,nm,np
   real :: xn
   real, dimension(:,:), allocatable :: elast,asymv,telast,tasymv
   real, dimension(:,:), allocatable :: resh,resl,resk,tresh,tresl,tresk
   real :: taurr,taurt,dmx,dmy

   real, dimension(:,:), allocatable :: beta0, cxy0, cxy
   real, dimension(:,:), allocatable :: topoxy, topoxy_m1, tinit
                                                 ! topoxy_m1: topogramy from the previous timestep (m1: minus one)
                                                 ! topoxy: topography at the currect timestep
                                                 ! tinit: initial topography within the TimeWindow

   complex, dimension(:,:), allocatable :: dS_converged
   integer :: nmelt, nfiles                                ! Number of melting episodes up to the current timestep
!   real :: starttime                                       ! start time of the simulation
!   integer :: nmelt, nfiles, iter, itersl, dtime           ! Number of melting episodes up to the current timestep
   integer :: iargc, nargs                                 ! Arguments read in from a bash script
   character(16) :: carg(20)                               ! Arguments from a bash script
   character(3) :: skip                                    ! variable used to skip lines in reading TPW file

   ! Exact hereditary-memory state.  The standalone legacy driver restarts the
   ! executable every physical timestep, so this state is checkpointed in
   ! outputfolder and reloaded on the next call.
   type(response_spectrum) :: exact_spec
   type(response_memory_state) :: exact_state
   complex(dp), allocatable :: exact_viscous(:), exact_radial(:), exact_current_load(:)
   real(dp) :: exact_sum_il(3,3), exact_sum_m(3)
   complex(dp) :: exact_tsl(0:2), exact_trr(0:2)
   ! Optional horizontal-displacement read-out (SLEVOLVE_HORIZ >= 1):
   ! Psi_lm = T_l [ l^E_l L_lm + sum_k (l_k/s_k) R_k,lm ] (+ degree-2 rotation),
   ! u_h = grad_h Psi.  Same modal states as the radial R; output only.
   complex(dp), allocatable :: exact_horiz(:)
   complex(dp) :: exact_thl(0:2)
   complex, dimension(:,:), allocatable :: hzlm
   real, dimension(:,:), allocatable :: hzxy
   logical :: pp_hz_needed
   character(len=1024) :: exact_statefile
   integer :: exact_idx
   integer, parameter :: exact_prev=1, exact_cur=2

   !--------------------------------------------------------------------
   ! everything below is state that the baseline re-creates from
   ! disk on every timestep.  In single-process mode it stays resident.
   !--------------------------------------------------------------------
   logical :: pp_arrays_ready = .false.   ! big arrays allocated once

   ! inner-loop acceleration state (SLEVOLVE_INNER)
   complex, allocatable :: extrap_dsm1(:,:), extrap_dsm2(:,:)      ! converged dS of steps n-1, n-2
   integer :: inner_nhist = 0                                   ! converged steps stored (max 2 used)
   real :: extrap_dtm1 = 0.0, extrap_dtm2 = 0.0                    ! time steps of steps n-1, n-2 (for the extrapolation)
   ! optional non-uniform time steps (times file, see sl_read_times): t_k in years, k = 0..N
   logical :: vt_on = .false., vt_checked = .false.
   real(8), allocatable :: vt_times(:)
   character(len=1024) :: vt_path = ''
   complex, allocatable :: inner_hf(:,:,:), inner_hg(:,:,:)      ! Anderson ring: residuals F_j, images G_j
   real(8) :: inner_gram(6,6)                                   ! <F_i,F_j> by ring slot
   integer :: inner_slot(6)                                     ! chronological list of ring slots
   integer :: inner_nh = 0                                      ! entries in the ring
   logical :: inner_lastmixed = .false.
   real(8) :: inner_lastfn = 0.0d0
   logical :: inner_frozen_it, inner_plain_rest, inner_did_reguess, inner_conv
   integer :: inner_ndiff, inner_ndiff2, inner_nrefreeze
   ! SLEVOLVE_REUSE=warm state (see warm_read / warm_write)
   logical :: warm_have = .false.                          ! stored state loaded for this step
   logical :: warm_infile = .false.                        ! reuse_<itersl-1>.bin present and valid
   logical :: warm_init = .false.
   complex, allocatable :: warm_ds(:,:)                    ! stored dS_n
   integer(1), allocatable :: warm_mask(:,:)               ! stored ocean mask of step n
   integer :: warm_nwarm = 0
   real(8) :: warm_tio = 0.0d0                             ! wall s spent in reuse I/O
   ! in-process outer iteration (driver run_outer_mem.f90)
   logical, public :: pp_mem_outer = .false.     ! set by the driver: passes loop inside this process
   logical :: pp_mem_have = .false.              ! pred_pres_topo / tinit_0_last / truetopo hold the previous pass
   logical, public :: pp_icecache = .true.       ! mem mode: keep the ice grids resident after first read
   real, allocatable :: pp_icec(:,:,:)
   logical, allocatable :: pp_icec_ok(:)
   complex, allocatable :: warm_sds(:,:)        ! mem-mode warm store: (ntri, nsteps) dS lower triangle
   integer(1), allocatable :: warm_smk(:,:,:)   ! mem-mode warm store: masks
   logical, allocatable :: warm_sok(:)
   real(8), public :: pp_t_init = 0.0d0, pp_t_ice = 0.0d0, pp_t_topo = 0.0d0, pp_t_quant = 0.0d0
   integer(8) :: inner_cnt(6) = 0   ! mixed, fallback, forced mask updates, frozen iters, cstar analyses skipped, steps
   logical :: pp_love_ready   = .false.   ! Earth model parsed once
   logical :: pp_t0lm_ready   = .false.   ! tinit constant => t0lm constant
   logical :: pp_state_ready  = .false.   ! hereditary state live in memory
   logical :: pp_final_free   = .false.   ! set to really deallocate
   logical :: pp_batch        = .true.    ! batch independent SH analyses
   real :: pp_sabs_cur = 0.0, pp_sabs_old = 0.0   ! cached sum(abs(load_inc_lm))
   real, dimension(:,:,:), allocatable :: pp_gbuf   ! batched-transform grids
   complex, dimension(:,:,:), allocatable :: pp_ubuf ! batched-transform coeffs
   ! transform-count reduction state
   logical :: opt_ready = .false.
   real, dimension(:), allocatable :: opt_g      ! (nlats) degree-0 quadrature row weights
   real, dimension(:), allocatable :: opt_k00    ! (nlats) synth(e00) per latitude row
   real, dimension(:,:), allocatable :: opt_rs   ! (nlats,4) latitude-row sums
   real :: opt_q00 = 0.0, opt_c00 = 0.0, opt_cst00 = 0.0, opt_slc00 = 0.0, opt_slcb00 = 0.0
   integer :: opt_nf = 1, opt_st = 0
   logical :: opt_rg_needed = .true.
   character(6) :: opt_stepstr = ''
   integer(8) :: opt_cnt(4) = 0_8   ! analysis fields, analysis sweeps, synthesis fields, synthesis sweeps

   contains

   !======================================================================================================================!
   subroutine set_planet
      ! Planetary values
      if (whichplanet == 'earth' .or. whichplanet == 'Earth' .or. whichplanet == 'EARTH') then
         call earth_init
      elseif (whichplanet == 'mars' .or. whichplanet == 'Mars' .or. whichplanet == 'MARS') then
         call mars_init
      else
         write(unit_num,*) 'The parameters for the planet you entered are not built in.'
         write(unit_num,*) 'Please check your spelling for the variable whichplanet in the user_specs_mod module.'
         write(unit_num,*) 'If you preferred, you could create a new subroutine for the new planet in the planets_mod module.'
         write(unit_num,*) 'Terminating: program sl_model'
         stop
      endif
   end subroutine set_planet
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_call_readnl

      call sl_readnl(inputfolder_ice, inputfolder, planetfolder, gridfolder, &
                        outputfolder, outputfolder_ice, folder_coupled, ext, fType_in, fType_out, &
                        planetmodel, icemodel, icemodel_out, timearray, &
                        topomodel, topo_initial, grid_lat, grid_lon, &
                        checkmarine, tpw, calcRG, input_times, &
                        initial_topo, iceVolume, coupling, patch_ice, &
                        L_sim, dt1, dt2, dt3, dt4, Ldt1, Ldt2, Ldt3, Ldt4, whichplanet)

      call sl_get_model_res(norder, nglv)

   end subroutine sl_call_readnl
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_timewindow(iter) !create a time window that points to which ice history to read in

      integer :: iter

      ! save internal timestep profiles to a big array 'int_dt'
      int_dt(1) = dt4
      int_dt(2) = dt3
      int_dt(3) = dt2
      int_dt(4) = dt1

      ! Ldt#: total length of time over which each dt# covers
      ! save into a big array
      Ldt(1) = Ldt4
      Ldt(2) = Ldt3
      Ldt(3) = Ldt2
      Ldt(4) = Ldt1

      if (sum(Ldt) > L_sim) then
          write(unit_num,*) 'Total length of simulation CANNOT be smaller than the total lengths of the internal time windows!!'
          write(unit_num,*) 'Make sure the sum(LdT) ≤ L_sim !!'
          write(unit_num,*) 'Terminating: program sl_model'
          stop
      endif

      ! Compute various variables
      L_TW = 0                    ! Initializing the length of the timewindow
      TW_nmelt = 0                 ! Initializing the total number of melting episodes within a FULL time window

      do i=1,size(int_dt)
         Rdt(i) = int_dt(i)/dt1        ! e.g. dt1/dt1, dt2/dt1, dt3/dt1, dt4/dt1
         Ndt(i) = Ldt(i)/int_dt(i)     ! the number of each internal timestep
         L_TW =  L_TW + Ldt(i)         ! total length of the time window
         TW_nmelt = TW_nmelt + Ndt(i)  ! total number of melting episodes in the full timewindow
      enddo

      ncalls = L_TW/dt1        ! total number of melting episodes (i.e. iterations) in case of uniform timesteps within the TW
      TW_nfiles = TW_nmelt + 1  ! total number of ice files within in the FULL time window
      Travel_total = (L_sim - L_TW)/dt1 !total number of marching steps taken by the TW

      if (Travel_total == 0) then
         write(unit_num,*) 'The length of a total simulation and the length of timewindow are the same.'
         write(unit_num,*)' TW will not march forward.'
      endif

      write(unit_num,*) 'Total length of time window in years (L_TW)', L_TW
      write(unit_num,*) 'Total number of melting episode within a Full TW (TW_nmelt)', TW_nmelt
      write(unit_num,*) 'Total number of icefiles in a FULL TW (TW_nfiles)', TW_nfiles

      if (iter.LT.ncalls) then !while the time window is growing, the time window has not started travelling when
         Travel = 0
         if (iter == 0) then
            nfiles = 1
         endif
      else
         Travel = iter - ncalls ! This is the number of marching steps the full TW has taken upto the current timestep
      endif

      ! allocate dimensions to arrays (allocated once and reused)
      if (allocated(mask)) then
         if (size(mask) /= ncalls+1) deallocate(mask,iceload,icefiles)
      endif
      if (.not.allocated(mask)) allocate (mask(ncalls+1),iceload(ncalls+1),icefiles(ncalls+1))
      mask(:) = 0
      iceload(:) = 0
      icefiles(:) = 0
      if (allocated(TIMEWINDOW)) then
         if (size(TIMEWINDOW) /= TW_nfiles) deallocate(TIMEWINDOW)
      endif
      if (.not.allocated(TIMEWINDOW)) allocate (TIMEWINDOW(TW_nfiles))
      TIMEWINDOW(:)=0

      ! Initialize and grow the time window when t>0 but Travel == 0
      ! initialize the mask array and icefile numbers for the timewindow
      do i = 1, ncalls+1
          iceload(i) = i-1 ! iceload number calculated from ice sheet model
          mask(i) = 0      ! set all elements of the mask vector to be zero
      end do
      mask(1) = 1          ! the first element of the mask will always be one

      ! fill in the mask array for the time window
      i = 1
      k = 1
      do while (i.LE.size(int_dt))
         do j = 1, Ndt(i)
            m =  k + Rdt(i)
            mask(m) = 1
            k = m
         enddo
         i=i+1
      enddo

      ! Multiply the mask created above to an array of iceloads
      !    & to find out which ice files to read in
      i=1
      do while (i.LE.ncalls+1)
          if (i.LE.iter+1) then
             j=1
             do while (j.LT.i+1)
                difference =ncalls+1-(i-j)
                icefiles(j) = iceload(j)*mask(difference)
                j=j+1
             end do
          elseif (i.GT.iter+1) then
             icefiles(i) = 0
          endif
          i=i+1
      end do

      ! Output the ice file numbers that are called in a simulation over the FIRST TIMEWINDOW
      k = 2
      do j=1,ncalls+1 !iter_end is equivalent to ncalls
          ice=icefiles(j)

          if (j.LE.1) then
             TIMEWINDOW(1) = 0
         endif
          if (j.GT.1 .AND. ice.GT.0) then
             TIMEWINDOW(k) = icefiles(j)
             k = k+1
             m = k
             do while (m.LT.TW_nmelt+1)
                TIMEWINDOW(m) = 0
                m = m + 1
             end do
          endif
      end do
      nmelt = k - 2    ! Update the number of metling episode  nmelt
      nfiles = nmelt+1 ! Number of icefiles read in at each time the sea-level model is called

      ! Find ice files to read in once the time window has fully grown and started travelling.
      if (iter==ncalls .OR. Travel.GT.0) then
         k=2
          do i=2,ncalls+1
             TIMEWINDOW(1) = 0
            masked_iceload = mask(i)*iceload(i)

             if (masked_iceload.GT.0) then
                TIMEWINDOW(k)=masked_iceload
                k=k+1
             endif
          end do

         ! TW is now marching forward
         do i=1,nfiles
            TIMEWINDOW(i) = TIMEWINDOW(i) + Travel
         enddo
         write(unit_num,*) 'Icefiles in the TW at current marching step:', TIMEWINDOW
         write(unit_num,*) 'The TW will stop marching when TRAVEL == Travel_total:'
         write(unit_num,*) ''
         write(unit_num,'(A,I4,A)') '   ', Travel_total - Travel, ' more TW steps to march!'
      endif

      !!!!!!!!!!!!!!!!!!!ATTENTION!!!!!!!!
      ! TIMEWINDOW = TIMEWINDOW+1
      ! TIMEWINDOW = -(TIMEWINDOW*100+starttime)
      ! CAUTION in modifying the ice file numbers

      !write(unit_num,*) 'Timewindow is growing, ice file numbers to read in :', TIMEWINDOW(1:nfiles)
      write(unit_num,*) 'number of melting episodes (nmelt):', nmelt
      write(unit_num,'(A,I4)') 'Total number of files in the current time window (nfiles) = ', nfiles
      write(unit_num,'(A,I4)') ' Number of marching steps the TW has taken =        ', Travel
      write(unit_num,*) ''

   end subroutine sl_timewindow
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_allocate_and_initialize_array

      ! the driver calls this every timestep; do the work once.
      ! compat bit 16 restores the baseline's per-step re-zeroing.
      call pp_read_options()
      if (pp_arrays_ready) then
         if (iand(pp_compat,16) == 0) return
         go to 9001
      endif
      pp_arrays_ready = .true.
      allocate (pp_gbuf(nglv,2*nglv,4), pp_ubuf(0:norder,0:norder,4))
      pp_gbuf = 0.0
      pp_ubuf = 0.0

      allocate (tinit_0(nglv,2*nglv), tinit_0_last(nglv,2*nglv), pred_pres_topo(nglv,2*nglv), &
                init_topo_corr(nglv,2*nglv), truetopo(nglv,2*nglv), glw_matrix(nglv,2*nglv), &
                slcxy_ocn(nglv,2*nglv), slcxy_ocnBeta(nglv,2*nglv), deltaslxy(nglv,2*nglv), &
                dslxy(nglv,2*nglv), icestarxy(nglv,2*nglv), beta(nglv,2*nglv), cstarxy(nglv,2*nglv), &
                cstar0(nglv,2*nglv), tOxy(nglv,2*nglv), rOxy(nglv,2*nglv), tTxy(nglv,2*nglv), &
                rrxy(nglv,2*nglv), drrxy_computed(nglv,2*nglv), beta0(nglv,2*nglv), cxy0(nglv,2*nglv), &
                cxy(nglv,2*nglv), topoxy(nglv,2*nglv), topoxy_m1(nglv,2*nglv), tinit(nglv,2*nglv), &
                latgrid(nglv), longrid(2*nglv))

      ! variables of size(npam,norder)
      allocate (rprime(npam,norder), r(npam,norder), s(npam,norder), rprimeT(npam,norder),&
                rT(npam,norder), resh(npam,norder), resl(npam,norder), resk(npam,norder),&
                tresh(npam,norder), tresl(npam,norder), tresk(npam,norder))

      ! variables of size (3,norder)
      allocate (elast(3,norder), asymv(3,norder), telast(3,norder), tasymv(3,norder))

      ! variables of size (norder)
      allocate (ke(norder), he(norder), kTE(norder), hTE(norder), legord(norder), &
                nmod(norder),nmodes(norder))

      ! variables of size (0:norder,0:norder)
      allocate (slclm_ocn(0:norder,0:norder), slclm_ocnBeta(0:norder,0:norder), clm(0:norder,0:norder), &
                cstarlm(0:norder,0:norder), oldcstarlm(0:norder,0:norder), tOlm(0:norder,0:norder), &
                rOlm(0:norder,0:norder), load_inc_lm(0:norder,0:norder), old_load_inc_lm(0:norder,0:norder), &
                icestarlm(0:norder,0:norder), dicestarlm(0:norder,0:norder), deltaicestarlm(0:norder,0:norder), &
                oldicestarlm(0:norder,0:norder), icestar0(0:norder,0:norder), t0lm(0:norder,0:norder), &
                oldt0lm(0:norder,0:norder), tTlm(0:norder,0:norder), oldtTlm(0:norder,0:norder), &
                sl_inc_lm(0:norder,0:norder), sl_tot_lm(0:norder,0:norder), viscous(0:norder,0:norder), &
                rrlm(0:norder,0:norder), dgglm(0:norder,0:norder), drrlm_computed(0:norder,0:norder), &
                dS_converged(0:norder,0:norder))

      ! Based on the calcualted number of files read in within the time window at each timestep,
      ! allocate arrays to the below variables.
      ! Exact memory: only previous/current rolling slots are required.
      ! TIMEWINDOW/nfiles still describe the physical chronology, but no field
      ! below grows with the number of elapsed timesteps.
      allocate (times(2), lovebetatt(1), lovebetattrr(1))
      allocate (lovebetarr(1,norder),lovebeta(1,norder))
      allocate (icexy(nglv,2*nglv,2),sl(nglv,2*nglv,2))
      allocate (dS(0:norder,0:norder,2),deltaS(0:norder,0:norder,2))
      allocate (dicestar(0:norder,0:norder,2), deltaicestar(0:norder,0:norder,2))
      allocate (rr(nglv,2*nglv,2),gg(nglv,2*nglv,2))
      allocate (dil(3,3,2), dlambda(0:2,0:2,2),deltalambda(0:2,0:2,2))
      allocate (dm(3,2))

      ! initialize the allocated arrays
9001  continue
      tinit_0 = 0.0
      tinit_0_last = 0.0
      pred_pres_topo = 0.0
      init_topo_corr = 0.0
      truetopo = 0.0
      glw_matrix = 0.0
      slcxy_ocn = 0.0
      slcxy_ocnBeta = 0.0
      deltaslxy = 0.0
      dslxy = 0.0
      icestarxy = 0.0
      beta = 0.0
      cstarxy = 0.0
      cstar0 = 0.0
      tOxy = 0.0
      rOxy = 0.0
      tTxy = 0.0
      rrxy = 0.0
      drrxy_computed = 0.0
      beta0 = 0.0
      cxy0 = 0.0
      cxy = 0.0
      topoxy = 0.0
      topoxy_m1 = 0.0
      tinit = 0.0
      latgrid = 0.0
      longrid = 0.0
      rprime = 0.0
      r = 0.0
      s = 0.0
      rprimeT = 0.0
      rT = 0.0
      resh = 0.0
      resl = 0.0
      resk = 0.0
      tresh = 0.0
      tresl = 0.0
      tresk = 0.0
      elast = 0.0
      asymv = 0.0
      telast = 0.0
      tasymv = 0.0
      ke = 0.0
      he = 0.0
      kTE = 0.0
      hTE = 0.0
      legord = 0.0
      nmod = 0.0
      nmodes = 0.0
      slclm_ocn = 0.0
      slclm_ocnBeta = 0.0
      clm = 0.0
      cstarlm = 0.0
      oldcstarlm = 0.0
      tOlm = 0.0
      rOlm = 0.0
      load_inc_lm = 0.0
      old_load_inc_lm = 0.0
      icestarlm = 0.0
      dicestarlm = 0.0
      deltaicestarlm = 0.0
      oldicestarlm = 0.0
      icestar0 = 0.0
      t0lm = 0.0
      oldt0lm = 0.0
      tTlm = 0.0
      oldtTlm = 0.0
      sl_inc_lm = 0.0
      sl_tot_lm = 0.0
      viscous = 0.0
      rrlm = 0.0
      dgglm = 0.0
      drrlm_computed = 0.0
      dS_converged = 0.0

      times = 0.0
      lovebetatt = 0.0
      lovebetattrr = 0.0
      lovebetarr = 0.0
      lovebeta = 0.0
      icexy = 0.0
      sl = 0.0
      dS = 0.0
      deltaS = 0.0
      dicestar = 0.0
      deltaicestar = 0.0
      rr = 0.0
      gg = 0.0
      dil = 0.0
      dlambda = 0.0
      deltalambda = 0.0
      dm = 0.0

   end subroutine sl_allocate_and_initialize_array
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_deallocate_array

      ! the arrays and the hereditary state stay resident for the
      ! whole run.  Only an explicit teardown frees them.
      if (.not. pp_final_free .and. iand(pp_compat,32) == 0) return
      if (iand(pp_compat,32) /= 0) then
         pp_arrays_ready = .false.
         pp_love_ready = .false.
         pp_state_ready = .false.
         pp_t0lm_ready = .false.
      endif
      deallocate(pp_gbuf, pp_ubuf)
      deallocate(mask, iceload, icefiles, TIMEWINDOW)
      deallocate(times, lovebetatt, lovebetattrr)
      deallocate(lovebetarr, lovebeta)
      deallocate(icexy, sl)
      deallocate(dS, deltaS)
      deallocate(dicestar, deltaicestar)
      deallocate(rr, gg)
      deallocate(dil, dlambda ,deltalambda)
      deallocate(dm)
      deallocate(tinit_0, tinit_0_last, pred_pres_topo, init_topo_corr, truetopo, glw_matrix, &
                 slcxy_ocn, slcxy_ocnBeta, deltaslxy, dslxy, icestarxy, beta, cstarxy, &
                 cstar0, tOxy, rOxy, tTxy, rrxy, drrxy_computed, beta0, cxy0, cxy, topoxy, topoxy_m1, tinit)
      deallocate(rprime, r, s, rprimeT, rT, resh, resl, resk, tresh, tresl, tresk)
      deallocate(elast, asymv, telast, tasymv)
      deallocate(ke, he, kTE, hTE, legord, nmod, nmodes)
      deallocate(slclm_ocn, slclm_ocnBeta, clm, cstarlm, oldcstarlm, tOlm, rOlm, load_inc_lm, old_load_inc_lm,&
                 icestarlm, dicestarlm, deltaicestarlm, oldicestarlm, icestar0, t0lm, oldt0lm,&
                 tTlm, oldtTlm, sl_inc_lm, sl_tot_lm, viscous, rrlm, dgglm, drrlm_computed, dS_converged)
      deallocate(latgrid, longrid)
      if (allocated(exact_viscous)) deallocate(exact_viscous)
      if (allocated(exact_radial)) deallocate(exact_radial)
      if (allocated(exact_current_load)) deallocate(exact_current_load)
      if (allocated(exact_horiz)) deallocate(exact_horiz)
      if (allocated(hzlm)) deallocate(hzlm)
      if (allocated(hzxy)) deallocate(hzxy)
      call exact_state%destroy()
      call exact_spec%destroy()

   end subroutine sl_deallocate_array
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_set_unit_num(unit_num_in)

      integer :: unit_num_in

         unit_num = unit_num_in

   end subroutine sl_set_unit_num
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_solver_checkpoint(itersl, dtime)

      integer :: itersl, dtime

      if (dtime /= dt1) then
         write(unit_num,*) 'dtime and dt1 should be equal to each other.'
         write(unit_num,*) 'Please check your set up for the variables'
         write(unit_num,*) 'Terminating: program sl_model'
         stop
      endif

      if (itersl.lt.1) then
          write(unit_num,*) 'itersl must be equal to or greather than 1'
          write(unit_num,*) 'itersl = 1: No topography correction'
          write(unit_num,*) 'itersl > 1: topography correction'
          write(unit_num,*) 'Terminating: program sl_model'
          stop
      endif

      if (coupling) then ! Set a default unit number for I/O
         write(unit_num,*) 'Sea level model is coupled to the ice sheet model. Unit number is:', unit_num
      else
         call sl_set_unit_num(6)! Set a default unit number for I/O
         write(unit_num,*) 'Sea level model is running standalone. Unit number is:', unit_num
      endif

      ! set up the planet profile
      call set_planet

      ! report the I/O mode and build only the Legendre tables the
      ! BLAS transform needs (saves the 1.1 GB reference pnm/hnm pair).
      call pp_read_options()
      call pp_report(unit_num)
#if defined(USE_BLAS_SHT) || defined(USE_EXACT_SHT)
      sht_lean = .true.
#endif
      call sht_env_init()
      write(unit_num,*) 'SLEvolve SHT backend: ', trim(sht_backend_name()), '  symmetric=', sht_sym, &
                        '  fftw=', sht_fftw, '  threads=', sht_nt()

      ! initalize spherical harmonics
      call spharmt_init(spheredat, 2*nglv, nglv, norder, radius) ! Initialize spheredat (for SH transform subroutines)

      call sl_read_times()

   end subroutine sl_solver_checkpoint
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_read_times
      ! Optional non-uniform time steps.  SLEVOLVE_TIMES=<path> selects a times file,
      ! SLEVOLVE_TIMES=off (or none / 0) disables it; when the variable is unset,
      ! <inputfolder>times (INPUT_FILES/others/times) is used if it exists.  The file holds
      ! one time t_k in years per line for the ice files k = 0..N (N = L_sim/dt1), strictly
      ! increasing; step k then has length t_k - t_(k-1) and output time stamps are t_k.
      ! Without a times file the build-time uniform step dt1 is used, exactly as before.
      character(len=1024) :: v
      integer :: st, u, ios, n, nneed
      logical :: ex
      real(8) :: x
      if (vt_checked) return
      vt_checked = .true.
      vt_on = .false.
      call get_environment_variable('SLEVOLVE_TIMES', v, status=st)
      if (st == 0) then
         if (trim(v) == 'off' .or. trim(v) == 'none' .or. trim(v) == '0' .or. len_trim(v) == 0) return
         vt_path = trim(v)
         inquire(file=trim(vt_path), exist=ex)
         if (.not. ex) then
            write(unit_num,*) 'SLEVOLVE_TIMES: file not found: ', trim(vt_path)
            error stop 'SLEVOLVE_TIMES file not found'
         endif
      else
         vt_path = trim(inputfolder)//'times'
         inquire(file=trim(vt_path), exist=ex)
         if (.not. ex) return
      endif
      nneed = L_sim/dt1 + 1
      open(newunit=u, file=trim(vt_path), status='old', action='read')
      n = 0
      do
         read(u,*,iostat=ios) x
         if (ios /= 0) exit
         n = n + 1
      enddo
      if (n < nneed) then
         write(unit_num,*) 'times file ', trim(vt_path), ': ', n, ' values, need ', nneed
         error stop 'times file too short'
      endif
      if (allocated(vt_times)) deallocate(vt_times)
      allocate(vt_times(0:nneed-1))
      rewind(u)
      do n = 0, nneed-1
         read(u,*) vt_times(n)
      enddo
      close(u)
      do n = 1, nneed-1
         if (.not. (vt_times(n) > vt_times(n-1))) error stop 'times file: times must be strictly increasing'
      enddo
      vt_on = .true.
      write(unit_num,'(A,A,A,I0,A,ES12.5,A,ES12.5,A,ES12.5,A)') ' non-uniform time steps from ', trim(vt_path), ': ', &
           nneed-1, ' steps, t0 = ', vt_times(0), ' yr, dt ', minval(vt_times(1:)-vt_times(:nneed-2)), &
           ' - ', maxval(vt_times(1:)-vt_times(:nneed-2)), ' yr'
   end subroutine sl_read_times
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================
   subroutine sl_solver_init(itersl, starttime, mali_iceload, mali_bedrock, mali_mask)

      integer :: itersl
      real :: starttime
      real, dimension(:,:), optional :: mali_iceload, mali_bedrock, mali_mask
      integer(8) :: pp_c0, pp_c1, pp_c2, pp_cr

      !================================================================================================
      !       Initialize the model at times(1) when there has been no melting episodes yet  NMELT = 0
      !================================================================================================

      write(unit_num,*) 'nmelt=0, INITIALIZING THE SEA LEVEL MODEL..'
      flush(unit_num)

      write(unit_num,*) 'Number of orders and degree of spherical harmonics is ', norder
      write(unit_num,*) 'The number of Gauss-Legendre nodes in latitude(nglv): ', nglv

      !initialize variables
      j = TIMEWINDOW(1) ! initial file number (i.e. 0)
      write(unit_num,'(A,I4)') 'initial file number:', j
      write(numstr,'(I4)') j
      numstr = trim(adjustl(numstr))
      !     k = -1*starttime
      !     write(numstr2,'(I6)') k
      !     write(unit_num,'(A,I6)') 'initial icefile suffix:', k
      !     numstr2 = trim(adjustl(numstr2))

      !====================== topography and ice load========================
      ! read in the initial iceload from the coupled ice input folder
      call system_clock(pp_c0, pp_cr)
      call pp_ice_get(icexy(:,:,1), j, numstr)
      !  Initialize topography (STEP 1)
      if (initial_topo) then
         write(unit_num,*) 'Reading in initial topo file'
         call pp_read_grid(tinit_0, topo_initial, inputfolder)
      else  ! if initial topo is unknown

         ! Present-day observed topography
         write(unit_num,*) 'Reading in ETOPO2 file'
         call system_clock(pp_c1)
         if (.not. (pp_mem_outer .and. pp_mem_have)) call read_sl(truetopo, topomodel, inputfolder)
         call system_clock(pp_c2); pp_t_topo = pp_t_topo + real(pp_c2-pp_c1,8)/real(pp_cr,8)

         ! if topography is not iteratively getting improved, or if at the first loop of outer-iteration
         if (itersl == 1) then
            write(unit_num,*) 'Initial topo is unknown, set the modern-observed topo to be the initial topo'
            ! truetopo(:,:)=truetopo(:,:)-icexy(:,:,TW_nfiles) ! If using ice topography
            tinit_0(:,:) = truetopo(:,:)! (eq. 48)
          else

             write(unit_num,*) 'Topographic correction is ON. Updating initial topography'
             ! read in predicted present topo from the previous outer-iteration 'itersl-1'
             write(iterstr,'(I2)') itersl-1
             iterstr = trim(adjustl(iterstr))

             if (pp_mem_outer .and. pp_mem_have) then
                write(unit_num,*) 'in-process outer iteration: pred_pres_topo and tgrid0 of pass ', &
                                  itersl-1, ' held in memory'
             else if (outer_set) then
                call pp_read_grid(pred_pres_topo, 'pred_pres_topo_', outputfolder, suffix=iterstr)
                call pp_read_grid(tinit_0_last, 'tgrid0_', outputfolder, suffix=iterstr)
             else
             call read_sl(pred_pres_topo, 'pred_pres_topo_', outputfolder, suffix=iterstr)

             ! read in tinit_0 from the previous outer-iteration 'itersl-1'
             call read_sl(tinit_0_last, 'tgrid0_', outputfolder, suffix=iterstr)
             endif

             ! compute topography correction for the initial topography
             init_topo_corr(:,:) = truetopo(:,:) - pred_pres_topo(:,:)
             if (outer_mode == 1 .and. itersl >= 3) then
                call outer_secant(itersl)
             else
             tinit_0(:,:) = tinit_0_last(:,:) + init_topo_corr(:,:)
             endif
          endif
       endif

       if (coupling) then  !if coupling the ICE SHEET - SEA LEVEL MODELs
          write(unit_num,*) 'Merge initial topography and iceload with the ISM data'

          ! merge intitial topography with bedrock provided by the ice sheet model.
          do j = 1,2*nglv
             do i = 1,nglv
                tinit_0(i,j) = mali_bedrock(i,j) + tinit_0(i,j)*(1 - mali_mask(i,j))
             enddo
          enddo
          ! merge intitial topography with bedrock provided by the ice sheet model.
          if (patch_ice) then
             do j = 1,2*nglv
                do i = 1,nglv
                   icexy(i,j,nfiles) = 0.0
                enddo
             enddo
          endif

          do j = 1,2*nglv
             do i = 1,nglv
                icexy(i,j,nfiles) = mali_iceload(i,j) + icexy(i,j,nfiles)*(1 - mali_mask(i,j))
             enddo
          enddo

          !write out the current ice load as a new file to the sea-level model ice folder
          call write_sl(icexy(:,:,exact_cur), icemodel_out, outputfolder_ice, suffix=numstr)
      endif ! end if (coupling)

      !write out the initial topo of the simulation, tinit_0
      ! the baseline re-reads tgrid0 (ES16.9E2 text) at every step, so
      ! its tinit_0 is the quantized value.  Reproduce that once here, then
      ! keep tinit_0 resident.
      call system_clock(pp_c1)
      if (pp_mem_outer) then
         ! same ES16.9E2 round trip, through an internal file (no 8 MB text file)
         call pp_quantize_es169(tinit_0)
      else
      call write_sl(tinit_0(:,:), 'tgrid', outputfolder, suffix=numstr)
      call read_sl(tinit_0(:,:), 'tgrid', outputfolder, suffix=numstr)
      endif
      call system_clock(pp_c2); pp_t_quant = pp_t_quant + real(pp_c2-pp_c1,8)/real(pp_cr,8)
      call pp_write_grid(tinit_0(:,:), 'tgrid', outputfolder, suffix=numstr)


      !========================== ocean function =========================
      ! Calculate an initial ocean function based on the present topography as a first guess
      do j = 1,2*nglv
         do i = 1,nglv
            if (tinit_0(i,j)<0) then
               cxy0(i,j) = 1
            else
               cxy0(i,j) = 0
            endif
         enddo
      enddo

      !  write out the initial ocean function as a file
      call pp_write_grid(cxy0(:,:), 'ocean', outputfolder, suffix=numstr)

      ! write out the ocean area
      call spat2spec(cxy0, clm, spheredat)
      open(unit = 201, file = trim(outputfolder)//'ocean_area', form = 'formatted', access ='sequential', &
      & status = 'replace')
      write(201,'(ES14.4E2)') real(clm(0,0))*4*pi*radius**2
      close(201)

      !========================== beta function =========================
      ! calculate initial beta
      do j = 1,2*nglv
         do i = 1,nglv
            if (icexy(i,j,1)==0) then
               beta0(i,j)=1
            else
               beta0(i,j)=0
            endif
         enddo
      enddo

      !  write out the initial beta function as a file
      call pp_write_grid(beta0(:,:), 'beta', outputfolder, suffix=numstr)

      ! write out the ocean area excluding the marine-based region
      cstarxy(:,:) = cxy0(:,:) * beta0(:,:)
      call spat2spec(cstarxy,cstarlm,spheredat)
      open(unit = 201, file = trim(outputfolder)//'oceanBeta_area', form = 'formatted', access ='sequential', &
      & status = 'replace')
      write(201,'(ES14.4E2)') real(cstarlm(0,0))*4*pi*radius**2
      close(201)

      !================== total ocean loading change =====================
      ! initialize the total ocean loading change and output as a file
      deltaS(:,:,1) = (0.0,0.0)
      call pp_write_cplx(deltaS(:,:,1), 'dS_converged', outputfolder, suffix=numstr, fmt='(ES16.9E2)')


      !========================== computing time =========================
      ! To write out how much time it took to compute sea-level change over one step
      ! Open a new file
      open(unit = 201, file = trim(outputfolder)//'elapsed_wall_time', form = 'formatted', access = 'sequential', &
      & status = 'replace')
      close(201)

      open(unit = 201, file = trim(outputfolder)//'elapsed_cpu_time', form = 'formatted', access = 'sequential', &
      & status = 'replace')
      close(201)


      !========================== time array =============================
      if (.not. input_times) then !if time array is not read in from a text file, make a new one
         ! write a new file
         open(unit = 201, file = trim(outputfolder)//trim(timearray), form = 'formatted', access = 'sequential', &
         & status = 'replace')
         if (vt_on) then
            write(201,'(ES24.16E3)') vt_times(0)
         else
            write(201,'(ES14.4E2)') starttime
         endif
         close(201)
      endif


      !=========================== TPW ====================================
      ! initialize the rotational components
      if (tpw) then
      !dil(:,:,1) = 0.0
      !dm(:,1) = 0.0
      !dlambda(:,:,1) = (0.0,0.0)

      il(:,:) = 0.0
      mm(:) = 0.0
      lambda(:,:) = 0.0

      ! write the values (0.0) for the first timestep
      open(unit = 201, file = trim(outputfolder)//'TPW'//trim(numstr), form = 'formatted', access = 'sequential', &
      & status = 'replace')
      write(201,'(9ES19.8E2/,3ES19.8E2/,18ES19.8E2)') il(:,:), mm(:), lambda(:,:)
      ! write(unit_num,'(9ES19.8E2/,3ES19.8E2/,18ES19.8E2)') dil(:,:,1), dm(:,1), dlambda(:,:,1)
      close(201)
      endif


       !=========================ice volume================================
      if (iceVolume) then
         if (checkmarine) then
            do j = 1,2*nglv
               do i = 1,nglv
                  if (tinit_0(i,j) > 0) then
                     ! If not marine...
                     icestarxy(i,j) = icexy(i,j,1)
                  elseif (icexy(i,j,1) > (abs(tinit_0(i,j)) * rhow / rhoi)) then
                     !...else if marine, but thick enough to be grounded
                     icestarxy(i,j) = icexy(i,j,1)
                  else
                     !...if floating ice
                     icestarxy(i,j) = 0
                  endif
               enddo
            enddo
                  ! Decompose ice field
            call spat2spec(icestarxy(:,:),icestarlm(:,:),spheredat)
         else ! If not checking for floating ice
            call spat2spec(icexy(:,:,1),icestarlm(:,:),spheredat) ! Decompose ice field
         endif

         ice_volume = icestarlm(0,0)*4*pi*radius**2

         open(unit = 201, file = trim(outputfolder)//'ice_volume', form = 'formatted', access ='sequential', &
         & status = 'replace')
         write(201,'(ES14.4E2)') ice_volume
         close(201)
      endif

      !=========================== GMSLE ==================================
      gmslc_ocn = 0.0
      open(unit = 201, file = trim(outputfolder)//'gmslc_ocnArea', form = 'formatted', access ='sequential', &
      & status = 'replace')
      write(201,'(ES14.4E2)') gmslc_ocn
      close(201)

      gmslc_ocnBeta = 0.0
      open(unit = 201, file = trim(outputfolder)//'gmslc_ocnBetaArea', form = 'formatted', access ='sequential', &
      & status = 'replace')
      write(201,'(ES14.4E2)') gmslc_ocnBeta
      close(201)

      !print out the number of iteration it takes for the inner convergence
      open(unit = 201, file = trim(outputfolder)//'numiter', form = 'formatted', access ='sequential', &
      & status = 'replace')
      close(201)

      !print out the number of melting episodes taken into account at the current timestep (i.e., nfiles-1)
      open(unit = 201, file = trim(outputfolder)//'nmelt', form = 'formatted', access ='sequential', &
      & status = 'replace')
      close(201)

      call system_clock(pp_c2); pp_t_init = pp_t_init + real(pp_c2-pp_c0,8)/real(pp_cr,8)
      write(unit_num,*) 'DONE INITIALIZATION.'
   !   call exit
   end subroutine sl_solver_init

   !========================================================================================================================
   subroutine sl_solver(itersl, iter, dtime, starttime, mali_iceload, mali_mask, slchange)
      ! Compute sea-level change associated with past ice loading changes

      real :: starttime
      integer :: iter, itersl, dtime
      real :: extrap_c1, extrap_c2
      real(8) :: dtr   ! length of this step in years: dtime, or t_iter - t_(iter-1) from the times file
      real, dimension(:,:), optional :: mali_iceload, mali_mask ! variables for coupled ISM-SLM simulations
      real, dimension(:,:), intent(out), optional :: slchange ! variable exchanged with the ISM

      !===========================================================
      !                   BEGIN TIMING & EXECUTION
      !___________________________________________________________
      call cpu_time(counti_cpu)
      call system_clock(count = counti, count_rate = countrate)  ! Total computation time
      call opt_setup()

      !-----------------------------------------------------------------
      ! the baseline re-allocates *and zeroes* every working array
      ! at the top of each timestep, and the solver depends on that in at
      ! least one place: the rotation block evaluates
      !     il(i,3) = ... rhow*deltaS(2,m,exact_cur) + ...
      ! on the FIRST inner iteration, before deltaS(:,:,exact_cur) has been
      ! assigned in that step -- so the baseline's first rotation guess uses
      ! a zero current ocean load.  With resident arrays it would instead see
      ! the previous timestep's converged load, which changes the converged
      ! solution (TPW by ~3e-3 relative at step 2) and occasionally the inner
      ! iteration count.  Re-zeroing here costs ~3 ms/step and keeps the result
      ! identical to the legacy code.
      !-----------------------------------------------------------------
      ! only `deltaS` is actually read stale (the rotation
      ! block above).  A SLEVOLVE_COMPAT=256 bisection found
      ! nothing else in the 220 MB list whose stale contents the baseline
      ! semantics depend on, and SLEVOLVE_ZERO=min was verified to give
      ! identical output over 50 steps.  SLEVOLVE_ZERO=full restores
      ! the full memset (~220 MB/step).
      if (pp_zero_full) then
      slcxy_ocn=0.0; slcxy_ocnBeta=0.0; deltaslxy=0.0; dslxy=0.0; icestarxy=0.0
      beta=0.0; cstarxy=0.0; cstar0=0.0; tOxy=0.0; rOxy=0.0; tTxy=0.0
      rrxy=0.0; drrxy_computed=0.0; topoxy=0.0; topoxy_m1=0.0
      slclm_ocn=0.0; slclm_ocnBeta=0.0; clm=0.0; cstarlm=0.0; oldcstarlm=0.0
      tOlm=0.0; rOlm=0.0; load_inc_lm=0.0; old_load_inc_lm=0.0; icestarlm=0.0; dicestarlm=0.0
      deltaicestarlm=0.0; oldicestarlm=0.0; icestar0=0.0; oldt0lm=0.0
      tTlm=0.0; oldtTlm=0.0; sl_inc_lm=0.0; sl_tot_lm=0.0; viscous=0.0
      rrlm=0.0; dgglm=0.0; drrlm_computed=0.0; dS_converged=0.0
      icexy(:,:,exact_prev)=0.0; sl=0.0
      dS=0.0; deltaS=0.0; dicestar=0.0; deltaicestar=0.0
      rr=0.0; gg=0.0; dil=0.0; dlambda=0.0; deltalambda=0.0; dm=0.0
      times=0.0
      else
      deltaS=0.0
      endif


      ! The exact-memory path is restricted to the
      ! standard/full-history chronology.  This prevents accidental mixing of
      ! a persistent all-history state with the legacy thinned/moving time window.
      if (nfiles /= iter + 1) then
         error stop 'exact-memory stage-2 requires standard/full-history TIMEWINDOW'
      endif
      do n = 1,nfiles
         if (TIMEWINDOW(n) /= n-1) then
            error stop 'exact-memory stage-2 requires consecutive history indices'
         endif
      enddo

      ! History-free input path.  The coupled MALI path can be
      ! added with the same rolling-state semantics, but this path
      ! is restricted to the standalone (uncoupled) configuration.
      if (coupling) error stop 'exact-memory stage-2 standalone path only'

      ! Read only the current ice field.  At iter=1 also read the time-zero
      ! field so the checkpoint can be initialized without any older history.
      if (iter == 1) then
         j = TIMEWINDOW(1)
         write(numstr,'(I6)') j
         numstr = trim(adjustl(numstr))
         call pp_ice_get(icexy(:,:,exact_prev), j, numstr)
      endif
      j = TIMEWINDOW(nfiles)
      write(numstr,'(I6)') j
      numstr = trim(adjustl(numstr))
      call pp_ice_get(icexy(:,:,exact_cur), j, numstr)

      ! Only previous/current physical times are needed by the exact-memory path.
      if (input_times) error stop 'exact-memory stage-2 input_times path not implemented'
      if (vt_on) then
         dtr = vt_times(iter) - vt_times(iter-1)
         times(exact_prev) = vt_times(TIMEWINDOW(nfiles-1))
         times(exact_cur)  = vt_times(TIMEWINDOW(nfiles))
      else
         dtr = real(dtime,8)
         times(exact_prev) = starttime + TIMEWINDOW(nfiles-1)*dt1
         times(exact_cur)  = starttime + TIMEWINDOW(nfiles)*dt1
      endif
      open(unit = 201, file = trim(outputfolder)//trim(timearray), form = 'formatted', access = 'sequential', &
      & status = 'old', position='append')
      if (vt_on) then
         write(201,'(ES24.16E3)') times(exact_cur)
      else
         write(201,'(ES14.4E2)') times(exact_cur)
      endif
      close(201)

      ! tinit_0, cxy0 and beta0 are written once at step 0 and never
      ! change, so they are kept resident instead of being re-read (3 x 8.4 MB
      ! of formatted text per timestep in the baseline).
      if (iand(pp_compat,1) /= 0) call pp_read_grid(tinit_0, 'tgrid0', outputfolder)
      j = TIMEWINDOW(1) ! first element of the time window as the initial file
      write(unit_num,'(A,I4)') 'file number of the first item in the TW:', j
      write(numstr,'(I4)') j
      numstr = trim(adjustl(numstr))

      if (iand(pp_compat,1) /= 0) then
         call pp_read_grid(cxy0(:,:), 'ocean', outputfolder, suffix=numstr)
         call pp_read_grid(beta0(:,:), 'beta', outputfolder, suffix=numstr)
      endif

      if (nmelt == 1) then
         cxy(:,:) = cxy0(:,:)
      endif

      ! Previous rotational totals are restored from the exact-memory checkpoint
      ! after the Love spectrum/state has been initialized below.

      ! find index of the previous timestep
      m = TIMEWINDOW(nfiles-1)
      write(numstr2,'(I4)') m
      numstr2 = trim(adjustl(numstr2))

      ! the baseline reads tgrid<prev> into topoxy_m1 here, but
      ! topoxy_m1 is only consumed on the (disabled) coupling path.  The
      ! 8.4 MB text read is dropped.

      ! condition for time window travel
      if (Travel.EQ.0) then
         !if the TW hasnt started travelling, initial topography of the TW is that of the total simulation
         tinit(:,:) = tinit_0(:,:)
      else
         !if The TW is moving, initial topography of the TW is read in
         j = TIMEWINDOW(1)
         write(numstr,'(I4)') j
         numstr = trim(adjustl(numstr))

        call read_sl(tinit(:,:), 'tgrid', outputfolder, suffix=numstr)
      endif
      !endif nmelt>0

      ! cxy already holds the converged ocean function of the previous
      ! timestep (the very field the baseline re-reads from ocean<prev>).
      if (nmelt.GT.1) then
         j = TIMEWINDOW(nfiles-1)
         write(unit_num,*) 'Previous ocean function held in memory (file number)', j
         if (iand(pp_compat,1) /= 0) then
            write(numstr,'(I4)') j
            numstr = trim(adjustl(numstr))
            call pp_read_grid(cxy(:,:), 'ocean', outputfolder, suffix=numstr)
         endif
      endif

      !-----------------------------------------------------------
      !  Read in Love numbers (collocation Love-number text file)
      !  the Earth model is time independent -- parse the 1.3 MB
      !  text file and build the response spectrum once per run.
      !-----------------------------------------------------------
      if (.not. pp_love_ready .or. iand(pp_compat,2) /= 0) then
      open(unit = 201, file = trim(planetfolder)//trim(planetmodel), status = 'old')
      ! Reader as in the legacy code of Han et al. (2022)
      read(201,*)
      do j = 1,norder
         read(201,*) legord(j), nmodes(j)
         nm = nmodes(j)
         xn = real(legord(j))
         ll = legord(j)
         nmod(ll) = nm
         read(201,*) (s(i,ll), i=1, nm)
         read(201,*) elast(1,ll), elast(2,ll), taurr, taurt, elast(3,ll)
         read(201,*) asymv(1,ll), asymv(2,ll), taurr, taurt, asymv(3,ll)
         read(201,*) (resh(i,ll), i=1,nm)
         read(201,*) (resl(i,ll), i=1,nm)
         read(201,*) (resk(i,ll), i=1,nm)
         if (xn .lt. 1.5) cycle
         read(201,*) telast(1,ll), telast(2,ll), dmx, dmy, telast(3,ll)
         read(201,*) tasymv(1,ll), tasymv(2,ll), dmx, dmy, tasymv(3,ll)
         read(201,*) (tresh(i,ll), i=1, nm)
         read(201,*) (tresl(i,ll), i=1, nm)
         read(201,*) (tresk(i,ll), i=1, nm)
      enddo
      close(201)

      ! Divide by l (numbers are multiplied by l in the Love-number file)
      do l = 1,norder
         resl(:,l) = resl(:,l) / real(l)
         resk(:,l) = resk(:,l) / real(l)
         tresl(:,l) = tresl(:,l) / real(l)
         tresk(:,l) = tresk(:,l) / real(l)
      enddo

      ! Assign Love number inputs to the letters used in Kendall et al
      do l = 1,norder
         he(l) = elast(1,l)
         ke(l) = elast(3,l) / real(l)
         hTE(l) = telast(1,l)
         kTE(l) = telast(3,l) / real(l)
      enddo
      rprime(:,:) = resk(:,:)
      r(:,:) = resh(:,:)
      rprimeT(:,:) = tresk(:,:)
      rT(:,:) = tresh(:,:)

      ! ------------------------------------------------------------------
      ! Exact hereditary state for all old loads.  The state is advanced
      ! exactly once per physical timestep and is frozen during the inner
      ! nonlinear shoreline/ocean iteration.
      exact_spec%lmax = norder
      exact_spec%max_modes = npam
      if (.not.allocated(exact_spec%nmod)) then
      allocate(exact_spec%nmod(0:norder))
      allocate(exact_spec%s(npam,0:norder))
      allocate(exact_spec%a_sl(npam,0:norder), exact_spec%a_rr(npam,0:norder))
      allocate(exact_spec%a_il(npam,0:norder), exact_spec%a_m(npam,0:norder))
      allocate(exact_spec%a_tsl(npam,0:norder), exact_spec%a_trr(npam,0:norder))
      endif
      exact_spec%nmod = 0
      exact_spec%s = 0.0_dp
      exact_spec%a_sl = 0.0_dp; exact_spec%a_rr = 0.0_dp
      exact_spec%a_il = 0.0_dp; exact_spec%a_m = 0.0_dp
      exact_spec%a_tsl = 0.0_dp; exact_spec%a_trr = 0.0_dp
      if (pp_horiz > 0) then
         if (.not.allocated(exact_spec%a_hl)) &
            allocate(exact_spec%a_hl(npam,0:norder), exact_spec%a_thl(npam,0:norder))
         exact_spec%a_hl = 0.0_dp; exact_spec%a_thl = 0.0_dp
         ! resl/tresl were divided by l above (file stores l*l_k)
         do l = 1,norder
            do k = 1,nmod(l)
               exact_spec%a_hl(k,l) = resl(k,l)/s(k,l)
               if (l == 2) exact_spec%a_thl(k,l) = tresl(k,l)/s(k,l)
            enddo
         enddo
      endif
      do l = 1,norder
         exact_spec%nmod(l) = nmod(l)
         do k = 1,nmod(l)
            exact_spec%s(k,l) = s(k,l)
            exact_spec%a_sl(k,l) = (rprime(k,l)-r(k,l))/s(k,l)
            exact_spec%a_rr(k,l) = r(k,l)/s(k,l)
            if (l == 2) then
               exact_spec%a_il(k,l) = rprime(k,l)/s(k,l)
               exact_spec%a_m(k,l) = rprimeT(k,l)/s(k,l)
               exact_spec%a_tsl(k,l) = (rprimeT(k,l)-rT(k,l))/s(k,l)
               exact_spec%a_trr(k,l) = rT(k,l)/s(k,l)
            endif
         enddo
      enddo
      if (.not.allocated(exact_state%p_load)) then
         call exact_state%init(exact_spec)
         allocate(exact_viscous(exact_state%ncoef), exact_radial(exact_state%ncoef), &
                  exact_current_load(exact_state%ncoef))
      endif
      write(iterstr,'(I2)') itersl
      iterstr = trim(adjustl(iterstr))
      exact_statefile = trim(outputfolder)//'exact_memory_state_'//trim(iterstr)//'.bin'
      pp_love_ready = .true.
      endif ! .not. pp_love_ready

      ! the 53 MB hereditary state stays in memory across timesteps;
      ! the baseline writes and reads it back on every single step.
      if (iter == 1) then
         call exact_state%zero()
         pp_state_ready = .true.
      else if (.not. pp_state_ready .or. iand(pp_compat,4) /= 0) then
         call exact_state%load_checkpoint(trim(exact_statefile),exact_spec)
         pp_state_ready = .true.
      endif
      call exact_state%advance_and_read(exact_spec,dtr,exact_viscous,exact_radial, &
                                        exact_sum_il,exact_sum_m,exact_tsl,exact_trr)
      viscous(:,:) = (0.0,0.0)
      do l = 1,norder
         do m = 0,l
            exact_idx = lm_index(l,m)
            viscous(l,m) = exact_viscous(exact_idx)
         enddo
      enddo

      ! Restore the O(1) cumulative states needed by the instantaneous
      ! (elastic/conservation/rotation) equations.  No historical dS/TPW
      ! files are replayed.
      deltaS(:,:,exact_prev) = (0.0,0.0)
      do l = 0,norder
         do m = 0,l
            exact_idx = lm_index(l,m)
            deltaS(l,m,exact_prev) = exact_state%p_ocean(exact_idx)
         enddo
      enddo
      oldil(:,:) = exact_state%p_il(:,:)
      oldm(:) = exact_state%p_m(:)
      deltalambda(:,:,exact_prev) = (0.0,0.0)
      if (norder >= 2) deltalambda(2,0:2,exact_prev) = exact_state%p_lambda(0:2)

      !===========================================================
      !                       CALCULATIONS
      !___________________________________________________________

      write(unit_num,*) ''
      write(unit_num,'(A,I4,A,EN15.4E2,A)') '  Timestep',iter,', from ',times(exact_prev),' years to '
      write(unit_num,'(A,EN15.4E2,A)') '                     ',times(exact_cur),' years'

      ! Decompose initial topography (STEP 2) (used to check convergence of outer loop)
      ! tinit is time independent while the time window is not
      ! travelling, so this transform is done once.
      if (Travel /= 0) pp_t0lm_ready = .false.
      if (.not. pp_t0lm_ready .or. iand(pp_compat,8) /= 0) then
         call opt_ana(tinit(:,:), t0lm)
         pp_t0lm_ready = .true.
      endif

      !=====================================================================================
      !                   1. BEGIN ICE PART
      !=====================================================================================

      ! Transform only the current grounded-ice field.  The previous absolute
      ! spectrum and cumulative ice change are carried in the checkpoint.
      ! the grounded-ice grid does not depend on anything
      ! the ocean part computes, so its analysis rides in the first-guess
      ! sweep below (slot 1 of pp_gbuf) instead of costing a sweep of its own.
      if (checkmarine) then
         !$omp parallel do schedule(static) default(shared) private(i,j)
         do j = 1,2*nglv
            do i = 1,nglv
               if (tinit(i,j) > 0) then
                  icestarxy(i,j) = icexy(i,j,exact_cur)
               elseif (icexy(i,j,exact_cur) > (abs(tinit(i,j))*rhow/rhoi)) then
                  icestarxy(i,j) = icexy(i,j,exact_cur)
               else
                  icestarxy(i,j) = 0
               endif
            enddo
         enddo
         if (opt_icebatch) then
            !$omp parallel do schedule(static) default(shared) private(j)
            do j = 1,2*nglv
               pp_gbuf(:,j,1) = icestarxy(:,j)
            enddo
            !$omp end parallel do
         else
            call opt_ana(icestarxy,icestarlm)
         endif
      else
         if (opt_icebatch) then
            !$omp parallel do schedule(static) default(shared) private(j)
            do j = 1,2*nglv
               pp_gbuf(:,j,1) = icexy(:,j,exact_cur)
            enddo
            !$omp end parallel do
         else
            call opt_ana(icexy(:,:,exact_cur),icestarlm)
         endif
      endif

      if (iter == 1) then
         if (checkmarine) then
            do j = 1,2*nglv
               do i = 1,nglv
                  if (tinit(i,j) > 0) then
                     icestarxy(i,j) = icexy(i,j,exact_prev)
                  elseif (icexy(i,j,exact_prev) > (abs(tinit(i,j))*rhow/rhoi)) then
                     icestarxy(i,j) = icexy(i,j,exact_prev)
                  else
                     icestarxy(i,j) = 0
                  endif
               enddo
            enddo
            call opt_ana(icestarxy,oldicestarlm)
         else
            call opt_ana(icexy(:,:,exact_prev),oldicestarlm)
         endif
      else
         oldicestarlm(:,:) = (0.0,0.0)
         do l = 0,norder
            do m = 0,l
               exact_idx = lm_index(l,m)
               oldicestarlm(l,m) = exact_state%prev_ice(exact_idx)
            enddo
         enddo
      endif

      ! the ice-increment spectra (dicestar, deltaicestar) are
      ! formed after the ocean-part analyses below, which is where the
      ! batched grounded-ice spectrum becomes available.  Nothing in the
      ! ocean grids depends on them, so the order change is arithmetic-free.

      !========================================================================================
      !                            BEGIN OCEAN PART
      !========================================================================================

      ! Calculate beta (STEP 3) (eq. 44)
      ! Calculate current beta based on iceload at the current timestep
      !$omp parallel do schedule(static) default(shared) private(i,j)
      do j = 1, 2 * nglv
         do i = 1, nglv
             if (icexy(i, j, exact_cur) < epsilon(0.0)) then
                beta(i,j) = 1
             else
                beta(i,j) = 0
             endif
         enddo
      enddo

      ! SLEVOLVE_REUSE=warm: converged ocean mask of step n from the
      ! previous outer iteration replaces the previous step's mask as first guess
      warm_have = .false.
      if (warm_mode > 0) then
         call warm_read(itersl, iter)
         if (warm_have) then
            !$omp parallel do schedule(static) default(shared) private(j)
            do j = 1,2*nglv
               cxy(:,j) = real(warm_mask(:,j))
            enddo
            !$omp end parallel do
         endif
      endif

      ! calculate initial (initial within the time window) cstar
      !$omp parallel do schedule(static) default(shared) private(j)
      do j = 1,2*nglv
         cstar0(:,j)  = cxy0(:,j) * beta0(:,j)
         cstarxy(:,j) = cxy(:,j) * beta(:,j)  ! First guess to the O.F using converged O.F the preivous timestep
      enddo
      !$omp end parallel do

      ! Calculate the first guess to topography correction
      ! cstar and tO are independent, so one sweep of the Legendre
      ! table serves both (SLEVOLVE_BATCH=0 restores two separate sweeps).
      ! with `rt` the tO spectrum is never needed (only the grid);
      ! with `icebatch` the grounded-ice field occupies slot 1.
      !$omp parallel do schedule(static) default(shared) private(j)
      do j = 1,2*nglv
         tOxy(:,j) = tinit(:,j) * (cstarxy(:,j) - cstar0(:,j)) ! (eq. 70)
      enddo
      !$omp end parallel do
      if (opt_icebatch) then
         opt_nf = 2
         if (.not. opt_rt) opt_nf = 3
         !$omp parallel do schedule(static) default(shared) private(j)
         do j = 1,2*nglv
            pp_gbuf(:,j,2) = cstarxy(:,j)
            if (.not. opt_rt) pp_gbuf(:,j,3) = tOxy(:,j)
         enddo
         !$omp end parallel do
         call opt_ana_multi(opt_nf)
         !$omp parallel do schedule(static) default(shared) private(m)
         do m = 0,norder
            icestarlm(:,m) = pp_ubuf(:,m,1)
            cstarlm(:,m)   = pp_ubuf(:,m,2)
            if (.not. opt_rt) tOlm(:,m) = pp_ubuf(:,m,3)
         enddo
         !$omp end parallel do
      else if (opt_rt) then
         call opt_ana(cstarxy,cstarlm)
      else if (pp_batch_opt) then
         !$omp parallel do schedule(static) default(shared) private(j)
         do j = 1,2*nglv
            pp_gbuf(:,j,1) = cstarxy(:,j)
            pp_gbuf(:,j,2) = tOxy(:,j)
         enddo
         !$omp end parallel do
         call opt_ana_multi(2)
         cstarlm(:,:) = pp_ubuf(:,:,1)
         tOlm(:,:)    = pp_ubuf(:,:,2)
      else
         call opt_ana(cstarxy,cstarlm) ! Decompose the current cstar
         call opt_ana(tOxy, tOlm) ! Decompose the topo correction term
      endif

      ! ice increments (moved here from the ice part, see above)
      !$omp parallel do schedule(static) default(shared) private(m)
      do m = 0,norder
         dicestarlm(:,m) = icestarlm(:,m) - oldicestarlm(:,m)
         deltaicestarlm(:,m) = dicestarlm(:,m)
      enddo
      !$omp end parallel do
      !$omp parallel do schedule(static) default(shared) private(l,m,exact_idx)
      do l = 0,norder
         do m = 0,l
            exact_idx = lm_index(l,m)
            deltaicestarlm(l,m) = exact_state%p_ice(exact_idx) + dicestarlm(l,m)
         enddo
      enddo
      !$omp end parallel do
      !$omp parallel do schedule(static) default(shared) private(m)
      do m = 0,norder
         dicestar(:,m,exact_cur) = dicestarlm(:,m)
         deltaicestar(:,m,exact_cur) = deltaicestarlm(:,m)
      enddo
      !$omp end parallel do


      ! Initial eustatic guess for the current increment.  The cumulative
      ! previous ocean load is already in deltaS(:,:,exact_prev).
      dS(:,:,exact_cur) = (cstarlm(:,:) / cstarlm(0,0)) * &
                          ((-rhoi / rhow) * dicestar(0,0,exact_cur))
      load_inc_lm(:,:) = dS(:,:,exact_cur)

      ! extrap: linear-in-time extrapolation of the converged
      ! increment; the cumulative current load that the first rotation
      ! evaluation reads is made consistent with it.  For uniform steps this
      ! is 2 dS(n-1) - dS(n-2); for variable steps the rate dS/dt is
      ! extrapolated linearly between step midpoints and scaled by dt(n).
      if (inner_mode > 0) then
         if (.not. allocated(extrap_dsm1)) then
            allocate(extrap_dsm1(0:norder,0:norder), extrap_dsm2(0:norder,0:norder))
            extrap_dsm1 = (0.0,0.0); extrap_dsm2 = (0.0,0.0)
         endif
         if (inner_nhist >= 2) then
            ! paper Sect. 2.5: dS0_j = c1 dS_{j-1} + c2 dS_{j-2} (modal_state_mod::extrap_coeffs)
            call extrap_coeffs(dtr, extrap_dtm1, extrap_dtm2, extrap_c1, extrap_c2)
            !$omp parallel do schedule(static) default(shared) private(m)
            do m = 0,norder
               load_inc_lm(:,m) = extrap_c1 * extrap_dsm1(:,m) + extrap_c2 * extrap_dsm2(:,m)
               deltaS(:,m,exact_cur) = deltaS(:,m,exact_prev) + load_inc_lm(:,m)
            enddo
            !$omp end parallel do
         endif
         inner_nh = 0; inner_lastmixed = .false.; inner_plain_rest = .false.; inner_nrefreeze = 0
      endif
      if (warm_have) then
         ! warm start: stored converged increment replaces the eustatic/extrap guess
         !$omp parallel do schedule(static) default(shared) private(m)
         do m = 0,norder
            load_inc_lm(:,m) = warm_ds(:,m)
            deltaS(:,m,exact_cur) = deltaS(:,m,exact_prev) + load_inc_lm(:,m)
         enddo
         !$omp end parallel do
         warm_nwarm = warm_nwarm + 1
      endif


      !========================================================================================
      !                   END OF THE OCEAN PART
      !========================================================================================


      !----------------------------------------------------------------------------------------
      !>>>>>>>>>>>>>>>>>>> Start of inner loop <<<<<<<<<<<<<<<<<<<
      !----------------------------------------------------------------------------------------
      ninner = 1
      pp_sabs_cur = sum(abs(load_inc_lm))   ! old_load_inc_lm == load_inc_lm at the head of iteration 1
      do ! Inner loop

         !-----\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    Rotation    \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/----!

         if (tpw) then ! If incorporating rotation

             ! MOI perturbations (Matsuyama et al, 2006) (only need (1,3),(2,3), and (3,3))

            il(3,3) = (-8 * pi * radius**4 / (3 * sqrt(5.0))) * real(rhow * deltaS(2,0,exact_cur) + rhoi * deltaicestar(2,0,exact_cur))
            il(1,3) = (8 * pi * radius**4 / (sqrt(30.0))) * real(rhow * deltaS(2,1,exact_cur) + rhoi * deltaicestar(2,1,exact_cur))
            il(2,3) = (-8 * pi * radius**4 / (sqrt(30.0))) &
                      & * real(-1.0 * ii * (rhow * deltaS(2,1,exact_cur) + rhoi * deltaicestar(2,1,exact_cur)))

            dil(:,:,exact_cur) = il(:,:) - oldil(:,:)

            ! Calculate m (rotation vector perturbation)
            !  - See Mitrovica, Wahr, Matsuyama, and Paulson (2005) eq. 2
            sum_il(:,:) = exact_sum_il(:,:)
            sum_m(:) = exact_sum_m(:)

            do i = 1,2
               mm(i) = (1 / (1 - kTE(2) / kf)) * &
                     & ( (1 / (moiC - moiA)) * ((1 + kE(2)) * il(i,3) + sum_il(i,3)) + (1 / kf) * sum_m(i) )
            enddo

            mm(3) = (1.0 / moiC) * ( (1 + kE(2)) * il(3,3) + sum_il(3,3) )

            dm(:,exact_cur) = mm(:) - oldm(:)

            ! Calculate lambda (rotational driving) from m (Milne and Mitrovica 1998)
            lambda(0,0) = (radius**2 * omega**2 / 3.0) * ((mm(1)**2 + mm(2)**2 + mm(3)**2) + 2.0 * mm(3))
            lambda(2,0) = (radius**2 * omega**2 / (6.0 * sqrt(5.0))) &
                        & * (mm(1)**2 + mm(2)**2 - 2.0 * mm(3)**2 - 4.0 * mm(3))
            lambda(2,1) = (radius**2 * omega**2 / sqrt(30.0)) * ((mm(1) - ii * mm(2)) * (1.0 + mm(3)))
            lambda(2,2) = (radius**2 * omega**2 / sqrt(5.0 * 24.0)) * (mm(2)**2 - mm(1)**2 + 2.0 * ii * mm(1) * mm(2))

            dlambda(:,:,exact_cur) = lambda(:,:) - deltalambda(:,:,exact_prev)

            ! Calculate effect on sea level (Love numbers) (Kendall)
            dsl_rot(:,:) = (0.0,0.0)
            ekhTE = 1 + kTE(2) - hTE(2)
            do m = 0,2
               viscoustt = exact_tsl(m)
               dsl_rot(2,m) = (ekhTE * (deltalambda(2,m,exact_prev) + dlambda(2,m,exact_cur)) + viscoustt) / gacc ! (eq. B28/B25)
            enddo

            if (calcRG) then ! For R calculations
               rr_rot(:,:) = (0.0,0.0)
               do m = 0,2
                  viscousttrr = exact_trr(m)
                  rr_rot(2,m) = (hTE(2) * (deltalambda(2,m,exact_prev) + dlambda(2,m,exact_cur)) + viscousttrr) / gacc
               enddo
            endif

         else ! If not incorporating rotation
            dsl_rot(:,:) = (0.0,0.0)
            rr_rot(:,:) = (0.0,0.0)
         endif ! TPW

         !-----/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\    End Rotation    /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\----!

         ! Old-load hereditary response was advanced once before the inner loop.
         ! `viscous` is intentionally frozen while the current ocean load converges.

         ! Calculate change in sea level from change in load (STEP 6)
         sl_inc_lm(:,:) = (0.0,0.0) ! No change on first timestep
         !$omp parallel do schedule(static) default(shared) private(l,m,ttl,ekhl,exact_idx,viscoustt)
         do l = 1,norder
            ttl = (4.0 * pi * (radius**3)) / (mass * (2 * real(l) + 1.0)) ! (eq. B16)
            ekhl = 1 + ke(l) - he(l) ! (eq. B13)
            do m = 0,l
               ! Total  = _____________________________elastic_____________________________ + _viscous_
               sl_inc_lm(l,m) = ttl * ekhl * (rhoi * deltaicestar(l,m,exact_cur) + rhow * deltaS(l,m,exact_prev) + rhow * load_inc_lm(l,m)) &
                             + ttl * viscous(l,m) ! (eq. B18)
            enddo
         enddo

         ! Add rotational effects (calculated above)
         sl_inc_lm(0:2,0:2) = sl_inc_lm(0:2,0:2) + dsl_rot(0:2,0:2)

         ! Convert dSL (total spatially heterogeneous change since time0) to spatial domain
         call opt_syn(dslxy, sl_inc_lm)

         if (opt_merge) then
            !--------------------------------------------------------------
            ! merge: the conservation term needs only the (0,0)
            ! coefficient of rO (or of rO - tO), which is the Gauss-Legendre
            ! quadrature of the grid mean -- so conserv, deltaslxy and the
            ! next ocean-function guess are all available BEFORE the rO
            ! analysis, and the rO analysis and the re-guess analysis
            ! (cstar [, tO]) share one sweep of the Legendre table.  If this
            ! iteration turns out to be the converged one, the re-guess
            ! spectra are simply discarded; the re-guess GRIDS (topoxy, cxy)
            ! are exactly what the post-loop block recomputes anyway.
            !--------------------------------------------------------------
            ! r0 (STEP 7) straight into batch slot 1, with its row sums
            call opt_rowsum_ro(dslxy, cstarxy, tOxy, pp_gbuf(:,:,1), opt_rs(:,1), &
                               spheredat%nlons, spheredat%nlats, opt_rt)
            opt_q00 = opt_quad(opt_rs(:,1))
            ! Compute conservation term (STEP 8) (eq. 78)
            if (opt_rt) then
               conserv = (1 / cstarlm(0,0)) * (-1.0 * (rhoi / rhow) * deltaicestar(0,0,exact_cur) - opt_q00)
            else
               conserv = (1 / cstarlm(0,0)) * (-1.0 * (rhoi / rhow) * deltaicestar(0,0,exact_cur) - opt_q00 + tOlm(0,0))
            endif
            ! total sea-level change on the grid: dsl + uniform conservation term
            call opt_addconst(deltaslxy, dslxy, conserv, spheredat%nlons, spheredat%nlats)
            ! anderson: from inner iteration 2 on the ocean function
            ! (cxy, cstarxy, tOxy, cstarlm) is frozen.  Only the number of cells
            ! whose ocean/land flag WOULD change is counted; if it is zero the
            ! frozen mask is the plain-iteration mask and nothing is lost.
            inner_frozen_it = (inner_mode == 2 .and. ninner >= 2 .and. .not. inner_plain_rest)
            inner_did_reguess = .not. inner_frozen_it
            inner_ndiff = 0
            if (inner_frozen_it) then
               inner_ndiff = 0
               !$omp parallel do schedule(static) default(shared) private(i,j) reduction(+:inner_ndiff)
               do j = 1,2*nglv
                  do i = 1,nglv
                     if ((tinit(i,j) - deltaslxy(i,j) >= 0.0) .neqv. (cxy(i,j) == 0.0)) inner_ndiff = inner_ndiff + 1
                  enddo
               enddo
               !$omp end parallel do
               inner_cnt(4) = inner_cnt(4) + 1
               opt_nf = 1                                   ! rO (- tO) only
               inner_cnt(5) = inner_cnt(5) + 1
            else
            ! speculative re-guess of topography / ocean function (eq. 39, 65, 70)
            !$omp parallel do schedule(static) default(shared) private(i,j)
            do j = 1,2*nglv
               do i = 1,nglv
                  topoxy(i,j) = tinit(i,j) - deltaslxy(i,j)
                  if (topoxy(i,j) >= 0.0) then
                     cxy(i,j) = 0.0
                  else
                     cxy(i,j) = 1.0
                  endif
                  cstarxy(i,j) = cxy(i,j) * beta(i,j)
                  tOxy(i,j) = tinit(i,j) * (cstarxy(i,j) - cstar0(i,j))
                  pp_gbuf(i,j,2) = cstarxy(i,j)
               enddo
               if (.not. opt_rt) pp_gbuf(:,j,3) = tOxy(:,j)
            enddo
            !$omp end parallel do
            opt_nf = 2
            if (.not. opt_rt) opt_nf = 3
            endif
            call opt_ana_multi(opt_nf)

            ! Compute change in ocean load again (STEP 9) (eq. 73)
            !$omp parallel do schedule(static) default(shared) private(m)
            do m = 0,norder
               old_load_inc_lm(:,m) = load_inc_lm(:,m) ! Save previous-iterate load to check convergence
            enddo
            !$omp end parallel do
            if (opt_rt) then
               !$omp parallel do schedule(static) default(shared) private(m)
               do m = 0,norder
                  load_inc_lm(:,m) = -1.0 * deltaS(:,m,exact_prev) + pp_ubuf(:,m,1) + conserv * cstarlm(:,m)
                  deltaS(:,m,exact_cur) = deltaS(:,m,exact_prev) + load_inc_lm(:,m)
               enddo
               !$omp end parallel do
            else
               !$omp parallel do schedule(static) default(shared) private(m)
               do m = 0,norder
                  load_inc_lm(:,m) = -1.0 * deltaS(:,m,exact_prev) + pp_ubuf(:,m,1) + conserv * cstarlm(:,m) - tOlm(:,m)
                  deltaS(:,m,exact_cur) = deltaS(:,m,exact_prev) + load_inc_lm(:,m)
               enddo
               !$omp end parallel do
            endif
         else
            ! Compute r0 (STEP 7)
            ! rt: rO and tO enter the solver only as (rO - tO), in
            ! conserv (degree 0) and in load_inc_lm (all degrees), and both are
            ! products with the SAME cstar -- one analysis of the difference
            ! replaces the rO analysis here and the tO analysis of every
            ! re-guess.
            if (opt_rt) then
               !$omp parallel do schedule(static) default(shared) private(j)
               do j = 1,2*nglv
                  rOxy(:,j) = dslxy(:,j) * cstarxy(:,j) - tOxy(:,j)
               enddo
               !$omp end parallel do
            else
               !$omp parallel do schedule(static) default(shared) private(j)
               do j = 1,2*nglv
                  rOxy(:,j) = dslxy(:,j) * cstarxy(:,j) ! (eq. 68)
               enddo
               !$omp end parallel do
            endif
            call opt_ana(rOxy, rOlm) ! Decompose r0 (or r0 - t0)

            ! Compute conservation term (STEP 8)
            if (opt_rt) then
               conserv = (1 / cstarlm(0,0)) * (-1.0 * (rhoi / rhow) * deltaicestar(0,0,exact_cur) - rOlm(0,0)) ! (eq. 78)
            else
               conserv = (1 / cstarlm(0,0)) * (-1.0 * (rhoi / rhow) * deltaicestar(0,0,exact_cur) - rOlm(0,0) + tOlm(0,0)) ! (eq. 78)
            endif

            ! Compute change in ocean load again (STEP 9)
            !$omp parallel do schedule(static) default(shared) private(m)
            do m = 0,norder
               old_load_inc_lm(:,m) = load_inc_lm(:,m) ! Save previous-iterate load to check convergence
            enddo
            !$omp end parallel do

            if (opt_rt) then
               !$omp parallel do schedule(static) default(shared) private(m)
               do m = 0,norder
                  load_inc_lm(:,m) = -1.0 * deltaS(:,m,exact_prev) + rOlm(:,m) + conserv * cstarlm(:,m) ! (eq. 73)
                  deltaS(:,m,exact_cur) = deltaS(:,m,exact_prev) + load_inc_lm(:,m)
               enddo
               !$omp end parallel do
            else
               !$omp parallel do schedule(static) default(shared) private(m)
               do m = 0,norder
                  load_inc_lm(:,m) = -1.0 * deltaS(:,m,exact_prev) + rOlm(:,m) + conserv * cstarlm(:,m) - tOlm(:,m) ! (eq. 73)
                  deltaS(:,m,exact_cur) = deltaS(:,m,exact_prev) + load_inc_lm(:,m)! Save total load changes, used in Love number calculation
               enddo
               !$omp end parallel do
            endif

            ! Update the total sea-level change up to the current time step
            if (opt_dsl) then
               ! dsl: sl_tot_lm differs from sl_inc_lm only by conserv in
               ! (0,0), and Y00 is a constant, so synth(sl_tot_lm) =
               ! synth(sl_inc_lm) + conserv*synth(e00) = dslxy + conserv*K00.
               call opt_addconst(deltaslxy, dslxy, conserv, spheredat%nlons, spheredat%nlats)
            else
               !$omp parallel do schedule(static) default(shared) private(m)
               do m = 0,norder
                  sl_tot_lm(:,m) = sl_inc_lm(:,m) ! Spatially heterogeneous component
               enddo
               !$omp end parallel do
               sl_tot_lm(0,0) = sl_tot_lm(0,0) + conserv ! Add uniform conservation term to (0,0)
               call opt_syn(deltaslxy, sl_tot_lm) ! Synthesize deltasl
            endif
         endif

         ! Calculate convergence criterion for inner loop
         ! the baseline spells out sum(abs(load_inc_lm)) and sum(abs(old_load_inc_lm))
         ! up to seven times in this one construct; each is a reduction over
         ! 263169 complex numbers with a libm hypot per element (2.7 % of the
         ! round-1 profile).  The two values are computed once here, and
         ! sum(abs(old_load_inc_lm)) is carried from the previous iteration because
         ! old_load_inc_lm was assigned `= load_inc_lm` (identical) just above.  Identical
         ! arithmetic, identical result.
         pp_sabs_old = pp_sabs_cur
         pp_sabs_cur = sum(abs(load_inc_lm))
         if ( abs(pp_sabs_cur - pp_sabs_old) < epsilon(0.0) .and. abs(pp_sabs_old) < epsilon(0.0)) then
             xi = 0 ! Otherwise xi = 0 / 0 = NaN at the first loop of the first timestep.
         elseif (abs(pp_sabs_old) < epsilon(0.0)) then
             xi = abs( (pp_sabs_cur - pp_sabs_old) / (epsilon(0.0) * 10) ) ! Avoid dividing by 0
         else
             xi = abs( (pp_sabs_cur - pp_sabs_old) / pp_sabs_old ) ! (eq. 83)
         endif

         ! If the ocean loading guess has not been converged,
         if (xi > epsilon1) then
            if (opt_merge .and. inner_mode == 2 .and. .not. inner_did_reguess) then
               continue                     ! ocean function frozen
            else if (opt_merge) then
               ! the re-guess was already analysed in this iteration's sweep
               !$omp parallel do schedule(static) default(shared) private(m)
               do m = 0,norder
                  cstarlm(:,m) = pp_ubuf(:,m,2)
                  if (.not. opt_rt) tOlm(:,m) = pp_ubuf(:,m,3)
               enddo
               !$omp end parallel do
            else
            ! new guess to the topography correction
            !$omp parallel do schedule(static) default(shared) private(j)
            do j = 1,2*nglv
               topoxy(:,j) = tinit(:,j) - deltaslxy(:,j) ! (eq. 39)
            enddo
            !$omp end parallel do

            ! new ocean function
            !$omp parallel do schedule(static) default(shared) private(i,j)
            do j = 1,2*nglv
               do i = 1,nglv
                  if (topoxy(i,j) >= 0.0) then
                     cxy(i,j) = 0.0
                  else
                     cxy(i,j) = 1.0
                  endif
               enddo
            enddo

            !new guess to ocean*beta function
            ! new guess to the topo correction
            !$omp parallel do schedule(static) default(shared) private(j)
            do j = 1,2*nglv
               cstarxy(:,j) = cxy(:,j) * beta(:,j) ! (eq. 65)
               tOxy(:,j) = tinit(:,j) * (cstarxy(:,j) - cstar0(:,j)) ! (eq. 70)
            enddo
            !$omp end parallel do
            if (opt_rt) then
               call opt_ana(cstarxy, cstarlm) ! tO is only needed on the grid
            else if (pp_batch_opt) then
               !$omp parallel do schedule(static) default(shared) private(j)
               do j = 1,2*nglv
                  pp_gbuf(:,j,1) = cstarxy(:,j)
                  pp_gbuf(:,j,2) = tOxy(:,j)
               enddo
               !$omp end parallel do
               call opt_ana_multi(2)
               cstarlm(:,:) = pp_ubuf(:,:,1)
               tOlm(:,:)    = pp_ubuf(:,:,2)
            else
               call opt_ana(cstarxy, cstarlm) ! Decompose cstar
               call opt_ana(tOxy, tOlm) ! Decompose tO
            endif
            endif
         endif

         inner_conv = (xi <= epsilon1)
         if (inner_log) then
            inner_ndiff2 = 0
            !$omp parallel do schedule(static) default(shared) private(i,j) reduction(+:inner_ndiff2)
            do j = 1,2*nglv
               do i = 1,nglv
                  if ((tinit(i,j) - deltaslxy(i,j) >= 0.0) .neqv. (cxy(i,j) == 0.0)) inner_ndiff2 = inner_ndiff2 + 1
               enddo
            enddo
            !$omp end parallel do
            write(unit_num,'(A,I6,I4,ES12.4,ES12.4,ES12.4,I8)') ' innerlog ', nfiles, ninner, xi, &
                 sqrt(sum(abs(load_inc_lm-old_load_inc_lm)**2)), sqrt(sum(abs(load_inc_lm)**2)), inner_ndiff2
         endif
         if (inner_mode == 2) then
            if (inner_refresh .and. inner_conv .and. .not. inner_did_reguess .and. inner_ndiff > 0) then
               ! converged with a stale (frozen) ocean function: apply the
               ! re-guess now (eq. 39, 65, 70, one extra analysis), restart the
               ! Anderson history, and keep iterating.  After 3 such updates
               ! the rest of the step is plain iteration.
               !$omp parallel do schedule(static) default(shared) private(i,j)
               do j = 1,2*nglv
                  do i = 1,nglv
                     topoxy(i,j) = tinit(i,j) - deltaslxy(i,j)
                     if (topoxy(i,j) >= 0.0) then
                        cxy(i,j) = 0.0
                     else
                        cxy(i,j) = 1.0
                     endif
                     cstarxy(i,j) = cxy(i,j) * beta(i,j)
                     tOxy(i,j) = tinit(i,j) * (cstarxy(i,j) - cstar0(i,j))
                  enddo
               enddo
               !$omp end parallel do
               call opt_ana(cstarxy, cstarlm)
               inner_cnt(3) = inner_cnt(3) + 1
               inner_cnt(5) = inner_cnt(5) - 1
               inner_nh = 0; inner_lastmixed = .false.
               inner_nrefreeze = inner_nrefreeze + 1
               if (inner_nrefreeze >= 3) inner_plain_rest = .true.
               inner_conv = .false.
            else if (.not. inner_conv .and. ninner >= 2 .and. .not. inner_plain_rest) then
               call inner_anderson()
               pp_sabs_cur = sum(abs(load_inc_lm))
            endif
         endif
         if (inner_conv) then ! If converged
            exit
         elseif (ninner == 9999) then ! If no convergence after a huge number of iterations
            write(unit_num,*)
            write(unit_num,'(A,I5,A)') 'WARNING: The inner loop failed to converge after the limit of ', ninner, ' iterations.'
            write(unit_num,'(A,ES15.3,A)') '         The variable xi finished with a value of ', xi, ', resulting from '
            write(unit_num,'(A,ES15.3)') '         sum(abs(old_load_inc_lm)) = ', sum(abs(old_load_inc_lm))
            write(unit_num,'(A,ES15.3)') '            sum(abs(load_inc_lm)) = ', sum(abs(load_inc_lm))
            write(unit_num,*)
            write(unit_num,'(A)') '!!!---- Program sl_model will now be terminated. ----!!!'
            call abort ! Terminate program
            !exit ! DEBUG line: Continue inner loop despite non-convergence. Normal operation: Enable the 2 lines above.
         endif
         ninner = ninner + 1

      enddo ! End inner loop
      !-----------------------------------------------------------
      !<<<<<<<<<<<<<<<<<<<< End of inner loop >>>>>>>>>>>>>>>>>>>>
      !-----------------------------------------------------------
      write(unit_num,'(A,I4,A)') '  ', ninner, ' inner-loop iterations'

      if (inner_mode > 0) then            ! extrap history
         !$omp parallel do schedule(static) default(shared) private(m)
         do m = 0,norder
            extrap_dsm2(:,m) = extrap_dsm1(:,m)
            extrap_dsm1(:,m) = load_inc_lm(:,m)
         enddo
         !$omp end parallel do
         extrap_dtm2 = extrap_dtm1; extrap_dtm1 = dtr
         inner_nhist = inner_nhist + 1
         inner_cnt(6) = inner_cnt(6) + 1
      endif

      ! lazyrg: R (and G = deltasl + R) are output-only fields.  When
      ! neither is written at this step (snapshot mode), their synthesis is
      ! skipped.  Outside snapshot mode pp_snap_want is always .true.
      write(opt_stepstr,'(I6)') TIMEWINDOW(nfiles)
      opt_stepstr = trim(adjustl(opt_stepstr))
      opt_rg_needed = calcRG
      if (calcRG .and. opt_lazyrg) then
         opt_rg_needed = pp_snap_want('R', trim(opt_stepstr), opt_st) .or. &
                         pp_snap_want('G', trim(opt_stepstr), opt_st)
      endif

      if (opt_rg_needed) then ! For R calculations
         if (nmelt == 0) then
            rrlm(:,:) = (0.0,0.0)! No change on first timestep
         else
            !$omp parallel do schedule(static) default(shared) private(l,m,ttl,exact_idx,viscousrr)
            do l = 1,norder
               ttl = (4.0 * pi * (radius**3)) / (mass * (2.0 * real(l) + 1.0)) ! (eq. B16)
               do m = 0,l
                  ! Exact old-load hereditary radial response.
                  exact_idx = lm_index(l,m)
                  viscousrr = exact_radial(exact_idx)
                  rrlm(l,m) = ttl * he(l) * (rhoi * deltaicestar(l,m,exact_cur) + rhow * deltaS(l,m,exact_prev) + rhow * load_inc_lm(l,m)) &
                           & + ttl * viscousrr
               enddo
            enddo
         endif
         rrlm(0:2,0:2) = rrlm(0:2,0:2) + rr_rot(0:2,0:2)
         call opt_syn(rrxy, rrlm)
         !$omp parallel do schedule(static) default(shared) private(j)
         do j = 1,2*nglv
            rr(:,j,exact_cur) = rrxy(:,j)
         enddo
         !$omp end parallel do
      endif


      ! Optional horizontal-displacement read-out (SLEVOLVE_HORIZ).  Same
      ! states and the same instantaneous load as rrlm above; weights l
      ! instead of h.  Nothing here feeds back into the solver.
      pp_hz_needed = .false.
      if (pp_horiz > 0) pp_hz_needed = pp_snap_want('horiz', trim(opt_stepstr), opt_st) .or. &
                                       (pp_horiz > 1 .and. pp_snap_want('Psi', trim(opt_stepstr), opt_st))
      if (pp_hz_needed) then
         if (.not.allocated(exact_horiz)) allocate(exact_horiz(exact_state%ncoef))
         if (.not.allocated(hzlm)) allocate(hzlm(0:norder,0:norder))
         call exact_state%read_horizontal(exact_spec, exact_horiz, exact_thl)
         hzlm(:,:) = (0.0,0.0)
         if (nmelt /= 0) then
            !$omp parallel do schedule(static) default(shared) private(l,m,ttl,exact_idx)
            do l = 1,norder
               ttl = (4.0 * pi * (radius**3)) / (mass * (2.0 * real(l) + 1.0)) ! (eq. B16)
               do m = 0,l
                  exact_idx = lm_index(l,m)
                  hzlm(l,m) = ttl * (elast(2,l)/real(l)) * (rhoi * deltaicestar(l,m,exact_cur) &
                           & + rhow * deltaS(l,m,exact_prev) + rhow * load_inc_lm(l,m)) &
                           & + ttl * exact_horiz(exact_idx)
               enddo
            enddo
            !$omp end parallel do
            if (tpw .and. norder >= 2) then
               do m = 0,2
                  hzlm(2,m) = hzlm(2,m) + ((telast(2,2)/2.0) * (deltalambda(2,m,exact_prev) &
                           & + dlambda(2,m,exact_cur)) + exact_thl(m)) / gacc
               enddo
            endif
         endif
         call pp_write_cplx(hzlm, 'horiz', outputfolder, suffix=trim(opt_stepstr), fmt='(ES16.9E2)')
         if (pp_horiz > 1) then
            if (pp_snap_want('Psi', trim(opt_stepstr), opt_st)) then
               if (.not.allocated(hzxy)) allocate(hzxy(nglv,2*nglv))
               call opt_syn(hzxy, hzlm)
               call pp_write_grid(hzxy, 'Psi', outputfolder, suffix=trim(opt_stepstr))
            endif
         endif
      endif

      ! Commit the converged current-step increments only now.  The current
      ! load has zero age at this instant, so it did not contribute to the
      ! hereditary response used above; it becomes memory for future steps.
      exact_current_load = (0.0_dp,0.0_dp)
      !$omp parallel do schedule(static) default(shared) private(l,m,exact_idx)
      do l = 1,norder
         do m = 0,l
            exact_idx = lm_index(l,m)
            exact_current_load(exact_idx) = rhoi*dicestar(l,m,exact_cur) + rhow*load_inc_lm(l,m)
         enddo
      enddo
      !$omp end parallel do
      if (tpw) then
         call exact_state%commit(exact_current_load,dil(:,:,exact_cur),dm(:,exact_cur),dlambda(2,0:2,exact_cur))
      else
         call exact_state%commit(exact_current_load)
      endif
      ! Persist only O(1)-in-time instantaneous/cumulative auxiliaries.
      !$omp parallel do schedule(static) default(shared) private(l,m,exact_idx)
      do l = 0,norder
         do m = 0,l
            exact_idx = lm_index(l,m)
            exact_state%p_ice(exact_idx) = deltaicestar(l,m,exact_cur)
            exact_state%p_ocean(exact_idx) = deltaS(l,m,exact_cur)
            exact_state%prev_ice(exact_idx) = icestarlm(l,m)
         enddo
      enddo
      if (iand(pp_compat,4) /= 0) then
         call exact_state%save_checkpoint(trim(exact_statefile))
      else if (pp_ckpt_every > 0) then
         if (mod(iter,pp_ckpt_every) == 0 .or. iter == L_sim/dt1) &
            call exact_state%save_checkpoint(trim(exact_statefile))
      endif

      ! Update topography fieldS, ocean functions (STEP 11)
      oldt0lm = t0lm ! Save old initial topography to check convergence

      ! Update the topography at the current time step
      !$omp parallel do schedule(static) default(shared) private(j)
      do j = 1,2*nglv
         topoxy(:,j) = tinit(:,j) - deltaslxy(:,j) ! (eq. 12)
      enddo
      !$omp end parallel do

      ! Update ocean function
      !$omp parallel do schedule(static) default(shared) private(i,j)
      do j = 1,2*nglv
         do i = 1,nglv
            if (topoxy(i,j) >= 0.0) then
               cxy(i,j) = 0.0
            else
               cxy(i,j) = 1.0
            endif
         enddo
      enddo

      if (warm_mode > 0) call warm_write(itersl, iter)

      ! Calculate global mean sea-level change where marine-based region is flooded
      ! Calculate global mean sea-level change where marine-based region is NOT flooded
      if (opt_diag) then
         ! diag: only the (0,0) coefficients of these four fields are
         ! ever used (two areas, two GMSL ratios), and no later step reads the
         ! spectra.  A (0,0) coefficient is the Gauss-Legendre quadrature of
         ! the grid mean, with exactly the SHT's normalisation (opt_g).
         call opt_rowsum_gmsl(deltaslxy, cxy, beta, opt_rs, spheredat%nlons, spheredat%nlats)
         opt_slc00  = opt_quad(opt_rs(:,1))
         opt_c00    = opt_quad(opt_rs(:,2))
         opt_slcb00 = opt_quad(opt_rs(:,3))
         opt_cst00  = opt_quad(opt_rs(:,4))
         gmslc_ocn = opt_slc00 / opt_c00
         gmslc_ocnBeta = opt_slcb00 / opt_cst00
      else
      ! these four analyses are mutually independent -> one sweep.
      !$omp parallel do schedule(static) default(shared) private(j)
      do j = 1,2*nglv
         slcxy_ocn(:,j) = deltaslxy(:,j) * cxy(:,j)
         slcxy_ocnBeta(:,j) = deltaslxy(:,j) * cxy(:,j) * beta(:,j)
         cstarxy(:,j) = cxy(:,j) * beta(:,j)
      enddo
      !$omp end parallel do
      if (pp_batch_opt) then
         !$omp parallel do schedule(static) default(shared) private(j)
         do j = 1,2*nglv
            pp_gbuf(:,j,1) = slcxy_ocn(:,j)
            pp_gbuf(:,j,2) = cxy(:,j)
            pp_gbuf(:,j,3) = slcxy_ocnBeta(:,j)
            pp_gbuf(:,j,4) = cstarxy(:,j)
         enddo
         !$omp end parallel do
         call opt_ana_multi(4)
         slclm_ocn(:,:)     = pp_ubuf(:,:,1)
         clm(:,:)           = pp_ubuf(:,:,2)
         slclm_ocnBeta(:,:) = pp_ubuf(:,:,3)
         cstarlm(:,:)       = pp_ubuf(:,:,4)
      else
         call opt_ana(slcxy_ocn, slclm_ocn)
         call opt_ana(cxy, clm)
         call opt_ana(slcxy_ocnBeta, slclm_ocnBeta)
         call opt_ana(cstarxy, cstarlm)
      endif
      gmslc_ocn = slclm_ocn(0,0) / clm(0,0)
      gmslc_ocnBeta = slclm_ocnBeta(0,0) / cstarlm(0,0)
      opt_c00 = real(clm(0,0))
      opt_cst00 = real(cstarlm(0,0))
      endif

      !=========================================================================================
      !                          OUTPUT
      !_________________________________________________________________________________________

      write(unit_num,'(A)') 'Writing output files...'
      j = TIMEWINDOW(nfiles)
      write(unit_num,*) 'FILENUMBER of new outputs : ',j
      write(numstr,'(I4)') j
      numstr = trim(adjustl(numstr))

      ! nmelt
      open(unit = 201, file = trim(outputfolder)//'nmelt', form = 'formatted', access ='sequential', &
      & status = 'old',position='append')
      write(201,'(I4)') nmelt
      close(201)

      ! number of iteration it takes for the inner convergence
      open(unit = 201, file = trim(outputfolder)//'numiter', form = 'formatted', access ='sequential', &
      & status = 'old', position='append')
      write(201,'(I5)') ninner
      close(201)

      ! converged rotation-related quantities
      if (tpw) then
         open(unit = 201, file = trim(outputfolder)//'TPW'//trim(numstr), form = 'formatted', &
         & access = 'sequential', status = 'replace')
         write(201,'(9ES19.8E2/,3ES19.8E2/,18ES19.8E2)') il(:,:), mm(:), lambda(:,:)
         close(201)
      endif

      ! topography at the current timestep
      call pp_write_grid(topoxy, 'tgrid', outputfolder, suffix=numstr)

      ! converged ocean function at the current timestep
      call pp_write_grid(cxy, 'ocean', outputfolder, suffix=numstr)

      ! converged beta function at the current timestpe
      call pp_write_grid(beta, 'beta', outputfolder, suffix=numstr)

      ! output the ocean areas
      open(unit = 201, file = trim(outputfolder)//'ocean_area', form = 'formatted', access ='sequential', &
      & status = 'old', position = 'append')
      write(201,'(ES14.4E2)') opt_c00*4*pi*radius**2
      close(201)

      open(unit = 201, file = trim(outputfolder)//'oceanBeta_area', form = 'formatted', access ='sequential', &
      & status = 'old', position = 'append')
      write(201,'(ES14.4E2)') opt_cst00*4*pi*radius**2
      close(201)

      ! output global mean sea level changes
      open(unit = 201, file = trim(outputfolder)//'gmslc_ocnArea', form = 'formatted', access ='sequential', &
      & status = 'old', position = 'append')
      write(201,'(ES14.4E2)') gmslc_ocn
      close(201)

      open(unit = 201, file = trim(outputfolder)//'gmslc_ocnBetaArea', form = 'formatted', access ='sequential', &
      & status = 'old', position = 'append')
      write(201,'(ES14.4E2)') gmslc_ocnBeta
      close(201)

      ! output converged total ocean loading changes
      call pp_write_cplx(deltaS(:,:,exact_cur), 'dS_converged', outputfolder, suffix=numstr, fmt='(ES14.7E2)')

      if (iceVolume) then
         ice_volume = icestarlm(0,0)*4*pi*radius**2 !multiply the (0,0) component of ice to the area of a sphere
         open(unit = 201, file = trim(outputfolder)//'ice_volume', form = 'formatted', access = 'sequential', &
         & status = 'old', position = 'append')
         write(201,'(ES14.4E2)') ice_volume
         close(201)
      endif

      current_time = iter * dt1     !time passed since the start of the simulation
      if (current_time == L_sim) then !if we are at the last time step of simulation
         write(iterstr,'(I2)') itersl
           iterstr = trim(adjustl(iterstr))
         write(unit_num,*) 'Last time step of the simulation! writing out files for next outer-iteration loop'

         ! write out the predicted present day topography into a file so it can be used in the next outer-iteration
         call pp_write_grid(topoxy, 'pred_pres_topo_', outputfolder, suffix=iterstr, force=outer_set)

         ! write out the initial topography of the simulation at the currect outer-loop into a file
         call pp_write_grid(tinit_0, 'tgrid0_', outputfolder, suffix=iterstr, force=outer_set)
         if (pp_mem_outer) then        ! hand the pass result to the next pass in memory
            pred_pres_topo(:,:) = topoxy(:,:)
            tinit_0_last(:,:) = tinit_0(:,:)
         endif
       endif

      if (opt_rg_needed) then
         call pp_write_grid(rr(:,:,exact_cur), 'R', outputfolder, suffix=numstr)

         ! Compute geoid displacement and write out
         gg(:,:,exact_cur) = deltaslxy(:,:)+rr(:,:,exact_cur)

         call pp_write_grid(gg(:,:,exact_cur), 'G', outputfolder, suffix=numstr)
       endif

      if (coupling) then
         ! topography change between the previous and the current timestep
         ! this is the information passed to the ice sheet model
         !call write_sl(topoxy_m1(:,:)-topoxy(:,:), 'bedrock', folder_coupled)
          slchange = topoxy_m1(:,:)-topoxy(:,:)
         !write out the current ice load as a new file
         call write_sl(icexy(:,:,exact_cur), icemodel_out, outputfolder_ice, suffix=numstr)
      endif !endif coupling

      call system_clock(countf) ! Total time
      call cpu_time(countf_cpu)

      ! Write out total compuatation time of sea level change over current timestep
      open(unit = 201, file = trim(outputfolder)//'elapsed_wall_time', form = 'formatted', access = 'sequential', &
      & status = 'old', position='append')
      write(201,'(ES14.4E2)') float(countf-counti)/float(countrate)
      close(201)

      open(unit = 201, file = trim(outputfolder)//'elapsed_cpu_time', form = 'formatted', access = 'sequential', &
      & status = 'old', position='append')
      write(201,'(ES14.4E2)') countf_cpu-counti_cpu
      close(201)


      write(unit_num,'(A,F7.2,A)') 'Done! Total time ', real(countf - counti) / real(countrate), ' seconds'
      write(unit_num,*) ''
      write(unit_num,*) ''

      if (Travel_total > 0 .and. Travel == Travel_total) then
         write(unit_num,*) ' THE TIMEWINDOW HAS REACHED THE END OF THE SIMULATION!'
         write(unit_num,*) ' GREAT JOB TW!'
      endif
      write(unit_num,*) ''
   end subroutine sl_solver


!----------------------------------------------------------------------------------------------------------------------!




   !=======================================================================================================================!
   ! transform-count reduction helpers
   !=======================================================================================================================!

   ! One-off setup: degree-0 quadrature row weights and the synthesis of a
   ! unit (0,0) coefficient, both measured/derived for the SHT actually built.
   subroutine opt_setup()
      integer :: nlo, nla, jj, ii2
      real, allocatable :: kg(:,:), tf(:,:)
      complex, allocatable :: e00(:,:)
      real :: q, a00, rsj
      integer :: nbad
      if (opt_ready) return
      opt_ready = .true.
      nlo = spheredat%nlons
      nla = spheredat%nlats
      if (mod(nlo,8) /= 0) stop 'opt_setup: nlons must be a multiple of 8'
      allocate(opt_g(nla), opt_k00(nla), opt_rs(nla,4))
      ! analysis(f)(0,0) = sum_j g_j * sum_i f(i,j) : legend's P_00 = sqrt(1/2),
      ! the Gauss weight w_j, the 1/nlon of the forward FFT, and the legacy
      ! 1/1.41421356 (spat2spec).  The BLAS backend applies exactly these
      ! factors (w_j*fsc*RSQRT2 in the interleave, P_00 in the DGEMM).
      do jj = 1,nla
         opt_g(jj) = (sqrt(0.5d0) * spheredat%weights(jj)) * SHT_RSQRT2 / real(nlo)
      enddo
      ! K00 = spec2spat(e00): measured through the backend itself.
      allocate(kg(nlo,nla), tf(nlo,nla), e00(0:norder,0:norder))
      e00 = (0.0,0.0)
      e00(0,0) = (1.0,0.0)
      call spec2spat(kg, e00, spheredat)
      nbad = 0
      do jj = 1,nla
         opt_k00(jj) = kg(1,jj)
         do ii2 = 1,nlo
            if (kg(ii2,jj) /= kg(1,jj)) nbad = nbad + 1
         enddo
      enddo
      write(unit_num,'(A,ES24.16,A,ES24.16,A,I8)') ' SLEvolve K00 = synth(e00): min ', minval(opt_k00), &
           '  max ', maxval(opt_k00), '  non-uniform along longitude: ', nbad
      ! self-test of the quadrature against a full analysis
      do jj = 1,nla
         do ii2 = 1,nlo
            tf(ii2,jj) = real(modulo(ii2*7919 + jj*104729, 1000)) / 1000.0 - 0.5 + &
                         2.0*spheredat%gaulats(jj)**3
         enddo
      enddo
      call spat2spec(tf, e00, spheredat)
      a00 = real(e00(0,0))
      q = 0.0
      do jj = 1,nla
         rsj = 0.0
         do ii2 = 1,nlo
            rsj = rsj + tf(ii2,jj)
         enddo
         q = q + opt_g(jj)*rsj
      enddo
      write(unit_num,'(A,ES24.16,A,ES24.16,A,ES10.2)') ' SLEvolve deg-0 quadrature self-test: analysis ', a00, &
           '  quadrature ', q, '  rel diff ', abs(q-a00)/abs(a00)
      deallocate(kg, tf, e00)
   end subroutine opt_setup

   ! counted transform wrappers (explicit-shape: sequence association as in spharmt)
   subroutine opt_ana(Z, U)
      real, intent(in) :: Z(nglv, 2*nglv)
      complex, intent(inout) :: U(0:norder, 0:norder)
      call spat2spec(Z, U, spheredat)
      opt_cnt(1) = opt_cnt(1) + 1
      opt_cnt(2) = opt_cnt(2) + 1
   end subroutine opt_ana

   subroutine opt_ana_multi(nf)
      integer, intent(in) :: nf
      call spat2spec_multi(pp_gbuf, pp_ubuf, spheredat, nf)
      opt_cnt(1) = opt_cnt(1) + nf
      opt_cnt(2) = opt_cnt(2) + 1
   end subroutine opt_ana_multi

   subroutine opt_syn(Z, U)
      real, intent(out) :: Z(nglv, 2*nglv)
      complex, intent(in) :: U(0:norder, 0:norder)
      call spec2spat(Z, U, spheredat)
      opt_cnt(3) = opt_cnt(3) + 1
      opt_cnt(4) = opt_cnt(4) + 1
   end subroutine opt_syn

   !=======================================================================================================================!
   ! Anderson mixing (type II, depth inner_aa_m) on the spectral ocean-load increment.
   ! On entry old_load_inc_lm = x_k (the iterate fed to this inner iteration) and load_inc_lm = g_k = G(x_k).
   ! F_k = g_k - x_k.  The ring holds (F_j, G_j) of the current step since the last reset;
   ! x_{k+1} = G_k - sum_a gamma_a (G_{a+1} - G_a), gamma = argmin |F_k - sum_a gamma_a (F_{a+1}-F_a)|.
   ! Safeguard: if the previous update was a mixed one and |F_k| > |F_{k-1}|, the ring is reset
   ! to (F_k, G_k) and the plain step x_{k+1} = G_k is taken.
   ! Inner products: Re sum conj(a) b over the (l,m) array (entries m > l are zero in all iterates).
   !=======================================================================================================================!
   subroutine inner_anderson()
      integer :: mm1, s, snew, a, b, n, info
      real(8) :: fn, d, amat(5,5), bvec(5), c(6), tr
      mm1 = inner_aa_m + 1
      if (.not. allocated(inner_hf)) then
         allocate(inner_hf(0:norder,0:norder,6), inner_hg(0:norder,0:norder,6))
      endif
      ! choose the ring slot: drop the oldest when full
      if (inner_nh == mm1) then
         snew = inner_slot(1)
         inner_slot(1:mm1-1) = inner_slot(2:mm1)
         inner_nh = inner_nh - 1
      else
         snew = 0
         do s = 1,6
            if (all(inner_slot(1:inner_nh) /= s)) then; snew = s; exit; endif
         enddo
      endif
      !$omp parallel do schedule(static) default(shared) private(m)
      do m = 0,norder
         inner_hf(:,m,snew) = load_inc_lm(:,m) - old_load_inc_lm(:,m)
         inner_hg(:,m,snew) = load_inc_lm(:,m)
      enddo
      !$omp end parallel do
      inner_nh = inner_nh + 1
      inner_slot(inner_nh) = snew
      do a = 1,inner_nh
         s = inner_slot(a)
         d = inner_dot(inner_hf(:,:,s), inner_hf(:,:,snew))
         inner_gram(s,snew) = d; inner_gram(snew,s) = d
      enddo
      fn = sqrt(inner_gram(snew,snew))
      if (inner_lastmixed .and. fn > inner_lastfn) then
         ! safeguard: the mixed step increased the residual -> plain step, restart
         inner_slot(1) = snew; inner_nh = 1
         inner_lastmixed = .false.; inner_lastfn = fn
         inner_cnt(2) = inner_cnt(2) + 1
         return                                     ! load_inc_lm = G_k (plain)
      endif
      inner_lastfn = fn
      n = inner_nh
      if (n < 2) then
         inner_lastmixed = .false.
         return
      endif
      ! normal equations for the n-1 differences
      do a = 1,n-1
         do b = 1,n-1
            amat(a,b) = inner_gram(inner_slot(a+1),inner_slot(b+1)) - inner_gram(inner_slot(a+1),inner_slot(b)) &
                      - inner_gram(inner_slot(a),inner_slot(b+1)) + inner_gram(inner_slot(a),inner_slot(b))
         enddo
         bvec(a) = inner_gram(inner_slot(a+1),inner_slot(n)) - inner_gram(inner_slot(a),inner_slot(n))
      enddo
      tr = 0.0d0
      do a = 1,n-1
         tr = tr + amat(a,a)
      enddo
      do a = 1,n-1
         amat(a,a) = amat(a,a) + 1.0d-12 * tr
      enddo
      call inner_solve(amat, bvec, n-1, info)
      if (info /= 0) then
         inner_slot(1) = snew; inner_nh = 1; inner_lastmixed = .false.
         return
      endif
      ! x_{k+1} = sum_j c_j G_j
      c = 0.0d0
      c(n) = 1.0d0
      do a = 1,n-1
         c(a+1) = c(a+1) - bvec(a)
         c(a)   = c(a)   + bvec(a)
      enddo
      !$omp parallel do schedule(static) default(shared) private(m,a)
      do m = 0,norder
         load_inc_lm(:,m) = (0.0,0.0)
         do a = 1,n
            load_inc_lm(:,m) = load_inc_lm(:,m) + c(a) * inner_hg(:,m,inner_slot(a))
         enddo
         deltaS(:,m,exact_cur) = deltaS(:,m,exact_prev) + load_inc_lm(:,m)
      enddo
      !$omp end parallel do
      inner_lastmixed = .true.
      inner_cnt(1) = inner_cnt(1) + 1
   end subroutine inner_anderson

   !=======================================================================================================================!
   ! outer secant (SLEVOLVE_OUTER=secant, outer iteration k >= 3).
   ! Outer map per cell: T0 -> P(T0) (predicted present topography); want P = Tobs.
   ! Picard:  T0_k = T0_{k-1} + (Tobs - P_{k-1})              (slope 1)
   ! Secant:  s = (P_{k-1} - P_{k-2}) / (T0_{k-1} - T0_{k-2}), clamped to [0.5, 1.5];
   !          T0_k = T0_{k-1} + (Tobs - P_{k-1}) / s
   !          Picard (s = 1) where |T0_{k-1} - T0_{k-2}| < 1e-2 m.
   ! On entry: pred_pres_topo = P_{k-1}, tinit_0_last = T0_{k-1}, truetopo = Tobs,
   !           init_topo_corr = Tobs - P_{k-1}.
   !=======================================================================================================================!
   subroutine outer_secant(itersl)
      integer, intent(in) :: itersl
      character(len=8) :: istr2
      real, allocatable :: p2(:,:), t2(:,:)
      real :: dt, sl
      integer :: nlo, nhi, nflat, nsec
      allocate(p2(nglv,2*nglv), t2(nglv,2*nglv))
      write(istr2,'(I2)') itersl-2
      istr2 = trim(adjustl(istr2))
      call pp_read_grid(p2, 'pred_pres_topo_', outputfolder, suffix=istr2)
      call pp_read_grid(t2, 'tgrid0_', outputfolder, suffix=istr2)
      nlo = 0; nhi = 0; nflat = 0; nsec = 0
      do j = 1,2*nglv
         do i = 1,nglv
            dt = tinit_0_last(i,j) - t2(i,j)
            if (abs(dt) < 1.0e-2) then
               sl = 1.0
               nflat = nflat + 1
            else
               sl = (pred_pres_topo(i,j) - p2(i,j)) / dt
               if (sl < 0.5) then
                  sl = 0.5; nlo = nlo + 1
               else if (sl > 1.5) then
                  sl = 1.5; nhi = nhi + 1
               else
                  nsec = nsec + 1
               endif
            endif
            tinit_0(i,j) = tinit_0_last(i,j) + init_topo_corr(i,j) / sl
         enddo
      enddo
      write(unit_num,'(A,I3,A,4I9)') ' SLEvolve outer secant it', itersl, &
           ': cells secant / clamped-low / clamped-high / picard(|dT0|<1cm) =', nsec, nlo, nhi, nflat
      deallocate(p2, t2)
   end subroutine outer_secant

   real(8) function inner_dot(x, y)
      complex, intent(in) :: x(0:norder,0:norder), y(0:norder,0:norder)
      real(8) :: acc
      integer :: mcol
      acc = 0.0d0
      !$omp parallel do schedule(static) default(shared) private(mcol) reduction(+:acc)
      do mcol = 0,norder
         acc = acc + sum(real(x(:,mcol))*real(y(:,mcol)) + aimag(x(:,mcol))*aimag(y(:,mcol)))
      enddo
      !$omp end parallel do
      inner_dot = acc
   end function inner_dot

   ! Gaussian elimination with partial pivoting, n <= 5; solution returned in b
   subroutine inner_solve(a, b, n, info)
      real(8), intent(inout) :: a(5,5), b(5)
      integer, intent(in) :: n
      integer, intent(out) :: info
      integer :: i, k, p
      real(8) :: t, rowt(5)
      info = 0
      do k = 1,n
         p = k
         do i = k+1,n
            if (abs(a(i,k)) > abs(a(p,k))) p = i
         enddo
         if (abs(a(p,k)) <= 0.0d0) then; info = 1; return; endif
         if (p /= k) then
            rowt(1:n) = a(k,1:n); a(k,1:n) = a(p,1:n); a(p,1:n) = rowt(1:n)
            t = b(k); b(k) = b(p); b(p) = t
         endif
         do i = k+1,n
            t = a(i,k) / a(k,k)
            a(i,k:n) = a(i,k:n) - t * a(k,k:n)
            b(i) = b(i) - t * b(k)
         enddo
      enddo
      do k = n,1,-1
         b(k) = (b(k) - sum(a(k,k+1:n)*b(k+1:n))) / a(k,k)
      enddo
   end subroutine inner_solve

   !---------------------------------------------------------------------
   ! SLEVOLVE_REUSE=warm file  OUTPUT_SLM/reuse_<ITERSL>.bin  (stream, little-endian)
   !   header 32 B: char(8) 'PPREUSE1', int32 norder, nlats, nlons, 0,
   !                int64 record length in bytes
   !   record for step n at offset 32 + (n-1)*reclen:
   !     int32 n, int32 0,
   !     complex(8) dS(l,m) for m = 0..norder, l = m..norder  (lower triangle)
   !     int8 ocean mask (nlats,nlons) column-major (Fortran order), 0/1
   ! dS(l,m>l) is identically 0 and is not stored.
   !---------------------------------------------------------------------
   subroutine warm_paths(itersl, fin, fout)
      integer, intent(in) :: itersl
      character(len=*), intent(out) :: fin, fout
      character(len=16) :: a
      write(a,'(I0)') itersl - 1; fin  = trim(outputfolder)//'reuse_'//trim(a)//'.bin'
      write(a,'(I0)') itersl;     fout = trim(outputfolder)//'reuse_'//trim(a)//'.bin'
   end subroutine warm_paths

   integer(8) function warm_reclen()
      warm_reclen = 8_8 + 16_8 * int(norder+1,8) * int(norder+2,8) / 2_8 &
                       + int(nglv,8) * int(2*nglv,8)
   end function warm_reclen

   subroutine warm_read(itersl, iter)
      integer, intent(in) :: itersl, iter
      character(len=1024) :: fin, fout
      character(len=8) :: magic
      integer :: u, ios, h(4), nstep, z, mm
      integer(8) :: rl, pos
      logical :: ex
      integer(8) :: c0, c1, cr
      call system_clock(c0, cr)
      if (.not. warm_init) then
         warm_init = .true.
         allocate(warm_ds(0:norder,0:norder), warm_mask(nglv,2*nglv))
         warm_ds = (0.0,0.0)
         call warm_paths(itersl, fin, fout)
         ex = .false.
         if (.not. pp_mem_outer) inquire(file=trim(fin), exist=ex)
         warm_infile = .false.
         if (ex .and. itersl > 1) then
            open(newunit=u, file=trim(fin), access='stream', form='unformatted', status='old', &
                 action='read', convert='little_endian', iostat=ios)
            if (ios == 0) then
               read(u, iostat=ios) magic, h, rl
               if (ios == 0 .and. magic == 'PPREUSE1' .and. h(1) == norder .and. h(2) == nglv &
                   .and. h(3) == 2*nglv .and. rl == warm_reclen()) warm_infile = .true.
               close(u)
            endif
            if (.not. warm_infile) write(unit_num,*) 'SLEvolve reuse: ignoring incompatible ', trim(fin)
         endif
         write(unit_num,*) 'SLEvolve reuse=warm: input ', trim(fin), ' present/valid=', warm_infile
      endif
      warm_have = .false.
      if (pp_mem_outer) then
         if (.not. allocated(warm_sok)) then
            allocate(warm_sds((norder+1)*(norder+2)/2, L_sim/dt1), &
                     warm_smk(nglv,2*nglv,L_sim/dt1), warm_sok(L_sim/dt1))
            warm_sok = .false.
         endif
         if (warm_sok(iter)) then
            pos = 0
            do mm = 0, norder
               warm_ds(mm:norder, mm) = warm_sds(pos+1:pos+norder-mm+1, iter)
               pos = pos + norder - mm + 1
            enddo
            warm_mask(:,:) = warm_smk(:,:,iter)
            warm_have = .true.
         endif
         call system_clock(c1)
         warm_tio = warm_tio + real(c1-c0,8)/real(cr,8)
         return
      endif
      if (.not. warm_infile) return
      call warm_paths(itersl, fin, fout)
      rl = warm_reclen()
      pos = 33_8 + int(iter-1,8) * rl
      open(newunit=u, file=trim(fin), access='stream', form='unformatted', status='old', &
           action='read', convert='little_endian', iostat=ios)
      if (ios /= 0) return
      read(u, pos=pos, iostat=ios) nstep, z
      if (ios == 0 .and. nstep == iter) then
         do mm = 0, norder
            read(u, iostat=ios) warm_ds(mm:norder, mm)
            if (ios /= 0) exit
         enddo
         if (ios == 0) read(u, iostat=ios) warm_mask
         warm_have = (ios == 0)
      endif
      close(u)
      call system_clock(c1)
      warm_tio = warm_tio + real(c1-c0,8)/real(cr,8)
   end subroutine warm_read

   subroutine warm_write(itersl, iter)
      integer, intent(in) :: itersl, iter
      character(len=1024) :: fin, fout
      integer :: u, mm, i, j
      integer(1), allocatable :: mk(:,:)
      integer(8) :: c0, c1, cr, pos
      real :: up
      call system_clock(c0, cr)
      if (pp_mem_outer) then
         pos = 0
         do mm = 0, norder
            warm_sds(pos+1:pos+norder-mm+1, iter) = load_inc_lm(mm:norder, mm)
            pos = pos + norder - mm + 1
         enddo
         !$omp parallel do schedule(static) default(shared) private(i,j)
         do j = 1, 2*nglv
            do i = 1, nglv
               warm_smk(i,j,iter) = int(merge(1, 0, cxy(i,j) /= 0.0), 1)
            enddo
         enddo
         !$omp end parallel do
         warm_sok(iter) = .true.
         call system_clock(c1)
         warm_tio = warm_tio + real(c1-c0,8)/real(cr,8)
         return
      endif
      call warm_paths(itersl, fin, fout)
      if (iter == 1) then
         open(newunit=u, file=trim(fout), access='stream', form='unformatted', status='replace', &
              action='write', convert='little_endian')
         write(u) 'PPREUSE1', norder, nglv, 2*nglv, 0, warm_reclen()
         close(u)
         up = 0.0
         do mm = 1, norder
            up = max(up, maxval(abs(load_inc_lm(0:mm-1, mm))))
         enddo
         write(unit_num,*) 'SLEvolve reuse=warm: writing ', trim(fout), '  max|dS(l<m)| =', up
      endif
      allocate(mk(nglv,2*nglv))
      do j = 1, 2*nglv
         do i = 1, nglv
            mk(i,j) = int(merge(1, 0, cxy(i,j) /= 0.0), 1)
         enddo
      enddo
      open(newunit=u, file=trim(fout), access='stream', form='unformatted', status='old', &
           action='write', position='append', convert='little_endian')
      write(u) iter, 0
      do mm = 0, norder
         write(u) load_inc_lm(mm:norder, mm)
      enddo
      write(u) mk
      close(u)
      deallocate(mk)
      call system_clock(c1)
      warm_tio = warm_tio + real(c1-c0,8)/real(cr,8)
   end subroutine warm_write

   !---------------------------------------------------------------------
   ! in-process outer iteration helpers (run_outer_mem.f90)
   !---------------------------------------------------------------------
   subroutine pp_ice_get(dst, idx, sfx)
      real, dimension(:,:), intent(inout) :: dst
      integer, intent(in) :: idx
      character(len=*), intent(in) :: sfx
      integer(8) :: c0, c1, cr
      call system_clock(c0, cr)
      if (pp_mem_outer .and. pp_icecache) then
         if (.not. allocated(pp_icec)) then
            allocate(pp_icec(nglv,2*nglv,0:L_sim/dt1), pp_icec_ok(0:L_sim/dt1))
            pp_icec_ok = .false.
         endif
         if (.not. pp_icec_ok(idx)) then
            call pp_read_grid(pp_icec(:,:,idx), icemodel, inputfolder_ice, suffix=sfx)
            pp_icec_ok(idx) = .true.
         endif
         dst(:,:) = pp_icec(:,:,idx)
      else
         call pp_read_grid(dst, icemodel, inputfolder_ice, suffix=sfx)
      endif
      call system_clock(c1); pp_t_ice = pp_t_ice + real(c1-c0,8)/real(cr,8)
   end subroutine pp_ice_get

   ! value -> ES16.9E2 text -> list-directed read, i.e. exactly what
   ! write_sl + read_sl do to tinit_0, without the file
   subroutine pp_quantize_es169(a)
      real, dimension(:,:), intent(inout) :: a
      character(len=16), allocatable :: buf(:)
      allocate(buf(size(a)))
      write(buf,'(ES16.9E2)') a      ! one record per value, as write_sl
      read(buf,*) a                  ! list-directed, as read_sl
      deallocate(buf)
   end subroutine pp_quantize_es169

   ! Between two in-process passes: restore the zero state a fresh process
   ! would have for every per-pass work array; keep the Earth model, the
   ! Legendre/FFTW plans, the hereditary-memory allocation (zeroed at step 1),
   ! the ice cache, the warm store and the pass result (pred_pres_topo,
   ! tinit_0_last, truetopo).
   subroutine sl_outer_reset()
      pp_mem_have = .true.
      pp_t0lm_ready = .false.        ! t0lm depends on tinit_0, which changes
      inner_nhist = 0                  ! extrap history is per pass
      if (allocated(extrap_dsm1)) then
         extrap_dsm1 = (0.0,0.0); extrap_dsm2 = (0.0,0.0)
      endif
      pp_gbuf = 0.0; pp_ubuf = 0.0
      tinit_0 = 0.0; init_topo_corr = 0.0
      slcxy_ocn = 0.0; slcxy_ocnBeta = 0.0; deltaslxy = 0.0; dslxy = 0.0; icestarxy = 0.0
      beta = 0.0; cstarxy = 0.0; cstar0 = 0.0; tOxy = 0.0; rOxy = 0.0; tTxy = 0.0
      rrxy = 0.0; drrxy_computed = 0.0; beta0 = 0.0; cxy0 = 0.0; cxy = 0.0
      topoxy = 0.0; topoxy_m1 = 0.0; tinit = 0.0
      slclm_ocn = 0.0; slclm_ocnBeta = 0.0; clm = 0.0; cstarlm = 0.0; oldcstarlm = 0.0
      tOlm = 0.0; rOlm = 0.0; load_inc_lm = 0.0; old_load_inc_lm = 0.0; icestarlm = 0.0; dicestarlm = 0.0
      deltaicestarlm = 0.0; oldicestarlm = 0.0; icestar0 = 0.0; t0lm = 0.0; oldt0lm = 0.0
      tTlm = 0.0; oldtTlm = 0.0; sl_inc_lm = 0.0; sl_tot_lm = 0.0; viscous = 0.0
      rrlm = 0.0; dgglm = 0.0; drrlm_computed = 0.0; dS_converged = 0.0
      times = 0.0; icexy = 0.0; sl = 0.0; dS = 0.0; deltaS = 0.0; dicestar = 0.0
      deltaicestar = 0.0; rr = 0.0; gg = 0.0; dil = 0.0; dlambda = 0.0; deltalambda = 0.0; dm = 0.0
      il = 0.0; mm = 0.0; lambda = 0.0
   end subroutine sl_outer_reset

   subroutine sl_print_counts(u)
      integer, intent(in) :: u
      write(u,'(A,I10,A,I10,A,I10,A,I10)') '== SLEvolve transforms: analyses ', opt_cnt(1), ' in sweeps ', opt_cnt(2), &
           ' ; syntheses ', opt_cnt(3), ' in sweeps ', opt_cnt(4)
      if (inner_mode /= 0) write(u,'(A,I2,A,6I9)') '== SLEvolve inner mode', inner_mode, &
           ' : mixed, fallbacks, forced-mask-updates, frozen-iters, cstar-analyses-skipped, steps =', inner_cnt
      write(u,'(A,4F9.3,A,F9.3)') '== SLEvolve setup s: init, ice-read, truetopo-read, tgrid0-quantize =', &
           pp_t_init, pp_t_ice, pp_t_topo, pp_t_quant, ' ; pp_read_grid total', pp_io_read
      if (warm_mode /= 0) write(u,'(A,I6,A,F8.3,A)') '== SLEvolve reuse=warm: warm-started steps', &
           warm_nwarm, ' ; reuse file I/O ', warm_tio, ' s'
   end subroutine sl_print_counts

   ! Weighted sum of latitude-row sums, serial in latitude order (so the
   ! result does not depend on the thread count).
   real function opt_quad(rs)
      real, intent(in) :: rs(:)
      integer :: jj
      opt_quad = 0.0
      do jj = 1,size(rs)
         opt_quad = opt_quad + opt_g(jj)*rs(jj)
      enddo
   end function opt_quad

   ! ro = dsl*cst [- tO]  and its latitude-row sums, one pass.  Arrays are the
   ! solver's grids viewed as (nlons, nlats), i.e. in the SHT's own order.
   subroutine opt_rowsum_ro(dsl, cst, tO, ro, rs, nlo, nla, sub_t)
      integer, intent(in) :: nlo, nla
      real, intent(in) :: dsl(nlo,nla), cst(nlo,nla), tO(nlo,nla)
      real, intent(out) :: ro(nlo,nla)
      real, intent(out) :: rs(nla)
      logical, intent(in) :: sub_t
      integer :: ii2, jj, kk
      real :: a(8)
      !$omp parallel do schedule(static) default(shared) private(ii2,jj,kk,a)
      do jj = 1,nla
         a = 0.0
         if (sub_t) then
            do ii2 = 1,nlo,8
               do kk = 0,7
                  ro(ii2+kk,jj) = dsl(ii2+kk,jj)*cst(ii2+kk,jj) - tO(ii2+kk,jj)
                  a(kk+1) = a(kk+1) + ro(ii2+kk,jj)
               enddo
            enddo
         else
            do ii2 = 1,nlo,8
               do kk = 0,7
                  ro(ii2+kk,jj) = dsl(ii2+kk,jj)*cst(ii2+kk,jj)
                  a(kk+1) = a(kk+1) + ro(ii2+kk,jj)
               enddo
            enddo
         endif
         rs(jj) = ((a(1)+a(5)) + (a(2)+a(6))) + ((a(3)+a(7)) + (a(4)+a(8)))
      enddo
      !$omp end parallel do
   end subroutine opt_rowsum_ro

   ! Row sums of  dsl*c,  c,  (dsl*c)*b,  c*b   (the four post-loop fields).
   subroutine opt_rowsum_gmsl(dsl, c, b, rs, nlo, nla)
      integer, intent(in) :: nlo, nla
      real, intent(in) :: dsl(nlo,nla), c(nlo,nla), b(nlo,nla)
      real, intent(out) :: rs(nla,4)
      integer :: ii2, jj
      real :: s1, s2, s3, s4, sl
      !$omp parallel do schedule(static) default(shared) private(ii2,jj,s1,s2,s3,s4,sl)
      do jj = 1,nla
         s1 = 0.0; s2 = 0.0; s3 = 0.0; s4 = 0.0
         do ii2 = 1,nlo
            sl = dsl(ii2,jj)*c(ii2,jj)
            s1 = s1 + sl
            s2 = s2 + c(ii2,jj)
            s3 = s3 + sl*b(ii2,jj)
            s4 = s4 + c(ii2,jj)*b(ii2,jj)
         enddo
         rs(jj,1) = s1; rs(jj,2) = s2; rs(jj,3) = s3; rs(jj,4) = s4
      enddo
      !$omp end parallel do
   end subroutine opt_rowsum_gmsl

   ! dst = src + c*K00  (K00 per latitude row; constant along longitude)
   subroutine opt_addconst(dst, src, c, nlo, nla)
      integer, intent(in) :: nlo, nla
      real, intent(out) :: dst(nlo,nla)
      real, intent(in) :: src(nlo,nla)
      real, intent(in) :: c
      integer :: jj
      real :: ck
      !$omp parallel do schedule(static) default(shared) private(jj,ck)
      do jj = 1,nla
         ck = c*opt_k00(jj)
         dst(:,jj) = src(:,jj) + ck
      enddo
      !$omp end parallel do
   end subroutine opt_addconst

end module sl_model_mod
