module spharmt
 implicit none
 private
 public :: sphere,spharmt_init,spharmt_destroy,spharm,spharm_ref,rfft,cosgrad,&
           specsmooth,getuv,getvrtdiv,sumnm,gaulw,legend,spat2spec,spec2spat, &
           spat2spec_multi,sht_backend_name,sht_set_threads, &
           spharm_blas_ana,spharm_blas_syn,spharm_exact_ana,spharm_exact_syn, &
           fft991,set99,fft991_mt,sht_env_init,spec2spat_multi,sht_nt, &
           spharm_blas_ana_u,spharm_blas_syn_u
 type sphere
   real    :: rsphere = 0.
   integer :: nlons = 0
   integer :: nlats = 0
   real, dimension(:), pointer :: trigs
   integer, dimension(:), pointer :: ifax
   integer :: ntrunc = 0
   real, dimension(:), pointer :: gaulats, weights, gwrc, lap, ilap
   real, dimension(:,:), pointer :: pnm, hnm
   integer, dimension(:), pointer :: indxm, indxn
   logical :: isinitialized=.FALSE.
   ! ---- additions: persistent work buffers and index tables --------
   ! (pointer components, so they stay mutable through an intent(in) dummy)
   real,    dimension(:), pointer, contiguous :: wrk1 => null()   ! nlats*(nlons+2) fft buffer
   real,    dimension(:), pointer, contiguous :: wrk2 => null()   ! nlats*(nlons+2) fft scratch
   real,    dimension(:), pointer, contiguous :: bbuf => null()   ! (2*nf*mwaves, nhs) analysis rhs / synthesis out
   real,    dimension(:), pointer, contiguous :: bbuf2 => null()  ! odd-parity partner of bbuf (symmetric mode)
   real,    dimension(:), pointer, contiguous :: ybuf => null()   ! (ntrunc+1)*2*nfmax legendre block
   integer, dimension(:), pointer, contiguous :: nmstrt_tab => null() ! offset of each m block (0-based)
   integer, dimension(:), pointer, contiguous :: s2s_idx => null()    ! packed index for spat2spec unpack
   ! Legendre tables laid out for the two fastest DGEMM shapes.
   !   pa(nmdim, nhs) : analysis   Y(lenm,ncol) = pa_m * B(ncol,nhs)^T   ('N','T')
   !   ps(nhs, nmdim) : synthesis  C(nhs,ncol)  = ps_m * Y(lenm,ncol)    ('N','N')
   ! In symmetric mode nhs = nlat/2 and the rows of each m block are
   ! permuted so that the even-parity degrees come first.
   real,    dimension(:,:), pointer, contiguous :: pa => null()
   real,    dimension(:,:), pointer, contiguous :: ps => null()
   ! mode 2 synthesis table, blocked over latitude:
   !   psb((nm-1)*SHT_JB+jj, blk) = pnm(nm, (blk-1)*SHT_JB+jj)
   ! so the SHT_JB latitudes needed for one coefficient are one cache line
   ! and consecutive coefficients are consecutive memory.
   real,    dimension(:,:), pointer, contiguous :: psb => null()
   ! mode 2 analysis table, blocked over coefficients:
   !   pab(slot + (j-1)*SHT_KB, bk) = pnm(nm,j)*weights(j)
   ! with nm = (bk-1)*SHT_KB + slot inside a zero-padded m block, so a whole
   ! coefficient block streams contiguously across all latitudes and the
   ! accumulators stay in registers.
   real,    dimension(:,:), pointer, contiguous :: pab => null()
   integer, dimension(:), pointer, contiguous :: nmpad_tab => null()
   integer :: nmpad = 0
   integer, dimension(:), pointer, contiguous :: parperm => null()  ! nmdim: k -> original n within m block
   integer :: nhs = 0
   integer :: nfmax = 0
   ! ---- FFTW plan cache (C pointers held as integer(8)) ---------
   ! plans are created for a given `howmany` row count; at most SHT_NPLAN
   ! distinct row counts occur (one per thread-chunk size).
   integer          :: nplan = 0          ! threads the plan table was built for
   integer          :: planh(128) = 0      ! rows per plan
   integer(8)       :: planf(128) = 0_8    ! r2c
   integer(8)       :: planb(128) = 0_8    ! c2r
 end type sphere
 ! Maximum number of fields a single batched transform may carry.
 integer, parameter, public :: SHT_NFMAX = 4
 ! latitudes carried simultaneously through the exact synthesis reduction
 integer, parameter, public :: SHT_JB = 8
 ! coefficient-tile length for the exact analysis (keeps both accumulators
 ! inside L1 while the Legendre table streams from memory)
 integer, parameter, public :: SHT_NB = 4096
 ! coefficients carried simultaneously through the exact analysis
 integer, parameter, public :: SHT_KB = 8
 ! Set before spharmt_init:
 !   sht_lean = .true.  -> do not keep the reference pnm/hnm tables (saves 1.1 GB)
 !   sht_sym  = .true.  -> exploit P_lm(-mu) = (-1)^(l+m) P_lm(mu) (halves the
 !                         Legendre table traffic).  Requires an even nlat; the
 !                         Gauss nodes are made exactly antisymmetric, which
 !                         perturbs the quadrature at the 1e-16 level.
 logical, public :: sht_lean = .false.
 ! Waive the Gauss-Legendre bandlimit guard (needed only to reproduce the
 ! historical norder=512 / nglv=512 runs).  Set from the environment variable
 ! SHT_UNGUARDED=1 by sht_env_init, called from the solver before spharmt_init.
 logical, public :: sht_unguarded = .false.
 ! FFT along longitude:  0 = Temperton fft991 (exact legacy results), 1 = FFTW3.
 ! fft991's forward transform already carries the 1/nlon normalisation;
 ! FFTW's r2c does not, so the analysis path multiplies by sht_fsc.  nlon is a
 ! power of two here, so this scaling is exact and changes no rounding.
 integer, public :: sht_fftw = 0
 real,    public :: sht_fsc  = 1.0
 ! The legacy spat2spec divides the coefficients by 1.41421356 and its spec2spat
 ! multiplies them by it.  In the direct-(l,m) path the analysis factor is
 ! folded into the Gauss weights and the synthesis factor into the coefficient
 ! pack, so neither costs a separate pass over the spectrum.
 real, parameter, public :: SHT_SQRT2  = 1.41421356
 real, parameter, public :: SHT_RSQRT2 = 1.0/1.41421356
#ifdef USE_EXACT_SHT
 logical, public :: sht_sym  = .false.
 integer, parameter, public :: sht_mode = 2   ! reordered kernels, exact results
#else
#ifdef USE_BLAS_SHT
 logical, public :: sht_sym  = .true.
 integer, parameter, public :: sht_mode = 1   ! BLAS (DGEMM) kernels
#else
 logical, public :: sht_sym  = .true.
 integer, parameter, public :: sht_mode = 0   ! original triple loop
#endif
#endif

 ! Optional internal profiling counters (seconds): 1 pack, 2 fft, 3 transpose,
 ! 4 gemm, 5 pack/unpack of spectral coefficients.
 real(8), public :: sht_times(8) = 0.0d0
 integer, public :: sht_counts(8) = 0
 logical, public :: sht_profile = .false.
 contains
 subroutine spharmt_init(sphere_dat,nlon,nlat,ntrunc,re)
   integer, intent(in) :: nlon,nlat,ntrunc
   real, intent(in) :: re
   type (sphere), intent(inout) :: sphere_dat
   real, dimension(:), allocatable :: pnm_tmp,hnm_tmp
   double precision, dimension(:), allocatable :: gaulats_tmp,weights_tmp,weights_tmp2
   integer :: nmdim,m,n,j,nmstrt0,lenm0,kk
   logical :: need_pnm,need_hnm,need_pa,need_ps,need_psb
   ! GUARD (from src_patched/spharmt_fft99bfix_guarded.f90): Gauss-Legendre
   ! quadrature with nlat nodes is exact only for degrees <= nlat-1.  the legacy
   ! original setting ntrunc=512 with nlat=512 loses the (512,0) mode entirely
   ! and doubles (512,1).  Refuse such configurations unless explicitly waived
   ! (SHT_UNGUARDED=1) so that the historical norder=512 runs can be reproduced.
   if (ntrunc > nlat-1 .and. .not. sht_unguarded) then
      print *, 'spharmt_init: ntrunc must be <= nlat-1 for exact GL analysis; got ntrunc=',ntrunc,' nlat=',nlat
      stop 'spharmt_init: spectral truncation exceeds Gauss-Legendre bandlimit'
   end if
   sphere_dat%nlons = nlon
   sphere_dat%nlats = nlat
   sphere_dat%ntrunc = ntrunc
   sphere_dat%rsphere = re
   nmdim = (ntrunc+1)*(ntrunc+2)/2
   allocate(gaulats_tmp(nlat))
   allocate(weights_tmp(nlat))
   call gaulw(gaulats_tmp,weights_tmp)
   allocate(sphere_dat%gwrc(nlat))
   sphere_dat%gwrc = weights_tmp/(dble(re)*(1.d0-gaulats_tmp**2))
   allocate(sphere_dat%gaulats(nlat))
   allocate(sphere_dat%weights(nlat))
   sphere_dat%gaulats = gaulats_tmp
   sphere_dat%weights = weights_tmp
   allocate(weights_tmp2(nlat)); weights_tmp2 = weights_tmp
   deallocate(weights_tmp)
   allocate(sphere_dat%indxm(nmdim))
   allocate(sphere_dat%indxn(nmdim))
   allocate(sphere_dat%lap(nmdim))
   allocate(sphere_dat%ilap(nmdim))
   sphere_dat%indxm = (/((m,n=m,ntrunc),m=0,ntrunc)/)
   sphere_dat%indxn = (/((n,n=m,ntrunc),m=0,ntrunc)/)
   sphere_dat%lap(:)=-real(sphere_dat%indxn(:))*real(sphere_dat%indxn(:)+1)/re**2
   sphere_dat%ilap(1) = 0.
   sphere_dat%ilap(2:nmdim) = 1./sphere_dat%lap(2:nmdim)
   allocate(sphere_dat%nmstrt_tab(ntrunc+1))
   nmdim = 0
   do m = 1, ntrunc+1
      sphere_dat%nmstrt_tab(m) = nmdim
      nmdim = nmdim + (ntrunc+1-m+1)
   enddo
   if (sht_sym) then
      if (mod(nlat,2) /= 0) then
         print *,'spharmt_init: sht_sym requires an even number of latitudes'; stop
      endif
      ! Make the Gauss nodes/weights exactly antisymmetric so that the parity
      ! relation of the Legendre functions holds exactly.
      do j = 1,nlat/2
         gaulats_tmp(nlat+1-j) = -gaulats_tmp(j)
         weights_tmp2(nlat+1-j) = weights_tmp2(j)
      enddo
      sphere_dat%gaulats = gaulats_tmp
      sphere_dat%weights = weights_tmp2
      sphere_dat%gwrc = weights_tmp2/(dble(re)*(1.d0-gaulats_tmp**2))
      sphere_dat%nhs = nlat/2
   else
      sphere_dat%nhs = nlat
   endif
   ! parity permutation within each m block: even-parity degrees (k odd) first
   allocate(sphere_dat%parperm(nmdim))
   nmstrt0 = 0
   do m = 0,ntrunc
      lenm0 = ntrunc+1-m
      if (sht_sym) then
         kk = 0
         do n = 1,lenm0,2
            kk = kk + 1; sphere_dat%parperm(nmstrt0+kk) = n
         enddo
         do n = 2,lenm0,2
            kk = kk + 1; sphere_dat%parperm(nmstrt0+kk) = n
         enddo
      else
         do kk = 1,lenm0
            sphere_dat%parperm(nmstrt0+kk) = kk
         enddo
      endif
      nmstrt0 = nmstrt0 + lenm0
   enddo
   ! Which Legendre tables does the selected backend actually need?
   !   mode 0 (reference)  pnm + hnm                        1.08 GB at L=512
   !   mode 1 (BLAS)       pa (analysis) + ps (synthesis)   0.54 GB symmetric
   !   mode 2 (exact)      pnm (synthesis) + pa = pnm*w     1.08 GB
   ! hnm is only used by getuv/cosgrad, which the solver never calls.
   need_pnm = (sht_mode == 0) .or. (.not. sht_lean)
   need_hnm = (sht_mode == 0) .or. (.not. sht_lean)
   need_pa  = (sht_mode == 1)
   need_ps  = (sht_mode == 1)
   need_psb = (sht_mode == 2)
   if (need_pnm) then; allocate(sphere_dat%pnm(nmdim,nlat)); else; allocate(sphere_dat%pnm(1,1)); endif
   if (need_hnm) then; allocate(sphere_dat%hnm(nmdim,nlat)); else; allocate(sphere_dat%hnm(1,1)); endif
   if (need_pa)  then; allocate(sphere_dat%pa(nmdim,sphere_dat%nhs)); else; allocate(sphere_dat%pa(1,1)); endif
   if (need_ps)  then; allocate(sphere_dat%ps(sphere_dat%nhs,nmdim)); else; allocate(sphere_dat%ps(1,1)); endif
   if (need_psb) then
      allocate(sphere_dat%nmpad_tab(ntrunc+1))
      nmstrt0 = 0
      do m = 0,ntrunc
         sphere_dat%nmpad_tab(m+1) = nmstrt0
         nmstrt0 = nmstrt0 + ((ntrunc+1-m + SHT_KB-1)/SHT_KB)*SHT_KB
      enddo
      sphere_dat%nmpad = nmstrt0
      allocate(sphere_dat%pab(SHT_KB*nlat, sphere_dat%nmpad/SHT_KB))
      sphere_dat%pab = 0.0
      if (mod(nlat,SHT_JB) /= 0) then
         print *,'spharmt_init: exact backend needs nlat divisible by', SHT_JB; stop
      endif
      allocate(sphere_dat%psb(SHT_JB*nmdim, nlat/SHT_JB))
   else
      allocate(sphere_dat%psb(1,1)); allocate(sphere_dat%pab(1,1))
      allocate(sphere_dat%nmpad_tab(1))
   endif
   allocate(pnm_tmp(nmdim))
   allocate(hnm_tmp(nmdim))
   do j=1,nlat
      call legend(gaulats_tmp(j),pnm_tmp,hnm_tmp)
      if (need_pnm) sphere_dat%pnm(:,j) = pnm_tmp(:)
      if (need_hnm) sphere_dat%hnm(:,j) = hnm_tmp(:)
      if (sht_mode == 2) then
         kk = (j-1)/SHT_JB + 1
         lenm0 = j - (kk-1)*SHT_JB
         do n = 1,nmdim
            sphere_dat%psb((n-1)*SHT_JB+lenm0, kk) = pnm_tmp(n)
         enddo
         do m = 0,ntrunc
            nmstrt0 = sphere_dat%nmstrt_tab(m+1)
            do n = 1,ntrunc+1-m
               kk = sphere_dat%nmpad_tab(m+1) + n - 1
               ! exactly the product the reference forms in its inner loop
               sphere_dat%pab(mod(kk,SHT_KB)+1 + (j-1)*SHT_KB, kk/SHT_KB+1) = &
                    pnm_tmp(nmstrt0+n)*sphere_dat%weights(j)
            enddo
         enddo
      else if (need_pa .and. j <= sphere_dat%nhs) then
         nmstrt0 = 0
         do m = 0,ntrunc
            lenm0 = ntrunc+1-m
            do kk = 1,lenm0
               sphere_dat%pa(nmstrt0+kk,j) = pnm_tmp(nmstrt0+sphere_dat%parperm(nmstrt0+kk))
               sphere_dat%ps(j,nmstrt0+kk) = pnm_tmp(nmstrt0+sphere_dat%parperm(nmstrt0+kk))
            enddo
            nmstrt0 = nmstrt0 + lenm0
         enddo
      endif
   enddo
   deallocate(gaulats_tmp)
   deallocate(weights_tmp2)
   deallocate(pnm_tmp)
   deallocate(hnm_tmp)
   allocate(sphere_dat%trigs((3*nlon/2)+1))
   allocate(sphere_dat%ifax(13))
   call set99(sphere_dat%trigs,sphere_dat%ifax,nlon)
   ! ---- persistent work buffers + index tables -------------------
   sphere_dat%nfmax = SHT_NFMAX
   allocate(sphere_dat%wrk1(nlat*(nlon+2)*SHT_NFMAX))
   allocate(sphere_dat%wrk2(nlat*(nlon+2)*SHT_NFMAX))
   allocate(sphere_dat%bbuf(sphere_dat%nhs*2*(ntrunc+1)*SHT_NFMAX))
   allocate(sphere_dat%bbuf2(sphere_dat%nhs*2*(ntrunc+1)*SHT_NFMAX))
   allocate(sphere_dat%ybuf(max(2*(ntrunc+1),2*nmdim)*SHT_NFMAX))
   sphere_dat%wrk1 = 0.0
   sphere_dat%wrk2 = 0.0
   sphere_dat%bbuf = 0.0
   sphere_dat%bbuf2 = 0.0
   sphere_dat%ybuf = 0.0
   ! spat2spec unpack: for (n,m) with m<=n, U(n,m)=Uvec(s2s_idx(n*(ntrunc+1)+m+1))
   allocate(sphere_dat%s2s_idx(0:(ntrunc+1)*(ntrunc+1)-1))
   sphere_dat%s2s_idx = 0
   do m = 0, ntrunc
      do n = m, ntrunc
         sphere_dat%s2s_idx(n*(ntrunc+1)+m) = sphere_dat%nmstrt_tab(m+1) + n - m + 1
      enddo
   enddo
   ! ---- FFTW plans (one row, in place; the row loop is threaded) --
#ifdef HAVE_FFTW
   if (sht_fftw /= 0) then
      block
        integer :: nn(1), it
        integer, parameter :: FFTW_MEASURE = 0, FFTW_ESTIMATE = 64, FFTW_UNALIGNED = 2
        nn(1) = nlon
        ! ONE BATCHED PLAN PER (nf, thread-chunk).
        !
        ! Two lessons paid for here.  (1) FFTW_MEASURE picks its algorithm by
        ! timing, and for this 1024-point r2c several candidates tie, so the
        ! choice differs from run to run -- the solver stopped being
        ! reproducible (~1 ulp in the deepest grid cell).  FFTW_ESTIMATE is a
        ! deterministic heuristic and restores bit-reproducibility.  (2) With
        ! ESTIMATE a row-at-a-time plan is 50 % slower than a batched one
        ! (1.04 vs 0.70 ms per 512x1024 r2c), while with MEASURE they tie --
        ! so plan the whole thread chunk at once.  Every concurrently
        ! executing chunk gets its own plan object, so no plan is ever
        ! executed by two threads at the same time.
        sphere_dat%nplan = min(max(sht_nt(),1), 16)
        do it = 1, SHT_NFMAX*sphere_dat%nplan
          block
            integer :: nf0, c0p, lot0, base0, rem0, h0
            nf0  = (it-1)/sphere_dat%nplan + 1
            c0p  = it - (nf0-1)*sphere_dat%nplan
            lot0 = nlat*nf0
            base0 = lot0/sphere_dat%nplan; rem0 = mod(lot0,sphere_dat%nplan)
            if (c0p <= rem0) then; h0 = base0+1; else; h0 = base0; endif
            sphere_dat%planh(it) = h0
            if (h0 <= 0) cycle
            call dfftw_plan_many_dft_r2c(sphere_dat%planf(it), 1, nn, h0, sphere_dat%wrk1, &
                 nn, 1, nlon+2, sphere_dat%wrk1, nn, 1, (nlon+2)/2, FFTW_ESTIMATE+FFTW_UNALIGNED)
            call dfftw_plan_many_dft_c2r(sphere_dat%planb(it), 1, nn, h0, sphere_dat%wrk1, &
                 nn, 1, (nlon+2)/2, sphere_dat%wrk1, nn, 1, nlon+2, FFTW_ESTIMATE+FFTW_UNALIGNED)
            if (sphere_dat%planf(it) == 0_8 .or. sphere_dat%planb(it) == 0_8) then
               print *,'spharmt_init: FFTW planning failed'; stop
            endif
          end block
        enddo
        sphere_dat%wrk1 = 0.0
      end block
   endif
#endif
   sht_fsc = 1.0
   if (sht_fftw /= 0) sht_fsc = 1.0/real(nlon)
   sphere_dat%isinitialized = .TRUE.
 end subroutine spharmt_init
 !--------------------------------------------------------------------------
 function sht_backend_name() result(nm)
   character(len=32) :: nm
#if defined(USE_EXACT_SHT)
   nm = 'exact'
#elif defined(USE_BLAS_SHT)
   nm = 'blas'
#else
   nm = 'reference'
#endif
 end function sht_backend_name
 subroutine sht_set_threads(n)
   integer, intent(in) :: n
   integer :: dummy
   dummy = n
 end subroutine sht_set_threads
 subroutine spharmt_destroy(sphere_dat)
   type (sphere), intent(inout) :: sphere_dat
   if (.not. sphere_dat%isinitialized) return
   deallocate(sphere_dat%gaulats)
   deallocate(sphere_dat%weights)
   deallocate(sphere_dat%gwrc)
   deallocate(sphere_dat%pnm)
   deallocate(sphere_dat%hnm)
   deallocate(sphere_dat%indxm)
   deallocate(sphere_dat%indxn)
   deallocate(sphere_dat%lap)
   deallocate(sphere_dat%ilap)
   deallocate(sphere_dat%trigs)
   deallocate(sphere_dat%ifax)
   if (associated(sphere_dat%wrk1)) deallocate(sphere_dat%wrk1)
   if (associated(sphere_dat%wrk2)) deallocate(sphere_dat%wrk2)
   if (associated(sphere_dat%bbuf)) deallocate(sphere_dat%bbuf)
   if (associated(sphere_dat%bbuf2)) deallocate(sphere_dat%bbuf2)
   if (associated(sphere_dat%pa)) deallocate(sphere_dat%pa)
   if (associated(sphere_dat%ps)) deallocate(sphere_dat%ps)
   if (associated(sphere_dat%psb)) deallocate(sphere_dat%psb)
   if (associated(sphere_dat%pab)) deallocate(sphere_dat%pab)
   if (associated(sphere_dat%nmpad_tab)) deallocate(sphere_dat%nmpad_tab)
   if (associated(sphere_dat%parperm)) deallocate(sphere_dat%parperm)
   if (associated(sphere_dat%ybuf)) deallocate(sphere_dat%ybuf)
   if (associated(sphere_dat%nmstrt_tab)) deallocate(sphere_dat%nmstrt_tab)
   if (associated(sphere_dat%s2s_idx)) deallocate(sphere_dat%s2s_idx)
   sphere_dat%isinitialized = .FALSE.
 end subroutine spharmt_destroy
 subroutine gaulw(sinlats, wts)
 double precision, intent (inout), dimension(:) :: sinlats, wts
 integer :: itermax
 integer :: i, iter, j, nlat, nprec
 double precision :: pi, pp, p1, p2, p3, z, z1, converg, ten
 ten = 10.d0
 converg = ten * EPSILON(ten)
 nprec = precision(converg)
 converg = .1**nprec
 itermax = 10
 pi = 4.0*ATAN(1.0)
 nlat = size(sinlats)
 if (size(sinlats) .ne. size(wts)) then
   print *, 'sinlats and wts must be same size in gaulw!'
   stop
 end if
 do i=1,nlat
   z = cos(pi*(i - 0.25)/(nlat + 0.5))
   do iter=1,itermax
     p1 = 1.0
     p2 = 0.0
     do j=1,nlat
        p3 = p2
        p2 = p1
        p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j
     end do
     pp = nlat*(z*p1 - p2)/(z*z - 1.0E+00)
     z1 = z
     z  = z1 - p1/pp
     if(ABS(z - z1) .LT. converg) go to 10
   end do
   print *, 'abscissas failed to converge in itermax iterations'
   print *, 'stopping in gaulw!'
   stop
10 continue
   sinlats (i)     = z
   wts (i)     = 2.0/((1.0 - z*z)*pp*pp)
 end do
 return
 end subroutine gaulw
 SUBROUTINE LEGEND(X,PMN,HMN)
      real, dimension(:), intent(inout) ::  PMN,  HMN
      double precision, intent(in) :: x
      integer :: m,n,nm,i,nmax,np1,nmstrt,j,nmdim,ntrunc
      double precision :: A, B, PROD, SINSQ, &
      eps,epsnmp1,EPSPMN, PMNJ, PMNJM1, PMNJM2
      SINSQ = 1.D0 - X * X
      A     = 1.D0
      B     = 0.D0
      PROD  = 1.D0
      nmdim = size(pmn)
      ntrunc = nint((-1.+sqrt(1+8*float(nmdim)))/2.)-1
      if (size(pmn) .ne. size(hmn)) then
         print *, 'pnm and hnm must be same size in subroutine legend!'
         stop
      end if
      NMSTRT = 0
      DO I = 1, NTRUNC+1
         M    = (I - 1)
         NMAX = NTRUNC+1-M
         IF (M .NE. 0) THEN
               A    = A + 2.D0
               B    = B + 2.D0
               PROD = PROD * SINSQ * A / B
         END IF
         PMNJM2   = SQRT(0.5D0 * PROD)
         NM = NMSTRT + 1
         PMN(NM) = PMNJM2
         PMNJM1   = SQRT( DBLE(2 * M + 3) ) * X * PMNJM2
         IF (NM .NE. NMDIM) PMN(NM+1) = PMNJM1
         NP1 = M + 1
         EPSNMP1 = SQRT( DBLE(NP1*NP1 - M*M) / DBLE(4*NP1*NP1 - 1) )
         EPSPMN   = X * PMNJM1 - EPSNMP1 * PMNJM2
         HMN(NM) = DBLE(M) * EPSNMP1 * PMNJM1
         IF (NM .NE. NMDIM) &
         HMN(NM+1) = DBLE(M+1) * EPSPMN  -  &
         DBLE(M+2) * EPSNMP1 * PMNJM2
         DO J    = 3, NMAX
            N = M + J - 1
            NM = NMSTRT + J
            EPS = SQRT( DBLE(N*N - M*M) / DBLE(4*N*N - 1) )
            PMNJ     = EPSPMN / EPS
            PMN(NM) = PMNJ
            EPSPMN   = X * PMNJ - EPS * PMNJM1
            HMN(NM) = DBLE(N) * EPSPMN -  &
            DBLE(N+1) * EPS * PMNJM1
            PMNJM2   = PMNJM1
            PMNJM1   = PMNJ
         ENDDO
         NMSTRT = NMSTRT + NMAX
     ENDDO
 END SUBROUTINE LEGEND
 subroutine rfft(sphere_dat, data, coeff, idir)
 type (sphere), intent(in) :: sphere_dat
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(inout) :: data
 real, dimension(sphere_dat%nlats*(sphere_dat%nlons+2)) :: wrk1
 real, dimension(sphere_dat%nlats*(sphere_dat%nlons+2)) :: wrk2
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats), intent(inout) :: coeff
 integer ::  nlons,nlats,ntrunc,mwaves,i,j,m,n
 integer, intent(in) :: idir
 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in rfft!'
    stop
 end if
 nlons = sphere_dat%nlons
 nlats = sphere_dat%nlats
 ntrunc = sphere_dat%ntrunc
 if (ntrunc .gt. nlons/2) then
    print *, 'ntrunc must be less than or equal to nlons in rfft'
    stop
 end if
 mwaves = ntrunc+1
 wrk2 = 0.0
 if (idir .eq. +1) then
   n = 0
   wrk1 = 0.
   do j=1,nlats
   do i=1,nlons+2
      n = n + 1
      wrk1(n) = 0.0
      if (i .le. nlons) then
         wrk1(n) = data(i,j)
      end if
   enddo
   enddo
   call fft991(wrk1,wrk2,sphere_dat%trigs,sphere_dat%ifax,1,nlons+2,nlons,nlats,-1)
   n = -1
   do j=1,nlats
   do m=1,(nlons/2)+1
      n = n + 2
      if (m .le. mwaves) then
        coeff(m,j) = cmplx(wrk1(n),wrk1(n+1))
      end if
   enddo
   enddo
 else if (idir .eq. -1) then
   wrk1 = 0.
   n = -1
   do j=1,nlats
   do m=1,(nlons/2)+1
      n = n + 2
      if (m .le. mwaves) then
         wrk1(n) = real(coeff(m,j))
         wrk1(n+1) = aimag(coeff(m,j))
      end if
   enddo
   enddo
   call fft991(wrk1,wrk2,sphere_dat%trigs,sphere_dat%ifax,1,nlons+2,nlons,nlats,1)
   n = 0
   do j=1,nlats
   do i=1,nlons+2
      n = n + 1
      if (i .le. nlons) then
         data(i,j) = wrk1(n)
      end if
   enddo
   enddo
 else
    write(6,*) ' idir must be +1 or -1 in RFFT!'
    write(6,*) ' execution terminated.'
    stop
 end if
 end subroutine rfft
 subroutine spharm_ref(sphere_dat, ugrid, anm, idir)
 type (sphere), intent(in) :: sphere_dat
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(inout) :: ugrid
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2), intent(inout) :: anm
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats) :: am
 integer ::  nlats,ntrunc,mwaves,nmstrt,nm,m,n,j
 integer, intent(in) :: idir
 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in spharm!'
    stop
 end if
 nlats = sphere_dat%nlats
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 IF (IDIR .eq. +1) then
   anm = 0.
   call rfft(sphere_dat, ugrid, am, 1)
   do j=1,nlats
     NMSTRT = 0
     DO m = 1, mwaves
       DO n = 1, mwaves-m+1
          NM = NMSTRT + n
          anm(NM)=anm(NM)+sphere_dat%pnm(nm,j)*sphere_dat%weights(j)*am(m,j)
       enddo
       nmstrt = nmstrt + mwaves-m+1
     enddo
   enddo
 else if (idir .eq. -1) then
 DO J = 1, NLATS
    NMSTRT = 0
    do m = 1, MWAVES
      am(m,j) = cmplx(0., 0.)
      DO n = 1, mwaves-m+1
         NM = NMSTRT + n
         am(m,j) = am(m,j)  +  anm(NM) * sphere_dat%pnm(NM,j)
      enddo
      NMSTRT = NMSTRT + mwaves-m+1
    enddo
 ENDDO
   call rfft(sphere_dat, ugrid, am, -1)
 else
   print *, 'error in spharm: idir must be -1 or +1!'
   print *, 'execution terminated in subroutine spharm'
 end if
 end subroutine spharm_ref
!=========================================================================
! BLAS-based spherical harmonic transform.
!
! Mathematically identical to spharm_ref; the Legendre sums that dominate
! the runtime are expressed as one DGEMM per zonal wavenumber m, and the
! FFT work arrays are persistent instead of being (re)created on the stack
! and zeroed on every call.  nf >= 1 fields are transformed in one sweep
! over the pnm table, so the 540 MB table is streamed once for nf fields.
!=========================================================================
 !=====================================================================
 ! analysis:  grid -> spectral, nf fields at once.
 !   step 1  pack the nf grids into the persistent FFT buffer
 !   step 2  one batched real FFT (Temperton fft991, lot = nlats*nf)
 !   step 3  interleave into bbuf(2*nf*mwaves, nhs), folding in the Gauss
 !           weights and (symmetric mode) the north/south sum and difference
 !   step 4  one DGEMM('N','T') per zonal wavenumber
 !=====================================================================
 subroutine spharm_blas_ana(sphere_dat, ugrid, anm, nf)
 type (sphere), intent(in) :: sphere_dat
 integer, intent(in) :: nf
 real, dimension(sphere_dat%nlons*sphere_dat%nlats*nf), intent(in) :: ugrid
 complex, dimension(((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2)*nf), intent(inout) :: anm
 integer :: nlats,nlons,ntrunc,mwaves,nmstrt,n,j,f,nl2,base,lenm,nmdim,nhs
 integer :: i2,off,goff,ldb,ne,no,c0,cs,basen,bases,k
 real :: wj, a1, a2
 real(8) :: tt
 nlats  = sphere_dat%nlats
 nlons  = sphere_dat%nlons
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 nmdim  = (ntrunc+1)*(ntrunc+2)/2
 nhs    = sphere_dat%nhs
 nl2    = nlons+2
 ldb    = 2*nf*mwaves
 call sht_check(sphere_dat,nf)
 call sht_tic(tt)
 !$omp parallel do collapse(2) schedule(static) default(shared) private(f,j,off,goff,base)
 do f = 1,nf
   do j = 1,nlats
     off  = (f-1)*nlats*nl2
     goff = (f-1)*nlons*nlats
     base = off + (j-1)*nl2
     sphere_dat%wrk1(base+1:base+nlons) = ugrid(goff+(j-1)*nlons+1 : goff+(j-1)*nlons+nlons)
     sphere_dat%wrk1(base+nlons+1) = 0.0
     sphere_dat%wrk1(base+nlons+2) = 0.0
   enddo
 enddo
 !$omp end parallel do
 call sht_toc(1,tt)
 call sht_fft(sphere_dat, sphere_dat%wrk1, sphere_dat%wrk2, nlats*nf, -1)
 call sht_toc(2,tt)
 ! ---- interleave (+ weights, + parity folding) ------------------------
 if (sht_sym) then
   !$omp parallel do schedule(static) default(shared) &
   !$omp   private(j,wj,c0,f,basen,bases,cs,i2,a1,a2,k)
   do j = 1,nhs
     wj = sphere_dat%weights(j)*sht_fsc
     c0 = (j-1)*ldb
     do f = 1,nf
       basen = (f-1)*nlats*nl2 + (j-1)*nl2
       bases = (f-1)*nlats*nl2 + (nlats-j)*nl2
       cs = c0 + 2*(f-1)
       do i2 = 1, 2*mwaves
         a1 = sphere_dat%wrk1(basen+i2)
         a2 = sphere_dat%wrk1(bases+i2)
         k = cs + ((i2+1)/2-1)*2*nf + (2-mod(i2,2))
         sphere_dat%bbuf(k)  = (a1+a2)*wj
         sphere_dat%bbuf2(k) = (a1-a2)*wj
       enddo
     enddo
   enddo
   !$omp end parallel do
 else
   !$omp parallel do schedule(static) default(shared) &
   !$omp   private(j,wj,c0,f,basen,cs,i2,k)
   do j = 1,nlats
     wj = sphere_dat%weights(j)*sht_fsc
     c0 = (j-1)*ldb
     do f = 1,nf
       basen = (f-1)*nlats*nl2 + (j-1)*nl2
       cs = c0 + 2*(f-1)
       do i2 = 1, 2*mwaves
         k = cs + ((i2+1)/2-1)*2*nf + (2-mod(i2,2))
         sphere_dat%bbuf(k) = sphere_dat%wrk1(basen+i2)*wj
       enddo
     enddo
   enddo
   !$omp end parallel do
 endif
 call sht_toc(3,tt)
 ! ---- Legendre analysis -------------------------------------------------
 ! One DGEMM per zonal wavenumber; the m loop is parallel (work per m is
 ! proportional to lenm, hence schedule(dynamic,1)).  Each m writes a
 ! disjoint slice of anm, so the result does not depend on the schedule.
 !$omp parallel do schedule(dynamic,1) default(shared) &
 !$omp   private(j,lenm,nmstrt,ne,no,f,n,off,goff)
 do j = 1,mwaves
   lenm   = mwaves-j+1
   nmstrt = sphere_dat%nmstrt_tab(j)
   block
     real :: yb(2*mwaves*nf)
     if (sht_sym) then
       ne = (lenm+1)/2
       no = lenm-ne
       call sht_gemm('N','T', ne, 2*nf, nhs, sphere_dat%pa, nmstrt+1, nmdim, &
                     sphere_dat%bbuf, (j-1)*2*nf+1, ldb, yb, 1, mwaves)
       if (no > 0) &
         call sht_gemm('N','T', no, 2*nf, nhs, sphere_dat%pa, nmstrt+ne+1, nmdim, &
                       sphere_dat%bbuf2, (j-1)*2*nf+1, ldb, yb, ne+1, mwaves)
     else
       call sht_gemm('N','T', lenm, 2*nf, nhs, sphere_dat%pa, nmstrt+1, nmdim, &
                     sphere_dat%bbuf, (j-1)*2*nf+1, ldb, yb, 1, mwaves)
     endif
     ! yb is (mwaves, 2*nf) with rows in parity-permuted order
     do f = 1,nf
       off  = (2*f-2)*mwaves
       goff = (2*f-1)*mwaves
       do n = 1,lenm
         anm((f-1)*nmdim+nmstrt+sphere_dat%parperm(nmstrt+n)) = &
              cmplx(yb(off+n), yb(goff+n))
       enddo
     enddo
   end block
 enddo
 !$omp end parallel do
 call sht_toc(4,tt)
 end subroutine spharm_blas_ana

 !=====================================================================
 ! the same analysis, writing straight into the solver's
 ! (0:ntrunc, 0:ntrunc) coefficient layout.
 !
 ! The legacy path built a packed complex vector Uvec of length
 ! (L+1)(L+2)/2, then ran a separate pass  U(n,m) = Uvec(idx)/1.41421356
 ! over it -- a 2.1 MB temporary, an extra read+write of the whole
 ! spectrum and 2*(L+1)(L+2)/2 floating-point *divisions* per analysis.
 ! Here the m-block of the DGEMM result is written directly into the
 ! column U(m:L, m, f) (contiguous!), and the 1/1.41421356 is folded into
 ! the Gauss weights, where it costs nothing at all.
 !=====================================================================
 subroutine spharm_blas_ana_u(sphere_dat, ugrid, U, nf, cscale)
 type (sphere), intent(in) :: sphere_dat
 integer, intent(in) :: nf
 real, intent(in) :: cscale
 real, dimension(sphere_dat%nlons*sphere_dat%nlats*nf), intent(in) :: ugrid
 complex, dimension(0:sphere_dat%ntrunc, 0:sphere_dat%ntrunc, nf), intent(inout) :: U
 integer :: nlats,nlons,ntrunc,mwaves,nmstrt,n,j,f,nl2,base,lenm,nmdim,nhs
 integer :: i2,off,goff,ldb,ne,no,c0,cs,basen,bases,k,mm
 real :: wj, a1, a2
 real(8) :: tt
 nlats  = sphere_dat%nlats
 nlons  = sphere_dat%nlons
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 nmdim  = (ntrunc+1)*(ntrunc+2)/2
 nhs    = sphere_dat%nhs
 nl2    = nlons+2
 ldb    = 2*nf*mwaves
 call sht_check(sphere_dat,nf)
 call sht_tic(tt)
 !$omp parallel do collapse(2) schedule(static) default(shared) private(f,j,off,goff,base)
 do f = 1,nf
   do j = 1,nlats
     off  = (f-1)*nlats*nl2
     goff = (f-1)*nlons*nlats
     base = off + (j-1)*nl2
     sphere_dat%wrk1(base+1:base+nlons) = ugrid(goff+(j-1)*nlons+1 : goff+(j-1)*nlons+nlons)
     sphere_dat%wrk1(base+nlons+1) = 0.0
     sphere_dat%wrk1(base+nlons+2) = 0.0
   enddo
 enddo
 !$omp end parallel do
 call sht_toc(1,tt)
 call sht_fft(sphere_dat, sphere_dat%wrk1, sphere_dat%wrk2, nlats*nf, -1)
 call sht_toc(2,tt)
 if (sht_sym) then
   !$omp parallel do schedule(static) default(shared) &
   !$omp   private(j,wj,c0,f,basen,bases,cs,i2,a1,a2,k)
   do j = 1,nhs
     wj = sphere_dat%weights(j)*sht_fsc*cscale
     c0 = (j-1)*ldb
     do f = 1,nf
       basen = (f-1)*nlats*nl2 + (j-1)*nl2
       bases = (f-1)*nlats*nl2 + (nlats-j)*nl2
       cs = c0 + 2*(f-1)
       do i2 = 1, 2*mwaves
         a1 = sphere_dat%wrk1(basen+i2)
         a2 = sphere_dat%wrk1(bases+i2)
         k = cs + ((i2+1)/2-1)*2*nf + (2-mod(i2,2))
         sphere_dat%bbuf(k)  = (a1+a2)*wj
         sphere_dat%bbuf2(k) = (a1-a2)*wj
       enddo
     enddo
   enddo
   !$omp end parallel do
 else
   !$omp parallel do schedule(static) default(shared) &
   !$omp   private(j,wj,c0,f,basen,cs,i2,k)
   do j = 1,nlats
     wj = sphere_dat%weights(j)*sht_fsc*cscale
     c0 = (j-1)*ldb
     do f = 1,nf
       basen = (f-1)*nlats*nl2 + (j-1)*nl2
       cs = c0 + 2*(f-1)
       do i2 = 1, 2*mwaves
         k = cs + ((i2+1)/2-1)*2*nf + (2-mod(i2,2))
         sphere_dat%bbuf(k) = sphere_dat%wrk1(basen+i2)*wj
       enddo
     enddo
   enddo
   !$omp end parallel do
 endif
 call sht_toc(3,tt)
 !$omp parallel do schedule(dynamic,1) default(shared) &
 !$omp   private(j,mm,lenm,nmstrt,ne,no,f,n,off,goff)
 do j = 1,mwaves
   mm     = j-1
   lenm   = mwaves-j+1
   nmstrt = sphere_dat%nmstrt_tab(j)
   block
     real :: yb(2*mwaves*nf)
     if (sht_sym) then
       ne = (lenm+1)/2
       no = lenm-ne
       call sht_gemm('N','T', ne, 2*nf, nhs, sphere_dat%pa, nmstrt+1, nmdim, &
                     sphere_dat%bbuf, (j-1)*2*nf+1, ldb, yb, 1, mwaves)
       if (no > 0) &
         call sht_gemm('N','T', no, 2*nf, nhs, sphere_dat%pa, nmstrt+ne+1, nmdim, &
                       sphere_dat%bbuf2, (j-1)*2*nf+1, ldb, yb, ne+1, mwaves)
     else
       call sht_gemm('N','T', lenm, 2*nf, nhs, sphere_dat%pa, nmstrt+1, nmdim, &
                     sphere_dat%bbuf, (j-1)*2*nf+1, ldb, yb, 1, mwaves)
     endif
     do f = 1,nf
       off  = (2*f-2)*mwaves
       goff = (2*f-1)*mwaves
       do n = 1,lenm
         U(mm+sphere_dat%parperm(nmstrt+n)-1, mm, f) = cmplx(yb(off+n), yb(goff+n))
       enddo
     enddo
   end block
 enddo
 !$omp end parallel do
 call sht_toc(4,tt)
 end subroutine spharm_blas_ana_u

 !=====================================================================
 ! synthesis reading straight from the (0:L,0:L) layout.
 ! The legacy path first built Uvec(i) = U(n,m)*1.41421356 in a separate
 ! pass; here that multiply happens in the pack loop that had to read the
 ! coefficients anyway, so the result is identical and one full pass
 ! over the 2.1 MB spectrum per synthesis disappears.
 !=====================================================================
 subroutine spharm_blas_syn_u(sphere_dat, ugrid, U, nf, cscale)
 type (sphere), intent(in) :: sphere_dat
 integer, intent(in) :: nf
 real, intent(in) :: cscale
 real, dimension(sphere_dat%nlons*sphere_dat%nlats*nf), intent(out) :: ugrid
 complex, dimension(0:sphere_dat%ntrunc, 0:sphere_dat%ntrunc, nf), intent(in) :: U
 integer :: nlats,nlons,ntrunc,mwaves,nmstrt,n,j,f,nl2,base,lenm,nmdim,nhs
 integer :: i2,off,goff,ldc,ne,no,basen,bases,cofs,m,mm
 real(8) :: tt
 real :: se, so
 nlats  = sphere_dat%nlats
 nlons  = sphere_dat%nlons
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 nmdim  = (ntrunc+1)*(ntrunc+2)/2
 nhs    = sphere_dat%nhs
 nl2    = nlons+2
 ldc    = nhs
 call sht_check(sphere_dat,nf)
 call sht_tic(tt)
 !$omp parallel do schedule(dynamic,1) default(shared) &
 !$omp   private(m,mm,lenm,nmstrt,f,n,off,goff,cofs,ne,no)
 do m = 1,mwaves
   mm     = m-1
   lenm   = mwaves-m+1
   nmstrt = sphere_dat%nmstrt_tab(m)
   block
     real :: yb(2*mwaves*nf)
     complex :: uu
     do f = 1,nf
       off  = (2*f-2)*mwaves
       goff = (2*f-1)*mwaves
       do n = 1,lenm
         uu = U(mm+sphere_dat%parperm(nmstrt+n)-1, mm, f)*cscale
         yb(off+n)  = real(uu)
         yb(goff+n) = aimag(uu)
       enddo
     enddo
     cofs = (m-1)*2*nf*nhs + 1
     if (sht_sym) then
       ne = (lenm+1)/2
       no = lenm-ne
       call sht_gemm('N','N', nhs, 2*nf, ne, sphere_dat%ps, (nmstrt)*nhs+1, nhs, &
                     yb, 1, mwaves, sphere_dat%bbuf, cofs, ldc)
       if (no > 0) then
         call sht_gemm('N','N', nhs, 2*nf, no, sphere_dat%ps, (nmstrt+ne)*nhs+1, nhs, &
                       yb, ne+1, mwaves, sphere_dat%bbuf2, cofs, ldc)
       else
         sphere_dat%bbuf2(cofs:cofs+2*nf*nhs-1) = 0.0
       endif
     else
       call sht_gemm('N','N', nhs, 2*nf, lenm, sphere_dat%ps, (nmstrt)*nhs+1, nhs, &
                     yb, 1, mwaves, sphere_dat%bbuf, cofs, ldc)
     endif
   end block
 enddo
 !$omp end parallel do
 call sht_toc(6,tt)
 if (sht_sym) then
   !$omp parallel do collapse(2) schedule(static) default(shared) &
   !$omp   private(f,j,basen,bases,i2,cofs,se,so)
   do f = 1,nf
     do j = 1,nhs
       basen = (f-1)*nlats*nl2 + (j-1)*nl2
       bases = (f-1)*nlats*nl2 + (nlats-j)*nl2
       do i2 = 1,2*mwaves
         cofs = ((i2+1)/2-1)*2*nf*nhs + (2*f-2+(2-mod(i2,2))-1)*nhs + j
         se = sphere_dat%bbuf(cofs)
         so = sphere_dat%bbuf2(cofs)
         sphere_dat%wrk1(basen+i2) = se + so
         sphere_dat%wrk1(bases+i2) = se - so
       enddo
       do i2 = 2*mwaves+1, nl2
         sphere_dat%wrk1(basen+i2) = 0.0
         sphere_dat%wrk1(bases+i2) = 0.0
       enddo
     enddo
   enddo
   !$omp end parallel do
 else
   !$omp parallel do collapse(2) schedule(static) default(shared) &
   !$omp   private(f,j,basen,i2,cofs)
   do f = 1,nf
     do j = 1,nlats
       basen = (f-1)*nlats*nl2 + (j-1)*nl2
       do i2 = 1,2*mwaves
         cofs = ((i2+1)/2-1)*2*nf*nhs + (2*f-2+(2-mod(i2,2))-1)*nhs + j
         sphere_dat%wrk1(basen+i2) = sphere_dat%bbuf(cofs)
       enddo
       do i2 = 2*mwaves+1, nl2
         sphere_dat%wrk1(basen+i2) = 0.0
       enddo
     enddo
   enddo
   !$omp end parallel do
 endif
 call sht_toc(7,tt)
 call sht_fft(sphere_dat, sphere_dat%wrk1, sphere_dat%wrk2, nlats*nf, 1)
 call sht_toc(8,tt)
 !$omp parallel do collapse(2) schedule(static) default(shared) private(f,j,off,goff,base)
 do f = 1,nf
   do j = 1,nlats
     off  = (f-1)*nlats*nl2
     goff = (f-1)*nlons*nlats
     base = off + (j-1)*nl2
     ugrid(goff+(j-1)*nlons+1 : goff+(j-1)*nlons+nlons) = sphere_dat%wrk1(base+1:base+nlons)
   enddo
 enddo
 !$omp end parallel do
 call sht_toc(5,tt)
 end subroutine spharm_blas_syn_u

 subroutine sht_tic(tt)
   real(8), intent(out) :: tt
   integer(8) :: c,r
   if (.not.sht_profile) then; tt=0.0d0; return; endif
   call system_clock(c,r); tt = real(c,8)/real(r,8)
 end subroutine sht_tic
 subroutine sht_toc(slot,tt)
   integer, intent(in) :: slot
   real(8), intent(inout) :: tt
   integer(8) :: c,r
   real(8) :: t2
   if (.not.sht_profile) return
   call system_clock(c,r); t2 = real(c,8)/real(r,8)
   sht_times(slot) = sht_times(slot) + (t2-tt)
   sht_counts(slot) = sht_counts(slot) + 1
   tt = t2
 end subroutine sht_toc

 !=====================================================================
 ! synthesis:  spectral -> grid, nf fields at once.
 !=====================================================================
 subroutine spharm_blas_syn(sphere_dat, ugrid, anm, nf)
 type (sphere), intent(in) :: sphere_dat
 integer, intent(in) :: nf
 real, dimension(sphere_dat%nlons*sphere_dat%nlats*nf), intent(out) :: ugrid
 complex, dimension(((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2)*nf), intent(in) :: anm
 integer :: nlats,nlons,ntrunc,mwaves,nmstrt,n,j,f,nl2,base,lenm,nmdim,nhs
 integer :: i2,off,goff,ldc,ne,no,basen,bases,cofs,m
 real(8) :: tt
 real :: se, so
 nlats  = sphere_dat%nlats
 nlons  = sphere_dat%nlons
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 nmdim  = (ntrunc+1)*(ntrunc+2)/2
 nhs    = sphere_dat%nhs
 nl2    = nlons+2
 ldc    = nhs
 call sht_check(sphere_dat,nf)
 call sht_tic(tt)
 !$omp parallel do schedule(dynamic,1) default(shared) &
 !$omp   private(m,lenm,nmstrt,f,n,off,goff,cofs,ne,no)
 do m = 1,mwaves
   lenm   = mwaves-m+1
   nmstrt = sphere_dat%nmstrt_tab(m)
   block
     real :: yb(2*mwaves*nf)
     do f = 1,nf
       off  = (2*f-2)*mwaves
       goff = (2*f-1)*mwaves
       do n = 1,lenm
         yb(off+n)  = real(anm((f-1)*nmdim+nmstrt+sphere_dat%parperm(nmstrt+n)))
         yb(goff+n) = aimag(anm((f-1)*nmdim+nmstrt+sphere_dat%parperm(nmstrt+n)))
       enddo
     enddo
     cofs = (m-1)*2*nf*nhs + 1
     if (sht_sym) then
       ne = (lenm+1)/2
       no = lenm-ne
       call sht_gemm('N','N', nhs, 2*nf, ne, sphere_dat%ps, (nmstrt)*nhs+1, nhs, &
                     yb, 1, mwaves, sphere_dat%bbuf, cofs, ldc)
       if (no > 0) then
         call sht_gemm('N','N', nhs, 2*nf, no, sphere_dat%ps, (nmstrt+ne)*nhs+1, nhs, &
                       yb, ne+1, mwaves, sphere_dat%bbuf2, cofs, ldc)
       else
         sphere_dat%bbuf2(cofs:cofs+2*nf*nhs-1) = 0.0
       endif
     else
       call sht_gemm('N','N', nhs, 2*nf, lenm, sphere_dat%ps, (nmstrt)*nhs+1, nhs, &
                     yb, 1, mwaves, sphere_dat%bbuf, cofs, ldc)
     endif
   end block
 enddo
 !$omp end parallel do
 call sht_toc(6,tt)
 ! ---- de-interleave into the fft buffer --------------------------------
 if (sht_sym) then
   !$omp parallel do collapse(2) schedule(static) default(shared) &
   !$omp   private(f,j,basen,bases,i2,cofs,se,so)
   do f = 1,nf
     do j = 1,nhs
       basen = (f-1)*nlats*nl2 + (j-1)*nl2
       bases = (f-1)*nlats*nl2 + (nlats-j)*nl2
       do i2 = 1,2*mwaves
         cofs = ((i2+1)/2-1)*2*nf*nhs + (2*f-2+(2-mod(i2,2))-1)*nhs + j
         se = sphere_dat%bbuf(cofs)
         so = sphere_dat%bbuf2(cofs)
         sphere_dat%wrk1(basen+i2) = se + so
         sphere_dat%wrk1(bases+i2) = se - so
       enddo
       do i2 = 2*mwaves+1, nl2
         sphere_dat%wrk1(basen+i2) = 0.0
         sphere_dat%wrk1(bases+i2) = 0.0
       enddo
     enddo
   enddo
   !$omp end parallel do
 else
   !$omp parallel do collapse(2) schedule(static) default(shared) &
   !$omp   private(f,j,basen,i2,cofs)
   do f = 1,nf
     do j = 1,nlats
       basen = (f-1)*nlats*nl2 + (j-1)*nl2
       do i2 = 1,2*mwaves
         cofs = ((i2+1)/2-1)*2*nf*nhs + (2*f-2+(2-mod(i2,2))-1)*nhs + j
         sphere_dat%wrk1(basen+i2) = sphere_dat%bbuf(cofs)
       enddo
       do i2 = 2*mwaves+1, nl2
         sphere_dat%wrk1(basen+i2) = 0.0
       enddo
     enddo
   enddo
   !$omp end parallel do
 endif
 call sht_toc(7,tt)
 call sht_fft(sphere_dat, sphere_dat%wrk1, sphere_dat%wrk2, nlats*nf, 1)
 call sht_toc(8,tt)
 !$omp parallel do collapse(2) schedule(static) default(shared) private(f,j,off,goff,base)
 do f = 1,nf
   do j = 1,nlats
     off  = (f-1)*nlats*nl2
     goff = (f-1)*nlons*nlats
     base = off + (j-1)*nl2
     ugrid(goff+(j-1)*nlons+1 : goff+(j-1)*nlons+nlons) = sphere_dat%wrk1(base+1:base+nlons)
   enddo
 enddo
 !$omp end parallel do
 call sht_toc(5,tt)
 end subroutine spharm_blas_syn


 !=====================================================================
 ! exact transform.
 !
 ! Same operations, same operand order and the same summation order as
 ! spharm_ref, so the result is identical; only the loop
 ! nesting and the placement of the accumulators change.
 !
 ! Reference analysis:  do j / do m / do n :
 !        anm(nm) = anm(nm) + pnm(nm,j)*weights(j)*am(m,j)
 ! i.e. anm(nm) accumulates over j in increasing j, and the operand is
 ! (pnm*weights) (a real) times the complex Fourier coefficient.  The
 ! pnm*weights product is precomputed once into the `pa` table (identical
 ! IEEE product), and the complex accumulator is split into two real
 ! arrays so the inner loop is a pair of FMA streams.
 !=====================================================================
 subroutine spharm_exact_ana(sphere_dat, ugrid, anm, nf)
 type (sphere), intent(in) :: sphere_dat
 integer, intent(in) :: nf
 real, dimension(sphere_dat%nlons*sphere_dat%nlats*nf), intent(in) :: ugrid
 complex, dimension(((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2)*nf), intent(inout) :: anm
 integer :: nlats,nlons,ntrunc,mwaves,nmstrt,n,j,f,nl2,base,lenm,nmdim
 integer :: i2,off,goff,ldb,c0,cs,fo,bo
 real(8) :: tt
 nlats  = sphere_dat%nlats
 nlons  = sphere_dat%nlons
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 nmdim  = (ntrunc+1)*(ntrunc+2)/2
 nl2    = nlons+2
 ldb    = 2*nf*mwaves
 call sht_check(sphere_dat,nf)
 call sht_tic(tt)
 !$omp parallel do collapse(2) schedule(static) default(shared) private(f,j,off,goff,base)
 do f = 1,nf
   do j = 1,nlats
     off  = (f-1)*nlats*nl2
     goff = (f-1)*nlons*nlats
     base = off + (j-1)*nl2
     sphere_dat%wrk1(base+1:base+nlons) = ugrid(goff+(j-1)*nlons+1 : goff+(j-1)*nlons+nlons)
     sphere_dat%wrk1(base+nlons+1) = 0.0
     sphere_dat%wrk1(base+nlons+2) = 0.0
   enddo
 enddo
 !$omp end parallel do
 call sht_toc(1,tt)
 call sht_fft(sphere_dat, sphere_dat%wrk1, sphere_dat%wrk2, nlats*nf, -1)
 call sht_toc(2,tt)
 ! interleave the Fourier coefficients: bbuf(2*nf*mwaves, nlats)
 !$omp parallel do schedule(static) default(shared) private(j,c0,f,base,cs,i2)
 do j = 1,nlats
   c0 = (j-1)*ldb
   do f = 1,nf
     base = (f-1)*nlats*nl2 + (j-1)*nl2
     cs = c0 + 2*(f-1)
     do i2 = 1, 2*mwaves
       sphere_dat%bbuf(cs + ((i2+1)/2-1)*2*nf + (2-mod(i2,2))) = sphere_dat%wrk1(base+i2)*sht_fsc
     enddo
   enddo
 enddo
 !$omp end parallel do
 call sht_toc(3,tt)
 call sht_ana_kernel(sphere_dat%pab, sphere_dat%nmpad, nmdim, nlats, mwaves, nf, ldb, &
                     sphere_dat%bbuf, sphere_dat%nmstrt_tab, sphere_dat%nmpad_tab, anm)
 call sht_toc(4,tt)
 end subroutine spharm_exact_ana

 subroutine sht_ana_kernel(pab, nmpad, nmdim, nlats, mwaves, nf, ldb, bbuf, &
                           nmstrt_tab, nmpad_tab, anm)
 ! Coefficient-block sweep.  SHT_KB consecutive coefficients of one zonal
 ! wavenumber are accumulated in registers while their Legendre rows stream
 ! contiguously across all latitudes; the 2*SHT_KB dependency chains hide
 ! the FMA latency.  Each coefficient is still summed over increasing
 ! latitude with the very same pnm*weights multiplier, hence identical
 ! to the reference triple loop.
 integer, intent(in) :: nmpad, nmdim, nlats, mwaves, nf, ldb
 real, intent(in) :: pab(SHT_KB*nlats, nmpad/SHT_KB)
 real, intent(in) :: bbuf(*)
 integer, intent(in) :: nmstrt_tab(*), nmpad_tab(*)
 complex, intent(inout) :: anm(*)
 integer :: m,n,j,f,lenm,nm0,bidx,bb,bk,kk,base,nleft
 real :: sr(SHT_KB), si(SHT_KB), amr, ami
 do f = 1,nf
   !$omp parallel do schedule(dynamic,1) default(shared) &
   !$omp   private(m,n,j,lenm,nm0,bidx,bb,bk,kk,base,nleft,sr,si,amr,ami)
   do m = 1,mwaves
     lenm = mwaves-m+1
     nm0  = nmstrt_tab(m)
     do bb = 1,(lenm+SHT_KB-1)/SHT_KB
       bk = nmpad_tab(m)/SHT_KB + bb
       do kk = 1,SHT_KB
         sr(kk) = 0.0
         si(kk) = 0.0
       enddo
       do j = 1,nlats
         bidx = (j-1)*ldb + (m-1)*2*nf + 2*(f-1)
         amr  = bbuf(bidx+1)
         ami  = bbuf(bidx+2)
         base = (j-1)*SHT_KB
         do kk = 1,SHT_KB
           sr(kk) = sr(kk) + pab(base+kk,bk)*amr
           si(kk) = si(kk) + pab(base+kk,bk)*ami
         enddo
       enddo
       nleft = min(SHT_KB, lenm-(bb-1)*SHT_KB)
       do kk = 1,nleft
         n = (bb-1)*SHT_KB + kk
         anm((f-1)*nmdim+nm0+n) = cmplx(sr(kk), si(kk))
       enddo
     enddo
   enddo
   !$omp end parallel do
 enddo
 end subroutine sht_ana_kernel

 !---------------------------------------------------------------------
 ! Reference synthesis:  do j / do m / do n :
 !        am(m,j) = am(m,j) + anm(nm)*pnm(nm,j)
 ! accumulating over n in increasing n, starting from (0,0).
 !---------------------------------------------------------------------
 subroutine spharm_exact_syn(sphere_dat, ugrid, anm, nf)
 type (sphere), intent(in) :: sphere_dat
 integer, intent(in) :: nf
 real, dimension(sphere_dat%nlons*sphere_dat%nlats*nf), intent(out) :: ugrid
 complex, dimension(((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2)*nf), intent(in) :: anm
 integer :: nlats,nlons,ntrunc,mwaves,n,j,f,nl2,base,nmdim
 integer :: i2,off,goff,ldb,c0,cs
 real(8) :: tt
 nlats  = sphere_dat%nlats
 nlons  = sphere_dat%nlons
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 nmdim  = (ntrunc+1)*(ntrunc+2)/2
 nl2    = nlons+2
 ldb    = 2*nf*mwaves
 call sht_check(sphere_dat,nf)
 call sht_tic(tt)
 call sht_syn_kernel(sphere_dat%psb, nmdim, nlats, mwaves, nf, ldb, &
                     sphere_dat%bbuf, sphere_dat%nmstrt_tab, sphere_dat%ybuf, anm)
 call sht_toc(6,tt)
 !$omp parallel do schedule(static) default(shared) private(j,c0,f,base,cs,i2)
 do j = 1,nlats
   c0 = (j-1)*ldb
   do f = 1,nf
     base = (f-1)*nlats*nl2 + (j-1)*nl2
     cs = c0 + 2*(f-1)
     do i2 = 1, 2*mwaves
       sphere_dat%wrk1(base+i2) = sphere_dat%bbuf(cs + ((i2+1)/2-1)*2*nf + (2-mod(i2,2)))
     enddo
     do i2 = 2*mwaves+1, nl2
       sphere_dat%wrk1(base+i2) = 0.0
     enddo
   enddo
 enddo
 !$omp end parallel do
 call sht_toc(7,tt)
 call sht_fft(sphere_dat, sphere_dat%wrk1, sphere_dat%wrk2, nlats*nf, 1)
 call sht_toc(8,tt)
 !$omp parallel do collapse(2) schedule(static) default(shared) private(f,j,off,goff,base)
 do f = 1,nf
   do j = 1,nlats
     off  = (f-1)*nlats*nl2
     goff = (f-1)*nlons*nlats
     base = off + (j-1)*nl2
     ugrid(goff+(j-1)*nlons+1 : goff+(j-1)*nlons+nlons) = sphere_dat%wrk1(base+1:base+nlons)
   enddo
 enddo
 !$omp end parallel do
 call sht_toc(5,tt)
 end subroutine spharm_exact_syn

 subroutine sht_syn_kernel(ps, nmdim, nlats, mwaves, nf, ldb, bbuf, nmstrt_tab, ybuf, anm)
 ! Synthesis is a reduction over n, so the exact summation order forbids
 ! splitting the accumulator.  Instead JB latitudes are carried through the
 ! n loop at once: JB independent dependency chains hide the FMA latency
 ! while every individual chain still adds in increasing n.  The Legendre
 ! table is stored transposed, ps(nlat, nmdim), so the JB values needed for
 ! one n sit in a single cache line.
 integer, intent(in) :: nmdim, nlats, mwaves, nf, ldb
 real, intent(in) :: ps(SHT_JB*nmdim, nlats/SHT_JB)
 real, intent(inout) :: bbuf(*)
 integer, intent(in) :: nmstrt_tab(*)
 real, intent(inout) :: ybuf(*)
 complex, intent(in) :: anm(*)
 integer, parameter :: JB = SHT_JB
 integer :: m,n,j,jb0,jj,f,lenm,nm0,bidx,yr0,yi0,blk
 real :: sr(JB), si(JB), ar, ai
 do f = 1,nf
   yr0 = (f-1)*2*nmdim
   yi0 = yr0 + nmdim
   do n = 1,nmdim
     ybuf(yr0+n) = real(anm((f-1)*nmdim+n))
     ybuf(yi0+n) = aimag(anm((f-1)*nmdim+n))
   enddo
 enddo
 if (mod(nlats,JB) /= 0) then
   print *,'sht_syn_kernel: nlats must be a multiple of', JB; stop
 endif
 !$omp parallel do schedule(static) default(shared) &
 !$omp   private(m,n,j,jb0,jj,f,lenm,nm0,bidx,yr0,yi0,blk,sr,si,ar,ai)
 do jb0 = 1,nlats,JB
   blk = (jb0-1)/JB + 1
   do f = 1,nf
     yr0 = (f-1)*2*nmdim
     yi0 = yr0 + nmdim
     do m = 1,mwaves
       lenm = mwaves-m+1
       nm0  = nmstrt_tab(m)
       do jj = 1,JB
         sr(jj) = 0.0
         si(jj) = 0.0
       enddo
       do n = 1,lenm
         ar = ybuf(yr0+nm0+n)
         ai = ybuf(yi0+nm0+n)
         do jj = 1,JB
           sr(jj) = sr(jj) + ar*ps((nm0+n-1)*JB+jj, blk)
           si(jj) = si(jj) + ai*ps((nm0+n-1)*JB+jj, blk)
         enddo
       enddo
       do jj = 1,JB
         bidx = (jb0+jj-2)*ldb + (m-1)*2*nf + 2*(f-1)
         bbuf(bidx+1) = sr(jj)
         bbuf(bidx+2) = si(jj)
       enddo
     enddo
   enddo
 enddo
 !$omp end parallel do
 end subroutine sht_syn_kernel
 !---- Temperton fft991 split over threads (each chunk is independent) ----
 subroutine fft991_mt(a, work, trigs, ifax, inc, jump, n, lot, isign)
   real :: a(*), work(*)
   real, intent(in) :: trigs(:)
   integer, intent(in) :: ifax(:)
   integer, intent(in) :: inc, jump, n, lot, isign
   integer :: nt, c, lo, hi, chunk, off
#ifdef _OPENMP
   integer :: omp_get_max_threads
   nt = omp_get_max_threads()
#else
   nt = 1
#endif
   if (nt <= 1 .or. lot < 2*nt) then
     call fft991(a, work, trigs, ifax, inc, jump, n, lot, isign)
     return
   endif
   chunk = (lot + nt - 1)/nt
   !$omp parallel do schedule(static) private(c,lo,hi,off)
   do c = 1, nt
     lo = (c-1)*chunk + 1
     hi = min(c*chunk, lot)
     if (lo <= hi) then
       off = (lo-1)*jump
       call fft991(a(off+1), work(off+1), trigs, ifax, inc, jump, n, hi-lo+1, isign)
     endif
   enddo
   !$omp end parallel do
 end subroutine fft991_mt

 !=======================================================================
 ! FFT layer.
 !
 ! `fft991` (Temperton) and FFTW3's real transform agree to the level of
 ! the algorithms' own rounding: for nlon = 1024, lot = 512 random rows the
 ! largest difference is 1.4e-17 forward (coefficients ~3e-2) and 5.6e-16
 ! backward (grid ~5e-1), i.e. pure round-off.  The layouts are identical:
 ! fft991 with inc=1, jump=nlon+2 stores (Re,Im) interleaved for
 ! m = 0..nlon/2, which is exactly FFTW's in-place r2c half-complex array.
 ! The only convention difference is the 1/nlon factor fft991 folds into the
 ! forward transform; it is applied through sht_fsc (exactly, nlon = 2^10).
 !=======================================================================
 integer function sht_nt() result(nt)
#ifdef _OPENMP
   integer :: omp_get_max_threads
   nt = omp_get_max_threads()
#else
   nt = 1
#endif
 end function sht_nt

 subroutine sht_fft(sd, a, work, lot, isign)
   type (sphere), intent(in) :: sd
   real :: a(*), work(*)
   integer, intent(in) :: lot, isign
   integer :: nt, base, rem, c, lo, h, off, ip, nl2, nfl
   nl2 = sd%nlons+2
#ifndef HAVE_FFTW
   call fft991_mt(a, work, sd%trigs, sd%ifax, 1, nl2, sd%nlons, lot, isign)
   return
#else
   if (sht_fftw == 0) then
     call fft991_mt(a, work, sd%trigs, sd%ifax, 1, nl2, sd%nlons, lot, isign)
     return
   endif
   nt = sd%nplan
   nfl = lot/sd%nlats
   if (nfl < 1 .or. nfl > SHT_NFMAX .or. nfl*sd%nlats /= lot) then
     print *,'sht_fft: unexpected lot',lot; stop
   endif
   base = lot/nt; rem = mod(lot,nt)
   !$omp parallel do schedule(static) num_threads(nt) default(shared) private(c,lo,h,off,ip) if(nt>1)
   do c = 1,nt
     if (c <= rem) then
        h = base+1; lo = (c-1)*(base+1)+1
     else
        h = base;   lo = rem*(base+1) + (c-1-rem)*base + 1
     endif
     if (h > 0) then
       off = (lo-1)*nl2
       ip  = (nfl-1)*nt + c
       if (isign == -1) then
         call dfftw_execute_dft_r2c(sd%planf(ip), a(off+1), a(off+1))
       else
         call dfftw_execute_dft_c2r(sd%planb(ip), a(off+1), a(off+1))
       endif
     endif
   enddo
   !$omp end parallel do
#endif
 end subroutine sht_fft

 !---- runtime configuration from the environment ------------------------
 subroutine sht_env_init()
   character(len=64) :: v
   integer :: st
#ifdef HAVE_FFTW
   sht_fftw = 1
#endif
   call get_environment_variable('SLEVOLVE_FFT', v, status=st)
   if (st == 0) then
     if (trim(v) == 'fftw') then
       sht_fftw = 1
     else if (trim(v) == 'fft991') then
       sht_fftw = 0
     else
       print *,'SLEVOLVE_FFT must be fftw or fft991'; stop
     endif
   endif
#ifndef HAVE_FFTW
   if (sht_fftw /= 0) then
     print *,'this build has no FFTW (rebuild with FFT=fftw)'; stop
   endif
#endif
   call get_environment_variable('SHT_UNGUARDED', v, status=st)
   if (st == 0) then
     if (trim(v) == '1') sht_unguarded = .true.
   endif
 end subroutine sht_env_init

 !---- single point of contact with BLAS ---------------------------------
 subroutine sht_gemm(ta,tb,m,n,k,a,ao,lda,b,bo,ldb,c,co,ldc)
   character(len=1), intent(in) :: ta,tb
   integer, intent(in) :: m,n,k,ao,lda,bo,ldb,co,ldc
   real, intent(in) :: a(*), b(*)
   real, intent(inout) :: c(*)
   call dgemm(ta,tb,m,n,k,1.0d0,a(ao),lda,b(bo),ldb,0.0d0,c(co),ldc)
 end subroutine sht_gemm
 subroutine sht_check(sphere_dat,nf)
 type (sphere), intent(in) :: sphere_dat
 integer, intent(in) :: nf
 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in spharm!'; stop
 end if
 if (nf < 1 .or. nf > sphere_dat%nfmax) then
    print *, 'spharm: nf out of range'; stop
 end if
 if (sphere_dat%ntrunc > sphere_dat%nlons/2) then
    print *, 'spharm(blas): requires ntrunc <= nlons/2'; stop
 end if
 end subroutine sht_check
 !---- dispatcher keeping the original single-field interface -------------
 subroutine spharm(sphere_dat, ugrid, anm, idir)
 type (sphere), intent(in) :: sphere_dat
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(inout) :: ugrid
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2), intent(inout) :: anm
 integer, intent(in) :: idir
#if defined(USE_EXACT_SHT)
 if (idir == 1) then
   call spharm_exact_ana(sphere_dat, ugrid, anm, 1)
 else
   call spharm_exact_syn(sphere_dat, ugrid, anm, 1)
 endif
#elif defined(USE_BLAS_SHT)
 if (idir == 1) then
   call spharm_blas_ana(sphere_dat, ugrid, anm, 1)
 else
   call spharm_blas_syn(sphere_dat, ugrid, anm, 1)
 endif
#else
 call spharm_ref(sphere_dat, ugrid, anm, idir)
#endif
 end subroutine spharm
 subroutine getuv(sphere_dat,vrtnm,divnm,ug,vg)
 type (sphere), intent(in) :: sphere_dat
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(out) ::  ug,vg
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2), intent(in) :: vrtnm,divnm
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats)  :: um,vm
 integer :: nlats,ntrunc,mwaves,m,j,n,nm,nmstrt
 real ::  rm
 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in getuv!'
    stop
 end if
 nlats = sphere_dat%nlats
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 do j=1,nlats
    nmstrt = 0
    do m=1,mwaves
       rm = m-1
       um(m,j) = cmplx(0.,0.)
       vm(m,j) = cmplx(0.,0.)
       do n=1,mwaves-m+1
          nm = nmstrt + n
          um(m,j) = um(m,j) + (sphere_dat%ilap(nm)/sphere_dat%rsphere)*( &
                    cmplx(0.,rm)*divnm(nm)*sphere_dat%pnm(nm,j) + &
                    vrtnm(nm)*sphere_dat%hnm(nm,j) )
          vm(m,j) = vm(m,j) + (sphere_dat%ilap(nm)/sphere_dat%rsphere)*( &
                    cmplx(0.,rm)*vrtnm(nm)*sphere_dat%pnm(nm,j) - &
                    divnm(nm)*sphere_dat%hnm(nm,j) )
       enddo
       nmstrt = nmstrt + mwaves-m+1
    enddo
 enddo
 call rfft(sphere_dat, ug, um, -1)
 call rfft(sphere_dat, vg, vm, -1)
 end subroutine getuv
 subroutine getvrtdiv(sphere_dat,vrtnm,divnm,ug,vg)
 type (sphere), intent(in) :: sphere_dat
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(inout) ::  ug,vg
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2), intent(in) :: vrtnm,divnm
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats) :: um,vm
 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in getvrtdiv!'
    stop
 end if
 call rfft(sphere_dat, ug, um, 1)
 call rfft(sphere_dat, vg, vm, 1)
 call sumnm(sphere_dat,um,vm,divnm,1,1)
 call sumnm(sphere_dat,vm,um,vrtnm,1,-1)
 end subroutine getvrtdiv
 subroutine sumnm(sphere_dat,am,bm,anm,isign1,isign2)
 type (sphere), intent(in) :: sphere_dat
 integer, intent(in) :: isign1,isign2
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2) :: anm
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats), intent(in) :: am,bm
 integer :: nlats,ntrunc,mwaves,j,m,n,nm,nmstrt
 real ::  sign1,sign2,rm
 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in sumnm!'
    stop
 end if
 sign1 = float(isign1)
 sign2 = float(isign2)
 if (isign2 .ne. 1 .and. isign2 .ne. -1) then
    print *, ' isign2 must be +1 or -1 in sumnm!'
    print *, ' execution terminated in sumnm'
 end if
 if (isign1 .ne. 1 .and. isign1 .ne. -1) then
    print *, ' isign1 must be +1 or -1 in sumnm!'
    print *, ' execution terminated in sumnm'
    stop
 end if
 nlats = sphere_dat%nlats
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 anm = 0.
 DO J=1,NLATS
    NMSTRT = 0
    DO m = 1, MWAVES
       RM = M-1
       DO n   = 1, mwaves-m+1
          NM = NMSTRT + n
          aNM(NM) = aNM(NM) + sign1*sphere_dat%GWrC(j)*(CMPLX(0.,RM) &
                    * sphere_dat%PNM(NM,J) * am(m,j) &
                    + sign2 * sphere_dat%HNM(NM,J) * bm(m,j))
       enddo
       nmstrt = nmstrt + mwaves - m +1
    enddo
 enddo
 end subroutine sumnm
 subroutine cosgrad(sphere_dat,divnm,ug,vg)
 type (sphere), intent(in) :: sphere_dat
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(out) ::  ug,vg
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2), intent(in) :: divnm
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats) :: um,vm
 integer :: nlats,ntrunc,mwaves,j,m,n,nm,nmstrt
 real :: rm
 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in cosgrad!'
    stop
 end if
 nlats = sphere_dat%nlats
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 do j=1,nlats
    nmstrt = 0
    do m=1,mwaves
       rm = (m-1)
       um(m,j) = cmplx(0.,0.)
       vm(m,j) = cmplx(0.,0.)
       do n=1,mwaves-m+1
          nm = nmstrt + n
          um(m,j) = um(m,j) + (1./sphere_dat%rsphere)* &
                    cmplx(0.,rm)*divnm(nm)*sphere_dat%pnm(nm,j)
          vm(m,j) = vm(m,j) - (1./sphere_dat%rsphere)* &
                    divnm(nm)*sphere_dat%hnm(nm,j)
       enddo
       nmstrt = nmstrt + mwaves - m +1
    enddo
 enddo
 call rfft(sphere_dat, ug, um, -1)
 call rfft(sphere_dat, vg, vm, -1)
 end subroutine cosgrad
 subroutine specsmooth(sphere_dat,datagrid,smooth)
 type (sphere), intent(in) :: sphere_dat
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(inout) :: datagrid
 real, dimension(sphere_dat%ntrunc+1), intent(in) :: smooth
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2) :: dataspec
 integer :: n,nm,nmdim
 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in specsmooth!'
    stop
 end if
 nmdim = (sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2
 call spharm(sphere_dat, datagrid, dataspec, 1)
 do nm=1,nmdim
    n = sphere_dat%indxn(nm)
    dataspec(nm) = dataspec(nm)*smooth(n+1)
 enddo
 call spharm(sphere_dat, datagrid, dataspec, -1)
 end subroutine specsmooth
 subroutine fft99 (a,work,trigs,ifax,inc,jump,n,lot,isign)
    integer, intent(in)    :: inc,jump,n,lot,isign
    integer, intent(in) :: ifax(:)
    real,    intent(in)    :: trigs(:)
    real,    intent(inout) :: a(*),work(*)
    integer :: nfax, nx, nh, ink, igo, ibase, jbase
    integer :: i, j, k, L, m, ia, la, ib
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=inc+1
      jbase=1
    do L=1,lot
      i=ibase
      j=jbase
      do m=1,n
        work(j)=a(i)
        i=i+inc
        j=j+1
      enddo
      ibase=ibase+jump
      jbase=jbase+nx
    enddo
      igo=60
      go to 40
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
   40 continue
      ia=inc+1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm (a(ia),a(ia+inc),work(1),work(2),trigs, &
                   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm (work(1),work(2),a(ia),a(ia+inc),trigs, &
                   2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
    if (isign.eq.-1) go to 130
    if (mod(nfax,2).ne.1) then
      ibase=1
      jbase=ia
      do L=1,lot
        i=ibase
        j=jbase
        do m=1,n
          a(j)=work(i)
          i=i+1
          j=j+inc
        enddo
        ibase=ibase+nx
        jbase=jbase+jump
      enddo
    endif
      ia=1
      ib=n*inc+1
      do L=1,lot
        a(ia)=a(ib)
        a(ib+inc)=a(ia+inc)
        ia=ia+jump
        ib=ib+jump
      enddo
      go to 140
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
  140 continue
 end subroutine fft99
 subroutine fft99a (a,work,trigs,inc,jump,n,lot)
    integer, intent(in)    :: inc,jump,n,lot
    real,    intent(in)    :: trigs(:)
    real,    intent(inout) :: a(*),work(*)
    integer :: nh, nx, ink, k, L
    integer :: ia, ib, ja, jb, iabase, ibbase, jabase, jbbase
    real    :: c, s
      nh=n/2
      nx=n+1
      ink=inc+inc
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
    do L=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
    enddo
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1
    do k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
      do L=1,lot
        work(ja)=(a(ia)+a(ib))- &
            (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
        work(jb)=(a(ia)+a(ib))+ &
            (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
        work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+ &
            (a(ia+inc)-a(ib+inc))
        work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))- &
            (a(ia+inc)-a(ib+inc))
        ia=ia+jump
        ib=ib+jump
        ja=ja+nx
        jb=jb+nx
      enddo
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
    enddo
    if (iabase.eq.ibbase) then
      ia=iabase
      ja=jabase
      do L=1,lot
        work(ja)=2.0*a(ia)
        work(ja+1)=-2.0*a(ia+inc)
        ia=ia+jump
        ja=ja+nx
      enddo
    endif
 end subroutine fft99a
 subroutine fft99b (work,a,trigs,inc,jump,n,lot)
    integer, intent(in)    :: inc,jump,n,lot
    real,    intent(in)    :: trigs(:)
    real,    intent(inout) :: a(*),work(*)
    integer :: nh, nx, ink, k, L
    integer :: ia, ib, ja, jb, iabase, ibbase, jabase, jbbase
    real    :: scale, c, s
      nh=n/2
      nx=n+1
      ink=inc+inc
      scale=1.0/real(n)
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
    do L=1,lot
      a(ja)=scale*(work(ia)+work(ib))
      a(jb)=scale*(work(ia)-work(ib))
      a(ja+inc)=0.0
      a(jb+inc)=0.0
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
    enddo
      scale=0.5*scale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1
    do k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
      do L=1,lot
        a(ja)=scale*((work(ia)+work(ib)) &
           +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
        a(jb)=scale*((work(ia)+work(ib)) &
           -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
        a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
            +(work(ib+1)-work(ia+1)))
        a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
            -(work(ib+1)-work(ia+1)))
        ia=ia+nx
        ib=ib+nx
        ja=ja+jump
        jb=jb+jump
      enddo
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
    enddo
    if (iabase.eq.ibbase) then
      ia=iabase
      ja=jabase
      scale=2.0*scale
      do L=1,lot
        a(ja)=scale*work(ia)
        a(ja+inc)=-scale*work(ia+1)
        ia=ia+nx
        ja=ja+jump
      enddo
    endif
 end subroutine fft99b
 subroutine fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
    integer, intent(in)    :: inc,jump,n,lot,isign
    integer, intent(in) :: ifax(:)
    real,    intent(in)    :: trigs(:)
    real,    intent(inout) :: a(*),work(*)
    integer :: nfax, nx, nh, ink, igo, ibase, jbase
    integer :: i, j, k, L, m, ia, la, ib
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=1
      jbase=1
    do L=1,lot
      i=ibase
      j=jbase
      do m=1,n
        work(j)=a(i)
        i=i+inc
        j=j+1
      enddo
      ibase=ibase+jump
      jbase=jbase+nx
    enddo
      igo=60
      go to 40
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
   40 continue
      ia=1
      la=1
    do k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
        call vpassm (a(ia),a(ia+inc),work(1),work(2),trigs, &
                     ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
        call vpassm (work(1),work(2),a(ia),a(ia+inc),trigs, &
                     2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
    enddo
    if (isign.eq.-1) go to 130
    if (mod(nfax,2).ne.1) then
      ibase=1
      jbase=1
      do L=1,lot
        i=ibase
        j=jbase
        do m=1,n
          a(j)=work(i)
          i=i+1
          j=j+inc
        enddo
        ibase=ibase+nx
        jbase=jbase+jump
      enddo
    endif
      ib=n*inc+1
      do L=1,lot
        a(ib)=0.0
        a(ib+inc)=0.0
        ib=ib+jump
      enddo
      go to 140
  130 continue
      call fft99b (work,a,trigs,inc,jump,n,lot)
  140 continue
 end subroutine fft991
 subroutine set99 (trigs, ifax, n)
    integer, intent(in)  :: n
    integer, intent(out) :: ifax(:)
    real,    intent(out) :: trigs(:)
    integer :: mode = 3
    integer :: i
      call fax (ifax, n, mode)
      i = ifax(1)
      if (ifax(i+1) .gt. 5 .or. n .le. 4) ifax(1) = -99
      if (ifax(1) .le. 0 ) then
        write(6,*) ' set99 -- invalid n'
        stop
      endif
      call fftrig (trigs, n, mode)
 end subroutine set99
 subroutine fax (ifax,n,mode)
    integer, intent(out) :: ifax(:)
    integer, intent(in)  :: n, mode
    integer :: nn, k, L, inc, nfax, ii, istop, i, item
      nn=n
      if (iabs(mode).eq.1) go to 10
      if (iabs(mode).eq.8) go to 10
      nn=n/2
      if ((nn+nn).eq.n) go to 10
      ifax(1)=-99
      return
   10 k=1
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
   50 L=5
      inc=2
   60 if (mod(nn,L).ne.0) go to 70
      k=k+1
      ifax(k)=L
      nn=nn/L
      if (nn.eq.1) go to 80
      go to 60
   70 L=L+inc
      inc=6-inc
      go to 60
   80 ifax(1)=k-1
      nfax=ifax(1)
      if (nfax.eq.1) go to 110
      do 100 ii=2,nfax
      istop=nfax+2-ii
      do 90 i=2,istop
      if (ifax(i+1).ge.ifax(i)) go to 90
      item=ifax(i)
      ifax(i)=ifax(i+1)
      ifax(i+1)=item
   90 continue
  100 continue
  110 continue
 end subroutine fax
 subroutine fftrig (trigs,n,mode)
    real,    intent(out) :: trigs(:)
    integer, intent(in)  :: n, mode
    real    :: del, angle,pi
    integer :: imode, nn, nh, i, L, la
      pi = 4.*atan(1.0)
      imode=iabs(mode)
      nn=n
      if (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/real(nn)
      L=nn+nn
      do i=1,L,2
        angle=0.5*real(i-1)*del
        trigs(i)=cos(angle)
        trigs(i+1)=sin(angle)
      enddo
      if (imode.eq.1) return
      if (imode.eq.8) return
      del=0.5*del
      nh=(nn+1)/2
      L=nh+nh
      la=nn+nn
      do i=1,L,2
        angle=0.5*real(i-1)*del
        trigs(la+i)=cos(angle)
        trigs(la+i+1)=sin(angle)
      enddo
      if (imode.le.3) return
      del=0.5*del
      la=la+nn
    if (mode.ne.5) then
      do i=2,nn
        angle=real(i-1)*del
        trigs(la+i)=2.0*sin(angle)
      enddo
      return
    endif
      del=0.5*del
      do i=2,n
        angle=real(i-1)*del
        trigs(la+i)=sin(angle)
      enddo
 end subroutine fftrig
 subroutine vpassm (a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
    integer, intent(in)  :: inc1, inc2, inc3, inc4, lot, n, ifac, la
    real,    intent(in)  :: a(n*(2*n+2)),b(n*(2*n+2)),trigs(n*(2*n+2))
    real,    intent(out) :: c(n*(2*n+2)),d(n*(2*n+2))
    real :: sin36=0.587785252292473
    real :: cos36=0.809016994374947
    real :: sin72=0.951056516295154
    real :: cos72=0.309016994374947
    real :: sin60=0.866025403784437
    integer :: i, j, k, L, m, iink, jink, jump, ibase, jbase, igo, ijk, la1
    integer :: ia, ja, ib, jb, kb, ic, jc, kc, id, jd, kd, ie, je, ke
    real    :: c1, s1, c2, s2, c3, s3, c4, s4
      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.gt.4) return
  select case (igo)
    case (1)
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      do 20 L=1,la
      i=ibase
      j=jbase
      do 15 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      do 30 L=1,la
      i=ibase
      j=jbase
      do 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 continue
      jbase=jbase+jump
   40 continue
    case (2)
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      do 60 L=1,la
      i=ibase
      j=jbase
      do 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      do 70 L=1,la
      i=ibase
      j=jbase
      do 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=                                                           &
          c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
         -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=                                                           &
          s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
         +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=                                                           &
          c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
         -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=                                                           &
          s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
         +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 continue
      jbase=jbase+jump
   80 continue
    case (3)
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      do 100 L=1,la
      i=ibase
      j=jbase
      do 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      do 110 L=1,la
      i=ibase
      j=jbase
      do 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=                                     &
          c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
         -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=                                     &
          s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
         +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=                                     &
          c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
         -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=                                     &
          s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
         +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=                                     &
          c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
         -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=                                     &
          s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
         +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 continue
      jbase=jbase+jump
  120 continue
    case (4)
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      do 140 L=1,la
      i=ibase
      j=jbase
      do 135 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
        -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
        +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
        +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
        -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
        -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
        +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
        +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
        -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  135 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      do 150 L=1,la
      i=ibase
      j=jbase
      do 145 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=                                                          &
          c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)=                                                          &
          s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)=                                                          &
          c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)=                                                          &
          s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)=                                                          &
          c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)=                                                          &
          s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)=                                                          &
          c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)=                                                          &
          s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      i=i+inc3
      j=j+inc4
  145 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 continue
      jbase=jbase+jump
  160 continue
  end select
 end subroutine vpassm
SUBROUTINE spat2spec(Z,U,spheredat)
IMPLICIT none
TYPE(sphere), INTENT(IN) :: spheredat
REAL, DIMENSION(spheredat%nlats, spheredat%nlons), INTENT(IN) :: Z
COMPLEX, DIMENSION(0:spheredat%ntrunc, 0:spheredat%ntrunc), INTENT(INOUT) :: U
#if defined(USE_BLAS_SHT)
CALL spharm_blas_ana_u(spheredat,Z,U,1,SHT_RSQRT2)
#else
BLOCK
COMPLEX, DIMENSION((spheredat%ntrunc+1)*(spheredat%ntrunc+2)/2) :: Uvec
#if defined(USE_EXACT_SHT)
CALL spharm_exact_ana(spheredat,Z,Uvec,1)
#else
BLOCK
REAL, DIMENSION(spheredat%nlats, spheredat%nlons) :: Zi
Zi(:,:) = Z(:,:)
CALL spharm_ref(spheredat,Zi,Uvec,1)
END BLOCK
#endif
CALL s2s_unpack(spheredat,Uvec,U)
END BLOCK
#endif
END SUBROUTINE spat2spec
!--- batched analysis of NF grids in one sweep over the pnm table
SUBROUTINE spat2spec_multi(Zmulti,Umulti,spheredat,nf)
IMPLICIT none
TYPE(sphere), INTENT(IN) :: spheredat
INTEGER, INTENT(IN) :: nf
REAL, DIMENSION(spheredat%nlats*spheredat%nlons*nf), INTENT(IN) :: Zmulti
COMPLEX, DIMENSION(0:spheredat%ntrunc, 0:spheredat%ntrunc, nf), INTENT(INOUT) :: Umulti
INTEGER :: f,nmdim
nmdim = (spheredat%ntrunc+1)*(spheredat%ntrunc+2)/2
#if defined(USE_BLAS_SHT)
CALL spharm_blas_ana_u(spheredat,Zmulti,Umulti,nf,SHT_RSQRT2)
#else
BLOCK
COMPLEX, DIMENSION(((spheredat%ntrunc+1)*(spheredat%ntrunc+2)/2)*nf) :: Uvec
#if defined(USE_EXACT_SHT)
CALL spharm_exact_ana(spheredat,Zmulti,Uvec,nf)
#else
BLOCK
REAL, DIMENSION(spheredat%nlats, spheredat%nlons) :: Zi
INTEGER :: ng
ng = spheredat%nlats*spheredat%nlons
DO f = 1,nf
  Zi = reshape(Zmulti((f-1)*ng+1:f*ng),(/spheredat%nlats,spheredat%nlons/))
  CALL spharm_ref(spheredat,Zi,Uvec((f-1)*nmdim+1:f*nmdim),1)
ENDDO
END BLOCK
#endif
DO f = 1,nf
  CALL s2s_unpack(spheredat,Uvec((f-1)*nmdim+1:f*nmdim),Umulti(:,:,f))
ENDDO
END BLOCK
#endif
END SUBROUTINE spat2spec_multi
SUBROUTINE s2s_unpack(spheredat,Uvec,U)
IMPLICIT none
TYPE(sphere), INTENT(IN) :: spheredat
COMPLEX, DIMENSION(*), INTENT(IN) :: Uvec
COMPLEX, DIMENSION(0:spheredat%ntrunc, 0:spheredat%ntrunc), INTENT(INOUT) :: U
INTEGER :: m,n,nt
nt = spheredat%ntrunc
DO m = 0,nt
  DO n = m,nt
    U(n,m) = Uvec(spheredat%s2s_idx(n*(nt+1)+m))/1.41421356
  ENDDO
ENDDO
END SUBROUTINE s2s_unpack
!--- batched synthesis of NF spectra in one sweep over the table
SUBROUTINE spec2spat_multi(Zmulti,Umulti,spheredat,nf)
IMPLICIT none
TYPE(sphere), INTENT(IN) :: spheredat
INTEGER, INTENT(IN) :: nf
REAL, DIMENSION(spheredat%nlats*spheredat%nlons*nf), INTENT(OUT) :: Zmulti
COMPLEX, DIMENSION(0:spheredat%ntrunc, 0:spheredat%ntrunc, nf), INTENT(IN) :: Umulti
INTEGER :: i,f,nmdim,ng
nmdim = (spheredat%ntrunc+1)*(spheredat%ntrunc+2)/2
#if defined(USE_BLAS_SHT)
CALL spharm_blas_syn_u(spheredat,Zmulti,Umulti,nf,SHT_SQRT2)
#else
BLOCK
COMPLEX, DIMENSION(((spheredat%ntrunc+1)*(spheredat%ntrunc+2)/2)*nf) :: Uvec
DO f = 1,nf
  DO i = 1,nmdim
    Uvec((f-1)*nmdim+i) = Umulti(spheredat%indxn(i),spheredat%indxm(i),f)*1.41421356
  ENDDO
ENDDO
#if defined(USE_EXACT_SHT)
CALL spharm_exact_syn(spheredat,Zmulti,Uvec,nf)
#else
BLOCK
REAL, DIMENSION(spheredat%nlats, spheredat%nlons) :: Zi
ng = spheredat%nlats*spheredat%nlons
DO f = 1,nf
  CALL spharm_ref(spheredat,Zi,Uvec((f-1)*nmdim+1:f*nmdim),-1)
  Zmulti((f-1)*ng+1:f*ng) = reshape(Zi,(/ng/))
ENDDO
END BLOCK
#endif
END BLOCK
#endif
END SUBROUTINE spec2spat_multi
SUBROUTINE spec2spat(Z,U,spheredat)
IMPLICIT none
TYPE(sphere), INTENT(IN) :: spheredat
REAL, DIMENSION(spheredat%nlats, spheredat%nlons), INTENT(OUT) :: Z
COMPLEX, DIMENSION(0:spheredat%ntrunc, 0:spheredat%ntrunc), INTENT(IN) :: U
INTEGER :: i
#if defined(USE_BLAS_SHT)
CALL spharm_blas_syn_u(spheredat,Z,U,1,SHT_SQRT2)
#else
BLOCK
COMPLEX, DIMENSION((spheredat%ntrunc+1)*(spheredat%ntrunc+2)/2) :: Uvec
DO i = 1,((spheredat%ntrunc+1)*(spheredat%ntrunc+2)/2)
  Uvec(i) = U(spheredat%indxn(i),spheredat%indxm(i))*1.41421356
ENDDO
#if defined(USE_EXACT_SHT)
CALL spharm_exact_syn(spheredat,Z,Uvec,1)
#else
CALL spharm_ref(spheredat,Z,Uvec,-1)
#endif
END BLOCK
#endif
END SUBROUTINE spec2spat
end module spharmt
