module modal_state_mod
  use iso_fortran_env, only: real64, int64
  implicit none
  private
  integer, parameter, public :: dp = real64

  public :: extrap_coeffs
  type, public :: response_spectrum
     integer :: lmax = 0
     integer :: max_modes = 0
     integer, allocatable :: nmod(:)
     real(dp), allocatable :: s(:,:)
     real(dp), allocatable :: a_sl(:,:), a_rr(:,:)
     real(dp), allocatable :: a_il(:,:), a_m(:,:)
     real(dp), allocatable :: a_tsl(:,:), a_trr(:,:)
     ! Optional horizontal-displacement weights (l residues / s and the
     ! degree-2 tidal l residues / s).  Allocated only when the solver is
     ! asked for the horizontal output; read by read_horizontal only.
     real(dp), allocatable :: a_hl(:,:), a_thl(:,:)
  contains
     procedure :: load => spectrum_load
     procedure :: destroy => spectrum_destroy
  end type response_spectrum

  type, public :: response_memory_state
     integer :: lmax = 0
     integer :: max_modes = 0
     integer :: ncoef = 0
     integer(int64) :: nstate = 0_int64
     integer, allocatable :: nmod(:)
     integer(int64), allocatable :: roffset(:)
     real(dp), allocatable :: decay(:,:)
     real(dp) :: cached_dt_years = -huge(1.0_dp)
     complex(dp), allocatable :: p_load(:)
     ! O(1) auxiliary restart state for the full solver.
     ! p_ice/p_ocean are cumulative spectral thickness changes since time0;
     ! prev_ice is the previous absolute grounded-ice spectrum.
     complex(dp), allocatable :: p_ice(:), p_ocean(:), prev_ice(:)
     complex(dp), allocatable :: r_load(:)
     real(dp) :: p_il(3,3) = 0.0_dp
     real(dp) :: p_m(3) = 0.0_dp
     complex(dp) :: p_lambda(0:2) = (0.0_dp,0.0_dp)
     real(dp), allocatable :: r_il(:,:,:)
     real(dp), allocatable :: r_m(:,:)
     complex(dp), allocatable :: r_lambda(:,:)
  contains
     procedure :: init => state_init
     procedure :: zero => state_zero
     procedure :: advance_and_read => state_advance_and_read
     procedure :: read_horizontal => state_read_horizontal
     procedure :: commit => state_commit
     procedure :: save_checkpoint => state_save_checkpoint
     procedure :: load_checkpoint => state_load_checkpoint
     procedure :: bytes => state_bytes
     procedure :: destroy => state_destroy
  end type response_memory_state

  public :: lm_index

contains

  ! Coefficients of the extrapolated inner start (variable steps):
  !   dS0_j = dt_j [ Sdot_{j-1} + (Sdot_{j-1} - Sdot_{j-2}) (dt_j + dt_{j-1}) / (dt_{j-1} + dt_{j-2}) ],
  !   Sdot_i = dS_i / dt_i,   i.e.  dS0_j = c1 dS_{j-1} + c2 dS_{j-2};
  ! equal steps give exactly c1 = 2, c2 = -1.  dt0 = dt_j, dtm1 = dt_{j-1}, dtm2 = dt_{j-2}.
  pure subroutine extrap_coeffs(dt0, dtm1, dtm2, c1, c2)
    real(dp), intent(in) :: dt0, dtm1, dtm2
    real(dp), intent(out) :: c1, c2
    if (dt0 == dtm1 .and. dtm1 == dtm2) then
       c1 = 2.0_dp; c2 = -1.0_dp
    else
       c1 = dt0 / dtm1 * (1.0_dp + (dt0 + dtm1) / (dtm1 + dtm2))
       c2 = -dt0 / dtm2 * ((dt0 + dtm1) / (dtm1 + dtm2))
    endif
  end subroutine extrap_coeffs

  pure integer function lm_index(l,m) result(idx)
    integer, intent(in) :: l,m
    idx = l*(l+1)/2 + m + 1
  end function lm_index

  subroutine spectrum_destroy(self)
    class(response_spectrum), intent(inout) :: self
    if (allocated(self%nmod)) deallocate(self%nmod)
    if (allocated(self%s)) deallocate(self%s)
    if (allocated(self%a_sl)) deallocate(self%a_sl)
    if (allocated(self%a_rr)) deallocate(self%a_rr)
    if (allocated(self%a_il)) deallocate(self%a_il)
    if (allocated(self%a_m)) deallocate(self%a_m)
    if (allocated(self%a_tsl)) deallocate(self%a_tsl)
    if (allocated(self%a_trr)) deallocate(self%a_trr)
    if (allocated(self%a_hl)) deallocate(self%a_hl)
    if (allocated(self%a_thl)) deallocate(self%a_thl)
    self%lmax=0; self%max_modes=0
  end subroutine spectrum_destroy

  subroutine spectrum_load(self, filename, requested_lmax)
    class(response_spectrum), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: requested_lmax
    integer :: u, file_lmax, file_max_modes, ll, lread, nm, kk, keep_lmax
    real(dp) :: sv, c1, c2, c3, c4, c5, c6

    open(newunit=u,file=filename,status='old',action='read')
    read(u,*) file_lmax, file_max_modes
    keep_lmax = file_lmax
    if (present(requested_lmax)) keep_lmax = min(requested_lmax,file_lmax)
    self%lmax = keep_lmax
    self%max_modes = file_max_modes
    allocate(self%nmod(0:keep_lmax))
    allocate(self%s(file_max_modes,0:keep_lmax))
    allocate(self%a_sl(file_max_modes,0:keep_lmax), self%a_rr(file_max_modes,0:keep_lmax))
    allocate(self%a_il(file_max_modes,0:keep_lmax), self%a_m(file_max_modes,0:keep_lmax))
    allocate(self%a_tsl(file_max_modes,0:keep_lmax), self%a_trr(file_max_modes,0:keep_lmax))
    self%nmod = 0
    self%s = 0.0_dp; self%a_sl=0.0_dp; self%a_rr=0.0_dp
    self%a_il=0.0_dp; self%a_m=0.0_dp
    self%a_tsl=0.0_dp; self%a_trr=0.0_dp

    do ll=1,file_lmax
       read(u,*) lread, nm
       if (lread /= ll) error stop 'spectrum degree ordering mismatch'
       if (ll <= keep_lmax) self%nmod(ll)=nm
       do kk=1,nm
          read(u,*) sv,c1,c2,c3,c4,c5,c6
          if (ll <= keep_lmax) then
             self%s(kk,ll)=sv
             self%a_sl(kk,ll)=c1
             self%a_rr(kk,ll)=c2
             self%a_il(kk,ll)=c3
             self%a_m(kk,ll)=c4
             self%a_tsl(kk,ll)=c5
             self%a_trr(kk,ll)=c6
          end if
       end do
    end do
    close(u)
  end subroutine spectrum_load

  subroutine state_init(self,spec)
    class(response_memory_state), intent(inout) :: self
    type(response_spectrum), intent(in) :: spec
    integer :: ll
    integer(int64) :: off

    self%lmax=spec%lmax; self%max_modes=spec%max_modes
    self%ncoef=(self%lmax+1)*(self%lmax+2)/2
    allocate(self%nmod(0:self%lmax), self%roffset(0:self%lmax))
    self%nmod=spec%nmod
    self%roffset=0_int64
    off=1_int64
    do ll=1,self%lmax
       self%roffset(ll)=off
       off = off + int(self%nmod(ll),int64)*int(ll+1,int64)
    end do
    self%nstate=off-1_int64
    allocate(self%decay(self%max_modes,0:self%lmax))
    allocate(self%p_load(self%ncoef), self%p_ice(self%ncoef), self%p_ocean(self%ncoef), &
             self%prev_ice(self%ncoef), self%r_load(self%nstate))
    allocate(self%r_il(self%max_modes,3,3), self%r_m(self%max_modes,3))
    allocate(self%r_lambda(self%max_modes,0:2))
    call self%zero()
  end subroutine state_init

  subroutine state_destroy(self)
    class(response_memory_state), intent(inout) :: self
    if (allocated(self%nmod)) deallocate(self%nmod)
    if (allocated(self%roffset)) deallocate(self%roffset)
    if (allocated(self%decay)) deallocate(self%decay)
    if (allocated(self%p_load)) deallocate(self%p_load)
    if (allocated(self%p_ice)) deallocate(self%p_ice)
    if (allocated(self%p_ocean)) deallocate(self%p_ocean)
    if (allocated(self%prev_ice)) deallocate(self%prev_ice)
    if (allocated(self%r_load)) deallocate(self%r_load)
    if (allocated(self%r_il)) deallocate(self%r_il)
    if (allocated(self%r_m)) deallocate(self%r_m)
    if (allocated(self%r_lambda)) deallocate(self%r_lambda)
    self%lmax=0; self%max_modes=0; self%ncoef=0; self%nstate=0_int64
    self%cached_dt_years=-huge(1.0_dp)
    self%p_il=0.0_dp; self%p_m=0.0_dp; self%p_lambda=(0.0_dp,0.0_dp)
  end subroutine state_destroy

  subroutine state_zero(self)
    class(response_memory_state), intent(inout) :: self
    self%cached_dt_years=-huge(1.0_dp)
    if (allocated(self%decay)) self%decay=0.0_dp
    if (allocated(self%p_load)) self%p_load=(0.0_dp,0.0_dp)
    if (allocated(self%p_ice)) self%p_ice=(0.0_dp,0.0_dp)
    if (allocated(self%p_ocean)) self%p_ocean=(0.0_dp,0.0_dp)
    if (allocated(self%prev_ice)) self%prev_ice=(0.0_dp,0.0_dp)
    if (allocated(self%r_load)) self%r_load=(0.0_dp,0.0_dp)
    self%p_il=0.0_dp; self%p_m=0.0_dp; self%p_lambda=(0.0_dp,0.0_dp)
    if (allocated(self%r_il)) self%r_il=0.0_dp
    if (allocated(self%r_m)) self%r_m=0.0_dp
    if (allocated(self%r_lambda)) self%r_lambda=(0.0_dp,0.0_dp)
  end subroutine state_zero

  subroutine update_decay(self,spec,dt_years)
    class(response_memory_state), intent(inout) :: self
    type(response_spectrum), intent(in) :: spec
    real(dp), intent(in) :: dt_years
    integer :: ll,kk
    if (abs(dt_years-self%cached_dt_years) <= epsilon(1.0_dp)*max(1.0_dp,abs(dt_years))) return
    self%decay=0.0_dp
    do ll=1,self%lmax
       do kk=1,self%nmod(ll)
          self%decay(kk,ll)=exp(-spec%s(kk,ll)*dt_years/1000.0_dp)
       end do
    end do
    self%cached_dt_years=dt_years
  end subroutine update_decay

  subroutine state_advance_and_read(self,spec,dt_years,viscous,radial,sum_il,sum_m,sum_lambda_sl,sum_lambda_rr)
    class(response_memory_state), intent(inout) :: self
    type(response_spectrum), intent(in) :: spec
    real(dp), intent(in) :: dt_years
    complex(dp), intent(out) :: viscous(:),radial(:)
    real(dp), intent(out) :: sum_il(3,3),sum_m(3)
    complex(dp), intent(out) :: sum_lambda_sl(0:2),sum_lambda_rr(0:2)
    integer :: ll,kk,ibase,rstart,rend
    real(dp) :: d,omd

    call update_decay(self,spec,dt_years)
    viscous=(0.0_dp,0.0_dp); radial=(0.0_dp,0.0_dp)

    ! Packed layout is [degree][pole][m=0:l], so each pole update is one
    ! contiguous vector operation.  R itself is the hereditary response state:
    ! R_new = d*R_old + (1-d)*DL_old   (DL = cumulative load, p_load).
    ! degrees are independent (each writes a disjoint slice of
    ! r_load / viscous / radial and accumulates over kk in the same order),
    ! so this loop threads identically.
    !$omp parallel do schedule(dynamic,4) default(shared) private(ll,kk,ibase,rstart,rend,d,omd)
    do ll=1,self%lmax
       ibase=lm_index(ll,0)
       do kk=1,self%nmod(ll)
          d=self%decay(kk,ll); omd=1.0_dp-d
          rstart=int(self%roffset(ll))+(kk-1)*(ll+1)
          rend=rstart+ll
          self%r_load(rstart:rend)=d*self%r_load(rstart:rend) + &
                                    omd*self%p_load(ibase:ibase+ll)
          viscous(ibase:ibase+ll)=viscous(ibase:ibase+ll) + &
                                  spec%a_sl(kk,ll)*self%r_load(rstart:rend)
          radial(ibase:ibase+ll)=radial(ibase:ibase+ll) + &
                                 spec%a_rr(kk,ll)*self%r_load(rstart:rend)
       end do
    end do
    !$omp end parallel do

    sum_il=0.0_dp; sum_m=0.0_dp
    sum_lambda_sl=(0.0_dp,0.0_dp); sum_lambda_rr=(0.0_dp,0.0_dp)
    if (self%lmax >= 2) then
       do kk=1,self%nmod(2)
          d=self%decay(kk,2); omd=1.0_dp-d
          self%r_il(kk,:,:)=d*self%r_il(kk,:,:)+omd*self%p_il(:,:)
          self%r_m(kk,:)=d*self%r_m(kk,:)+omd*self%p_m(:)
          self%r_lambda(kk,:)=d*self%r_lambda(kk,:)+omd*self%p_lambda(:)
          sum_il=sum_il+spec%a_il(kk,2)*self%r_il(kk,:,:)
          sum_m=sum_m+spec%a_m(kk,2)*self%r_m(kk,:)
          sum_lambda_sl=sum_lambda_sl+spec%a_tsl(kk,2)*self%r_lambda(kk,:)
          sum_lambda_rr=sum_lambda_rr+spec%a_trr(kk,2)*self%r_lambda(kk,:)
       end do
    end if
  end subroutine state_advance_and_read

  ! Horizontal-displacement read-out of the CURRENT modal states.  Call it
  ! between advance_and_read and commit of the same step: it reads the same
  ! states R_k that advance_and_read just used for the radial sum and only
  ! applies different weights (a_hl = l residues / s).  The state is not
  ! modified, so the solver path is unchanged whether or not this is called.
  subroutine state_read_horizontal(self,spec,horiz,sum_lambda_hl)
    class(response_memory_state), intent(in) :: self
    type(response_spectrum), intent(in) :: spec
    complex(dp), intent(out) :: horiz(:)
    complex(dp), intent(out) :: sum_lambda_hl(0:2)
    integer :: ll,kk,ibase,rstart,rend

    horiz=(0.0_dp,0.0_dp)
    sum_lambda_hl=(0.0_dp,0.0_dp)
    if (.not.allocated(spec%a_hl)) error stop 'read_horizontal: horizontal weights not set'
    !$omp parallel do schedule(dynamic,4) default(shared) private(ll,kk,ibase,rstart,rend)
    do ll=1,self%lmax
       ibase=lm_index(ll,0)
       do kk=1,self%nmod(ll)
          rstart=int(self%roffset(ll))+(kk-1)*(ll+1)
          rend=rstart+ll
          horiz(ibase:ibase+ll)=horiz(ibase:ibase+ll) + &
                                spec%a_hl(kk,ll)*self%r_load(rstart:rend)
       end do
    end do
    !$omp end parallel do
    if (self%lmax >= 2) then
       do kk=1,self%nmod(2)
          sum_lambda_hl=sum_lambda_hl+spec%a_thl(kk,2)*self%r_lambda(kk,:)
       end do
    end if
  end subroutine state_read_horizontal

  subroutine state_commit(self,current_load,dil,dm,dlambda)
    class(response_memory_state), intent(inout) :: self
    complex(dp), intent(in) :: current_load(:)
    real(dp), intent(in), optional :: dil(3,3),dm(3)
    complex(dp), intent(in), optional :: dlambda(0:2)
    self%p_load=self%p_load+current_load
    if (present(dil)) self%p_il=self%p_il+dil
    if (present(dm)) self%p_m=self%p_m+dm
    if (present(dlambda)) self%p_lambda=self%p_lambda+dlambda
  end subroutine state_commit

  subroutine state_save_checkpoint(self,filename)
    class(response_memory_state), intent(in) :: self
    character(len=*), intent(in) :: filename
    integer :: u
    open(newunit=u,file=filename,status='replace',access='stream',form='unformatted',action='write')
    write(u) int(z'45584D53'), self%lmax, self%max_modes, self%ncoef, self%nstate
    write(u) self%nmod, self%roffset, self%cached_dt_years
    write(u) self%decay
    write(u) self%p_load,self%p_ice,self%p_ocean,self%prev_ice,self%r_load
    write(u) self%p_il,self%p_m,self%p_lambda
    write(u) self%r_il,self%r_m,self%r_lambda
    close(u)
  end subroutine state_save_checkpoint

  subroutine state_load_checkpoint(self,filename,spec)
    class(response_memory_state), intent(inout) :: self
    character(len=*), intent(in) :: filename
    type(response_spectrum), intent(in) :: spec
    integer :: u,magic,lm,mm,nc
    integer(int64) :: ns
    if (.not.allocated(self%p_load)) call self%init(spec)
    open(newunit=u,file=filename,status='old',access='stream',form='unformatted',action='read')
    read(u) magic,lm,mm,nc,ns
    if (magic /= int(z'45584D53')) error stop 'bad stage-2 response-memory checkpoint magic'
    if (lm/=self%lmax .or. mm/=self%max_modes .or. nc/=self%ncoef .or. ns/=self%nstate) &
       error stop 'response-memory checkpoint shape mismatch'
    read(u) self%nmod,self%roffset,self%cached_dt_years
    read(u) self%decay
    read(u) self%p_load,self%p_ice,self%p_ocean,self%prev_ice,self%r_load
    read(u) self%p_il,self%p_m,self%p_lambda
    read(u) self%r_il,self%r_m,self%r_lambda
    close(u)
  end subroutine state_load_checkpoint

  integer(int64) function state_bytes(self) result(nbytes)
    class(response_memory_state), intent(in) :: self
    nbytes=0_int64
    if (allocated(self%p_load)) nbytes=nbytes+16_int64*size(self%p_load,kind=int64)
    if (allocated(self%p_ice)) nbytes=nbytes+16_int64*size(self%p_ice,kind=int64)
    if (allocated(self%p_ocean)) nbytes=nbytes+16_int64*size(self%p_ocean,kind=int64)
    if (allocated(self%prev_ice)) nbytes=nbytes+16_int64*size(self%prev_ice,kind=int64)
    if (allocated(self%r_load)) nbytes=nbytes+16_int64*size(self%r_load,kind=int64)
    if (allocated(self%r_il)) nbytes=nbytes+8_int64*size(self%r_il,kind=int64)
    if (allocated(self%r_m)) nbytes=nbytes+8_int64*size(self%r_m,kind=int64)
    if (allocated(self%r_lambda)) nbytes=nbytes+16_int64*size(self%r_lambda,kind=int64)
    nbytes=nbytes+8_int64*9_int64+8_int64*3_int64+16_int64*3_int64
  end function state_bytes

end module modal_state_mod
