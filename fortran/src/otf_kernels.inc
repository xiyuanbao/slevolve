!=======================================================================
! On-the-fly Legendre backend (included into module spharmt by
! spharmt_otf.f90; active when built with -DUSE_BLAS_SHT -DUSE_OTF_SHT).
!
! No P_lm table.  Per zonal wavenumber m the normalised associated Legendre
! functions of the legacy `legend`,
!     P_mm        = sqrt(0.5*prod_m)                (tabulated, nhsp*(L+1) values)
!     P_m+1,m     = sqrt(2m+3) mu P_mm
!     P_n,m       = ( mu P_n-1,m - eps_n-1 P_n-2,m ) / eps_n
! are regenerated in the monic (rescaled) form
!     P_n = c_n Q_n,  c_n = prod_{i=m+1..n} a_i,  a_m+1 = sqrt(2m+3), a_i = 1/eps_i (i >= m+2)
!     Q_m = P_mm,  Q_m+1 = mu Q_m,  Q_n = mu Q_n-1 + be_n Q_n-2,  be_n = -eps_n-1/eps_n/(a_n a_n-1)
! so one degree step costs one multiply and one FMA per latitude; c_n is
! folded into the coefficients (synthesis) or the result (analysis), O(L)
! per m.  c_n <= 2^L ~ 1e154, so Q never overflows; it can only underflow
! where P itself is < ~1e-150, i.e. far below the polar-skip threshold.
!
! Kernel shape: the latitude loop is innermost (vectorised by gfortran,
! 2 doubles per NEON register) and FOUR degrees are advanced per pass over
! the latitudes, so each pass loads mu, Q_n-1, Q_n-2 once for four
! degrees and the partial sums / dot products stay in registers.
! Even degrees pair with the N+S, odd degrees with the N-S Fourier
! coefficients (equatorial fold).  Latitudes nearer the pole than
! otf_jb0(m) are skipped (polar optimisation, threshold sht_otf_thr).
!=======================================================================
 subroutine otf_setup(sd, gl)
   type (sphere), intent(inout) :: sd
   double precision, intent(in) :: gl(:)
   integer :: L, nhs, nhsp, j, m, n, jf, nmstrt, nmdim
   double precision :: sinsq, a, b, prod, epsn, epsn1, an, an1, bn, cn, t
   L = sd%ntrunc; nhs = sd%nhs; nhsp = sd%nhsp
   nmdim = (L+1)*(L+2)/2
   allocate(sd%otf_x(nhsp), sd%otf_pmm(nhsp,0:L), sd%otf_jb0(0:L))
   allocate(sd%otf_be(nmdim+8), sd%otf_c(nmdim+8))
   sd%otf_x = 0.0; sd%otf_pmm = 0.0; sd%otf_be = 0.0; sd%otf_c = 0.0
   do j = 1, nhs
      sd%otf_x(j) = gl(j)
      sinsq = 1.d0 - gl(j)*gl(j)
      a = 1.d0; b = 0.d0; prod = 1.d0
      do m = 0, L
         if (m /= 0) then
            a = a + 2.d0; b = b + 2.d0
            prod = prod*sinsq*a/b
         endif
         sd%otf_pmm(j,m) = sqrt(0.5d0*prod)
      enddo
   enddo
   an1 = 1.d0
   do m = 0, L
      nmstrt = sd%nmstrt_tab(m+1)
      sd%otf_c(nmstrt+1) = 1.d0
      if (m < L) then
         an1 = sqrt(dble(2*m+3))
         sd%otf_c(nmstrt+2) = an1
      endif
      cn = an1
      do n = m+2, L
         epsn  = sqrt(dble(n*n-m*m)/dble(4*n*n-1))
         epsn1 = sqrt(dble((n-1)*(n-1)-m*m)/dble(4*(n-1)*(n-1)-1))
         an = 1.d0/epsn
         bn = -epsn1/epsn
         sd%otf_be(nmstrt+n-m+1) = bn/(an*an1)
         cn = cn*an
         sd%otf_c(nmstrt+n-m+1) = cn
         an1 = an
      enddo
   enddo
   ! polar skip: first latitude where some P_lm (l = m..L) exceeds the
   ! threshold, scanned with legend's own recurrence
   !$omp parallel do schedule(dynamic,4) default(shared) &
   !$omp   private(m,j,n,jf,a,b,t,prod,epsn,epsn1)
   do m = 0, L
      jf = nhs+1
      if (sht_otf_thr <= 0.0) jf = 1
      do j = 1, nhs
         if (jf <= nhs) exit
         a = sd%otf_pmm(j,m); prod = abs(a)
         if (m < L) then
            b = sqrt(dble(2*m+3))*sd%otf_x(j)*a; prod = max(prod,abs(b))
            epsn1 = sqrt(dble((m+1)*(m+1)-m*m)/dble(4*(m+1)*(m+1)-1))
            do n = m+2, L
               epsn = sqrt(dble(n*n-m*m)/dble(4*n*n-1))
               t = (sd%otf_x(j)*b - epsn1*a)/epsn
               a = b; b = t; epsn1 = epsn
               prod = max(prod, abs(b))
            enddo
         endif
         if (prod > sht_otf_thr) jf = j
      enddo
      sd%otf_jb0(m) = min(((jf-1)/2)*2 + 1, nhsp+1)
   enddo
   !$omp end parallel do
 end subroutine otf_setup

 !---- analysis: Fourier coefficients in wrk1 -> U(0:L,0:L,nf) ----------
 subroutine otf_ana_stage(sd, U, nf, cscale, tt)
   type (sphere), intent(in) :: sd
   integer, intent(in) :: nf
   real, intent(in) :: cscale
   complex, dimension(0:sd%ntrunc, 0:sd%ntrunc, nf), intent(inout) :: U
   real(8), intent(inout) :: tt
   integer :: nlats, nl2, nhs, nhsp, mwaves, nc, j, f, i2, m, ri, c, basen, bases, k, mg
   real :: wj, a1, a2
   nlats = sd%nlats; nl2 = sd%nlons+2; nhs = sd%nhs; nhsp = sd%nhsp
   mwaves = sd%ntrunc+1; nc = 2*nf
   ! E(j, c, m) = w_j (N+S),  O(j, c, m) = w_j (N-S);  c = 2(f-1)+{1 re, 2 im}
   ! m in groups of 4: one 64-byte line of a wrk1 row holds 4 (re,im) pairs,
   ! and the writes are unit-stride in j for every (m, c) stream.
   !$omp parallel do schedule(static) default(shared) &
   !$omp   private(mg,j,wj,f,basen,bases,i2,m,ri,c,k,a1,a2)
   do mg = 0, mwaves-1, 4
     do f = 1, nf
       do j = 1, nhs
         wj = sd%weights(j)*sht_fsc*cscale
         basen = (f-1)*nlats*nl2 + (j-1)*nl2
         bases = (f-1)*nlats*nl2 + (nlats-j)*nl2
         do i2 = 2*mg+1, min(2*mg+8, 2*mwaves)
           m = (i2+1)/2-1; ri = 2-mod(i2,2); c = 2*(f-1)+ri
           k = j + nhsp*((c-1) + nc*m)
           a1 = sd%wrk1(basen+i2); a2 = sd%wrk1(bases+i2)
           sd%bbuf(k)  = (a1+a2)*wj
           sd%bbuf2(k) = (a1-a2)*wj
         enddo
       enddo
       do j = nhs+1, nhsp
         do m = mg, min(mg+3, mwaves-1)
           do ri = 1, 2
             k = j + nhsp*((2*(f-1)+ri-1) + nc*m)
             sd%bbuf(k) = 0.0; sd%bbuf2(k) = 0.0
           enddo
         enddo
       enddo
     enddo
   enddo
   !$omp end parallel do
   call sht_toc(3,tt)
   !$omp parallel do schedule(dynamic,1) default(shared) private(m)
   do m = 0, mwaves-1
     call otf_ana_m(sd, m, nc, sd%bbuf(1+nhsp*nc*m:), sd%bbuf2(1+nhsp*nc*m:), U(:,m,:), nf)
   enddo
   !$omp end parallel do
   call sht_toc(4,tt)
 end subroutine otf_ana_stage

 subroutine otf_ana_m(sd, m, nc, E, O, Um, nf)
   type (sphere), intent(in) :: sd
   integer, intent(in) :: m, nc, nf
   real, intent(in) :: E(sd%nhsp, nc), O(sd%nhsp, nc)
   complex, intent(inout) :: Um(0:sd%ntrunc, nf)
   real :: ar(0:sd%ntrunc-m+4, 4), q1(sd%nhsp)
   integer :: L, lenm, j0, nh, n, k, f, nm0
   L = sd%ntrunc; lenm = L-m+1; nm0 = sd%nmstrt_tab(m+1)+1
   j0 = sd%otf_jb0(m); nh = sd%nhsp
   if (j0 > nh) then
     do f = 1, nf
       Um(m:L, f) = (0.0, 0.0)
     enddo
     return
   endif
   n = nh-j0+1
   if (nf >= sht_otf_gemm) then
     call otf_ana_gemm(sd, m, nc, E, O, Um, nf, j0, n)
     return
   endif
   q1(1:n) = sd%otf_x(j0:nh)*sd%otf_pmm(j0:nh, m)
   f = 1
   do while (f <= nf)
     if (sht_otf_pair .and. f+1 <= nf) then
       call otf_ana_k2(lenm, n, sd%otf_x(j0:nh), sd%otf_pmm(j0:nh,m), q1, sd%otf_be(nm0:nm0+lenm+3), &
            E(j0,2*f-1), E(j0,2*f), O(j0,2*f-1), O(j0,2*f), &
            E(j0,2*f+1), E(j0,2*f+2), O(j0,2*f+1), O(j0,2*f+2), ar(0,1), ar(0,2), ar(0,3), ar(0,4))
       do k = 0, lenm-1
         Um(m+k, f)   = cmplx(ar(k,1)*sd%otf_c(nm0+k), ar(k,2)*sd%otf_c(nm0+k))
         Um(m+k, f+1) = cmplx(ar(k,3)*sd%otf_c(nm0+k), ar(k,4)*sd%otf_c(nm0+k))
       enddo
       f = f+2
     else
       call otf_ana_k1(lenm, n, sd%otf_x(j0:nh), sd%otf_pmm(j0:nh,m), q1, sd%otf_be(nm0:nm0+lenm+3), &
            E(j0,2*f-1), E(j0,2*f), O(j0,2*f-1), O(j0,2*f), ar(0,1), ar(0,2))
       do k = 0, lenm-1
         Um(m+k, f) = cmplx(ar(k,1)*sd%otf_c(nm0+k), ar(k,2)*sd%otf_c(nm0+k))
       enddo
       f = f+1
     endif
   enddo
 end subroutine otf_ana_m

 ! one field: ar/ai(k) = sum_j Q_{m+k}(j) * (E or O)(j)
 subroutine otf_ana_k1(lenm, n, x, q0, q1i, be, Er, Ei, Or, Oi, ar, ai)
   integer, intent(in) :: lenm, n
   real, intent(in) :: x(n), q0(n), q1i(n), be(0:lenm+3), Er(n), Ei(n), Or(n), Oi(n)
   real, intent(out) :: ar(0:lenm+3), ai(0:lenm+3)
   real :: p1(n), p2(n), a, b, c, d, s1, s2, s3, s4, t1, t2, t3, t4
   integer :: j, k
   s1 = 0; t1 = 0; s2 = 0; t2 = 0
   !$omp simd reduction(+:s1,s2,t1,t2)
   do j = 1, n
     p2(j) = q0(j); p1(j) = q1i(j)
     s1 = s1 + q0(j)*Er(j); t1 = t1 + q0(j)*Ei(j)
     s2 = s2 + q1i(j)*Or(j); t2 = t2 + q1i(j)*Oi(j)
   enddo
   ar(0) = s1; ai(0) = t1; ar(1) = s2; ai(1) = t2
   do k = 2, lenm-1, 4
     s1 = 0; s2 = 0; s3 = 0; s4 = 0; t1 = 0; t2 = 0; t3 = 0; t4 = 0
     !$omp simd private(a,b,c,d) reduction(+:s1,s2,s3,s4,t1,t2,t3,t4)
     do j = 1, n
       a = x(j)*p1(j) + be(k)*p2(j)
       b = x(j)*a + be(k+1)*p1(j)
       c = x(j)*b + be(k+2)*a
       d = x(j)*c + be(k+3)*b
       s1 = s1 + a*Er(j); t1 = t1 + a*Ei(j)
       s2 = s2 + b*Or(j); t2 = t2 + b*Oi(j)
       s3 = s3 + c*Er(j); t3 = t3 + c*Ei(j)
       s4 = s4 + d*Or(j); t4 = t4 + d*Oi(j)
       p2(j) = c; p1(j) = d
     enddo
     ar(k) = s1; ar(k+1) = s2; ar(k+2) = s3; ar(k+3) = s4
     ai(k) = t1; ai(k+1) = t2; ai(k+2) = t3; ai(k+3) = t4
   enddo
 end subroutine otf_ana_k1

 ! two fields sharing one recurrence (two degrees per pass keeps the
 ! 16 accumulators of four degrees from spilling)
 subroutine otf_ana_k2(lenm, n, x, q0, q1i, be, Er, Ei, Or, Oi, Fr, Fi, Pr, Pi, ar, ai, br, bi)
   integer, intent(in) :: lenm, n
   real, intent(in) :: x(n), q0(n), q1i(n), be(0:lenm+3)
   real, intent(in) :: Er(n), Ei(n), Or(n), Oi(n), Fr(n), Fi(n), Pr(n), Pi(n)
   real, intent(out) :: ar(0:lenm+3), ai(0:lenm+3), br(0:lenm+3), bi(0:lenm+3)
   real :: p1(n), p2(n), a, b, s1, s2, t1, t2, u1, u2, v1, v2
   integer :: j, k
   s1 = 0; t1 = 0; s2 = 0; t2 = 0; u1 = 0; v1 = 0; u2 = 0; v2 = 0
   !$omp simd reduction(+:s1,s2,t1,t2,u1,u2,v1,v2)
   do j = 1, n
     p2(j) = q0(j); p1(j) = q1i(j)
     s1 = s1 + q0(j)*Er(j); t1 = t1 + q0(j)*Ei(j)
     s2 = s2 + q1i(j)*Or(j); t2 = t2 + q1i(j)*Oi(j)
     u1 = u1 + q0(j)*Fr(j); v1 = v1 + q0(j)*Fi(j)
     u2 = u2 + q1i(j)*Pr(j); v2 = v2 + q1i(j)*Pi(j)
   enddo
   ar(0) = s1; ai(0) = t1; ar(1) = s2; ai(1) = t2
   br(0) = u1; bi(0) = v1; br(1) = u2; bi(1) = v2
   do k = 2, lenm-1, 2
     s1 = 0; t1 = 0; s2 = 0; t2 = 0; u1 = 0; v1 = 0; u2 = 0; v2 = 0
     !$omp simd private(a,b) reduction(+:s1,s2,t1,t2,u1,u2,v1,v2)
     do j = 1, n
       a = x(j)*p1(j) + be(k)*p2(j)
       b = x(j)*a + be(k+1)*p1(j)
       s1 = s1 + a*Er(j); t1 = t1 + a*Ei(j)
       u1 = u1 + a*Fr(j); v1 = v1 + a*Fi(j)
       s2 = s2 + b*Or(j); t2 = t2 + b*Oi(j)
       u2 = u2 + b*Pr(j); v2 = v2 + b*Pi(j)
       p2(j) = a; p1(j) = b
     enddo
     ar(k) = s1; ai(k) = t1; ar(k+1) = s2; ai(k+1) = t2
     br(k) = u1; bi(k) = v1; br(k+1) = u2; bi(k+1) = v2
   enddo
 end subroutine otf_ana_k2

 !---- synthesis: U(0:L,0:L,nf) -> grid ---------------------------------
 subroutine otf_syn_stage(sd, ugrid, U, nf, cscale, tt)
   type (sphere), intent(in) :: sd
   integer, intent(in) :: nf
   real, intent(in) :: cscale
   real, dimension(sd%nlons*sd%nlats*nf), intent(out) :: ugrid
   complex, dimension(0:sd%ntrunc, 0:sd%ntrunc, nf), intent(in) :: U
   real(8), intent(inout) :: tt
   integer :: nlats, nlons, nl2, nhs, nhsp, mwaves, nc, j, f, i2, m, ri, c, basen, bases, k, off, goff, base, mg
   real :: se, so
   nlats = sd%nlats; nlons = sd%nlons; nl2 = nlons+2; nhs = sd%nhs; nhsp = sd%nhsp
   mwaves = sd%ntrunc+1; nc = 2*nf
   !$omp parallel do schedule(dynamic,1) default(shared) private(m)
   do m = 0, mwaves-1
     call otf_syn_m(sd, m, nc, sd%bbuf(1+nhsp*nc*m:), sd%bbuf2(1+nhsp*nc*m:), U(:,m,:), nf, cscale)
   enddo
   !$omp end parallel do
   call sht_toc(6,tt)
   !$omp parallel do schedule(static) default(shared) &
   !$omp   private(mg,f,j,basen,bases,i2,m,ri,c,k,se,so)
   do mg = 0, nl2/2-1, 4
     do f = 1, nf
       do j = 1, nhs
         basen = (f-1)*nlats*nl2 + (j-1)*nl2
         bases = (f-1)*nlats*nl2 + (nlats-j)*nl2
         do i2 = 2*mg+1, min(2*mg+8, nl2)
           if (i2 <= 2*mwaves) then
             m = (i2+1)/2-1; ri = 2-mod(i2,2); c = 2*(f-1)+ri
             k = j + nhsp*((c-1) + nc*m)
             se = sd%bbuf(k); so = sd%bbuf2(k)
             sd%wrk1(basen+i2) = se + so
             sd%wrk1(bases+i2) = se - so
           else
             sd%wrk1(basen+i2) = 0.0
             sd%wrk1(bases+i2) = 0.0
           endif
         enddo
       enddo
     enddo
   enddo
   !$omp end parallel do
   call sht_toc(7,tt)
   call sht_fft(sd, sd%wrk1, sd%wrk2, nlats*nf, 1)
   call sht_toc(8,tt)
   !$omp parallel do collapse(2) schedule(static) default(shared) private(f,j,off,goff,base)
   do f = 1, nf
     do j = 1, nlats
       off  = (f-1)*nlats*nl2
       goff = (f-1)*nlons*nlats
       base = off + (j-1)*nl2
       ugrid(goff+(j-1)*nlons+1 : goff+(j-1)*nlons+nlons) = sd%wrk1(base+1:base+nlons)
     enddo
   enddo
   !$omp end parallel do
   call sht_toc(5,tt)
 end subroutine otf_syn_stage

 subroutine otf_syn_m(sd, m, nc, E, O, Um, nf, cscale)
   type (sphere), intent(in) :: sd
   integer, intent(in) :: m, nc, nf
   real, intent(out) :: E(sd%nhsp, nc), O(sd%nhsp, nc)
   complex, intent(in) :: Um(0:sd%ntrunc, nf)
   real, intent(in) :: cscale
   real :: y(0:sd%ntrunc-m+4, 4), q1(sd%nhsp)
   integer :: L, lenm, j0, nh, n, k, f, nm0
   complex :: uu
   L = sd%ntrunc; lenm = L-m+1; nm0 = sd%nmstrt_tab(m+1)+1
   j0 = sd%otf_jb0(m); nh = sd%nhsp
   if (j0 > 1) then
     E(1:min(j0-1,nh), :) = 0.0; O(1:min(j0-1,nh), :) = 0.0
   endif
   if (j0 > nh) return
   n = nh-j0+1
   if (nf >= sht_otf_gemm) then
     call otf_syn_gemm(sd, m, nc, E, O, Um, nf, cscale, j0, n)
     return
   endif
   q1(1:n) = sd%otf_x(j0:nh)*sd%otf_pmm(j0:nh, m)
   f = 1
   do while (f <= nf)
     y = 0.0
     do k = 0, lenm-1
       uu = Um(m+k, f)*(cscale*sd%otf_c(nm0+k))
       y(k,1) = real(uu); y(k,2) = aimag(uu)
       if (sht_otf_pair .and. f+1 <= nf) then
         uu = Um(m+k, f+1)*(cscale*sd%otf_c(nm0+k))
         y(k,3) = real(uu); y(k,4) = aimag(uu)
       endif
     enddo
     if (sht_otf_pair .and. f+1 <= nf) then
       call otf_syn_k2(lenm, n, sd%otf_x(j0:nh), sd%otf_pmm(j0:nh,m), q1, sd%otf_be(nm0:nm0+lenm+3), &
            y(0,1), y(0,2), y(0,3), y(0,4), &
            E(j0,2*f-1), E(j0,2*f), O(j0,2*f-1), O(j0,2*f), &
            E(j0,2*f+1), E(j0,2*f+2), O(j0,2*f+1), O(j0,2*f+2))
       f = f+2
     else
       call otf_syn_k1(lenm, n, sd%otf_x(j0:nh), sd%otf_pmm(j0:nh,m), q1, sd%otf_be(nm0:nm0+lenm+3), &
            y(0,1), y(0,2), E(j0,2*f-1), E(j0,2*f), O(j0,2*f-1), O(j0,2*f))
       f = f+1
     endif
   enddo
 end subroutine otf_syn_m

 subroutine otf_syn_k1(lenm, n, x, q0, q1i, be, yr, yi, er, ei, orr, oi)
   integer, intent(in) :: lenm, n
   real, intent(in) :: x(n), q0(n), q1i(n), be(0:lenm+3), yr(0:lenm+3), yi(0:lenm+3)
   real, intent(out) :: er(n), ei(n), orr(n), oi(n)
   real :: p1(n), p2(n), a, b, c, d
   integer :: j, k
   do j = 1, n
     p2(j) = q0(j); p1(j) = q1i(j)
     er(j) = q0(j)*yr(0); ei(j) = q0(j)*yi(0)
     orr(j) = q1i(j)*yr(1); oi(j) = q1i(j)*yi(1)
   enddo
   do k = 2, lenm-1, 4
     !$omp simd private(a,b,c,d)
     do j = 1, n
       a = x(j)*p1(j) + be(k)*p2(j)
       b = x(j)*a + be(k+1)*p1(j)
       c = x(j)*b + be(k+2)*a
       d = x(j)*c + be(k+3)*b
       er(j)  = er(j)  + a*yr(k)   + c*yr(k+2)
       ei(j)  = ei(j)  + a*yi(k)   + c*yi(k+2)
       orr(j) = orr(j) + b*yr(k+1) + d*yr(k+3)
       oi(j)  = oi(j)  + b*yi(k+1) + d*yi(k+3)
       p2(j) = c; p1(j) = d
     enddo
   enddo
 end subroutine otf_syn_k1

 subroutine otf_syn_k2(lenm, n, x, q0, q1i, be, yr, yi, zr, zi, er, ei, orr, oi, fr, fi, pr, pi)
   integer, intent(in) :: lenm, n
   real, intent(in) :: x(n), q0(n), q1i(n), be(0:lenm+3)
   real, intent(in) :: yr(0:lenm+3), yi(0:lenm+3), zr(0:lenm+3), zi(0:lenm+3)
   real, intent(out) :: er(n), ei(n), orr(n), oi(n), fr(n), fi(n), pr(n), pi(n)
   real :: p1(n), p2(n), a, b, c, d
   integer :: j, k
   do j = 1, n
     p2(j) = q0(j); p1(j) = q1i(j)
     er(j) = q0(j)*yr(0); ei(j) = q0(j)*yi(0)
     orr(j) = q1i(j)*yr(1); oi(j) = q1i(j)*yi(1)
     fr(j) = q0(j)*zr(0); fi(j) = q0(j)*zi(0)
     pr(j) = q1i(j)*zr(1); pi(j) = q1i(j)*zi(1)
   enddo
   do k = 2, lenm-1, 4
     !$omp simd private(a,b,c,d)
     do j = 1, n
       a = x(j)*p1(j) + be(k)*p2(j)
       b = x(j)*a + be(k+1)*p1(j)
       c = x(j)*b + be(k+2)*a
       d = x(j)*c + be(k+3)*b
       er(j)  = er(j)  + a*yr(k)   + c*yr(k+2)
       ei(j)  = ei(j)  + a*yi(k)   + c*yi(k+2)
       orr(j) = orr(j) + b*yr(k+1) + d*yr(k+3)
       oi(j)  = oi(j)  + b*yi(k+1) + d*yi(k+3)
       fr(j)  = fr(j)  + a*zr(k)   + c*zr(k+2)
       fi(j)  = fi(j)  + a*zi(k)   + c*zi(k+2)
       pr(j)  = pr(j)  + b*zr(k+1) + d*zr(k+3)
       pi(j)  = pi(j)  + b*zi(k+1) + d*zi(k+3)
       p2(j) = c; p1(j) = d
     enddo
   enddo
 end subroutine otf_syn_k2

 !---- hybrid: regenerate P_lm of one m into an L2-resident block, then DGEMM
 ! Qe(n, ne) holds c_n-scaled P of the even degrees m, m+2, ..., Qo(n, no)
 ! the odd ones, i.e. exactly the table backend's pa/ps block for this m
 ! restricted to the non-skipped latitudes.
 subroutine otf_gen(lenm, n, x, q0, be, cc, Qe, Qo)
   integer, intent(in) :: lenm, n
   real, intent(in) :: x(n), q0(n), be(0:lenm+3), cc(0:lenm+3)
   real, intent(out) :: Qe(n, (lenm+1)/2+2), Qo(n, lenm/2+2)
   real :: p1(n), p2(n), a, b, c, d, c0, c1
   integer :: j, k, ke
   c1 = cc(1)
   do j = 1, n
     p2(j) = q0(j); p1(j) = x(j)*q0(j)
     Qe(j,1) = q0(j); Qo(j,1) = c1*p1(j)
   enddo
   ke = 2
   do k = 2, lenm-1, 4
     c0 = cc(k); c1 = cc(k+1); a = cc(k+2); b = cc(k+3)
     !$omp simd private(c,d)
     do j = 1, n
       c = x(j)*p1(j) + be(k)*p2(j)
       d = x(j)*c + be(k+1)*p1(j)
       Qe(j,ke) = c0*c; Qo(j,ke) = c1*d
       p2(j) = x(j)*d + be(k+2)*c
       p1(j) = x(j)*p2(j) + be(k+3)*d
       Qe(j,ke+1) = a*p2(j); Qo(j,ke+1) = b*p1(j)
     enddo
     ke = ke+2
   enddo
 end subroutine otf_gen

 subroutine otf_qbuf_get(nwords)
   integer, intent(in) :: nwords
   if (allocated(otf_qbuf)) then
     if (size(otf_qbuf) >= nwords) return
     deallocate(otf_qbuf)
   endif
   allocate(otf_qbuf(nwords))
 end subroutine otf_qbuf_get

 subroutine otf_ana_gemm(sd, m, nc, E, O, Um, nf, j0, n)
   type (sphere), intent(in) :: sd
   integer, intent(in) :: m, nc, nf, j0, n
   real, intent(in) :: E(sd%nhsp, nc), O(sd%nhsp, nc)
   complex, intent(inout) :: Um(0:sd%ntrunc, nf)
   real :: ye((sd%ntrunc-m+2)/2+2, 8), yo((sd%ntrunc-m+2)/2+2, 8)
   integer :: L, lenm, ne, no, nm0, f, k, lde
   L = sd%ntrunc; lenm = L-m+1; nm0 = sd%nmstrt_tab(m+1)+1
   ne = (lenm+1)/2; no = lenm/2; lde = (L-m+2)/2+2
   call otf_qbuf_get(n*((lenm+1)/2+2 + lenm/2+2))
   call otf_gen(lenm, n, sd%otf_x(j0:j0+n-1), sd%otf_pmm(j0:j0+n-1,m), sd%otf_be(nm0:nm0+lenm+3), &
                sd%otf_c(nm0:nm0+lenm+3), otf_qbuf(1:), otf_qbuf(n*((lenm+1)/2+2)+1:))
   call dgemm('T','N', ne, nc, n, 1.0d0, otf_qbuf(1), n, E(j0,1), sd%nhsp, 0.0d0, ye, lde)
   if (no > 0) call dgemm('T','N', no, nc, n, 1.0d0, otf_qbuf(n*((lenm+1)/2+2)+1), n, &
                          O(j0,1), sd%nhsp, 0.0d0, yo, lde)
   do f = 1, nf
     do k = 1, ne
       Um(m+2*k-2, f) = cmplx(ye(k,2*f-1), ye(k,2*f))
     enddo
     do k = 1, no
       Um(m+2*k-1, f) = cmplx(yo(k,2*f-1), yo(k,2*f))
     enddo
   enddo
 end subroutine otf_ana_gemm

 subroutine otf_syn_gemm(sd, m, nc, E, O, Um, nf, cscale, j0, n)
   type (sphere), intent(in) :: sd
   integer, intent(in) :: m, nc, nf, j0, n
   real, intent(out) :: E(sd%nhsp, nc), O(sd%nhsp, nc)
   complex, intent(in) :: Um(0:sd%ntrunc, nf)
   real, intent(in) :: cscale
   real :: ye((sd%ntrunc-m+2)/2+2, 8), yo((sd%ntrunc-m+2)/2+2, 8)
   integer :: L, lenm, ne, no, nm0, f, k, lde
   complex :: uu
   L = sd%ntrunc; lenm = L-m+1; nm0 = sd%nmstrt_tab(m+1)+1
   ne = (lenm+1)/2; no = lenm/2; lde = (L-m+2)/2+2
   call otf_qbuf_get(n*((lenm+1)/2+2 + lenm/2+2))
   call otf_gen(lenm, n, sd%otf_x(j0:j0+n-1), sd%otf_pmm(j0:j0+n-1,m), sd%otf_be(nm0:nm0+lenm+3), &
                sd%otf_c(nm0:nm0+lenm+3), otf_qbuf(1:), otf_qbuf(n*((lenm+1)/2+2)+1:))
   do f = 1, nf
     do k = 1, ne
       uu = Um(m+2*k-2, f)*cscale; ye(k,2*f-1) = real(uu); ye(k,2*f) = aimag(uu)
     enddo
     do k = 1, no
       uu = Um(m+2*k-1, f)*cscale; yo(k,2*f-1) = real(uu); yo(k,2*f) = aimag(uu)
     enddo
   enddo
   call dgemm('N','N', n, nc, ne, 1.0d0, otf_qbuf(1), n, ye, lde, 0.0d0, E(j0,1), sd%nhsp)
   if (no > 0) then
     call dgemm('N','N', n, nc, no, 1.0d0, otf_qbuf(n*((lenm+1)/2+2)+1), n, yo, lde, 0.0d0, O(j0,1), sd%nhsp)
   else
     O(j0:j0+n-1, :) = 0.0
   endif
 end subroutine otf_syn_gemm
