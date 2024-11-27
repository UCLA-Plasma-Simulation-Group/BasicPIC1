!-----------------------------------------------------------------------
! VMXT1 calculates maximum of absolute value of particle velocities and
!       energy, and scales to next highest power of dv
! VDIST1 calculates 1 component velocity distribution, velocity moments,
!        and entropy.
!-----------------------------------------------------------------------
      subroutine VMXT1(part,fvmx,idimp,nop)
! for 1d code, this subroutine calculates maximum of absolute value of
! particle velocities and energy, and scales to next highest power of dv
! input: all except fvmx, output: fvmx
! part(2,n) = velocity vx of particle n 
! fvmx = maximum of absolute values of velocities and energy
! idimp = size of phase space = 2
! nop = number of particles
      implicit none
      integer idimp, nop
      real part, fvmx
      dimension part(idimp,nop), fvmx(2)
! local data
      integer i, j, iv
      real dv, alg2, vx, v2, at1, at2
      double precision svmx, svm2
      dv = 2.0
      alg2 = 1.0/alog(dv)
      svmx = 0.0d0
      svm2 = 0.0d0
! loop over particles
      do 10 j = 1, nop
      vx = abs(part(2,j))
      v2 = vx*vx
      svmx = amax1(svmx,vx)
      svm2 = amax1(svm2,v2)
   10 continue
      fvmx(1) = svmx
      fvmx(2) = svm2
      do 20 i = 1, 2
      at1 = fvmx(i)
      iv = alog(at1)*alg2
      if (at1.ge.1.0) iv = iv + 1
      at2 = dv**iv
      if (at1.gt.at2) at2 = 2.0*at2
      fvmx(i) = at2
   20 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine VDIST1(part,fv,fvm,idimp,nop,nmv,nmvf)
! for 1d code, this subroutine calculates 1d velocity distribution,
! velocity moments, and entropy
! input: all except fvm, output: fv, fvm
! part(2,n) = velocity vx of particle n
! fv = distribution function, number of particles in each velocity range
! maximum velocity (used for scaling) is contained in last of element fv
! vdrift is contained in fvm(1)
! vth is contained in fvm(2)
! entropy is contained in fvm(3), defined to be:
! s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
! is uniform in space
! idimp = size of phase space = 2
! nop = number of particles
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nop, nmv, nmvf
      real part, fv, fvm
      dimension part(idimp,nop), fv(nmvf), fvm(3)
! local data
      integer j, nmv21, nvx
      real anmv, svx, vx
      double precision sumvx, sumvx2, anp
! velocity scaling
      nmv21 = 2*nmv + 1
      anmv = real(nmv)
      svx = anmv/fv(nmv21+1)
! zero out distribution
      do 10 j = 1, nmv21
      fv(j) = 0.0
   10 continue
! count particles in each velocity region
      anmv = anmv + 1.5
      sumvx = 0.0d0
      sumvx2 = 0.0d0
      do 20 j = 1, nop
      vx = part(2,j)
      nvx = vx*svx + anmv
      sumvx = sumvx + vx
      sumvx2 = sumvx2 + vx*vx
      if ((nvx.ge.1).and.(nvx.le.nmv21)) fv(nvx) = fv(nvx) + 1.0
   20 continue
! calculate velocity moments
      anp = 0.0d0
      if (nop.gt.0) anp = 1.0d0/dble(nop)
      sumvx = sumvx*anp
      fvm(1) = sumvx
      fvm(2) = dsqrt(sumvx2*anp - sumvx**2)
! calculate entropy
      sumvx = 0.0d0
      sumvx2 = 0.0d0
      do 30 j = 1, nmv21
      if (fv(j).gt.0.0) then
         sumvx = sumvx + fv(j)
         sumvx2 = sumvx2 + fv(j)*dlog(dble(fv(j)*svx))
      endif
   30 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      fvm(3) = sumvx
      return
      end
