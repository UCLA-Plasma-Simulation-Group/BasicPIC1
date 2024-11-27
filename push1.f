! Fortran Library for Basic 1D Electrostatic PIC Code
! written by Viktor K. Decyk, UCLA
!-----------------------------------------------------------------------
      subroutine DISTR1(part,vtx,vdx,npx,idimp,nop,nx,ipbc)
! for 1d code, this subroutine calculates initial particle co-ordinate
! and velocity, with uniform density and maxwellian velocity with drift
! part(1,n) = position x of particle n
! part(2,n) = velocity vx of particle n
! vtx = thermal velocity of particles in x direction
! vdx = drift velocity of particles x direction
! npx = number of particles distributed in x direction
! idimp = size of phase space = 2
! nop = number of particles
! nx = system length in x direction
! ipbc = particle boundary condition = (0,1,2) =
! (none,2d periodic,2d reflecting)
! ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, idimp, nop, nx, ipbc
      real part, vtx, vdx
      dimension part(idimp,nop)
! local data
      integer j
      real edgelx, at1, sum1
      double precision dsum1
      double precision ranorm
! set boundary values
      edgelx = 0.0
      at1 = real(nx)/real(npx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
! uniform density profile
      do 10 j = 1, npx
      part(1,j) = edgelx + at1*(real(j) - .5)
   10 continue
! maxwellian velocity distribution
      do 20 j = 1, npx
      part(2,j) = vtx*ranorm()
   20 continue
! add correct drift
      dsum1 = 0.0d0
      do 30 j = 1, npx
      dsum1 = dsum1 + part(2,j)
   30 continue
      sum1 = dsum1
      sum1 = sum1/real(npx) - vdx
      do 40 j = 1, npx
      part(2,j) = part(2,j) - sum1
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine GPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
! for 1d code, this subroutine updates particle co-ordinate and velocity
! using leap-frog scheme in time and first-order linear interpolation
! in space, with various boundary conditions.
! scalar version using guard cells
! 16 flops/particle, 4 loads, 2 stores
! input: all, output: part, ek
! equations used are:
! v(t+dt/2) = v(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
! and x(t+dt) = x(t) + v(t+dt/2)*dt
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! part(1,n) = position x of particle n
! part(2,n) = velocity vx of particle n
! fx(j) = force/charge at grid point j, that is convolution of electric
! field over particle shape
! qbm = particle charge/mass
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .125*sum((v(t+dt/2)+v(t-dt/2))**2)
! idimp = size of phase space = 2
! nop = number of particles
! nx = system length in x direction
! nxv = first dimension of charge array, must be >= nx+1
! ipbc = particle boundary condition = (0,1,2) =
! (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fx, qbm, dt, ek
      dimension part(idimp,nop), fx(nxv)
! local data
      integer j, nn
      real qtm, edgelx, edgerx, dx
      double precision sum1
      qtm = qbm*dt
      sum1 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop
! find interpolation weights
      nn = part(1,j)
      dx = part(1,j) - real(nn)
      nn = nn + 1
! find acceleration
      dx = (1.0 - dx)*fx(nn) + dx*fx(nn+1)
! new velocity
      dx = part(2,j) + qtm*dx
! average kinetic energy
      sum1 = sum1 + (part(2,j) + dx)**2
      part(2,j) = dx
! new position
      dx = part(1,j) + dx*dt
! periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
! reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(2,j) = -part(2,j)
         endif
      endif
! set new position
      part(1,j) = dx
   10 continue
! normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MGPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nvp,nx,nxv,ipbc)
! for 1d code, this subroutine updates particle co-ordinate and velocity
! using leap-frog scheme in time and first-order linear interpolation
! in space, with various boundary conditions.
! OpenMP version using guard cells and particle decomposition
! 16 flops/particle, 4 loads, 2 stores
! input: all, output: part, ek
! equations used are:
! v(t+dt/2) = v(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
! and x(t+dt) = x(t) + v(t+dt/2)*dt
! fx(x(t)) is approximated by interpolation from the nearest grid points
! fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
! where n = nearest grid point and dx = x-n
! part(1,n) = position x of particle n
! part(2,n) = velocity vx of particle n
! fx(j) = force/charge at grid point j, that is convolution of electric
! field over particle shape
! qbm = particle charge/mass
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .125*sum((v(t+dt/2)+v(t-dt/2))**2)
! idimp = size of phase space = 2
! nop = number of particles
! nvp = number of real or virtual processors
! nx = system length in x direction
! nxv = first dimension of charge array, must be >= nx+1
! ipbc = particle boundary condition = (0,1,2) =
! (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nvp, nx, nxv, ipbc
      real part, fx, qbm, dt, ek
      dimension part(idimp,nop), fx(nxv)
! local data
      integer j, k, npp, joff, nps, nn
      real qtm, edgelx, edgerx, x, vx, dx
      double precision sum1, sum2
      npp = (nop - 1)/nvp + 1
      qtm = qbm*dt
      sum2 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
!$OMP PARALLEL DO                                                       &
!$OMP& PRIVATE(j,k,joff,nps,nn,x,vx,dx,sum1) REDUCTION(+:sum2)
      do 20 k = 1, nvp
      joff = npp*(k - 1)
      nps = min(npp,max(0,nop-joff))
      sum1 = 0.0d0
      do 10 j = 1, nps
! find interpolation weights
      x = part(1,j+joff)
      vx = part(2,j+joff)
      nn = x
      dx = x - real(nn)
      nn = nn + 1
! find acceleration
      dx = (1.0 - dx)*fx(nn) + dx*fx(nn+1)
! new velocity
      dx = vx + qtm*dx
! average kinetic energy
      sum1 = sum1 + (vx + dx)**2
      part(2,j+joff) = dx
! new position
      dx = x + dx*dt
! periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
! reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(2,j+joff) = -part(2,j+joff)
         endif
      endif
! set new position
      part(1,j+joff) = dx
   10 continue
      sum2 = sum2 + sum1
   20 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine GPOST1L(part,q,qm,nop,idimp,nxv)
! for 1d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! scalar version using guard cells
! 7 flops/particle, 3 loads, 3 stores
! input: all, output: q
! charge density is approximated by values at the nearest grid points
! q(n)=qm*(1.-dx) and q(n+1)=qm*dx
! where n = nearest grid point and dx = x-n
! part(1,n) = position x of particle n
! q(j) = charge density at grid point j
! qm = charge on particle, in units of e
! nop = number of particles
! idimp = size of phase space = 2
! nxv = first dimension of charge array, must be >= nx+1
      implicit none
      integer nop, idimp, nxv
      real part, q, qm
      dimension part(idimp,nop), q(nxv)
! local data
      integer j, nn
      real dx
! find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) 
      dx = qm*(part(1,j) - real(nn))
      nn = nn + 1
! deposit charge
      q(nn) = q(nn) + (qm - dx)
      q(nn+1) = q(nn+1) + dx
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MGPOST1L(part,q,qp,qm,nop,nvp,idimp,nxv)
! for 1d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! OpenMP version using guard cells and particle decomposition
! 7 flops/particle, 3 loads, 3 stores
! input: all, output: q, qp
! charge density is approximated by values at the nearest grid points
! q(n)=qm*(1.-dx) and q(n+1)=qm*dx
! where n = nearest grid point and dx = x-n
! part(1,n) = position x of particle n
! q(j) = charge density at grid point j
! qp(j,k) = charge density at grid point j and partition k
! qm = charge on particle, in units of e
! nop = number of particles
! nvp = number of real or virtual processors
! idimp = size of phase space = 2
! nxv = first dimension of charge array, must be >= nx+1
      implicit none
      integer nop, nvp, idimp, nxv
      real part, q, qp, qm
      dimension part(idimp,nop), q(nxv), qp(nxv,nvp)
! local data
      integer j, k, npp, joff, nps, nn
      real x, dx
      npp = (nop - 1)/nvp + 1
!$OMP PARALLEL DO PRIVATE(j,k,joff,nps,nn,x,dx)
      do 30 k = 1, nvp
      joff = npp*(k - 1)
      nps = min(npp,max(0,nop-joff))
      do 10 j = 1, nxv
      qp(j,k) = 0.0
   10 continue
! find interpolation weights
      do 20 j = 1, nps
      x = part(1,j+joff)
      nn = x 
      dx = qm*(x - real(nn))
      nn = nn + 1
! deposit charge
      qp(nn,k) = qp(nn,k) + (qm - dx)
      qp(nn+1,k) = qp(nn+1,k) + dx
   20 continue
   30 continue
!$OMP END PARALLEL DO
! sum up the partitions
      do 50 k = 1, nvp
      do 40 j = 1, nxv
      q(j) = q(j) + qp(j,k)
   40 continue
   50 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine CGUARD1L(fx,nx,nxe)
! replicate extended periodic field fx
! linear interpolation
! nx = system length in x direction
! nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real fx
      integer nx, nxe
      dimension fx(nxe)
! copy edge of extended field
      fx(nx+1) = fx(1)
      return
      end
!-----------------------------------------------------------------------
      subroutine AGUARD1L(q,nx,nxe)
! accumulate extended periodic field
! linear interpolation
! nx = system length in x direction
! nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real q
      integer nx, nxe
      dimension q(nxe)
! accumulate edge of extended field
      q(1) = q(1) + q(nx+1)
      q(nx+1) = 0.0
      return
      end
!-----------------------------------------------------------------------
      subroutine POIS1(q,fx,isign,ffc,ax,affp,we,nx)
! this subroutine solves 1d poisson's equation in fourier space for
! force/charge (or convolution of electric field over particle shape)
! with periodic boundary conditions.
! for isign = 0, input: isign,ax,affp,nx, output: ffc
! for isign  /= 0, input: q,ffc,isign,nx, output: fx,we
! approximate flop count is: 6*nx
! equation used is:
! fx(k) = -sqrt(-1)*k*g(k)*s(k)*q(k), where k = 2pi*j/nx, j=fourier mode,
! g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2), except for
! fx(k=0) = fx(k=pi) = 0.
! cmplx(q(2*j-1),q(2*j)) = complex charge density for fourier mode j-1
! cmplx(fx(2*j-1),fx(2*j)) = complex force/charge for fourier mode j-1
! if isign = 0, form factor array is prepared
! if isign is not equal to 0, force/charge is calculated
! ffc(2*j) = finite-size particle shape factor s for fourier mode j-1
! ffc(2*j-1) = potential green's function g for fourier mode j-1
! ax = half-width of particle in x direction
! affp = normalization constant = nx/np, where np = number of particles
! electric field energy is also calculated, using
! we = nx*sum((affp/k**2)*|q(k)*s(k)|**2)
! nx = system length in x direction
      implicit none
      integer isign, nx
      real ax, affp, we
      real q, fx, ffc
      dimension q(nx), fx(nx), ffc(nx)
! local data
      integer j, nxh
      real dnx, dkx, at1, at2
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      if (isign.ne.0) go to 20
! prepare form factor array
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      ffc(2*j) = exp(-.5*(dkx*ax)**2)
      ffc(2*j-1) = affp*ffc(2*j)/(dkx*dkx)
   10 continue
      ffc(1) = affp
      ffc(2) = 1.0
      return
! calculate force/charge and sum field energy
   20 wp = 0.0d0
! mode numbers 0 < kx < nx/2
      do 30 j = 2, nxh
      at1 = ffc(2*j-1)*ffc(2*j)
      at2 = dnx*real(j - 1)*at1
      fx(2*j-1) = at2*q(2*j)
      fx(2*j) = -at2*q(2*j-1)
      wp = wp + at1*(q(2*j-1)**2 + q(2*j)**2)
   30 continue
! mode number kx = 0
      fx(1) = 0.
      fx(2) = 0.
      we = real(nx)*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine WFFT1RINIT(mixup,sct,indx,nxhd)
! this subroutine calculates tables needed by a one dimensional
! real to complex fast fourier transform and its inverse.
! input: indx, nxhd
! output: mixup, sct
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! indx = exponent which determines length in x direction,
! where nx=2**indx
! nxhd = nx/2
! written by viktor k. decyk, ucla
      implicit none
      integer indx, nxhd
      integer mixup
      complex sct
      dimension mixup(nxhd), sct(nxhd)
! local data
      integer indx1, nx, nxh
      integer j, k, lb, ll, jb, it
      real dnx, arg
      indx1 = indx - 1
      nx = 2**indx
      nxh = nx/2
! bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxh
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
! sine/cosine table for the angles 2*n*pi/nx
      dnx = 6.28318530717959/real(nx)
      do 30 j = 1, nxh
      arg = dnx*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
! this subroutine performs a one dimensional real to complex fast
! fourier transform and its inverse, using complex arithmetic
! for isign = (-1,1), input: all except t, output: f, t
! for isign = -1, approximate flop count: N*(5*log2(N) + 10)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = nx/2
! f = input and output data
! t = complex scratch array
! indx = power of 2 which determines length of transform, nx = 2**indx
! if isign = -1, an inverse fourier transform is performed
! f(n) = (1/nx)*sum(f(j)*exp(-sqrt(-1)*2pi*n*j/nx))
! if isign = 1, a forward fourier transform is performed
! f(j) = sum(f(n)*exp(sqrt(-1)*2pi*n*j/nx))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! fourier coefficients are stored as follows:
! f(1) = real part of mode 0, f(2) = real part of mode nx/2
! f(2*j-1),f(2*j) = real,imaginary part of mode j-1, 0 < j < nx/2
! written by viktor k. decyk, ucla
! scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sct
      dimension f(nxd), t(nxhd), mixup(nxhd), sct(nxhd)
! local data
      integer indx1, nx, nxh, nxhh, nxh2, j, k, l, j1
      integer nxs, nxs2, km, km2, k1, k2
      real ani
      complex t1, t2
      if (isign.eq.0) return
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      indx1 = indx - 1
      if (isign.gt.0) go to 70
! inverse fourier transform
! bit-reverse array elements to complex temporary
      do 10 j = 1, nxh
      j1 = mixup(j)
      t(j) = cmplx(f(2*j1-1),f(2*j1))
   10 continue
! transform
      do 40 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 30 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 20 j = 1, nxs
      t1 = sct(1+km2*(j-1))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
   20 continue
   30 continue
   40 continue
! unscramble coefficients and normalize
      ani = 1./real(2*nx)
      do 50 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sct(j)),-real(sct(j)))
      t(j) = ani*(t1 + t2)
      t(nxh2-j) = ani*conjg(t1 - t2)
   50 continue
      ani = 2.*ani
      t(nxhh+1) = ani*conjg(t(nxhh+1))
      t(1) = ani*cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1))
     1)
! move to real destination
      do 60 j = 1, nxh
      f(2*j-1) = real(t(j))
      f(2*j) = aimag(t(j))
   60 continue
      return
! forward fourier transform
! move to complex temporary
   70 do 80 j = 1, nxh
      t(j) = cmplx(f(2*j-1),f(2*j))
   80 continue
! scramble coefficients
      do 90 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sct(j)),real(sct(j)))
      t(j) = t1 + t2
      t(nxh2-j) = conjg(t1 - t2)
   90 continue
      t(nxhh+1) = 2.*conjg(t(nxhh+1))
      t(1) = cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1)))
! bit-reverse array elements to real destination
      do 100 j = 1, nxh
      j1 = mixup(j)
      f(2*j-1) = real(t(j1))
      f(2*j) = aimag(t(j1))
  100 continue
! move back to complex temporary
      do 110 j = 1, nxh
      t(j) = cmplx(f(2*j-1),f(2*j))
  110 continue
! transform
      do 140 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 130 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 120 j = 1, nxs
      t1 = conjg(sct(1+km2*(j-1)))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
  120 continue
  130 continue
  140 continue
! move to real destination
      do 150 j = 1, nxh
      f(2*j-1) = real(t(j))
      f(2*j) = aimag(t(j))
  150 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine EADDEXT1(fxe,amodex,freq,time,trmp,toff,el0,er0,nx,nxe)
! add external traveling wave field to electric field
! fxe(x) = fxe(x) + e1*sin(k0*x + freq*time) + e2*cos(k0*x - freq*time)
! where
!     e1 = el0*(time/trmp), e2 = er0*(time/trmp), if time < trmp
!     e1 = el0,             e2 = er0,             if trmp < time < toff
!     e1 = el0*((toff+trmp-time)/trmp), e2 = er0*((toff+trmp-time)/trmp)
!                                            if toff < time < toff+trmp
!     e1 = 0,               e2 = 0,               if time > toff+trmp
! if toff < 0 => toff = + infinity
      implicit none
      integer nx, nxe
      real amodex, freq, time, trmp, toff, el0, er0
      real fxe
      dimension fxe(nxe)
! local data
      integer j
      real tr, at, ft, dkx, xk
      if ((el0==0.0).and.(er0==0.0)) return
      tr = toff + trmp
! ramp up
      if (time < trmp) then
         at = time/trmp
      else if ((toff >= 0.0).and.(time > toff)) then
! ramp down
         if (time < tr) then
            at = (tr - time)/trmp
! shutdown
         else
            at = 0.0
         endif
! constant amplitude
      else
         at = 1.0
      endif
      ft = freq*time
      dkx = (6.28318530717959*amodex)/real(nx)
      if (at > 0.0) then
         do 10 j = 1, nx
         xk = dkx*real(j - 1)
         fxe(j) = fxe(j) + at*(er0*cos(xk - ft) + el0*sin(xk + ft))
   10    continue
      endif
      return
      end
!-----------------------------------------------------------------------
      function ranorm()
! this program calculates a random number y from a gaussian distribution
! with zero mean and unit variance, according to the method of
! mueller and box:
!    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
!    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
! where x is a random number uniformly distributed on (0,1).
! written for the ibm by viktor k. decyk, ucla
      implicit none
      integer iflg,isc,i1,r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
!-----------------------------------------------------------------------
      function randum()
! this is a version of the random number generator dprandom due to
! c. bingham and the yale computer center, producing numbers
! in the interval (0,1).  written for the sun by viktor k. decyk, ucla
      implicit none
      integer isc,i1,r1,r2
      double precision randum,h1l,h1u,r0,r3,asc,bsc
      save r1,r2,h1l,h1u
      data r1,r2 /1271199957,1013501921/
      data h1l,h1u /65533.0d0,32767.0d0/
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      randum = (dble(r1) + dble(r2)*asc)*asc
      return
      end
!-----------------------------------------------------------------------
      subroutine NEXTRAN1(nextrand,np)
! for 1d code, this subroutine skips over nextrand groups of random
! numbers in order to initialize different random ensembles
! nextrand = (0,N) = generate (default,Nth block) of random numbers
! np = number of particles in distribution
      implicit none
      integer nextrand, np
! local data
      integer j, n
      double precision d
      double precision ranorm
      n = np*nextrand
      do 10 j = 1, n
      d = ranorm()
   10 continue
      return
      end
