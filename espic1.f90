!-----------------------------------------------------------------------
! Basic 1D Electrostatic PIC code
! written by Viktor K. Decyk, UCLA
      program espic1
      use wpush1
      use wdiag1
      use omplib
      implicit none
! indx = exponent which determines grid points in x direction:
! nx = 2**indx.
      integer, parameter :: indx =   9
! npx = number of electrons distributed in x direction.
      integer, parameter :: npx = 184320
! npxb = number of beam electrons in x direction
      integer, parameter :: npxb = 18432

! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 100.0, dt = 0.2, qme = -1.0
! vtx = thermal velocity of electrons in x direction
! vx0 = drift velocity of electrons in x direction.
      real, parameter :: vtx = 1.0, vx0 = 0.0
! vtdx = thermal velocity of beam electrons in x direction
! vdx = drift velocity of beam electrons in x direction
      real, parameter :: vtdx = 0.5, vdx = 5.0
! ax = smoothed particle size in x direction
      real :: ax = .912871
! idimp = number of particle coordinates = 2
      integer :: idimp = 2
!
! Multi-tasking Parameter:
! nvp = number of shared memory threads (0=default)
      integer :: nvp = 0
!
! External Electrostatic Traveling Wave Driver
! amodex = wave number which determines k0 = 2*pi*amodex/NX
! freq = frequency of external wave driver
! trmp = ramp-up time for external wave driver
! toff = shut-off time for external wave driver
      real :: amodex = 0.0, freq = 0.0, trmp = 0.0, toff = 0.0
! el0/er0 = external pump amplitude for left-going/right-going wave
      real :: el0 = 0.0, er0 = 0.0
!
! diagnostics parameters
! nte = number of time steps between electric field display
      integer :: nte = 10
! ntv = number of time steps between velocity-space diagnostic
! nmv = number of segments in v for velocity distribution
      integer :: ntv = 10, nmv = 40
! nts = number of time steps between phase space display
      integer :: nts = 10
! ntw = number of time steps between energy display
      integer :: ntw = 10
!
! wke/we = particle kinetic/electric field/total energy
      real, dimension(1) :: wke = 0.0, we = 0.0, wtt = 0.0
!
! declare scalars for standard code
      integer :: n
      integer :: np, nx, nxh, nxe
      integer :: ntime, nloop, isign
      real :: qbme, affp, ts
      real, dimension(1) :: wt0
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), allocatable :: part
! qe = electron charge density with guard cells
      real, dimension(:), allocatable :: qe
! fxe = smoothed electric field with guard cells
      real, dimension(:), allocatable :: fxe
! ffc = form factor array for poisson solver
      complex, dimension(:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime, ltime
      real :: tinit = 0.0, tloop = 0.0
      real :: tdpost = 0.0, tpush = 0.0
      double precision :: dtime
!
! declare data for diagnostics
! velocity diagnostics
! fv = global electron velocity distribution functions
! fvm = electron vdrift, vth, entropy for global distribution
      real, dimension(:), allocatable :: fv, fvm
! fvmx/fvms = scaled maximum of absolute values of electron velocities
      real, dimension(:), allocatable :: fvmx
! phase space diagnostics
      real, dimension(:), allocatable :: fvms
! wt = energy time history array
      integer :: mtw, itw
      real :: dtw
      real, dimension(:,:), allocatable :: wt
!
! start timing initialization
      call dtimer(dtime,itime,-1)
!
! nvp = number of shared memory threads (0=default)
!     write (*,*) 'enter number of nodes:'
!     read (5,*) nvp
! initialize for shared memory parallel processing
      call INIT_OMP(nvp)
      nvp = GETNTHSIZE()
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx = number of grid points in x direction
      np = npx + npxb; nx = 2**indx; nxh = nx/2
      nxe = nx + 2
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = int(tend/dt + .0001) + 1; ntime = 0
      qbme = qme
      affp = real(nx)/real(np)
!
! allocate data for standard code
      allocate(part(idimp,np))
      allocate(qe(nxe),fxe(nxe))
      allocate(ffc(nxh),mixup(nxh),sct(nxh))
!
! allocate data for diagnostics
! velocity diagnostics
      if (ntv > 0) allocate(fv(2*nmv+2),fvm(3),fvmx(2))
! phase space diagnostics
      if (nts > 0) allocate(fvms(2))
! energy diagnostics
      if (ntw > 0) then
         mtw = int((nloop - 1)/ntw + 1); itw = 0; dtw = dt*float(ntw)
         allocate(wt(mtw,4))
      endif
!
! prepare fft tables
      call fft1_init(mixup,sct,indx)
! calculate form factors
      call pois1_init(ffc,ax,affp,nx)
! initialize electrons
! uniform density in space, maxwellian in velocity space: updates part
      if (npx > 0) call wdistr1(part,vtx,vx0,npx,nx)
! initialize beam electrons
! uniform density in space, maxwellian in velocity space: updates part
      if (npxb > 0) call wdistr1(part(:,npx+1:np),vtdx,vdx,npxb,nx)
!
! initialization time
      call dtimer(dtime,itime,1)
      tinit = tinit + real(dtime)
! start timing loop
      call dtimer(dtime,ltime,-1)
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop
      ntime = n - 1
!     write (*,*) 'ntime = ', ntime
!
! deposit charge with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call wgpost1(part,qe,qme)
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
!
! add guard cells with standard procedure: updates qe
      call waguard1(qe,nx)
!
! transform charge to fourier space with standard procedure:
! updates qe, fxe
      isign = -1
      call wfft1r(qe,isign,mixup,sct,indx)
!
! calculate force/charge in fourier space with standard procedure:
! updates fxe, we
      call wpois1(qe,fxe,ffc,we,nx)
!
! transform force to real space with standard procedure: updates fxe, qe
      isign = 1
      call wfft1r(fxe,isign,mixup,sct,indx)
!
! add external traveling wave field: updates fxe
      ts = abs(dt)*real(ntime)
      call meaddext1(fxe,amodex,freq,ts,trmp,toff,el0,er0,nx)
!
! copy guard cells with standard procedure: updates fxe
      call wcguard1(fxe,nx)
!
! velocity diagnostic
      if (ntv > 0) then
         if (ntime==ntv*(ntime/ntv)) then
! find maximum of electron velocity and energy
            call wvmxt1(part,fvmx)
            fv(2*nmv+2) = fvmx(1)
! calculate 1d velocity distribution
            call wvdist1(part,fv,fvm,nmv)
         endif
      endif
!
! phase space diagnostic
      if (nts > 0) then
         if (ntime==nts*(ntime/nts)) then
! find maximum of electron velocity and energy
            call wvmxt1(part,fvms)
         endif
      endif
!
! push particles: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
! use parallel version
      if (nvp > 1) then
         call wgmpush1(part,fxe,qbme,dt,wke,nvp,nx)
! use serial version
      else
         call wgpush1(part,fxe,qbme,dt,wke,nx)
      endif
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
!
      if (ntime==0) then
         write (*,*) 'Initial Field, Kinetic and Total Energies:'
         wt0 = wke + we
         write (*,'(3e14.7)') we, wke, wt0
      endif
!
! energy diagnostic: updates wt
      if (ntw > 0) then
         if (ntime==ntw*(ntime/ntw)) then
            itw = itw + 1
            wt(itw,:) = (/we,wke,0.0,we/)
         endif
      endif
!
      enddo
      ntime = ntime + 1
!
! loop time
      call dtimer(dtime,ltime,1)
      tloop = tloop + real(dtime)
!
! * * * end main iteration loop * * *
!
      write (*,*) 'ntime = ', ntime
      write (*,*) 'Final Field, Kinetic and Total Energies:'
      wtt = wke + we
      write (*,'(3e14.7)') we, wke, wtt
      write (*,*)
!
      write (*,*) 'Energy Conservation(%) = ', 100.0*(wtt - wt0)/wt0
!
      write (*,*)
      write (*,*) 'initialization time = ', tinit
      write (*,*) 'deposit time = ', tdpost
      write (*,*) 'push time = ', tpush
      time = tdpost + tpush
      write (*,*) 'total particle time = ', time
      write (*,*) 'total time = ', tloop
      write (*,*)
!
      wtt = 1.0e+09/(real(nloop)*real(np))
      write (*,*) 'Push Time (nsec) = ', tpush*wtt
      write (*,*) 'Deposit Time (nsec) = ', tdpost*wtt
      write (*,*) 'Total Particle Time (nsec) = ', time*wtt
!
      stop
      end program
