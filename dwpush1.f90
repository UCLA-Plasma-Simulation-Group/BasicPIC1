!-----------------------------------------------------------------------
!
      module wpush1
!
! Fortran90 wrappers to 1d PIC library push1.f
! Using double precision reals
! wdistr1 calculates initial particle co-ordinates with uniform plasma
!         maxwellian velocities
!         calls DISTR1
! wgpush1 push particles with leap-frog
!         calls GPUSH1L
! wgmpush1 push particles with leap-frog using OpenMP
!          calls MGPUSH1L
! wgpost1 deposit charge
!         calls GPOST1L
! wcguard1 copy guard cells for periodic 1d scalar data
!          calls CGUARD1L
! waguard1 add guard cells for periodic 1d scalar data
!          calls AGUARD1L
! pois1_init calculates table needed by 1d poisson solver
!            calls POIS1
! wpois1 poisson solver for periodic 1d electric field
!        calls POIS1
! fft1_init calculates tables needed by 1d FFTs
!           calls WFFT1RINIT
! wfft1r wrapper function for scalar 1d real/complex FFT
!        calls FFT1RXX
! meaddext1 add external traveling wave field to electric field for
!           1d code
!           calls EADDEXT1
! fprecision determines if default reals are actually doubles
! written by viktor k. decyk, ucla
! copyright 2021, regents of the university of california
! update: december 18, 2023
!
      use push1_h
      implicit none
!
! t = scratch array for mfft1rn
      double complex, dimension(:), allocatable :: t
      integer :: szt = 0
      save
!
      private :: t, szt
!
!
      contains
!
!-----------------------------------------------------------------------
      subroutine wdistr1(part,vtx,vdx,npx,nx)
! calculates initial particle co-ordinates with uniform plasma and
! maxwellian velocities
      implicit none
      integer, intent(in) :: npx, nx
      double precision, intent(in) :: vtx, vdx
      double precision, dimension(:,:), intent(inout) :: part
! local data
      integer :: ipbc = 1
      integer :: idimp, nop
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
! call low level procedure
      call DISTR1(part,vtx,vdx,npx,idimp,nop,nx,ipbc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wgpush1(part,fx,qbm,dt,ek,nx)
! push particles with leap-frog
      implicit none
      integer, intent(in) :: nx
      double precision, intent(in) :: qbm, dt
      double precision, dimension(1), intent(inout) :: ek
      double precision, dimension(:,:), intent(inout) :: part
      double precision, dimension(:), intent(in) :: fx
! local data
      integer :: ipbc = 1
      integer :: idimp, nop, nxv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(fx,1)
! call low level procedure
      call GPUSH1L(part,fx,qbm,dt,ek(1),idimp,nop,nx,nxv,ipbc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wgmpush1(part,fx,qbm,dt,ek,nvp,nx)
! push particles with leap-frog
      implicit none
      integer, intent(in) :: nvp, nx
      double precision, intent(in) :: qbm, dt
      double precision, dimension(1), intent(inout) :: ek
      double precision, dimension(:,:), intent(inout) :: part
      double precision, dimension(:), intent(in) :: fx
! local data
      integer :: ipbc = 1
      integer :: idimp, nop, nxv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(fx,1)
! call low level procedure
      call MGPUSH1L(part,fx,qbm,dt,ek(1),idimp,nop,nvp,nx,nxv,ipbc)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wgpost1(part,q,qm)
! deposit charge
      implicit none
      double precision, intent(in) :: qm
      double precision, dimension(:,:), intent(in) :: part
      double precision, dimension(:), intent(inout) :: q
! local data
      integer :: idimp, nop, nxv
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nxv = size(q,1)
! call low level procedure
      call GPOST1L(part,q,qm,nop,idimp,nxv)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wcguard1(fx,nx)
! copy guard cells for periodic 1d scalar data
      implicit none
      integer, intent(in) :: nx
      double precision, dimension(:), intent(inout) :: fx
! local data
      integer :: nxe
! extract dimensions
      nxe = size(fx,1)
! call low level procedure
      call CGUARD1L(fx,nx,nxe)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine waguard1(q,nx)
! add guard cells for periodic 1d scalar data
      implicit none
      integer, intent(in) :: nx
      double precision, dimension(:), intent(inout) :: q
! local data
      integer :: nxe
! extract dimensions
      nxe = size(q,1)
! call low level procedure
      call AGUARD1L(q,nx,nxe)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine pois1_init(ffc,ax,affp,nx)
! calculates table needed by 1d poisson solver
      implicit none
      integer, intent(in) :: nx
      double precision, intent(in) :: ax, affp
      double complex, dimension(:), intent(inout) :: ffc
! local data
      integer :: isign = 0
      double precision, dimension(1) :: we
      double precision, dimension(1)  :: q, fx
! call low level procedure
      call POIS1(q,fx,isign,ffc,ax,affp,we(1),nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wpois1(q,fx,ffc,we,nx)
! poisson solver for periodic 1d electric field
      implicit none
      integer, intent(in) :: nx
      double precision, dimension(1), intent(inout) :: we
      double precision, dimension(:), intent(in)  :: q
      double precision, dimension(:), intent(inout) :: fx
      double complex, dimension(:), intent(inout) :: ffc
! local data
      integer :: isign = -1
      double precision :: ax, affp
! call low level procedure
      call POIS1(q,fx,isign,ffc,ax,affp,we(1),nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fft1_init(mixup,sct,indx)
! calculates tables needed by 1d FFTs
      implicit none
      integer, intent(in) :: indx
      integer, dimension(:), intent(inout) :: mixup
      double complex, dimension(:), intent(inout) :: sct
! local data
      integer :: nxhd
! extract dimensions
      nxhd = size(mixup,1)
! call low level procedure
      call WFFT1RINIT(mixup,sct,indx,nxhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wfft1r(f,isign,mixup,sct,indx)
! wrapper function for scalar 1d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx
      double precision, dimension(:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      double complex, dimension(:), intent(in) :: sct
! local data
      integer :: nxd, nxhd
! extract dimensions
      nxd = size(f,1)
      nxhd = size(mixup,1)
! check if required size of buffer t has increased
      if (szt < nxhd) then
         if (szt /= 0) deallocate(t)
! allocate new buffers
         allocate(t(nxhd))
         szt = nxhd
      endif
! call low level procedure
      call FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine meaddext1(fxe,amodex,freq,time,trmp,toff,el0,er0,nx)
! add external traveling wave field to electric field for 1-2/2d code
      implicit none
      integer, intent(in) :: nx
      double precision, intent(in) :: amodex, freq, time, trmp, toff
      double precision, intent(in) :: el0, er0
      double precision, dimension(:), intent(inout) :: fxe
! local data
      integer :: nxe
! extract dimensions
      nxe = size(fxe,1)
! call low level procedure
      call EADDEXT1(fxe,amodex,freq,time,trmp,toff,el0,er0,nx,nxe)
      end subroutine
!
!-----------------------------------------------------------------------
      integer function fprecision()
! function determines if default reals are actually doubles
      implicit none
      double precision :: prec
! ndprec = (0,1) = (no,yes) = (normal,autodouble) precision used
      if (digits(prec) > 24) then
         fprecision = 1
      else
         fprecision = 0
      endif
      end function fprecision
!
      end module
