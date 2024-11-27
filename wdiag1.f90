!-----------------------------------------------------------------------
!
      module wdiag1
!
! Fortran90 wrappers to 1d PIC library diag1.f
! wvmxt1 calculates maximum of absolute value of particle velocities and
!        energy
!        calls VMXT1
! wvdist1 calculates 1d velocity distribution, velocity moments, and
!         entropy
!         calls VDIST1
! written by viktor k. decyk, ucla
! copyright 2023, regents of the university of california
! update: october 11, 2023
!
      use diag1_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine wvmxt1(part,fvmx)
! calculates maximum of absolute value of particle velocities and energy
      implicit none
      real, dimension(:,:), intent(in) :: part
      real, dimension(2), intent(inout) :: fvmx
! local data
      integer :: idimp, nop
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
! call low level procedure
      call VMXT1(part,fvmx,idimp,nop)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wvdist1(part,fv,fvm,nmv)
! calculates 1d velocity distribution, velocity moments, and entropy
      implicit none
      integer, intent(in) :: nmv
      real, dimension(:,:), intent(in) :: part
      real, dimension(:), intent(inout) :: fv
      real, dimension(3), intent(inout) :: fvm
! local data
      integer :: idimp, nop, nmvf
! extract dimensions
      idimp = size(part,1); nop = size(part,2)
      nmvf = size(fv,1)
! call low level procedure
      call VDIST1(part,fv,fvm,idimp,nop,nmv,nmvf)
      end subroutine
!
      end module
