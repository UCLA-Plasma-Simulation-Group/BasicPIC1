!-----------------------------------------------------------------------
! Interface file for diag1.f
      module diag1_h
      implicit none
!
      interface
         subroutine VMXT1(part,fvmx,idimp,nop)
         implicit none
         integer, intent(in) :: idimp, nop
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(2), intent(inout) :: fvmx
         end subroutine
      end interface
!
      interface
         subroutine VDIST1(part,fv,fvm,idimp,nop,nmv,nmvf)
         implicit none
         integer, intent(in) :: idimp, nop, nmv, nmvf
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(nmvf), intent(inout) :: fv
         real, dimension(3), intent(inout) :: fvm
         end subroutine
      end interface
!
      end module
