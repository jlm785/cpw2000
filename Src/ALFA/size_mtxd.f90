!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

!>     calculates the dimension of the hamiltonian matrix

       subroutine size_mtxd(emax,rkpt,adot,ng,kgv,mtxd)
       
!      Written December 20, 2013 from part of the matrix_diag_kb code
!      Modified, documentation, January 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.94

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input
       real(REAL64), intent(in)           ::  emax                       !  largest kinetic energy included in hamiltonian diagonal. (hartree).
       real(REAL64), intent(in)           ::  rkpt(3)                    !  j-th component in lattice coordinates of the k-point
       real(REAL64), intent(in)           ::  adot(3,3)                  !  metric in direct space
       integer, intent(in)                ::  ng                         !  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,ng)                  !  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

!      output
       integer, intent(out)               ::  mtxd                       !  dimension of the hamiltonian

!      allocatable local arrays
       real(REAL64), allocatable          ::  ekin(:) 
       integer, allocatable               ::  irow(:)

!      local variables
       real(REAL64)      ::  qk(3), vcell, bdot(3,3)

!      counters
       integer    ::  i

!      constants
       real(REAL64), parameter ::  EPS = 0.0000000001_REAL64

       call adot_to_bdot(adot,vcell,bdot)
       
       allocate(ekin(ng))
       allocate(irow(ng))

!      calculate the kinetic energies

       do i=1,ng
         qk(1) = rkpt(1) + kgv(1,i)
         qk(2) = rkpt(2) + kgv(2,i)
         qk(3) = rkpt(3) + kgv(3,i)
         ekin(i) = (qk(1)*bdot(1,1) + qk(2)*bdot(2,1) +                  &
     &              qk(3)*bdot(3,1))*qk(1) +                             &
     &             (qk(1)*bdot(1,2) + qk(2)*bdot(2,2) +                  &
     &              qk(3)*bdot(3,2))*qk(2) +                             &
     &             (qk(1)*bdot(1,3) + qk(2)*bdot(2,3) +                  &
     &              qk(3)*bdot(3,3))*qk(3)
         ekin(i) = ekin(i)/2
       enddo


!      sorts if it is not the gamma point

       if(ekin(1) < EPS) then
         do i=1,ng
           irow(i) = i
         enddo
       else
         call sort(ng,ekin,irow)
       endif

!      find mtxd

       mtxd = 1
       do i=1,ng

         if(ekin(irow(i)) > emax) exit

         mtxd = i
       enddo

       if(mtxd == ng) then
         write(6,*)
         write(6,'("   stopped in size_mtxd   something is wrong...",    &
     &      " the size of g-space is",i8,"  emax = ",f10.4)') ng,emax

         stop

       endif
       
       deallocate(ekin)
       deallocate(irow)

       return
       end subroutine size_mtxd
