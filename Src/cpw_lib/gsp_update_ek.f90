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

!>     updates the kinetic energies for fixed G-space structure and fixed k+G

       subroutine gsp_update_ek(ns,mstar,ek,kgv,adot,                    &
     &     mxdgve,mxdnst)

!      Written 6 August 2019. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars

       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
       integer, intent(in)                ::  mstar(mxdnst)              !<  number of g-vectors in the j-th star
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space

!      output

       real(REAL64), intent(out)          ::  ek(mxdnst)                 !<  kinetic energy (hartree) of g-vectors in star j

!      local variables

       real(REAL64)  ::  vcell, bdot(3,3)

!      counters

       integer       ::  n, iprot

       call adot_to_bdot(adot,vcell,bdot)

       iprot = 1
       do n = 1,ns
         ek(n) = kgv(1,iprot)*bdot(1,1)*kgv(1,iprot) +                   &
     &           kgv(1,iprot)*bdot(1,2)*kgv(2,iprot) +                   &
     &           kgv(1,iprot)*bdot(1,3)*kgv(3,iprot) +                   &
     &           kgv(2,iprot)*bdot(2,1)*kgv(1,iprot) +                   &
     &           kgv(2,iprot)*bdot(2,2)*kgv(2,iprot) +                   &
     &           kgv(2,iprot)*bdot(2,3)*kgv(3,iprot) +                   &
     &           kgv(3,iprot)*bdot(3,1)*kgv(1,iprot) +                   &
     &           kgv(3,iprot)*bdot(3,2)*kgv(2,iprot) +                   &
     &           kgv(3,iprot)*bdot(3,3)*kgv(3,iprot)
         ek(n) = ek(n)/2
         iprot = iprot + mstar(n)
       enddo

       return
       end subroutine gsp_update_ek
