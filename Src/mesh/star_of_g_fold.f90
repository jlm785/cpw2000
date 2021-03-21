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

!>     Given a scalar field (density) represented in the G-space, it constructs
!>     a field that has the full symmetry of the crystal
!>     and represents it on the prototype G-vector.

       subroutine star_of_g_fold(denk,denu,ladd,                         &
     & ng,phase,conj,ns,inds,mstar,                                      &
     & mxdgve,mxdnst)

!      Written September 30, 2015 from charge_by_fft
!      Modified 12 December 2019.  Documentation
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars

       logical, intent(in)                ::  ladd                       !<  If true adds to previous denk, otherwise initializes denk

       integer, intent(in)                ::  ng                         !<  size of g-space
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase
       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
       integer, intent(in)                ::  inds(mxdgve)               !<  star to which g-vector n belongs
       integer, intent(in)                ::  mstar(mxdnst)              !<  number of g-vectors in the j-th star
       
       complex(REAL64), intent(in)        ::  denu(mxdgve)               !<  unsymmetrized density G-vectors

!      output

       complex(REAL64), intent(out)       ::  denk(mxdnst)               !<  symmetrized density in prototype G-vector

!      counters

       integer         ::  i
 
!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)


       if(.not. ladd) then
         do i=1,ns
           denk(i) = C_ZERO
         enddo
       endif


!      CONJ SHOULD BE CONVERTED TO INTEGER OR LOGICAL

       do i = 1,ng
         if(conj(i) > 0) then
           denk(inds(i)) = denk(inds(i)) +                           &
     &         denu(i)*phase(i) / mstar(inds(i))
         else
           denk(inds(i)) = denk(inds(i)) +                           &
     &         conjg(denu(i))*phase(i) / mstar(inds(i))
         endif
       enddo

       return
       end subroutine star_of_g_fold
