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

!>     Given a scalar field (density) represented in the prototype G-vector, 
!>     unfolds it to the full G-space.

       subroutine star_of_g_unfold(deng,denk,ladd,                       &
     & ng,phase,conj,inds,                                               &
     & mxdgve,mxdnst)

!      Written October 6, 2015. JLM
!      Modified 12 December 2019. JLM
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
       integer, intent(in)                ::  inds(mxdgve)               !<  star to which g-vector n belongs

       complex(REAL64), intent(in)        ::  denk(mxdnst)               !<  symmetrized density in prototype G-vector

!      output
       
       complex(REAL64), intent(out)       ::  deng(mxdgve)               !<  density in G-vectors

!      counters

       integer         ::  i
 
!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)


       if(.not. ladd) then
         do i=1,ng
           deng(i) = C_ZERO
         enddo
       endif


!      CONJ SHOULD BE CONVERTED TO INTEGER OR LOGICAL

       do i = 1,ng
         if(conj(i) > 0) then
           deng(i) = deng(i) + denk(inds(i))*conjg(phase(i))
         else
           deng(i) = deng(i) + conjg(denk(inds(i)))*phase(i)
         endif
       enddo

       return
       end subroutine star_of_g_unfold
