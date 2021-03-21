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

!>     Unfolds a charge density or other quantity represented on 
!>     prototype G-vectors on a uniform mesh.

       subroutine mesh_unfold(den,chd,id,n1,n2,n3,lwrap,                 &
     & ng,kgv,phase,conj,inds,                                           &
     & mxdgve,mxdnst,mxdfft)

!      Written September 30, 2015 from v_hartree_xc
!      Modified 12 December 2019.  Documentation.  JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.62


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars
       integer, intent(in)                ::  mxdfft                     !<  array dimension for chd
       
       complex(REAL64), intent(in)        ::  den(mxdnst)                !<  density or other quantity in prototype G-vector
       logical, intent(in)                ::  lwrap                      !<  indicates if it should wrap around wrong results will be obtained with inconsistent choice.
       integer, intent(in)                ::  id, n1, n2, n3             !<  dimensions of mesh
       
       integer, intent(in)                ::  ng                         !<  size of g-space
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  G-vectors in reciprocal lattice coordinates 
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase
       integer, intent(in)                ::  inds(mxdgve)               !<  star to which g-vector n belongs

!      output

       complex(REAL64), intent(out)       ::  chd(mxdfft)                !<  density or other quantity on regular mesh in G-space

!      local varaibles

       integer         ::  k1, k2, k3, kd, iadd

!      counters

       integer         ::  i
 
!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)


!      initialize mesh array

       do i=1,id*n2*n3
         chd(i) = C_ZERO
       enddo

       if(lwrap) then

!        wraps around the charge density (aliasing in FFT)

         do i=1,ng

           k1 = kgv(1,i)
           kd = n1*(k1/n1)
           if (k1 < 0) kd = kd - n1
           k1 = k1 - kd

           k2 = kgv(2,i)
           kd = n2*(k2/n2)
           if (k2 < 0) kd = kd - n2
           k2 = k2 - kd

           k3 = kgv(3,i)
           kd = n3*(k3/n3)
           if (k3 .lt. 0) kd = kd - n3
           k3 = k3 - kd

           iadd = (k3*n2 + k2)*id + k1 + 1
           
           if(conj(i) > ZERO) then
             chd(iadd) = chd(iadd) + den(inds(i))*conjg(phase(i))
           else
             chd(iadd) = chd(iadd) + conjg(den(inds(i)))*phase(i)
           endif

         enddo

       else

!        does not wrap (this should be the normal case)

         do i=1,ng
           k1 = kgv(1,i)
           if (k1 < 0) k1 = n1 + k1
           k2 = kgv(2,i)
           if (k2 < 0) k2 = n2 + k2
           k3 = kgv(3,i)
           if (k3 < 0) k3 = n3 + k3
           
           iadd = (k3*n2 + k2)*id + k1 + 1
           
           if(conj(i) > ZERO) then
             chd(iadd) = den(inds(i))*conjg(phase(i))
           else
             chd(iadd) = conjg(den(inds(i)))*phase(i)
           endif

         enddo
       endif

       return
       end subroutine mesh_unfold
