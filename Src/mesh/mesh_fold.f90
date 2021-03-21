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

!>     Folds a charge density or other quantity represented on 
!>     a uniform mesh into the corresponding G-vectors.

       subroutine mesh_fold(den,chd,id,n1,n2,n3,                         &
     & ng,kgv,                                                           &
     & mxdgve,mxdfft)

!      Folds a charge density or other quantity represented on 
!      a uniform mesh into the corresponding G-vectors.

!      Written September 30, 2015 from v_hartree_xc
!      Modified 12 December 2019, documentation.  JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdfft                     !<  array dimension for chd
       
       integer, intent(in)                ::  id, n1, n2, n3             !<  dimensions of mesh
       complex(REAL64), intent(in)        ::  chd(mxdfft)                !<  density or other quantity on regular mesh in G-space
       
       integer, intent(in)                ::  ng                         !<  size of g-space
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  G-vectors in reciprocal lattice coordinates 

!      output

       complex(REAL64), intent(out)       ::  den(mxdgve)                !<  density or other quantity in prototype G-vector

!      local varaibles

       integer         ::  k1, k2, k3, iadd
       integer         ::  nn1, nn2, nn3

!      counters

       integer         ::  i
 
!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

       nn1 = (n1-1) / 2
       nn2 = (n2-1) / 2
       nn3 = (n3-1) / 2

!      initialize array

       do i = 1,ng
         den(i) = C_ZERO
       enddo

       do i=1,ng

         k1 = kgv(1,i)
         if (iabs(k1) <= nn1) then
           if (k1 < 0) k1 = n1 + k1
           k2 = kgv(2,i)
           if (iabs(k2) <= nn2) then
             if (k2 < 0) k2 = n2 + k2
             k3 = kgv(3,i)
             if (iabs(k3) <= nn3) then
               if (k3 < 0) k3 = n3 + k3
               iadd = (k3*n2 + k2)*id + k1 + 1
               den(i) = chd(iadd)
             endif
           endif
         endif
       enddo

       return
       end subroutine mesh_fold
