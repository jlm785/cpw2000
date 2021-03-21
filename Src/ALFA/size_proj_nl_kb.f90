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

!>     calculates the number of Kleinman-Bylander projectors

       subroutine size_proj_nl_kb(ntype,natom,nkb,nanl,nanlso,mxdtyp)

!      written December 20, 2013. jlm
!      adapted from proj_nl_kb.
!      modified February 10, 2014. jlm
!      Modified, documentation, January 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.94

       implicit none

!       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i

       integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)       !<  KB pseudo.  normalization for atom k, ang. mom. l

!      output

       integer, intent(out)               ::  nanl                       !<  number of Kleinman-Bylander projectorswithout spin-orbit.
       integer, intent(out)               ::  nanlso                     !<  number of Kleinman-Bylander projectors with spin-orbit.
       
!      counters

       integer    ::  k

       nanl = 0
       do k=1,ntype
         if(nkb(0,0,k) /= 0) nanl = nanl + natom(k)
         if(nkb(1,0,k) /= 0) nanl = nanl + 3*natom(k)
         if(nkb(2,0,k) /= 0) nanl = nanl + 5*natom(k)
         if(nkb(3,0,k) /= 0) nanl = nanl + 7*natom(k)
       enddo

       nanlso = 0
       do k=1,ntype
         if(nkb(1,-1,k) /= 0) nanlso = nanlso + 2*natom(k)
         if(nkb(2,-1,k) /= 0) nanlso = nanlso + 4*natom(k)
         if(nkb(3,-1,k) /= 0) nanlso = nanlso + 6*natom(k)
         if(nkb(0, 1,k) /= 0) nanlso = nanlso + 2*natom(k)
         if(nkb(1, 1,k) /= 0) nanlso = nanlso + 4*natom(k)
         if(nkb(2, 1,k) /= 0) nanlso = nanlso + 6*natom(k)
         if(nkb(3, 1,k) /= 0) nanlso = nanlso + 8*natom(k)
       enddo

       return

       end subroutine size_proj_nl_kb
