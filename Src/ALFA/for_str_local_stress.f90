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

!>     Computes the local-pseudopotential and hartree contribution to the stress,
!>     plus the correction due to the partial core charge in xc.

       subroutine for_str_local_stress(ealpha,strhl,                     &
     & ng,kgv,ns,mstar,ek,                                               &
     & vion,vhar,vxc,den,                                                &
     & adot,                                                             & 
     & dvql,ddc,                                                         &
     & mxdgve,mxdnst)

!      Adapted from Sverre Froyen plane wave program
!      Written January 19 1988. jlm
!      Modified March 9, 1999. jlm
!      Modified June 23,1999. jlm
!      Modified f90, December 2016. JLM
!      Modified, documentation, January 2020. JLM
!      copyright inesc-mn/Jose Luis Martins

!      version 4.94

       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars

       real(REAL64), intent(in)           ::  ealpha                     !<  alpha term. 
       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
       integer, intent(in)                ::  mstar(mxdnst)              !<  number of G-vectors in the j-th star
       real(REAL64), intent(in)           ::  ek(mxdnst)                 !<  kinetic energy (hartree) of g-vectors in star j

       complex(REAL64), intent(in)        ::  vion(mxdnst)               !<  ionic potential for the prototype G-vector in star j
       complex(REAL64), intent(in)        ::  den(mxdnst)                !<  total charge density for the prototype G-vector
       complex(REAL64), intent(in)        ::  vhar(mxdnst)               !<  Hartree potential for the prototype G-vector
       complex(REAL64), intent(in)        ::  vxc(mxdnst)                !<  exchange+correlation potential for the prototype G-vector

       complex(REAL64), intent(in)        ::  dvql(mxdnst)               !<  derivative of the local pseudopotential for the prototype g-vector in star j
       complex(REAL64), intent(in)        ::  ddc(mxdnst)                !<  derivative of the core charge for the prototype g-vector in star j

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in real space

!      output:

       real(REAL64), intent(out)          ::  strhl(3,3)                 !<  unsymmetrized stress tensor from hartree and local potential in covariant lattice coordinates (hartree,bohr)

!      local variables

       real(REAL64)           ::  bdot(3,3), vcell
       real(REAL64)           ::  ssum, vld, dvld, vhd, dvxd
       real(REAL64)           ::  exp1, exp2
       integer                ::  kadd

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64
       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64

!      counters

       integer       ::  i, j, k, l, ii 


       call adot_to_bdot(adot,vcell,bdot)

!      initialize ssum array

       do i = 1,3
       do j = 1,3
         strhl(i,j) = ZERO
       enddo
       enddo

!      stress contribution - local parts

!      loop over stars

       kadd = ng + 1
       ssum = ZERO

       do ii = 2,ns
         i = ns - ii + 2

!        vlocal * den

         vld = real(vion(i)*conjg(den(i)),REAL64)

!        derivative of v local * den

         dvld = real(dvql(i)*conjg(den(i)),REAL64)

!        vhartree * den / 2

         vhd = real(vhar(i)*conjg(den(i)),REAL64) / 2

!        correction for partial core charge in xc
!        vxc * d(dcor)/d(g)

         dvxd = real(vxc(i)*conjg(ddc(i)),REAL64)

!        diagonal (exp1) and nondiag (exp2) contributions

         exp1 = (vld + vhd) * mstar(i)
         exp2 = 2 * (dvld - vhd/(2*ek(i)) + dvxd)

!        add diagonal terms 

         ssum = ssum + exp1

!        loop over g vectors in star and add nondiag terms

         do j = 1,mstar(i)
           kadd = kadd - 1
           do k = 1,3
           do l = 1,3
             strhl(l,k) = strhl(l,k)                                     &
     &              + exp2 * (kgv(l,kadd)*kgv(k,kadd))
           enddo
           enddo
         enddo
       enddo

!      ealpha contribution + ssum

       do j = 1,3
       do i = 1,3
         strhl(i,j) = strhl(i,j)*4*PI*PI                                 &
     &                   + (ealpha + ssum) * adot(i,j)
       enddo
       enddo

       return
       end subroutine for_str_local_stress
