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

!>     Computes the local potential contribution
!>     to the Hellman-Feynman forces. It uses
!>     the selfconsistent charge density and potentials

       subroutine for_str_local_force(floc,                              &
     & ng,kgv,phase,conj,ns,mstar,                                       &
     & vxc,den,                                                          &
     & ntype,natom,rat,                                                  &
     & vql,dnc,                                                          &
     & mxdgve,mxdnst,mxdtyp,mxdatm)

!      Adapted from Sverre Froyen plane wave program
!      Written January 15 1988. jlm
!      Modified March 12 and 22 1999. jlm
!      Modified March 17 2000. jlm (bug found by maosheng)
!      Modified December 2016, f90. JLM
!      Modified, documentation, January 2020. JLM
!      copyright inesc-mn/Jose Luis Martins

!      version 4.94

       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type
       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  real part of the phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase

       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
       integer, intent(in)                ::  mstar(mxdnst)              !<  number of G-vectors in the j-th star

       complex(REAL64), intent(in)        ::  vxc(mxdnst)                !<  exchange+correlation potential for the prototype G-vector
       complex(REAL64), intent(in)        ::  den(mxdnst)                !<  total charge density for the prototype G-vector

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  lattice coordinates of atom j of type i
       
       real(REAL64), intent(in)           ::  vql(mxdtyp,mxdnst)         !<  local pseudopotential for atom type i and prototype g-vector in star j       real*8 floc(3,mxdatm,mxdtyp)
       real(REAL64), intent(in)           ::  dnc(mxdtyp,mxdnst)         !<  core charge for atom type i and prototype g-vector in star j

!      output

       real(REAL64), intent(out)          ::  floc(3,mxdatm,mxdtyp)      !<  unsymmetrized local force (in covariant lattice coordinates)

!      local variables

       real(REAL64)           ::  fi, exp1
       integer                ::  kadd, jj
!       real(REAL64)           ::  vdsr,vdsi,vdgr,vdgi
       complex(REAL64)        ::  vds, vdg

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       complex(REAL64), parameter :: C_I = cmplx(ZERO,UM,REAL64)

!      counters

       integer       ::  i, j, k, l




!      initialize floc

       do i = 1,ntype
         do j = 1,natom(i)
         do k=1,3
           floc(k,j,i) = ZERO
         enddo
         enddo
       enddo

!      force contribution - local parts

!      loop over atomic types

       do i = 1,ntype

!        loop over stars

         kadd = ng + 1
         do jj = 2,ns
           j = ns - jj + 2

           vds =  vql(i,j) * conjg(den(j)) + dnc(i,j) * conjg(vxc(j))

!          loop over g vectors in star

           do k = 1,mstar(j)
             kadd = kadd - 1
             vdg = vds * phase(kadd)
             if(conj(kadd) < ZERO) vdg = conjg(vdg)

!            loop over atoms of same type

             do l = 1,natom(i)
               fi = kgv(1,kadd) * rat(1,l,i) + kgv(2,kadd) * rat(2,l,i)  &
     &                                       + kgv(3,kadd) * rat(3,l,i)
               fi = 2*PI*fi
               exp1 = real(C_I*exp(-C_I*fi)*vdg,REAL64)

!              add to forces

               floc(1,l,i) = floc(1,l,i) + kgv(1,kadd) * exp1
               floc(2,l,i) = floc(2,l,i) + kgv(2,kadd) * exp1
               floc(3,l,i) = floc(3,l,i) + kgv(3,kadd) * exp1
             enddo
           enddo
         enddo
       enddo

!      bdot/adotm1

       do i = 1,ntype
         do j = 1,natom(i)
           do k=1,3
             floc(k,j,i) = floc(k,j,i) * 2*PI
           enddo
         enddo
       enddo

       return
       end subroutine for_str_local_force
