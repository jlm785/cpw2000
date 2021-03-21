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

!>     Calculates the kinetic energy of the wave-functions

       subroutine kinetic_energy(neig,mtxd,ekpg,psi,ekpsi,               &
     & mxddim,mxdbnd)

!      written february 18 1990. jlm
!      version 4.0. 15 october 93. jlm
!      modified (f90) 14 January 2014. jlm
!      icmplx removed, 18 October 2015. JLM
!      copyright INESC-MN/Jose Luis Martins


!      version 4.94

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands
       integer, intent(in)                ::  neig                       !<  number of eigenvectors (requested on input, modified by degeneracies on output)
       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian
       real(REAL64), intent(in)           ::  ekpg(mxddim)               !<  kinetic energy (hartree) of k+g-vector of row/column i

       complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)         !<  component j of eigenvector i

!      output

       real(REAL64), intent(out)          ::  ekpsi(mxdbnd)              !<  kinetic energy of eigenvector i. (Hartree)

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64

!      counters

       integer    ::  i, j


       do i=1,neig
         ekpsi(i) = ZERO
         do j=1,mtxd
           ekpsi(i) = ekpsi(i) + ekpg(j)*                              &
     &        real(psi(j,i)*conjg(psi(j,i)),REAL64)
         enddo
       enddo

       return
       end
