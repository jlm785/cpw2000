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

!>     Computes and adds contributions to the
!>     stress from the kinetic energy.

       subroutine for_str_kinetic_stress(strkin,                         &
     & mtxd,rkpt,neig,occp,                                              &
     & isort,psi,                                                        &
     & kgv,                                                              &
     & mxdgve,mxddim,mxdbnd)

!      Adapted from Sverre Froyen plane wave program
!      Written January 22 1988. jlm
!      Modified March 7 and 22 1999. jlm
!      Modified, documentation, January 2020. JLM
!      copyright inesc-mn/Jose Luis Martins

!      version 4.94

       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands

       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian
       real(REAL64), intent(in)           ::  rkpt(3)                    !<  component in lattice coordinates of the k-point
       integer, intent(in)                ::  neig                       !<  number of eigenvectors (requested on input, modified by degeneracies on output)
       real(REAL64), intent(in)           ::  occp(mxdbnd)               !<  fractional ocupation of level j

       integer, intent(in)                ::  isort(mxddim)              !<  g-vector associated with row/column i of hamiltonian
       complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)         !<  component j of vector i

       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

!      output

       real(REAL64), intent(out)          ::  strkin(3,3)                !<  unsymmetrized electron-kinetic energy contribution to the stress tensor in covariant lattice coordinates (hartree,bohr) 

!      local variables

       real(REAL64)           ::  qki(3), zsumr
       integer                ::  jmin, jmax, iadd

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64

!      counters

       integer       ::  i, j, k


       if (neig <= 0) return

       do j = 1,3
       do k = 1,3
         strkin(k,j) = ZERO
       enddo
       enddo

!      find min and max band with nonzero occupation

       jmin = 0
       jmax = 0
       do i = 1,neig
         if (occp(i) /= zero .and. jmin == 0) jmin = i
         if (occp(i) /= zero .and. jmin /= 0) jmax = i
       enddo

       if (jmin == 0) return

!      start loop over g(i)

       do i = 1,mtxd
         iadd = isort(i)

!        compute |g+k|

         qki(1) = rkpt(1) + kgv(1,iadd)*UM
         qki(2) = rkpt(2) + kgv(2,iadd)*UM
         qki(3) = rkpt(3) + kgv(3,iadd)*UM
      
!        do eigenvector multiplication and band sum

         zsumr = ZERO
         do k = jmin,jmax
           zsumr = zsumr + occp(k)*                                      &
     &                       real(psi(i,k)*conjg(psi(i,k)),REAL64)
         enddo
         zsumr = 2 * zsumr

!        stress contribution from kinetic energy

         do j = 1,3
         do k = 1,3
           strkin(k,j) = strkin(k,j) + zsumr * qki(k) * qki(j)
         enddo
         enddo

       enddo

       do j=1,3
       do k=1,3
         strkin(k,j) = strkin(k,j) * 2*PI*PI
       enddo
       enddo

       return
       end subroutine for_str_kinetic_stress
