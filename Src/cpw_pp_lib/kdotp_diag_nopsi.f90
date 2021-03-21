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

!>     Calculates the variational k.p (Langreth-Kohn) energies and eigenvectors.
!>     Works for both spin-orbit and no spin-orbit cases

       subroutine kdotp_diag_nopsi(rkpt,neig,ei,                         &
     & rk0,h0,dh0drk,d2h0drk2,                                           &
     & mxdbnd)

!      Written 7 February 2014. JLM
!      modified, vkb dimensions, March 31, 2014. JLM
!      modified, does not calculate psi, May 4, 2014. JLM
!      modified, documentation, February 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.95

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands

       real(REAL64), intent(in)           ::  rkpt(3)                    !<  component in lattice coordinates of the k-point
       integer, intent(in)                ::  neig                       !<  number of eigenvectors (requested on input, modified by degeneracies on output)

       real(REAL64), intent(in)           ::  rk0(3)                     !<  component in lattice coordinates of the reference k-point
       complex(REAL64), intent(in)        ::  h0(mxdbnd,mxdbnd)          !<  <Psi|p|Psi>
       complex(REAL64), intent(in)     ::  dh0drk(mxdbnd,mxdbnd,3)       !<  d <Psi|V_NL|Psi> d k
       complex(REAL64), intent(in)     ::  d2h0drk2(mxdbnd,mxdbnd,3,3)   !<  d^2 <Psi|V_NL|Psi> d k^2

!      output

       real(REAL64), intent(out)          ::  ei(mxdbnd)                 !<  eigenvalues (Hartree)

!      allocatable arrays

       complex(REAL64), allocatable       ::  vec(:,:)                   !  overlap matrix S
       complex(REAL64), allocatable       ::  hred(:,:)                  !  reduced hamiltonian

!      local

       integer            ::  info

!      counters

       integer    ::  i, j, m, n

       allocate(hred(neig,neig))
       allocate(vec(neig,neig))
       
       do i=1,neig
       do j=1,neig
         hred(j,i) = h0(j,i)
         do n=1,3
           hred(j,i) =  hred(j,i) + dh0drk(j,i,n)*(rkpt(n)-rk0(n))
         enddo
         do n=1,3
         do m=1,3
           hred(j,i) =  hred(j,i) + (rkpt(m)-rk0(m))*d2h0drk2(j,i,m,n)   &
     &                                      *(rkpt(n)-rk0(n)) / 2
         enddo
         enddo
       enddo
       enddo


       call diag_c16(neig,hred,ei,vec,neig,info)

  
       if(info /= 0) stop


       deallocate(vec)
       deallocate(hred)

       return

       end subroutine kdotp_diag_nopsi
