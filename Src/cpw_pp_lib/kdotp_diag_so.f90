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
!>     Spin-orbit version

       subroutine kdotp_diag_so(emax,rkpt,neig,                          &
     & rk0,mtxd0,isort0,psi0,hso0,dhso0drk,d2hso0drk2,                   &
     & psi,ei,mtxd,isort,qmod,ekpg,                                      &
     & ng,kgv,adot,                                                      &
     & mxdgve,mxddim,mxdbnd)

!      Written 7 February 2014. JLM
!      modified, vkb dimensions, March 31, 2014. jlm
!      Adapted for spin-orbit, April 9, 2014. jlm
!      modified, documentation, February 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.95

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands

       real(REAL64), intent(in)           ::  emax                       !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
       real(REAL64), intent(in)           ::  rkpt(3)                    !<  component in lattice coordinates of the k-point
       integer, intent(in)                ::  neig                       !<  number of eigenvectors (requested on input, modified by degeneracies on output)

       integer, intent(in)                ::  mtxd0                      !<  dimension of the reference hamiltonian
       integer, intent(in)                ::  isort0(mxddim)             !<  g-vector associated with row/column i of reference hamiltonian
       complex(REAL64), intent(in)        ::  psi0(2*mxddim,2*mxdbnd)        !<  component j of reference eigenvector

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space

       real(REAL64), intent(in)           ::  rk0(3)                     !<  component in lattice coordinates of the reference k-point
       
       complex(REAL64), intent(in)        ::  hso0(2*mxdbnd,2*mxdbnd)    !<  <Psi|H|Psi> with spin-orbit
       complex(REAL64), intent(in)   ::  dhso0drk(2*mxdbnd,2*mxdbnd,3)       !<  d <Psi|H|Psi> d k
       complex(REAL64), intent(in)   ::  d2hso0drk2(2*mxdbnd,2*mxdbnd,3,3)   !<  d^2 <Psi|H|Psi> d k^2

!      output

       complex(REAL64), intent(out)       ::  psi(2*mxddim,2*mxdbnd)     !<  component j of eigenvector i
       real(REAL64), intent(out)          ::  ei(2*mxdbnd)               !<  eigenvalues (Hartree)
       integer, intent(out)               ::  mtxd                       !<  dimension of the hamiltonian
       integer, intent(out)               ::  isort(mxddim)              !<  g-vector associated with row/column i of hamiltonian
       real(REAL64), intent(out)          ::  qmod(mxddim)               !<  length of k+g-vector of row/column i
       real(REAL64), intent(out)          ::  ekpg(mxddim)               !<  kinetic energy (hartree) of k+g-vector of row/column i

!      allocatable arrays

       complex(REAL64), allocatable       ::  psi_lk(:,:)                !  component j of eigenvector i
       complex(REAL64), allocatable       ::  vec(:,:)                   !  overlap matrix S
       complex(REAL64), allocatable       ::  hred(:,:)                  !  reduced hamiltonian

!      other variables

       real(REAL64)      ::  vcell, bdot(3,3)

       integer            ::  info

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       complex(REAL64), parameter  :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
       complex(REAL64), parameter  :: C_UM = cmplx(UM,ZERO,REAL64)

!      counters

       integer    ::  i, j, m, n


       call adot_to_bdot(adot,vcell,bdot)
       
       call hamilt_struct(emax,rkpt,mtxd,isort,qmod,ekpg,.FALSE.,        &
     & ng,kgv,adot,                                                      &
     & mxdgve,mxddim)

         
       allocate(psi_lk(2*mxddim,2*mxdbnd))

       call psi_convert_so(neig,mtxd0,isort0,psi0,mtxd,isort,psi_lk,     &
     & mxddim,mxdbnd)

       allocate(hred(2*neig,2*neig))
       allocate(vec(2*neig,2*neig))
       
       do i=1,2*neig
       do j=1,2*neig
         hred(j,i) = hso0(j,i)
         do n=1,3
           hred(j,i) =  hred(j,i) + dhso0drk(j,i,n)*(rkpt(n)-rk0(n))
         enddo
         do n=1,3
         do m=1,3
           hred(j,i) =  hred(j,i) + (rkpt(m)-rk0(m))*d2hso0drk2(j,i,m,n) &
     &                                      *(rkpt(n)-rk0(n)) / 2
         enddo
         enddo
       enddo
       enddo


       call diag_c16(2*neig,hred,ei,vec,2*neig,info)

  
       if(info /= 0) stop


       call zgemm('n','n',2*mtxd,2*neig,2*neig,C_UM,psi_lk,              &
     &          2*mxddim,vec,2*neig,C_ZERO,psi,2*mxddim)

       deallocate(vec)
       deallocate(hred)

       deallocate(psi_lk)

       return

       end subroutine kdotp_diag_so
