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

!>     calculates the matrix <Psi|V_NL|Psi> for a separable non-local
!>     pseudopotential V_NL for neig wavevectors.  complex version
!>     spin wave-functions are the same
!>     appropriate for spin-perturbation

       subroutine psi_vnl_sp_psi(mtxd,neig,nanlso,psi,vnl,               &
     & anlp,anlm,xnlkb,                                                  &
     & mxddim,mxdbnd,mxdaso)

!      written June 2020 from ps_vnl_psi. jlm
!      copyright INESC-MN/Jose Luis Martins

!      version 4.95

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdaso                     !<  array dimension of number of projectors
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands

       integer, intent(in)                ::  mtxd                       !<  wavefunction dimension
       integer, intent(in)                ::  neig                       !<  wavefunction dimension

       complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)         !<  wavevectors

       integer, intent(in)                ::  nanlso                     !<  number of projectors
       complex(REAL64), intent(in)        ::  anlp(mxddim,mxdaso)        !<  Kleinman-Bylander projectors s = +1/2
       complex(REAL64), intent(in)        ::  anlm(mxddim,mxdaso)        !<  Kleinman-Bylander projectors s = -1/2
       real(REAL64), intent(in)           ::  xnlkb(mxdaso)              !<  Kleinman-Bylander normalization

!      output

       complex(REAL64), intent(out)       ::  vnl(2*mxdbnd,2*mxdbnd)     !<  <Psi|V_NL|Psi>

!      local variables

       complex(REAL64),allocatable      ::  dhdp(:,:)
       complex(REAL64),allocatable      ::  dhdm(:,:)
       complex(REAL64),allocatable      ::  xdhdp(:,:)
       complex(REAL64),allocatable      ::  xdhdm(:,:)

       complex(REAL64),allocatable      ::  vnlhalf(:,:)

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
       complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

!      counters

       integer   ::   i,n


       allocate(dhdp(nanlso,neig))
       allocate(dhdm(nanlso,neig))
       allocate(xdhdp(nanlso,neig))
       allocate(xdhdm(nanlso,neig))

       allocate(vnlhalf(neig,neig))

!      dhd = < anl | psi >

       call zgemm('c','n',nanlso,neig,mtxd,C_UM,anlp,mxddim,psi,         &
     &                mxddim,C_ZERO,dhdp,nanlso)

       call zgemm('c','n',nanlso,neig,mtxd,C_UM,anlm,mxddim,psi,         &
     &                mxddim,C_ZERO,dhdm,nanlso)

!      xdhd := Diag(xnl) dhd

       do n=1,neig
         do i=1,nanlso
           xdhdp(i,n) = xnlkb(i)*dhdp(i,n)
         enddo
          do i=1,nanlso
           xdhdm(i,n) = xnlkb(i)*dhdm(i,n)
         enddo
       enddo

!      <Psi|V_NL|Psi> = < psi | anl > Diag(xnl)  < anl | psi >

       call zgemm('c','n',neig,neig,nanlso,C_UM,dhdp,nanlso,xdhdp,nanlso,    &
     &                C_ZERO,vnlhalf,mxdbnd)

       do i = 1,neig
       do n = 1,neig
         vnl(n,i) = vnlhalf(n,i)
       enddo
       enddo

       call zgemm('c','n',neig,neig,nanlso,C_UM,dhdm,nanlso,xdhdm,nanlso,    &
     &                C_ZERO,vnlhalf,mxdbnd)

       do i = 1,neig
       do n = 1,neig
         vnl(n+neig,i+neig) = vnlhalf(n,i)
       enddo
       enddo

       call zgemm('c','n',neig,neig,nanlso,C_UM,dhdp,nanlso,xdhdm,nanlso,    &
     &                C_ZERO,vnlhalf,mxdbnd)

       do i = 1,neig
       do n = 1,neig
         vnl(n,i+neig) = vnlhalf(n,i)
       enddo
       enddo

       do i = 1,neig
       do n = 1,neig
         vnl(n+neig,i) = conjg(vnl(i,n+neig))
       enddo
       enddo

       deallocate(dhdp)
       deallocate(xdhdp)
       deallocate(dhdm)
       deallocate(xdhdm)

       deallocate(vnlhalf)

       return

       end subroutine psi_vnl_sp_psi
