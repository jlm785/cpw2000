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

!>     Given a basis |bas> and the product of the hamiltonian on the basis,
!>     Calculates the energies and eigenvectors |psibas> of the wave-function
!>     in the basis |bas>
!>     If requested calculates also |psi> and H |psi>

       subroutine hk_from_bas_diag(mtxd, nbasorb, neig, ndeg,            &
     & lortho, lbas, lpsi, lhpsi,                                        &
     & bas, hbas, ei, psibas, psi, hpsi,                                 &
     & mxddim, mxdorb, mxdbnd)

!      Adapted from h_kb_dia_ao.f90, 2 June 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.98

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands
       integer, intent(in)                ::  mxdorb                     !<  array dimension for number of orbitals

       integer, intent(in)                ::  mtxd                       !<  dimension of the PW hamiltonian
       integer, intent(in)                ::  nbasorb                    !<  dimension of the atomic orbital hamiltonian

       integer, intent(in)                ::  neig                       !<  number of eigenvectors

       logical, intent(in)                ::  lortho                     !<  indicates if bas is othogonal
       logical, intent(in)                ::  lbas                       !<  indicates if |psi> in the |Bas> basis should be calculated
       logical, intent(in)                ::  lpsi                       !<  indicates if |psi> in PW basis should be calculated
       logical, intent(in)                ::  lhpsi                      !<  indicates if H |psi> in PW basis should be calculated

       complex(REAL64), intent(in)        ::  bas(mxddim,mxdorb)         !<  | bas >  atomic-orbital wave-functions
       complex(REAL64), intent(in)        ::  hbas(mxddim,mxdorb)        !<  H | bas >

!      output

       integer, intent(out)               ::  ndeg                       !<  number of eigenvectors including degeneracy
       real(REAL64), intent(out)          ::  ei(mxdbnd)                 !<  eigenvalues (Hartree)

       complex(REAL64), intent(out)       ::  psibas(mxdorb,mxdbnd)      !<  | psi > expanded in |bas>
       complex(REAL64), intent(out)       ::  psi(mxddim,mxdbnd)         !<  | psi > component j of eigenvector i
       complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)        !<  H | psi >

!      local allocatable arrays

       real(REAL64), allocatable          ::  eg(:)                     !  eigenvalues
       complex(REAL64), allocatable       ::  vec(:,:)                  !  wavefunctions in atomic basis

       complex(REAL64), allocatable       ::  hamsm(:,:)                !  atomic orbital (small) hamiltonian
       complex(REAL64), allocatable       ::  ovsm(:,:)                 !  overlap matrix

!      local variables

       real(REAL64)  ::  epsa, diff

       integer       ::  info

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
       complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

       real(REAL64), parameter :: EPS = 0.0001_REAL64

!      counters

       integer       ::  i

       allocate(hamsm(mxdorb,mxdorb))

       call zgemm('c','n',nbasorb,nbasorb,mtxd,C_UM,hbas,mxddim,         &
     &          bas,mxddim,C_ZERO,hamsm,mxdorb)

       allocate(vec(mxdorb,mxdorb))
       allocate(eg(mxdorb))

       if(lortho) then
         call diag_c16(nbasorb,hamsm,eg,vec,mxdorb,info)

  
         if(info /= 0) stop

       else
         allocate(ovsm(mxdorb,mxdorb))
         call zgemm('c','n',nbasorb,nbasorb,mtxd,C_UM,bas,mxddim,         &
     &          bas,mxddim,C_ZERO,ovsm,mxdorb)

         call diag_c16_gen(nbasorb,hamsm,ovsm,eg,vec,mxdorb,info)
 
  
         if(info /= 0) stop

        deallocate(ovsm)
       endif

       deallocate(hamsm)

!      recalculates ndeg taking into account degeneracies

       epsa = max(EPS,abs((eg(neig)-eg(1))/(50*neig*UM)))
       ndeg = neig
       do i = neig+1,nbasorb
         diff = abs(eg(i) - eg(i-1))

         if(diff > epsa) exit

         ndeg = i
       enddo

       if(ndeg > mxdbnd) ndeg = mxdbnd

       do i = 1,ndeg
         ei(i) = eg(i)
       enddo

!      wave-function in |bas>

       if(lbas) then
         do i = 1,ndeg
           call zcopy(nbasorb,vec(:,i),1,psibas(:,i),1)
         enddo
       endif

!      |psi> wave-function in plane-waves

       if(lpsi) then
         call zgemm('n','n',mtxd,ndeg,nbasorb,C_UM,bas,mxddim,           &
     &            vec,mxdorb,C_ZERO,psi,mxddim)
       endif

!      H |psi> in plane-waves

       if(lhpsi) then
         call zgemm('n','n',mtxd,ndeg,nbasorb,C_UM,hbas,mxddim,          &
     &            vec,mxdorb,C_ZERO,hpsi,mxddim)
       endif

       deallocate(vec)
       deallocate(eg)

       return
       end subroutine hk_from_bas_diag
