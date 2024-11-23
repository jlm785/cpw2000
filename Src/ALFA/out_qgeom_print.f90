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

!>  Prints the topological tensors for a given k-vector
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         27 December 2023, 23 November 2024
!>  \copyright    GNU Public License v2

subroutine out_qgeom_print(ioreplay, nlevel, levdeg, leveigs,            &
      ei, adot, efermi,                                                  &
      tfqg, tgammamf, td2hdk2,                                           &
      mxdbnd, mxdlev, mxddeg)

! written 29 December 2023. JLM
! major reworking, April 2024. JLM
! Removed double double counting in orbital magnetization. 23 November 2024. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdlev                          !<  array dimension for number of levels
  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

  integer, intent(in)                ::  nlevel                          !<  number of energy levels
  integer, intent(in)                ::  levdeg(mxdlev)                  !<  degeneracy of level
  integer, intent(in)                ::  leveigs(mxdlev,mxddeg)          !<  points to degenerate level

  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalue E

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  efermi                          !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

  complex(REAL64), intent(in)        ::  tfqg(3,3,mxddeg,mxddeg,mxdlev)      !<  Tensor F quantum geomtric, curvature and metric, (lattice coordinates)
  complex(REAL64), intent(in)        ::  tgammamf(3,3,mxddeg,mxddeg,mxdlev)  !<  Tensor Gamma, orbital magnetization and contribution to effective mass (lattice coordinates)
  complex(REAL64), intent(in)        ::  td2hdk2(3,3,mxddeg,mxddeg,mxdlev)   !<  <psi| d^2 H / d k^2 |psi> (lattice coordinates)

! local allocatable arrays

  real(REAL64), allocatable          ::  tbcurv(:,:,:,:)
  real(REAL64), allocatable          ::  tmag(:,:,:,:)
  COMPLEX(REAL64), ALLOCATABLE       ::  TMASS(:,:,:,:)
  real(REAL64), allocatable          ::  rtmass(:,:,:,:)
  real(REAL64), allocatable          ::  tqmetric(:,:,:,:)

! local variables

  integer          ::  nl, nq

! constants

  real(REAL64), parameter     ::  HARTREE = 27.21138386_REAL64

! counters

  integer    ::  i, j, n, k, m, jrepeat
  integer    ::  nk, mk


  allocate(tbcurv(mxddeg,mxddeg,3,3))
  allocate(tqmetric(mxddeg,mxddeg,3,3))
  allocate(tmag(mxddeg,mxddeg,3,3))
  allocate(tmass(mxddeg,mxddeg,3,3))
  allocate(rtmass(mxddeg,mxddeg,3,3))

  write(6,*)
  write(6,*) ' E-level     states        energy (eV)   Energy-E_F (eV)'
  write(6,*)
  do n = 1,nlevel

    if(levdeg(n) == 1) then
      write(6,'(i5,5x,i5,7x,f13.3,5x,f13.3)') n, leveigs(n,1),           &
          ei(leveigs(n,1))*HARTREE, (ei(leveigs(n,1))-efermi)*HARTREE
    elseif(levdeg(n) == 2) then
      write(6,'(i5,5x,i5,", ",i3,2x,f13.3,5x,f13.3)') n,                 &
          leveigs(n,1), leveigs(n,2),                                    &
          ei(leveigs(n,1))*HARTREE, (ei(leveigs(n,1))-efermi)*HARTREE
    else
      write(6,'(i5,5x,i5," -",i3,2x,f13.3,5x,f13.3)') n,                 &
          leveigs(n,1), leveigs(n,levdeg(n)),                            &
          ei(leveigs(n,1))*HARTREE, (ei(leveigs(n,1))-efermi)*HARTREE
    endif

  enddo
  write(6,*)

  do jrepeat = 1,1000

    write(6,*)
    write(6,*) '  Which energy level you want to analyze'
    write(6,*) '  enter 0 to finish analysis'
    write(6,*)

    read(5,*) nl
    write(ioreplay,*) nl, '         chosen level'

    if(nl > nlevel .or. nl < 0) nl = 0

    if(nl == 0) exit

    write(6,*)
    write(6,*) '  Which quantity you want to print'
    write(6,*)
    write(6,*) '  Berry curvature             enter 1'
    write(6,*) '  Quantum metric              enter 2'
    write(6,*) '  Orbital magnetization       enter 3'
    write(6,*) '  Effective mass tensor       enter 4'
    write(6,*)

    read(5,*) nq
    write(ioreplay,*) nq, '         chosen quantity'

    if(nq == 1) then

      do nk = 1,levdeg(nl)
      do mk = 1,levdeg(nl)
      do i = 1,3
      do j = 1,3
        tbcurv(nk,mk,i,j) = -2*dimag(tfqg(i,j,nk,mk,nl))
      enddo
      enddo
      enddo
      enddo

      write(6,*)
      write(6,'("   Berry curvature tensor for level with energy ",      &
           &   f13.3)') ei(leveigs(nl,1))*HARTREE
      write(6,*)

      call berry_tensor_print(levdeg(nl), adot, tbcurv, 'A', mxddeg)

    elseif(nq == 2) then

      do nk = 1,levdeg(nl)
      do mk = 1,levdeg(nl)
      do i = 1,3
      do j = 1,3
        tqmetric(nk,mk,i,j) = real(tfqg(i,j,nk,mk,nl),REAL64)
      enddo
      enddo
      enddo
      enddo

      write(6,*)
      write(6,'("Quantum metric tensor for level with energy ",          &
           &   f13.3)') ei(leveigs(nl,1))*HARTREE
      write(6,*)

      call berry_tensor_print(levdeg(nl), adot, tqmetric, 'S', mxddeg)

    elseif(nq == 3) then

      do nk = 1,levdeg(nl)
      do mk = 1,levdeg(nl)
      n = leveigs(nl,nk)
      m = leveigs(nl,mk)
      do i = 1,3
      do j = 1,3
        tmag(nk,mk,i,j) = dimag(tgammamf(i,j,nk,mk,nl))
      enddo
      enddo
      enddo
      enddo

      write(6,*)
      write(6,'("Orbital magnetization for level with energy ",          &
           &   f13.3)') ei(leveigs(nl,1))*HARTREE
      write(6,*)

      call berry_tensor_print(levdeg(nl), adot, tmag, 'A', mxddeg)

    elseif(nq == 4) then

      do nk = 1,levdeg(nl)
      do mk = 1,levdeg(nl)
      do i = 1,3
      do j = 1,3
        tmass(nk,mk,i,j) =                                                   &
                - tgammamf(i,j,nk,mk,nl) - conjg(tgammamf(i,j,mk,nk,nl))     &
                + (td2hdk2(i,j,nk,mk,nl) + conjg(td2hdk2(i,j,mk,nk,nl))) / 2
      enddo
      enddo
      enddo
      enddo

      rtmass(:,:,:,:) = real(tmass(:,:,:,:), REAL64)

      write(6,*)
      write(6,'("Real part mass tensor for level with energy ",          &
           &   f13.3)') ei(leveigs(nl,1))*HARTREE
      write(6,*)

      call berry_tensor_print(levdeg(nl), adot, rtmass, 'S', mxddeg)

      write(6,*)
      write(6,'("Imaginary part of that mass tensor")')
      write(6,*)
      do n = 1,levdeg(nl)
        do j = 1,3
          write(6,'(8(3f10.3,4x))') ((aimag(tmass(n,m,k,j)), k = 1,3), m = 1,levdeg(nl))
        enddo
        write(6,*)
      enddo
      write(6,*)

    endif

  enddo

  deallocate(tbcurv)
  deallocate(tqmetric)
  deallocate(tmag)
  deallocate(tmass)
  deallocate(rtmass)

  return

end subroutine out_qgeom_print
