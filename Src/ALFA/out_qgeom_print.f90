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
!>  \version      5.10
!>  \date         27 December 2023.
!>  \copyright    GNU Public License v2

subroutine out_qgeom_print(ioreplay, nlevel, levdeg, leveigs,            &
      ei, adot, efermi,                                                  &
      tbcurv, tmag, tmass, tqmetric,                                     &
      mxdbnd, mxdlev, mxddeg)

! written 29 December 2023. JLM

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

  real(REAL64), intent(in)           ::  tbcurv(mxddeg,mxddeg,3,3,mxdlev)    !<  Tensor of Berry curvature  (lattice coordinates)
  real(REAL64), intent(in)           ::  tmag(mxddeg,mxddeg,3,3,mxdlev)  !<  Tensor associated with orbital magnetization  (lattice coordinates)

  complex(REAL64), intent(in)        ::  tmass(mxddeg,mxddeg,3,3,mxdlev) !<  Tensor associated with effective mass, d^2 E_i /d k_1 d k_2 (lattice coordinates)
  real(REAL64), intent(in)           ::  tqmetric(mxddeg,mxddeg,3,3,mxdlev)  !<  Tensor of the quantum metic (lattice coordinates)

! local variables

  integer          ::  nl, nq
  real(REAL64)     ::  pcar(3)

! constants

  real(REAL64), parameter     ::  HARTREE = 27.21138386_REAL64

! counters

  integer    ::  j, n, k, m, jrepeat


  write(6,*)
  write(6,*) ' E-level     states        energy (eV)   Energy-E_F (eV)'
  write(6,*)
  do n = 1,nlevel

    if(levdeg(n) == 1) then
      write(6,'(i5,5x,i5,6x,f13.3,5x,f13.3)') n, leveigs(n,1),           &
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

      write(6,*)
      write(6,'("Berry curvature vector for level with energy ",         &
           &   f13.3)') ei(leveigs(nl,1))*HARTREE
      write(6,*)
      do n = 1,levdeg(nl)
        do j = 1,3
          write(6,'(8(3f10.3,4x))') ((tbcurv(n,m,k,j,nl), k = 1,3), m = 1,levdeg(nl))
        enddo
        write(6,*)
      enddo
!       do n = 1,levdeg(nl)
!         call berry_pseudovec_lat2car(adot, bcurv(:,leveigs(nl,n)), pcar)
!         write(6,'(3f14.5)') (pcar(k), k = 1,3)
!       enddo
      write(6,*)

    elseif(nq == 2) then

      write(6,*)
      write(6,'("Quantum metric tensor for level with energy ",          &
           &   f13.3)') ei(leveigs(nl,1))*HARTREE
      write(6,*)
      do n = 1,levdeg(nl)
        do j = 1,3
          write(6,'(8(3f10.3,4x))') ((tqmetric(n,m,k,j,nl), k = 1,3), m = 1,levdeg(nl))
        enddo
        write(6,*)
      enddo
!       do n = 1,levdeg(nl)
!         do m = 1,3
!           write(6,'(3f14.5)') (qmetric(m,k,leveigs(nl,n)), k = 1,3)
!         enddo
!         write(6,*)
!       enddo
      write(6,*)

    elseif(nq == 3) then

      write(6,*)
      write(6,'("Orbital magnetization for level with energy ",          &
           &   f13.3)') ei(leveigs(nl,1))*HARTREE
      write(6,*)
      do n = 1,levdeg(nl)
        do j = 1,3
          write(6,'(8(3f10.3,4x))') ((tmag(n,m,k,j,nl), k = 1,3), m = 1,levdeg(nl))
        enddo
        write(6,*)
      enddo
      write(6,*)

    elseif(nq == 4) then

      write(6,*)
      write(6,'("Real part mass tensor for level with energy ",          &
           &   f13.3)') ei(leveigs(nl,1))*HARTREE
      write(6,*)
      do n = 1,levdeg(nl)
        do j = 1,3
          write(6,'(8(3f10.3,4x))') ((real(tmass(n,m,k,j,nl)), k = 1,3), m = 1,levdeg(nl))
        enddo
        write(6,*)
      enddo

      write(6,*)
      write(6,'("Imaginary part of that mass tensor")')
      write(6,*)
      do n = 1,levdeg(nl)
        do j = 1,3
          write(6,'(8(3f10.3,4x))') ((aimag(tmass(n,m,k,j,nl)), k = 1,3), m = 1,levdeg(nl))
        enddo
        write(6,*)
      enddo
      write(6,*)

    endif

  enddo

  return

end subroutine out_qgeom_print
