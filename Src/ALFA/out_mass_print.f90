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

!>  Prints the results of the calculation of the effective mass
!>
!>  \author       Jose Luis Martins
!>  \version      5.08
!>  \date         9 November 2023.
!>  \copyright    GNU Public License v2

subroutine out_mass_print(nlevel, levdeg, leveigs,                       &
     ei, deidxk, d2eidxk2,                                               &
     mxdbnd, mxdlev, mxddeg)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdlev                          !<  array dimension for number of levels
  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels

  integer, intent(in)                ::  nlevel                          !<  number of energy levels
  integer, intent(in)                ::  levdeg(mxdlev)                  !<  degeneragy of level
  integer, intent(in)                ::  leveigs(mxdlev,mxddeg)          !<  points to degenerate level

  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalue E
  real(REAL64), intent(in)           ::  deidxk(mxdbnd)                  !<  d E / d xk  (lattice coordinates)
  real(REAL64), intent(in)           ::  d2eidxk2(mxdbnd)                !<  d^2 E / d xk^2  (lattice coordinates)

! constants

  real(REAL64), parameter     ::  UM = 1.0_REAL64
  real(REAL64), parameter     ::  HARTREE = 27.21138386_REAL64

! counters

  integer       ::  n, nl, nk



  write(6,*)
  write(6,*)
  write(6,*) '  first derivatives of energies'
  write(6,*)

  do nl = 1,nlevel

    if(levdeg(nl) == 1) then
      n = leveigs(nl,1)
      write(6,'(i5,f14.6)') n, deidxk(n)
    else

      do nk = 1,levdeg(nl)
        n = leveigs(nl,nk)
        write(6,'(i5,f14.6)') n, deidxk(n)
      enddo

    endif
  enddo
  write(6,*)


  write(6,*)
  write(6,*)
  write(6,*) '  second derivatives of energies'
  write(6,*)

  do nl = 1,nlevel

    if(levdeg(nl) == 1) then
      n = leveigs(nl,1)
      write(6,'(i5,f14.6)') n, d2eidxk2(n)
    else

      do nk = 1,levdeg(nl)
        n = leveigs(nl,nk)
        write(6,'(i5,f14.6)') n, d2eidxk2(n)
      enddo

    endif
  enddo
  write(6,*)


  write(6,*)
  write(6,*)
  write(6,*) '  Effective masses'
  write(6,*)
  write(6,*) '  n       energy(eV)        Effective mass'
  write(6,*)

  do nl = 1,nlevel

    if(levdeg(nl) == 1) then
      n = leveigs(nl,1)
      write(6,'(i5,f12.3,5x,f14.6)') n, ei(n)*HARTREE, UM/d2eidxk2(n)
    else

      do nk = 1,levdeg(nl)
        n = leveigs(nl,nk)
        write(6,'(i5,f12.3,5x,f14.6)') n, ei(n)*HARTREE, UM/d2eidxk2(n)
      enddo

    endif
  enddo

  return

end subroutine out_mass_print
