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
!>  \version      5.11
!>  \date         8 - 14 November 2023.
!>  \copyright    GNU Public License v2

subroutine out_mass_print(neig, ei, deidxk, d2eidxk2,                    &
     mxdbnd)

! Written November 2024. JLM
! Modified print statement. 5 March 2024. JLM




  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  neig                            !<  wavefunction dimension
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalue E
  real(REAL64), intent(in)           ::  deidxk(mxdbnd)                  !<  d E / d xk  (lattice coordinates)
  real(REAL64), intent(in)           ::  d2eidxk2(mxdbnd)                !<  d^2 E / d xk^2  (lattice coordinates)

! allocatable arrays

  integer, allocatable               ::  levdeg(:)                       !  degeneracy of energy level
  integer, allocatable               ::  leveigs(:,:)                    !  states belonging to level

! local variables

  integer           ::  mxdlev                                           !  array dimension for number of levels
  integer           ::  mxddeg                                           !  array dimension for number of levels
  integer           ::  neigtmp                                          !  temporary value of neig.  neig should not be changed here

  integer           ::  nlevel              !  number of energy levels (not counting degeneracies)
  integer           ::  maxdeg              !  maximum number of degeneracies

! constants

  real(REAL64), parameter     ::  UM = 1.0_REAL64
  real(REAL64), parameter     ::  TOL = 1.0E-4_REAL64
  real(REAL64), parameter     ::  HARTREE = 27.21138386_REAL64

! counters

  integer       ::  n, nl, nk

! energy levels.  first finds dimensions

  allocate(levdeg(1))
  allocate(leveigs(1,1))

  neigtmp = neig

  call berry_degeneracy(.TRUE., neig, neigtmp, ei, TOL,                  &
         nlevel, maxdeg, levdeg, leveigs,                                &
         mxdbnd, 1, 1)

  mxdlev = nlevel
  mxddeg = maxdeg

  deallocate(levdeg)
  deallocate(leveigs)

  allocate(levdeg(mxdlev))
  allocate(leveigs(mxdlev,mxddeg))

! fills the information

  call berry_degeneracy(.FALSE., neig, neigtmp, ei, TOL,                 &
         nlevel, maxdeg, levdeg, leveigs,                                &
         mxdbnd, mxdlev, mxddeg)


  write(6,*)
  write(6,*)
  write(6,*) '  Effective masses and derivatives of eigenvalues'
  write(6,*) '  (atomic units unless otherwise stated)'
  write(6,*)
  write(6,*) '  n     energy(eV)      Effective mass      d^2 E / d k^2     d E / d k'
  write(6,*)

  do nl = 1,nlevel

    do nk = 1,levdeg(nl)
      n = leveigs(nl,nk)
      write(6,'(i5,f12.3,5x,f14.6,5x,f12.3,5x,f14.6)') n, ei(n)*HARTREE, &
             UM/d2eidxk2(n), d2eidxk2(n), deidxk(n)
    enddo
    write(6,*)
  enddo

  deallocate(levdeg)
  deallocate(leveigs)

  return

end subroutine out_mass_print
