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

!>  Asks the user information to determine
!>  the range of initial and final states for
!>  the calculation of the oscillator sterengths
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         14 May 2025.
!>  \copyright    GNU Public License v2

subroutine out_band_oscillator_range(ioreplay, lso,                      &
            neig, ei, ztot,                                              &
            ninitbeg, ninitend, nfinalbeg, nfinalend, lpair, lexcit,     &
            mxdbnd)

! Written July 2, 2014. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations
  logical, intent(in)                ::  lso                             !<  true if with spin-orbit

  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalues (Hartree)
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)

! output

  integer, intent(out)               ::  ninitbeg, ninitend              !<  begin and end of initial state index
  integer, intent(out)               ::  nfinalbeg, nfinalend            !<  begin and end of final state index
  logical, intent(out)               ::  lpair                           !<  prints the oscillator strengths for pairs of bands
  logical, intent(out)               ::  lexcit                          !<  prints the oscillator strengths by excitation energies

! local variables

  integer           ::  nvbm, nmini
  character(len=1)  ::  yesno, yesno2, yesno3



  if(lso) then
    nvbm = nint(ztot)
  else
    nvbm = nint(0.5*ztot)
  endif

  if(nvbm > neig .or. nvbm > mxdbnd) then
    write(6,*)
    write(6,*) '   Stopped in out_band_oscillator_range'
    write(6,*) '   Inconsistent nvbm, neig, mxdbnd',nvbm, neig, mxdbnd
    write(6,*)

    stop

  endif

  write(6,*)
  write(6,*) '   Do you want the oscillator strengths'
  write(6,*) '   between occupied and empty states'
  write(6,*) '   in case the material is an insulator/semiconductor (y/n)'
  write(6,*)

  read(5,*) yesno
  write(ioreplay,*) yesno,'   material is insulator'

  if(yesno == 'y' .or. yesno == 'Y') then

    ninitbeg = 1
    ninitend = nvbm
    nfinalbeg = ninitend + 1
    nfinalend = neig
    lpair = .FALSE.
    lexcit = .TRUE.

  else

    write(6,*)
    write(6,*) '   Do you want the oscillator strengths'
    write(6,*) '   between states near one of the band edges,'
    write(6,*) '   in case you have superlattices minibands (y/n)'
    write(6,*)

    read(5,*) yesno2
    write(ioreplay,*) yesno2,'   material has minibands'

    if(yesno2 == 'y' .or. yesno2 == 'Y') then

      write(6,*)
      write(6,*) '   How many minibands you want to analyze?'
      write(6,*)
      read(5,*) nmini
      write(ioreplay,*) nmini,'   how many minibands'
      write(6,*)
      write(6,*) '   Do you want to analyze the minibands'
      write(6,*) '   near the valence band maximum?'
      write(6,*) '   "no" will default to conduction band minimum. (y/n)'
      write(6,*)
      read(5,*) yesno3
      write(ioreplay,*) yesno3,'   valence or conduction'
      write(6,*)

      if(yesno3 == 'y' .or. yesno2 == 'Y') then
        ninitbeg = nvbm - nmini + 1
        if(ninitbeg < 1) ninitbeg = 1
        ninitend = nvbm
        nfinalbeg = ninitbeg
        nfinalend = ninitend
      else
        ninitbeg = nvbm + 1
        ninitend = nvbm + nmini
        if(ninitend > neig) ninitend = neig
        nfinalbeg = ninitbeg
        nfinalend = ninitend
      endif
      lpair = .TRUE.
      lexcit = .FALSE.

    else

      write(6,*)
      write(6,*) '   Dropping into the general case'
      write(6,*) '   Top of valence band is around', nvbm
      write(6,*)
      write(6,*) '   Enter the band indices of the beginiing and end of the initial states'
      write(6,*)
      read(5,*) ninitbeg, ninitend
      write(ioreplay,*) ninitbeg, ninitend,'   initial states'
      if(ninitbeg < 1) ninitbeg = 1
      if(ninitend > neig) ninitend = neig
      write(6,*) '   Enter the band indices of the beginiing and end of the final states'
      write(6,*)
      read(5,*) nfinalbeg, nfinalend
      write(ioreplay,*) nfinalbeg, nfinalend,'   final states'
      if(nfinalbeg < 1) nfinalbeg = 1
      if(nfinalend > neig) nfinalend = neig
      write(6,*)

      lpair = .TRUE.
      lexcit = .FALSE.

    endif

  endif

  return

end subroutine out_band_oscillator_range
