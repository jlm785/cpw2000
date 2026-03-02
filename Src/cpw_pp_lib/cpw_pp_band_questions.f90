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

!>  Sets up the defults for PW energy cutoff,
!>  eigenvalue precision, Fermi energy estimate
!>  and dual approximation.
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         2 March 2026.
!>  \copyright    GNU Public License v2

subroutine cpw_pp_band_questions(ioreplay,                               &
     emax_in, epspsi_in, efermi_in,                                      &
     emax, flgdal, epspsi, efermi)

! written 2 March 2026 from cpw_pp_band_dos_opt cpw_pp_prepare. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

  real(REAL64), intent(in)           ::  emax_in                         !<  default value of plane wave cutoff
  real(REAL64), intent(in)           ::  epspsi_in                       !<  requested precision of the eigenvectors
  real(REAL64), intent(in)           ::  efermi_in                       !<  default value of Fermi energy or highest occupied state (T=0)

! output

  real(REAL64), intent(out)          ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4), intent(out)      ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  real(REAL64), intent(out)          ::  epspsi                          !<  requested precision of the eigenvectors
  real(REAL64), intent(out)          ::  efermi                          !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

! local variables


  logical            ::  ldefault                                        !  use default values of parameters

  real(REAL64)       ::  xprec

  character(len=1)   ::  yesno, yesno2

! constants

  real(REAL64), parameter  :: UM = 1.0_REAL64
  real(REAL64), parameter  :: HARTREE = 27.21138386_REAL64


  write(6,*)
  write(6,*) '  Do you want the default values for PW energy cutoff, '
  write(6,*) '  eigenvalue precision, fermi energy estimate and use'
  write(6,*) '  the dual approximation (y/n)?'
  read(5,*) yesno
  write(ioreplay,*) yesno,'   default parameters'

  ldefault = .FALSE.
  if(yesno == 'y' .or. yesno == 'Y') then
    ldefault = .TRUE.
    write(6,*)
    write(6,*) '  Will use default values.'
  endif

  if(ldefault) then

!   use the defaults

    emax = emax_in
    write(6,*)
    write(6,'("  Using a PW cutoff of ",f10.3," Hartree")') emax
    write(6,*)

    flgdal = 'DUAL'
    write(6,*)
    write(6,*) '  Will use DUAL approximation'
    write(6,*)

    epspsi = epspsi_in
    write(6,*)
    write(6,'("  The eigenvalue precision wil be: ",g14.5)') epspsi
    write(6,*)

    efermi = efermi_in
    write(6,*)
    write(6,'("   Will use the crude estimate for Fermi energy: ",g14.5, "eV")')   &
          efermi*HARTREE
    write(6,*)

  else

!   ask for new values

    write(6,*)
    write(6,'("  The original calculation used a maximum energy",          &
        & " PW cutoff of",f10.3," Hartree")') emax_in
    write(6,*) '  Enter desired maximum energy in Hartree '
    read(5,*) emax
    write(ioreplay,*) emax,'   emax'

    write(6,*)
    write(6,*) '  Do you want to use the dual approximation (y/n)?'
    read(5,*) yesno
    write(ioreplay,*) yesno,'   dual'


    flgdal = '    '
    if(yesno == 'y' .or. yesno == 'Y') then
      flgdal = 'DUAL'
      write(6,*)
      write(6,*) '  Will use DUAL approximation'
      write(6,*)
    else
      write(6,*)
      write(6,*) '  Will not use dual approximation'
      write(6,*)
    endif

    epspsi = epspsi_in
    write(6,*)
    write(6,*) '  Do you want to modify eigenvalue precision (y/n)?'
    write(6,'("  Current value of epspsi is: ",g14.5)') epspsi

    read(5,*) yesno2
    write(ioreplay,*) yesno2,'   modify eigenvalue precision'

    if(yesno2 == 'y' .or. yesno2 == 'Y') then
      write(6,*)
      write(6,*) '  Enter number of decimals for eigenvalue precision (eV)'
      read(5,*) xprec
      write(ioreplay,*) xprec,'   decimals in precision (eV)'
      if(xprec < 3.0) then
        xprec = 3*UM
        write(6,*)
        write(6,*) '  low number, will use three decimals'
        write(6,*)
      endif
      if(xprec > 8.0) then
        xprec = 8*UM
        write(6,*)
        write(6,*) '  high number, will use eight decimals'
        write(6,*)
      endif
      epspsi = UM / ((10*UM)**xprec)
      epspsi = epspsi / HARTREE
    endif


    efermi = efermi_in
    write(6,*)
    write(6,*) '  Do you want to enter a precise Fermi energy'
    write(6,*) '  or top of the valence band value (y/n)?'
    write(6,*)
    write(6,'("   Current estimate is: ",g14.5)') efermi*HARTREE
    write(6,*) '  It is an estimate from the SCF, may be improved'
    write(6,*) '  if you have a value from a band structure or DOS.'

    read(5,*) yesno2
    write(ioreplay,*) yesno2,'   modify Fermi energy'

    if(yesno2 == 'y' .or. yesno2 == 'Y') then
      write(6,*)
      write(6,*) '  Enter your better value (eV)'
      read(5,*) efermi
      write(ioreplay,*) efermi,'   new Fermi energy (eV)'
      if(abs(efermi_in*HARTREE - efermi) > 2.0) then
        write(6,*)
        write(6,*) '  WARNING  WARNING  WARNING:  value looks suspicious'
        write(6,*) '  eigenvalues may not be what you expect'
        write(6,*) '  but all other quantities are unaffected'
        write(6,*)
      endif
      efermi = efermi / HARTREE
    endif

  endif


  return

end subroutine cpw_pp_band_questions

