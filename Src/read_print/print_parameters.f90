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

!>     prints the information about the parameters used in the calculation
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         20 september 2002. 25 November 2025.
!>  \copyright    GNU Public License v2

subroutine print_parameters(flgcal, flgdal, flgscf,                      &
    emax, author, tblaha, epscv, epscvao, epspsi,                        &
    tempk, teleck, tempinik, nstep, tstep, beta, iseed,                  &
    press, strext, celmas, flgkplusg, epskplusg)

! written 20 september 2002. jlm
! modified 6 february 2008. jlm
! modified 18 february 2008. jlm
! Modified f90, tblaha, 24 Outubro 2015. JLM
! Modified, documentation, kplusg,August 10 2019. JLM
! Modified, indentation, types of correlation, len=* in author. 12 January 2024. JLM
! Modified, write statement continuation, 22 February 2024. JLM
! Modified, prints information about more functionals. 25 November 2025. JLM


  implicit none

  integer, parameter          :: real64 = selected_real_kind(12)

! intput

  character(len=6), intent(in)       ::  flgcal                          !<  type of calculation
  character(len=4), intent(in)       ::  flgdal                          !<  whether the dual approximation is used
  character(len=6), intent(in)       ::  flgscf                          !<  type of self consistent field and diagonalization

  real(REAL64), intent(in)           ::  emax                            !<  kinetic energy cutoff of plane wave expansion (Hartree).
  character(len=*), intent(in)       ::  author                          !<  type of xc wanted (ca=pz , pw92 , pbe,...)
  real(REAL64), intent(in)           ::  tblaha                          !<  Tran-Blaha constant

  real(REAL64), intent(in)           ::  epscv                           !<  convergence criteria for potential self consistency
  real(REAL64), intent(in)           ::  epscvao                         !<  convergence criteria for potential self consistency in atomic orbitals
  real(REAL64), intent(in)           ::  epspsi                          !<  convergence criteria for iterative diagonalization

  real(REAL64), intent(in)           ::  tempk                           !<  ionic temperature (in Kelvin)
  real(REAL64), intent(in)           ::  teleck                          !<  electronic temperature (in Kelvin)
  real(REAL64), intent(in)           ::  tempinik                        !<  initial ionic temperature (in Kelvin)
  integer, intent(in)                ::  nstep                           !<  number of MD or minimization steps
  real(REAL64), intent(in)           ::  tstep                           !<  time step for MD (in atomic units)
  real(REAL64), intent(in)           ::  beta                            !<  friction coefficient/mass (in a.u.)
  integer, intent(in)                ::  iseed                           !<  initial seed for random number generator

  real(REAL64), intent(in)           ::  press                           !<  external pressure
  real(REAL64), intent(in)           ::  strext(3,3)                     !<  external applied stress
  real(real64), intent(in)           ::  celmas                          !<  fictitious cell mass

  logical, intent(in)                ::  flgkplusg                       !<  finish cell minimization with fixed k+G
  real(REAL64), intent(in)           ::  epskplusg                       !<  criteria for switching to fixed k+G

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks

! parameters

  real(REAL64), parameter :: AUTOGPA = 29421.58_REAL64

! counters

  integer   ::  i, j


  write(6,*)

  call xc_author_print(author)

  write(6,*)

  if(flgdal == 'DUAL') then
    write(6,*) 'dual aproximation is used'
    write(6,*)
  endif

  write(6,'("    The energy cutoff for wave-function kinetic ",          &
        &   "energy is",f12.3," Hartree")') emax
  write(6,*)

  write(6,'("    SCF is converged when the difference in ",              &
        &   "potentials is less then ",f14.8)') epscv
  if(flgscf == 'AO    ' .or. flgscf == 'AOJC  ' .or.                     &
     flgscf == 'AOJCPW') then
    write(6,'("   For atomic orbitals the SCF parameter is",f14.8)') epscvao
  endif
  write(6,'("    Iterative diagonalizationis converged when the",        &
        &   "error in |h psi - e psi| is less then ",f14.8)') epspsi
  write(6,*)

  write(6,'("    The temperature for electron Fermi distribution",       &
        &   " is",f12.3,"Kelvin")') teleck
  write(6,*)

  write(6,*)
  if(flgscf == 'AO    ') write(6,*) '    Atomic orbital basis ',         &
           'calculation'
  if(flgscf == 'AOJC  ') write(6,*) '    Atomic orbital plus',           &
     'single relaxation in plane-wave basis calculation'
  if(flgscf == 'AOJCPW') write(6,*) '    Plane-wave basis ',             &
     'calculation with initial potential from atomic basis'
  if(flgscf == '    PW') write(6,*) '    Plane-wave basis ',             &
     'calculation'
  write(6,*)

  if(flgcal == 'MICRO ' .or. flgcal == 'LANG  ' .or.                     &
     flgcal == 'VCSLNG' .or. flgcal == 'VCSMIC' .or.                     &
     flgcal == 'EPILNG') then
    write(6,'("    The initial ion temperature is ",f13.3,               &
          &    " Kelvin")') tempinik
    write(6,'("    The time step is ",f12.3," au ",                      &
          &   "      ( * 2.4 10**-17 s)")') tstep
    write(6,'("    The number of MD steps is : ",i10)') nstep
  endif

  if(flgcal == 'LANG  ' .or. flgcal == 'VCSLNG' .or.                     &
     flgcal == 'VCSLNG') then
    write(6,'("    The thermal bath temperature is ",f12.3," Kelvin")') tempk
    write(6,'("    The relaxation time is",f10.2," time steps")')        &
           1.0/(beta*tstep)
    write(6,'("    The random number seed is : ",i10)') iseed
  endif

  if(flgcal == 'VCSLNG' .or. flgcal == 'VCSLBF' .or.                     &
     flgcal == 'EPILBF' .or. flgcal == 'VCSMIC' .or.                     &
     flgcal == 'EPILNG') then
    write(6,'("    The applied pressure is ",f12.3," GPa.  ",            &
          &   "The stress tensor is",3f10.2,/,65x,3f10.2,/,65x,3f10.2)') &
              press*AUTOGPA, ((strext(i,j)*AUTOGPA,i=1,3),j=1,3)
  endif
  if(flgcal == 'VCSLNG' .or. flgcal == 'VCSMIC' .or.                     &
     flgcal == 'EPILNG') then
    write(6,'("    The fictitious cell mass is",f12.3)') celmas
  endif

  if(flgcal == 'VCSLBF' .or. flgcal == 'LBFSYM' .or.                     &
     flgcal == 'EPILBF') then
    write(6,'("    Maximum number of steps is : ",i10)') nstep
  endif
  write(6,*)

  if(flgcal == 'VCSLBF' .or. flgcal == 'EPILBF') then
    if(flgkplusg) then
      write(6,'("    Switches to constant k+G when error is ",            &
            &   "smaller then: ",f16.8)') epskplusg
    endif
  endif
  write(6,*)


  return

end subroutine print_parameters
