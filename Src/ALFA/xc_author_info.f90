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

!>  Sets flags according to the exchange and correlation functionals.
!>  Concentrates dispersed conditional statements for ease
!>  of adding new functionals.
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         22 November 2025, 6 October 2026.
!>  \copyright    GNU Public License v2

subroutine xc_author_info(author, lxcgrad, lxclap, lxctau,               &
       lxctb09, lxccalc)

! Written 22 November 2025. JLM
! meta-GGA flags from xc_author_tau. 6 October 2026. JLM+claude

  implicit none

! input

  character(len=4), intent(in)       ::  author                          !<  type of xc wanted (ca=pz , pw92 , pbe,...)

!output

  logical, intent(out)               ::  lxcgrad                         !<  gradient of charge density should be calculated
  logical, intent(out)               ::  lxclap                          !<  laplacian of charge density should be calculated
  logical, intent(out)               ::  lxctau                          !<  Kinetic energy density should be calculated
  logical, intent(out)               ::  lxctb09                         !<  Tran-Blaha constant is present
  logical, intent(out)               ::  lxccalc                         !<  xc energy is calculted

! local variables

  character(len=4)                   ::  xcbase                          !  functional used in xc_mgga
  character(len=4)                   ::  tausrc                          !  source of tau

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks

  lxcgrad = .FALSE.
  lxclap = .FALSE.
  lxctau = .FALSE.
  lxctb09 = .FALSE.
  lxccalc = .TRUE.

  if(chrsameinfo(author, 'PBE')) then
       lxcgrad = .TRUE.
  endif

  if(chrsameinfo(author, 'TBL') .or. chrsameinfo(author, 'TB09') .or.    &
     chrsameinfo(author, 'BR') ) then
       lxcgrad = .TRUE.
       lxclap = .TRUE.
       lxctau = .TRUE.
       lxccalc = .FALSE.
  endif

  if(chrsameinfo(author, 'TBL') .or. chrsameinfo(author, 'TB09') ) then
       lxctb09 = .TRUE.
  endif

  call xc_author_tau(author, xcbase, tausrc)

  if(tausrc /= 'NONE') then
       lxcgrad = .TRUE.
       if(tausrc == 'PSI ') lxctau = .TRUE.
  endif

  return

end subroutine xc_author_info


!>  Sets flags according to the families of exchange and correlation functionals.


subroutine xc_author_family(author, lxclda, lxcgga, lxcmgga, lxcmggavxc)

! Written 22 November 2025. JLM
! meta-GGA flag from xc_author_tau. 6 October 2026. JLM+claude

  implicit none

! input

  character(len=4), intent(in)       ::  author                          !<  type of xc wanted (ca=pz , pw92 , pbe,...)

!output

  logical, intent(out)               ::  lxclda                          !<  LDA functionals
  logical, intent(out)               ::  lxcgga                          !<  GGA functionals
  logical, intent(out)               ::  lxcmgga                         !<  meta-GGA functionals with energy and potential
  logical, intent(out)               ::  lxcmggavxc                      !<  meta-GGA functionals with only the potential

! local variables

  character(len=4)                   ::  xcbase                          !  functional used in xc_mgga
  character(len=4)                   ::  tausrc                          !  source of tau

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks


! set to false.  Only one can be true...

  lxclda = .FALSE.
  lxcgga = .FALSE.
  lxcmgga = .FALSE.
  lxcmggavxc = .FALSE.

  if( chrsameinfo(author, 'PZ' ) .or. chrsameinfo(author, 'CA' ) .or.    &
      chrsameinfo(author, 'PW92' ) .or. chrsameinfo(author, 'VWN' ) .or. &
      chrsameinfo(author, 'WI' ) ) then
        lxclda = .TRUE.
  endif

  if(chrsameinfo(author, 'PBE' )) then
       lxcgga = .TRUE.
  endif

  call xc_author_tau(author, xcbase, tausrc)

  if(tausrc /= 'NONE') then
       lxcmgga = .TRUE.
  endif

  if(chrsameinfo(author, 'TBL') .or. chrsameinfo(author, 'TB09') .or.    &
      chrsameinfo(author, 'BR') ) then
       lxcmggavxc = .TRUE.
  endif

  return

end subroutine xc_author_family



!>  prints the information about exchange and correlation functionals.


subroutine xc_author_print(author)

! Written 25 November 2025. JLM

  implicit none

! input

  character(len=4), intent(in)       ::  author                          !<  type of xc wanted (ca=pz , pw92 , pbe,...)

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks

  if( chrsameinfo(author, 'CA' ) .or. chrsameinfo(author, 'PZ' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the local ",            &
      &   "density aproximation using Ceperley and Alder correlation")')
    write(6,'("  (as parametrized by Perdew and Zunger)")')
  elseif( chrsameinfo(author, 'PW92' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the local ",            &
      &   "density aproximation using Ceperley and Alder correlation")')
    write(6,'("  (as parametrized by Perdew and Wang (1992) )")')
  elseif( chrsameinfo(author, 'VWN' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the local ",            &
      &   "density aproximation using Ceperley and Alder correlation")')
    write(6,'("  (as parametrized by  Vosko, Wilk and Nusair)")')
  elseif( chrsameinfo(author, 'WI' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the local ",            &
      &   "density aproximation using Wigner correlation")')
  elseif( chrsameinfo(author, 'PBE' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the generalized",       &
      &   " gradient aproximation as parametrized by Perdew, Burke ",    &
      &   "and Ernzerhof")')
  elseif( chrsameinfo(author, 'TBL' ) .or. chrsameinfo(author, 'TB09' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the modified",          &
      &   " Becke-Johnson meta-GGA aproximation of Tran-Blaha ")')
  elseif( chrsameinfo(author, 'BR' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the meta-GGA",          &
      &   " of Becke and Roussel ")')
  elseif( chrsameinfo(author, 'TBL' ) .or. chrsameinfo(author, 'TB09' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the modified",          &
      &   " Becke-Johnson meta-GGA aproximation of Tran-Blaha ")')
  elseif( chrsameinfo(author, 'LAK' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the meta-GGA",          &
      &   " of Lebeda, Aschebrock and Kummel:  LAK")')
  elseif( chrsameinfo(author, 'TASK' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the meta-GGA",          &
      &   " of Aschebrock and Kummel:   TASK")')
  elseif( chrsameinfo(author, 'R2SC' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated in the meta-GGA",          &
      &   " r2SCAN of Furness et al.:   R2SC")')
  elseif( chrsameinfo(author, 'R2TF' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated with the meta-GGA",        &
      &   " r2SCAN with the Thomas-Fermi tau:   R2TF")')
  elseif( chrsameinfo(author, 'R2TW' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated with the meta-GGA",        &
      &   " r2SCAN with the Thomas-Fermi-von Weizsacker tau:   R2TW")')
  elseif( chrsameinfo(author, 'TATF' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated with the meta-GGA",        &
      &   " TASK with the Thomas-Fermi tau:   TATF")')
  elseif( chrsameinfo(author, 'TATW' ) ) then
    write(6,*)
    write(6,'("  The potential was calculated with the meta-GGA",        &
      &   " TASK with the Thomas-Fermi-von Weizsacker tau:   TATW")')
  else
    write(6,*)
    write(6,'("  The XC flag is:   ",a4)') author
  endif

  return

end subroutine xc_author_print


!>  For meta-GGA functionals gives the functional used in xc_mgga
!>  and the source of the kinetic energy density tau.
!>
!>    tausrc = 'PSI '  tau from the wave-functions (generalized Kohn-Sham)
!>    tausrc = 'TF  '  Thomas-Fermi tau of the density (deorbitalized, Kohn-Sham)
!>    tausrc = 'TFVW'  Thomas-Fermi + von Weizsacker tau of the density (deorbitalized)
!>    tausrc = 'NONE'  not a meta-GGA with energy and potential
!>
!>  New variants of a functional need only one line here
!>  (and a line in xc_author_print).
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         6 October 2026.
!>  \copyright    GNU Public License v2

subroutine xc_author_tau(author, xcbase, tausrc)

! Written 6 October 2026. JLM+claude

  implicit none

! input

  character(len=4), intent(in)       ::  author                          !<  type of xc wanted (ca=pz , pw92 , pbe,...)

! output

  character(len=4), intent(out)      ::  xcbase                          !<  meta-GGA functional used in xc_mgga
  character(len=4), intent(out)      ::  tausrc                          !<  source of tau: 'PSI ', 'TF  ', 'TFVW', or 'NONE'

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks


  xcbase = author
  tausrc = 'NONE'

  if(chrsameinfo(author, 'LAK')) then
    xcbase = 'LAK '
    tausrc = 'PSI '
  elseif(chrsameinfo(author, 'TASK')) then
    xcbase = 'TASK'
    tausrc = 'PSI '
  elseif(chrsameinfo(author, 'R2SC')) then
    xcbase = 'R2SC'
    tausrc = 'PSI '
  elseif(chrsameinfo(author, 'R2TF')) then
    xcbase = 'R2SC'
    tausrc = 'TF  '
  elseif(chrsameinfo(author, 'R2TW')) then
    xcbase = 'R2SC'
    tausrc = 'TFVW'
  elseif(chrsameinfo(author, 'TATF')) then
    xcbase = 'TASK'
    tausrc = 'TF  '
  elseif(chrsameinfo(author, 'TATW')) then
    xcbase = 'TASK'
    tausrc = 'TFVW'
  endif

  return

end subroutine xc_author_tau
