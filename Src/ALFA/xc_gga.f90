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

!>  Computes the exchange correlation energy and potential
!>  in the generalized gradient approximation (GGA).
!>  non-spin-polarized
!>  Lengths in Bohr, energies in Hartrees.
!>  Adapted from a package by L.C. Balbas and J.M. Soler, Dec'96. Version 0.5.
!>
!>  \author       L.C. Balbas, J.M. Soler, J.L. Martins and J.M. Pacheco
!>  \version      5.12
!>  \date         23 february 1999, 12 November 2025.
!>  \copyright    GNU Public License v2

subroutine xc_gga( author, rho, grho,                                    &
                   epsx, epsc, dexdr, decdr, dexdgr, decdgr )


! Written/adapted 23 february 1999. jlm
! Modified, documentation, December 2019. JLM
! Modified, indentation, separated PBE from main code. 12 November 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author                          !<  type of xc wanted (ca=pz , pw92 , vwn, wi)
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdr                           !<  derivative with respect to rho of rho*epsx (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdr                           !<  derivative with respect to rho of rho*epsc (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdgr                          !<  derivative with respect to grho of rho*epsx (hartree/bohr^2)
  real(REAL64), intent(out)          ::  decdgr                          !<  derivative with respect to grho of rho*epsc (hartree/bohr^2)

! parameters

  real(REAL64), parameter  :: ZERO = 0.0_REAL64
  real(REAL64), parameter  :: EPS = 1.0E-24_REAL64


  if(rho < EPS) then

    epsx = ZERO
    epsc = ZERO
    dexdr = ZERO
    decdr = ZERO
    dexdgr = ZERO
    decdgr = ZERO

  else


    if (author == 'pbe' .or. author == 'PBE') then

!     Perdew, Burke and Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)
!     Phys. Rev. Lett. 78, 1396 (1997)

      call xc_gga_c_pbe(rho, grho, epsc, decdr, decdgr)

    else

      write(6,'("   STOPPED in xc_gga:   unknown correlation")')
      write(6,*) '     ',author

      stop

    endif

!   E x c h a n g e

    call xc_gga_x_pbe(rho, grho, epsx, dexdr, dexdgr)

  endif

  return

end subroutine xc_gga


!>  GGA exchange accoding to PBE

subroutine xc_gga_x_pbe(rho, grho, epsx, dexdr, dexdgr)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdr                           !<  derivative with respect to rho of rho*epsx (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdgr                          !<  derivative with respect to grho of rho*epsx (hartree/bohr^2)

! local variables

  real(REAL64)     ::  rs                                                !  Wigner-Seitz radius
  real(REAL64)     ::  drsdr                                             !  d rs / d rho
  real(REAL64)     ::  grloc                                             !  local value of gradient
  real(REAL64)     ::  a0                                                !  exchange constant
  real(REAL64)     ::  exlda                                             !  LDA exchange
  real(REAL64)     ::  xkf                                               !  Fermi wave-vector
  real(REAL64)     ::  fx, dfxds                                         !  enhancement factor

  real(REAL64)     ::  s, dsdr, dsdgr                                    !  function s(rho,grho) and serivatives

! parameters

  real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! PBE parameters

  real(REAL64), parameter  ::  BETA = 0.066725_REAL64
  real(REAL64), parameter  ::  KAPPA = 0.8040_REAL64
  real(REAL64), parameter  ::  MU = BETA * PI*PI / 3


  rs = (3/(4*PI*rho))**(UM/3)
  drsdr = - rs / (3*rho)

  grloc = grho
  if(grho < ZERO) grloc = ZERO

  a0 = (4 / (9*PI))**(UM/3)
  exlda = -((3*UM) / (4*UM)) / (PI*a0*rs)

  xkf = (3*PI*PI * rho)**(UM/3)
  s = grloc / (2 * xkf * rho)
  dsdr = -4*s / (3*rho)
  dsdgr = UM / (2*xkf*rho)

  fx = UM + KAPPA - KAPPA / (UM + MU*s*s/KAPPA)

  epsx = exlda * fx

  dfxds = 2*MU*s /  ((UM + MU*s*s/KAPPA)*(UM + MU*s*s/KAPPA))

  dexdr = ((4*UM) / (3*UM))*exlda*fx + rho*exlda*dfxds*dsdr
  dexdgr = rho*exlda*dfxds*dsdgr

  return

end subroutine xc_gga_x_pbe


!>  GGA correlation according to PBE

subroutine xc_gga_c_pbe(rho, grho, epsc, decdr, decdgr)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)

! output

  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdr                           !<  derivative with respect to rho of rho*epsc (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdgr                          !<  derivative with respect to grho of rho*epsc (hartree/bohr^2)

! local variables

  real(REAL64)     ::  grloc                                             !  local value of the gradient

  real(REAL64)     ::  ecunif, vcunif                                    !  LDA correlation (PW92)
  real(REAL64)     ::  xkf, xks                                          !  Fermi wave-vector, Thomas-Fermi screening
  real(REAL64)     ::  t, dtdr, dtdgr                                    !  function t(rho,grho) and derivatives
  real(REAL64)     ::  expec
  real(REAL64)     ::  a, dadr                                           !  function a(rho)
  real(REAL64)     ::  at2,dat2dr                                        !  a*t*t
  real(REAL64)     ::  qat2, dqat2dat2                                   !  function q(at2) and derivative
  real(REAL64)     ::  arglog, dargdr, dargdt                            !  argument of logarithm and derivatives
  real(REAL64)     ::  h, dhdr, dhdt                                     !  correction to LDA and derivatives
  real(REAL64)     ::  decudr

! parameters

  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64


! PBE parameters

  real(REAL64), parameter  ::  BETA = 0.066725_REAL64                    !  values from original paper
  real(REAL64), parameter  ::  GAM =  0.031091_REAL64                            !  values from original paper
!  real(REAL64), parameter  ::  BETA = 0.0667245506031492_REAL64         !  values from libxc
!  real(REAL64), parameter  ::  GAM =  0.0310906908696549_REAL64         !  values from libxc


  call xc_lda_c_pw92(rho, ecunif, vcunif)
  decudr = (vcunif - ecunif) / rho

  grloc = grho
  if(grho < ZERO) grloc = ZERO

  xkf = (3*PI*PI * rho)**(UM/3)
  xks = sqrt(4*xkf / PI )

  t = grloc / (2 * xks * rho)
  dtdr = -7*t / (6*rho)
  dtdgr = UM / (2*xks*rho)

  expec = exp(-ecunif / GAM)
  a = (BETA / GAM) * (UM / (expec - UM))
  dadr = +(BETA / (GAM*GAM))*                                 &
     ( expec / ((expec - UM)*(expec - UM)) )*decudr

  at2 = a*t*t
  dat2dr = t*t*dadr + 2*a*t*dtdr

  qat2 = (UM + at2)/(UM + at2 + at2*at2)
  dqat2dat2 = - at2 * (2*UM + at2) /                              &
                ((UM + at2 + at2*at2)*(UM +at2 + at2*at2))

  arglog = UM + (BETA/GAM) * t*t * qat2
  dargdr = 2*(BETA/GAM)*t*qat2*dtdr +                           &
              (BETA/GAM)*t*t*dqat2dat2*dat2dr
  dargdt = 2*(BETA/GAM)*t*qat2 +                                &
           2*(BETA/GAM)*t*t*dqat2dat2*a*t

  h = GAM * log(arglog)
  dhdr = (GAM/arglog)*dargdr
  dhdt = (GAM/arglog)*dargdt

  epsc = ecunif + h

  decdr = vcunif + h + rho*dhdr

  decdgr = rho*dhdt*dtdgr

  return

end subroutine xc_gga_c_pbe
