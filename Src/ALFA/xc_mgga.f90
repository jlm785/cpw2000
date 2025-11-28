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
!>  in the meta generalized gradient approximation (GGA).
!>  non-spin-polarized
!>  Lengths in Bohr, energies in Hartrees.
!>
!>  \author       J.L. Martins
!>  \version      5.12
!>  \date         25 November 2025.
!>  \copyright    GNU Public License v2

subroutine xc_mgga( author, rho, grho, tau,                              &
                   epsx, epsc, dexdr, decdr, dexdgr, decdgr,             &
                   dexdtau, decdtau )


! Written 21 November 2025 from xc_gga . jlm
! Added other functionals. Removed (for now) lap_rho. 25 November 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author                          !<  type of xc wanted (ca=pz , pw92 , vwn, wi)
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)
!  real(REAL64), intent(in)           ::  lap_rho                         !<  Laplacian of charge density
  real(REAL64), intent(in)           ::  tau                             !<  kinetic energy density (hartree/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdr                           !<  derivative with respect to rho of rho*epsx (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdr                           !<  derivative with respect to rho of rho*epsc (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdgr                          !<  derivative with respect to grho of rho*epsx (hartree/bohr^2)
  real(REAL64), intent(out)          ::  decdgr                          !<  derivative with respect to grho of rho*epsc (hartree/bohr^2)
  real(REAL64), intent(out)          ::  dexdtau                         !<  derivative with respect to tau of rho*epsx (1/bohr^3)
  real(REAL64), intent(out)          ::  decdtau                         !<  derivative with respect to tau of rho*epsc (1/bohr^3)

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks

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
    dexdtau = ZERO
    decdtau = ZERO

  else


    if (chrsameinfo(author, 'LAK' ) ) then

!     T. Lebeda, T. Aschebrock, and S. Kummel,
!     Phys. Rev. Lett. 133, 136402 (2024)

      call xc_mgga_x_lak(rho, grho, tau, epsx, dexdr, dexdgr, dexdtau)

      call xc_mgga_c_lak(rho, grho, tau, epsc, decdr, decdgr, decdtau)

    elseif (chrsameinfo(author, 'TASK' ) ) then

!     T. Aschebrock, and S. Kummel,
!     Phys.Rev.Research. 1, 033082 (2019)

      call xc_mgga_x_task(rho, grho, tau, epsx, dexdr, dexdgr, dexdtau)

      call xc_lda_c_pw92(rho, epsc, decdr)
      decdgr = ZERO
      decdtau = ZERO

    elseif (chrsameinfo(author, 'R2SC' ) ) then

!     Furness, Kaplan, Ning, Perdew and Sun,
!     J.Phys.Chem.Lett. 11, 8208 (2020)

      call xc_mgga_x_r2scan(rho, grho, tau, epsx, dexdr, dexdgr, dexdtau)

      call xc_mgga_c_r2scan(rho, grho, tau, epsc, decdr, decdgr, decdtau)


    else

      write(6,'("   STOPPED in xc_gga:   unknown correlation")')
      write(6,*) '     ',author

      stop

    endif


  endif

  return

end subroutine xc_mgga





!>  Computes the exchange energy and potential
!>  in the meta generalized gradient approximation (mGGA).
!>  T. Lebeda, T. Aschebrock, and S. Kummel, Phys. Rev. Lett. 133, 136402 (2024)
!>  Lengths in Bohr, energies in Hartrees.
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         12 November 2025.
!>  \copyright    GNU Public License v2

subroutine xc_mgga_x_lak(rho, grho, tau, epsx, dexdr, dexdgr, dexdtau )


! Written 12 November 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)
  real(REAL64), intent(in)           ::  tau                             !<  kinetic energy density (hartree/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdr                           !<  derivative with respect to rho of rho*epsx (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdgr                          !<  derivative with respect to grho of rho*epsx (hartree/bohr^2)
  real(REAL64), intent(out)          ::  dexdtau                         !<  derivative with respect to tau of rho*epsx (dimensionless)

! local variables

  real(REAL64)     ::  rs                                                !  Wigner-Seitz radius
  real(REAL64)     ::  a0                                                !  exchange constant
  real(REAL64)     ::  exlda, d_exlda_dr                                 !  LDA exchange energy
  real(REAL64)     ::  xkf, d_xkf_dr                                     !  Fermi wave-vector and derivative
  real(REAL64)     ::  grloc                                             !  local value of grho
  real(REAL64)     ::  fx, d_fx_dgr, d_fx_dr, d_fx_ds, d_fx_dtau         !  enancement factor, and derivative

  real(REAL64)     ::  s, d_s_dr, d_s_dgr                                !  s parameter and derivatives

  real(REAL64)     ::  tausingle, d_tausingle_dr, d_tausingle_dgr        !  single orbital limit of tau
  real(REAL64)     ::  tauunif, d_tauunif_dr                             !  uniform limit of tau
  real(REAL64)     ::  alpha, d_alpha_dr, d_alpha_dgr, d_alpha_dtau      !  alpha parameter and derivatives

  real(REAL64)     ::  fxalpha, d_fxalpha_dalpha
  real(REAL64)     ::  d_fxalpha_dr, d_fxalpha_dgr, d_fxalpha_dtau
  real(REAL64)     ::  gxs, d_gxs_ds
  real(REAL64)     ::  hxge4, d_hxge4_ds
  real(REAL64)     ::  xkx, d_xkx_ds
  real(REAL64)     ::  h1x, d_h1x_ds
  real(REAL64)     ::  gnum, d_gnum_ds

  real(REAL64)     ::  arg, d_arg_dalpha, d_arg_ds

! parameters

  real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  AKF = (9*PI/4)**(UM/(3*UM))

! LAK_x parameters

  real(REAL64), parameter  ::  H0X = 1.174_REAL64
  real(REAL64), parameter  ::  BX = 4.9479_REAL64
  real(REAL64), parameter  ::  AX = 1.1_REAL64
  real(REAL64), parameter  ::  ANUM = 5.0_REAL64

! LAK_x composed constants

  real(REAL64), parameter  ::  XMUAX = -( 97*UM + 3*H0X + sqrt(9*H0X*H0X + 74166*H0X-64175*UM) ) / (1200*UM)
  real(REAL64), parameter  ::  XNUA = (73*UM-50*XMUAX) / (5000*UM)
  real(REAL64), parameter  ::  XMUSX = (10*UM + 60*XMUAX) / (81*UM)
  real(REAL64), parameter  ::  XNUS = -(1606*UM - 50*XMUAX) / (18225*UM)

  real(REAL64), parameter  ::  C1 = XMUAX / (H0X - UM)
  real(REAL64), parameter  ::  C2 = (XMUAX + XNUA) / (H0X - UM)

! Corrections to fxalpha expansion

  real(REAL64), parameter  ::  EPS = 1.0E-4_REAL64
  real(REAL64), parameter  ::  CFX1 = 4/(PI*PI*C1)
  real(REAL64), parameter  ::  CFX2 = -(32/(PI*PI*PI*C1*C1))*(1-0.0148812_REAL64)
  real(REAL64), parameter  ::  CFX3 = -0.6963_REAL64

! Avoids calculating irrelevant exponentials compared to unity

  real(REAL64), parameter  ::  ARGMAX = 50.0_REAL64



  rs = (3/(4*PI*rho))**(UM/3)

  grloc = grho
  if(grho < ZERO) grloc = ZERO

  a0 = (4 / (9*PI))**(UM/3)
  exlda = -((3*UM) / (4*UM)) / (PI*a0*rs)
  d_exlda_dr = exlda / (3*rho)


  xkf = AKF / rs                                                         !  (3*PI*PI * rho)**(UM/3)
  d_xkf_dr = (UM/3) * xkf / rho                                          !  (UM/3)*(3*PI*PI * rho)**(UM/3) / rho

  s = grloc / (2 * xkf * rho)
  d_s_dr = -s/rho - (s/xkf)*d_xkf_dr
  d_s_dgr = UM / (2 * xkf * rho)

  tausingle = grloc*grloc / (8*rho)                                      !   xkf*xkf * rho * s*s /2
  d_tausingle_dr = -tausingle/rho
  d_tausingle_dgr = 2* grloc / (8*rho)
  tauunif = (3*UM / 10) * xkf*xkf * rho
  d_tauunif_dr = (3*UM / 10) * xkf*xkf + 2*(3*UM / 10) * xkf*d_xkf_dr * rho
  alpha = (tau - tausingle) / tauunif
  d_alpha_dr = -d_tausingle_dr/tauunif - alpha*d_tauunif_dr/tauunif
  d_alpha_dgr = -d_tausingle_dgr/tauunif
  d_alpha_dtau = UM/tauunif

  if(abs(alpha) > EPS) then
    arg = (PI/2) * ( C1*(alpha-UM)/alpha + C2*(alpha-UM)*(alpha-UM) )
    d_arg_dalpha = (PI/2) * ( C1/(alpha*alpha) + 2*C2*(alpha-UM) )
    fxalpha = (2/PI) * atan(arg)
    d_fxalpha_dalpha = (2/PI) * (UM/(UM+arg*arg)) * d_arg_dalpha
  else
    if(alpha > 0) then
      fxalpha =  1
    else
      fxalpha =  -1
    endif
    fxalpha = fxalpha + alpha*(CFX1 + alpha*(CFX2 + alpha*CFX3))
    d_fxalpha_dalpha = CFX1 + alpha*(2*CFX2 + alpha*3*CFX3)
  endif
  d_fxalpha_dr = d_fxalpha_dalpha*d_alpha_dr
  d_fxalpha_dgr = d_fxalpha_dalpha*d_alpha_dgr
  d_fxalpha_dtau = d_fxalpha_dalpha*d_alpha_dtau

  if(s > (BX/ARGMAX)*(BX/ARGMAX)) then
    gxs = UM - exp(-BX/sqrt(s))
    d_gxs_ds = (gxs-UM) * BX / (2*sqrt(s)*s)
  else
    gxs = UM
    d_gxs_ds = ZERO
  endif

  hxge4 = UM + XMUSX*s*s + XNUS*s*s*s*s + H0X*(UM - gxs)
  d_hxge4_ds = 2*XMUSX*s + 4*XNUS*s*s*s - H0X*d_gxs_ds

  arg = ((s*s)/(AX*AX))*(UM+s*s)
  d_arg_ds = ((2*s)/(AX*AX))*(UM+s*s) + ((s*s)/(AX*AX))*2*s
  if(arg < UM/ARGMAX) then
    xkx = ZERO
    d_xkx_ds = ZERO
  else
    xkx = exp(-UM/arg)
    d_xkx_ds = (UM / (arg*arg)) * xkx * d_arg_ds
  endif

  h1x = hxge4 + xkx*(AX-hxge4)
  d_h1x_ds = d_hxge4_ds + d_xkx_ds*(AX-hxge4) - xkx*d_hxge4_ds


  if(s > ANUM/sqrt(ARGMAX)) then
    arg = ANUM*ANUM/(s*s)
    gnum = UM - exp(-arg)
    d_gnum_ds = -exp(-arg) * 2*arg / s
  else
    gnum = UM
    d_gnum_ds = ZERO
  endif

  fx = H0X*gxs + (UM - fxalpha)*(h1x - H0X)*gnum
  d_fx_ds = H0X*d_gxs_ds + (UM - fxalpha)*d_h1x_ds*gnum + (UM - fxalpha)*(h1x - H0X)*d_gnum_ds
  d_fx_dr = d_fx_ds*d_s_dr - d_fxalpha_dr*(h1x - H0X)*gnum
  d_fx_dgr = d_fx_ds*d_s_dgr  - d_fxalpha_dgr*(h1x - H0X)*gnum
  d_fx_dtau = - d_fxalpha_dtau*(h1x - H0X)*gnum

  epsx = exlda * fx
  dexdr = epsx + rho*(d_exlda_dr*fx + exlda*d_fx_dr)
  dexdgr = exlda * d_fx_dgr * rho
  dexdtau = exlda * d_fx_dtau * rho

  return

end subroutine xc_mgga_x_lak




!>  Computes the exchange energy and potential
!>  in the meta generalized gradient approximation (mGGA).
!>  T. Lebeda, T. Aschebrock, and S. Kummel, Phys. Rev. Lett. 133, 136402 (2024)
!>  Lengths in Bohr, energies in Hartrees.
!>  uses comments on mgga_c_lak.mpl of libxc from git November 2025
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         20 November 2025.
!>  \copyright    GNU Public License v2

subroutine xc_mgga_c_lak( rho, grho, tau, epsc, decdr, decdgr, decdtau )


! Written 20 November 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)
  real(REAL64), intent(in)           ::  tau                             !<  kinetic energy density (hartree/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdr                           !<  derivative with respect to rho of rho*epsc (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdgr                          !<  derivative with respect to grho of rho*epsc (hartree/bohr^2)
  real(REAL64), intent(out)          ::  decdtau                         !<  derivative with respect to tau of rho*epsc (dimensionless)

! local variables

  real(REAL64)     ::  rs                                                !  Wigner-Seitz radius
  real(REAL64)     ::  xkf, d_xkf_dr                                     !  Fermi wave-vector and derivative
  real(REAL64)     ::  grloc                                             !  local value of grho

  real(REAL64)     ::  s, d_s_dr, d_s_dgr                                !  s parameter and derivatives
  real(REAL64)     ::  d_rs_dr                                           !  d rs / d rho

  real(REAL64)     ::  tausingle, d_tausingle_dr, d_tausingle_dgr        !  single orbital limit of tau
  real(REAL64)     ::  tauunif, d_tauunif_dr                             !  uniform limit of tau
  real(REAL64)     ::  alpha, d_alpha_dr, d_alpha_dgr, d_alpha_dtau      !  alpha parameter and derivatives

  real(REAL64)     ::  gnum, d_gnum_ds

  real(REAL64)     ::  eclda, vclda, d_eclda_dr
  real(REAL64)     ::  eclda0, d_eclda0_drs

  real(REAL64)     ::  cmua, xmuac, d_cmua_drs, d_xmuac_drs
  real(REAL64)     ::  betalfa, d_betalfa_drs
  real(REAL64)     ::  tt, d_tt_drs, d_tt_ds
  real(REAL64)     ::  xmusc, d_xmusc_drs
  real(REAL64)     ::  betat, d_betat_drs
  real(REAL64)     ::  w1, d_w1_drs
  real(REAL64)     ::  afunc, d_afunc_drs
  real(REAL64)     ::  g1, g2, g3
  real(REAL64)     ::  d_g1_drs, d_g2_drs, d_g3_drs
  real(REAL64)     ::  d_g1_ds, d_g2_ds, d_g3_ds
  real(REAL64)     ::  h1, d_h1_drs, d_h1_ds
  real(REAL64)     ::  w0, d_w0_drs
  real(REAL64)     ::  h0, d_h0_drs, d_h0_ds
  real(REAL64)     ::  ec0, d_ec0_dr, d_ec0_dgr
  real(REAL64)     ::  ec1, d_ec1_dr,d_ec1_dgr
  real(REAL64)     ::  abar, d_abar_dr, d_abar_dgr, d_abar_dtau, d_abar_dalpha
  real(REAL64)     ::  fc, d_fc_dr, d_fc_dgr, d_fc_dtau
  real(REAL64)     ::  fcge2, d_fcge2_dr

  real(REAL64)     ::  xnumer, xdenom1,  xdenom2
  real(REAL64)     ::  d_xnumer_drs, d_xdenom1_drs, d_xdenom2_drs
  real(REAL64)     ::  arg, d_arg_ds

! parameters

  real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  AKF = (9*PI/4)**(UM/(3*UM))

! LAK_c parameters

  real(REAL64), parameter  ::  H0X = 1.174_REAL64
  real(REAL64), parameter  ::  B1C = 0.0468_REAL64
  real(REAL64), parameter  ::  B2C = 0.205601_REAL64
  real(REAL64), parameter  ::  B3C = 2.85_REAL64
  real(REAL64), parameter  ::  CHI0 = 1.55344_REAL64
  real(REAL64), parameter  ::  AC = 10.0_REAL64
  real(REAL64), parameter  ::  ANUM = 5.0_REAL64
  real(REAL64), parameter  ::  AX = -(3*UM/4)*(3/PI)**(UM/3)
  real(REAL64), parameter  ::  AXC = (3*UM/(4*PI))**(UM/3)*AX
  real(REAL64), parameter  ::  GAM = (UM-log(2*UM)) / (PI*PI)
  real(REAL64), parameter  ::  CTP = (3*PI*PI/(16*UM))**(2*UM/3)

! LAK_c composed constants

  real(REAL64), parameter  ::  XMUAX = -( 97*UM + 3*H0X + sqrt(9*H0X*H0X + 74166*H0X-64175*UM) ) / (1200*UM)   ! XMUAX = -0.209897

  real(REAL64), parameter  ::  CS0 = -16*PI*(3*PI*PI)**(UM/3) * 2.568_REAL64 / (3000*10*UM / 81*UM)     !  CS0 = -1.0782

! Avoids calculating irrelevant exponentials compared to unity

  real(REAL64), parameter  ::  ARGMAX = 50.0_REAL64

! Uses expansion for arctan

  real(REAL64), parameter  ::  ATANMAX = 500.0_REAL64

! Avoids division by zero

  real(REAL64), parameter  ::  DBL_MIN = 10*tiny(UM), DBL_EPSILON = 4*epsilon(UM)


  rs = (3/(4*PI*rho))**(UM/3)
  d_rs_dr = - rs / (3*rho)

  grloc = grho
  if(grho < ZERO) grloc = ZERO

  call xc_lda_c_pw92(rho, eclda, vclda)
  d_eclda_dr = (vclda-eclda) / rho

  xkf = AKF / rs                                                         !  (3*PI*PI * rho)**(UM/3)
  d_xkf_dr = (UM/3) * xkf / rho                                          !  (UM/3)*(3*PI*PI * rho)**(UM/3) / rho

  s = grloc / (2 * xkf * rho)
  d_s_dr = -s/rho - (s/xkf)*d_xkf_dr
  d_s_dgr = UM / (2 * xkf * rho)

  tausingle = grloc*grloc / (8*rho)
  d_tausingle_dr = -tausingle/rho
  d_tausingle_dgr = 2* grloc / (8*rho)
  tauunif = (3*UM / 10) * xkf*xkf * rho
  d_tauunif_dr = (3*UM / 10) * xkf*xkf + 2*(3*UM / 10) * xkf*d_xkf_dr * rho
  alpha = (tau - tausingle) / tauunif
  d_alpha_dr = -d_tausingle_dr/tauunif - alpha*d_tauunif_dr/tauunif
  d_alpha_dgr = -d_tausingle_dgr/tauunif
  d_alpha_dtau = UM/tauunif

  xnumer = UM + 0.1_REAL64*rs**0.65_REAL64
  d_xnumer_drs = 0.065_REAL64 / rs**0.35_REAL64
  xdenom1 = (UM + 0.065_REAL64*rs**0.9_REAL64)
  d_xdenom1_drs = 0.0585_REAL64 / rs**0.1_REAL64
  xdenom2 = (UM+0.03_REAL64*rs**1.2_REAL64)
  d_xdenom2_drs = 0.036_REAL64*rs**0.2_REAL64

  cmua = CS0*xnumer / ( xdenom1*xdenom2 )
  d_cmua_drs = CS0*d_xnumer_drs / ( xdenom1*xdenom2 ) -                        &
             CS0*xnumer*d_xdenom1_drs / ( xdenom1*xdenom1*xdenom2 ) -          &
             CS0*xnumer*d_xdenom2_drs / ( xdenom1*xdenom2*xdenom2 )

  xmuac = -XMUAX*( UM + cmua/2 )
  betalfa = -AXC*xmuac
  d_xmuac_drs = -XMUAX*d_cmua_drs / 2
  d_betalfa_drs = -AXC*d_xmuac_drs

  tt = CTP*s*s/rs
  d_tt_drs = -CTP*s*s/(rs*rs)
  d_tt_ds = 2*CTP*s/rs
  xmusc = (10*UM/81*UM) * ( cmua*(UM-3*XMUAX) - (UM+6*XMUAX) )
  d_xmusc_drs = (10*UM/81*UM)*(UM-3*XMUAX)*d_cmua_drs
  betat = (AXC/CTP) * xmusc

  d_betat_drs = (AXC/CTP) * d_xmusc_drs

  w1 = exp(-eclda/GAM) - UM
  d_w1_drs = -exp(-eclda/GAM)*(d_eclda_dr / d_rs_dr) / GAM

  afunc = betat / (GAM*w1)
  d_afunc_drs = d_betat_drs / (GAM*w1) - afunc*d_w1_drs / w1

  g2 = UM / (UM + afunc*afunc*tt*tt)
  d_g2_drs = -2*(afunc*tt * g2*g2)*(d_afunc_drs*tt + d_tt_drs*afunc)
  d_g2_ds = -2*d_tt_ds*afunc*afunc*tt * g2*g2
  if(betat > ZERO) then
    if(UM + 4*afunc*tt > DBL_EPSILON) then
      g1 = UM / sqrt(sqrt(UM + 4*afunc*tt))
      d_g1_drs = -g1 *(d_afunc_drs*tt+afunc*d_tt_drs) / (UM + 4*afunc*tt)
      d_g1_ds = -g1 *afunc*d_tt_ds / (UM + 4*afunc*tt)
    else
      g1 = UM / DBL_EPSILON
      d_g1_drs = ZERO
      d_g1_ds = ZERO
    endif
    if(abs(UM+AC*afunc*tt) > DBL_EPSILON) then
      g3 = UM / (UM+AC*afunc*tt)
      d_g3_drs = -g3*g3 * AC * (d_afunc_drs*tt + afunc*d_tt_drs)
      d_g3_ds = -g3*g3 * AC * afunc*d_tt_ds
    else
      g3 = UM / DBL_EPSILON
      d_g3_drs = ZERO
      d_g3_ds = ZERO
    endif
  else
    if(UM - 4*afunc*tt > DBL_EPSILON) then
      g1 = UM / sqrt(sqrt(UM - 4*afunc*tt))
      d_g1_drs = g1 *(d_afunc_drs*tt+afunc*d_tt_drs) / (UM - 4*afunc*tt)
      d_g1_ds = g1 *afunc*d_tt_ds / (UM - 4*afunc*tt)
    else
      g1 = UM / DBL_EPSILON
      d_g1_drs = ZERO
      d_g1_ds = ZERO
    endif
    if(abs(UM-(w1+B3C)*afunc*tt) > DBL_EPSILON) then
      g3 =-UM / (UM-(w1+B3C)*afunc*tt)
      d_g3_drs = -g3*g3 * ((w1+B3C)*(d_afunc_drs*tt + afunc*d_tt_drs) + d_w1_drs*afunc*tt)
      d_g3_ds = -g3*g3 * (w1+B3C)*afunc*d_tt_ds
    else
      g3 = UM / DBL_EPSILON
      d_g3_drs = ZERO
      d_g3_ds = ZERO
    endif
  endif

  if(UM + w1*(UM-g1)*(UM-g2+g3) > DBL_MIN) then
    h1 = GAM*log( UM + w1*(UM-g1)*(UM-g2+g3) )
    d_h1_drs = GAM*( d_w1_drs*(UM-g1)*(UM-g2+g3) - w1*d_g1_drs*(UM-g2+g3)        &
                     - w1*(UM-g1)*(d_g2_drs-d_g3_drs) )  /                       &
                     ( UM + w1*(UM-g1)*(UM-g2+g3) )

    d_h1_ds = -GAM*w1*(d_g1_ds*(UM-g2+g3) + (UM-g1)*(d_g2_ds-d_g3_ds) ) /        &
                     ( UM + w1*(UM-g1)*(UM-g2+g3) )
  else
    h1 = GAM*log(DBL_MIN)
    d_h1_drs = ZERO
    d_h1_ds = ZERO
  endif

  eclda0 = -B1C / (UM + B2C*rs)
  d_eclda0_drs = eclda0*eclda0*B2C/B1C
  w0 = exp(-eclda0/B1C) - UM
  d_w0_drs = -exp(-eclda0/B1C)*d_eclda0_drs/B1C

  arg = UM + 4*CHI0*s*s
  d_arg_ds = 8*CHI0*s
  h0 = B1C * log( UM + w0*( UM - UM / (sqrt(sqrt(arg))) ) )
  d_h0_drs = B1C*d_w0_drs*( UM - UM / (sqrt(sqrt(arg))))         /             &
              ( UM + w0*( UM - UM / (sqrt(sqrt(arg))) ) )
  d_h0_ds = B1C*(w0/4)*(UM / (sqrt(sqrt(arg))))*(UM/arg)*d_arg_ds   /          &
              ( UM + w0*( UM - UM / (sqrt(sqrt(arg))) ) )

  ec0 = eclda0 + h0
  d_ec0_dr = d_eclda0_drs*d_rs_dr + d_h0_drs*d_rs_dr + d_h0_ds*d_s_dr
  d_ec0_dgr = d_h0_ds*d_s_dgr

  ec1 = eclda + h1
  d_ec1_dr = d_eclda_dr + d_h1_drs*d_rs_dr + d_h1_ds*d_s_dr
  d_ec1_dgr = d_h1_ds*d_s_dgr

  if(eclda - eclda0 < -DBL_EPSILON) then                                        !  from maple file
    fcge2 = betalfa / (eclda - eclda0)
!    fcge2 = betalfa / (eclda - ec0)
    d_fcge2_dr = d_betalfa_drs*d_rs_dr / (eclda - eclda0) -                    &
             fcge2*d_eclda_dr / (eclda - eclda0) +                             &
             fcge2*d_eclda0_drs*d_rs_dr / (eclda - eclda0)
  else
    fcge2 = -betalfa / DBL_EPSILON
    d_fcge2_dr = -d_betalfa_drs*d_rs_dr / DBL_EPSILON
  endif

  if(abs(alpha*rs) > DBL_EPSILON) then
    abar = (alpha - UM) / (alpha*rs)
    d_abar_dalpha = UM / (alpha*rs) - abar*rs / (alpha*rs)
  else
    abar = (alpha - UM) / DBL_EPSILON
    d_abar_dalpha = UM / DBL_EPSILON
  endif
  d_abar_dr = d_abar_dalpha*d_alpha_dr - abar*d_rs_dr / rs
  d_abar_dgr = d_abar_dalpha*d_alpha_dgr
  d_abar_dtau = d_abar_dalpha*d_alpha_dtau

  arg = (PI/2)*fcge2*abar
  if(arg < ATANMAX) then
    fc = (2/PI) * atan( arg )
  else
    if(arg > ZERO) then
      fc = UM
    else
      fc =-UM
    endif
    fc = fc -(2/PI)*(UM/arg*(UM - UM/(3*arg*arg)))
  endif
  d_fc_dr = (UM/(UM+arg*arg)) * (d_fcge2_dr*abar + fcge2*d_abar_dr)
  d_fc_dgr = (UM/(UM+arg*arg)) * fcge2*d_abar_dgr
  d_fc_dtau = (UM/(UM+arg*arg)) * fcge2*d_abar_dtau


  if(s > ANUM/sqrt(ARGMAX)) then
    arg = ANUM*ANUM/(s*s)
    gnum = UM - exp(-arg)
    d_gnum_ds = -exp(-arg) * 2*arg / s
  else
    gnum = UM
    d_gnum_ds = ZERO
  endif

  epsc = ec0 + (UM-fc)*(ec1-ec0)*gnum
  decdr = epsc + rho*( d_ec0_dr - d_fc_dr*(ec1-ec0)*gnum                 &
                    + (UM-fc)*(d_ec1_dr-d_ec0_dr)*gnum                   &
                    + (UM-fc)*(ec1-ec0)*d_gnum_ds*d_s_dr )
  decdgr = rho * ( d_ec0_dgr - d_fc_dgr*(ec1-ec0)*gnum                   &
                 + (UM-fc)*d_ec1_dgr*gnum - (UM-fc)*d_ec0_dgr*gnum       &
                 + (UM-fc)*(ec1-ec0)*d_gnum_ds*d_s_dgr )
  decdtau =-d_fc_dtau*(ec1-ec0)*gnum*rho

  return

end subroutine xc_mgga_c_lak



!>  T. Aschebrock and S. Kummel, Phys.Rev.Research. 1, 033082 (2019)

subroutine xc_mgga_x_task(rho, grho, tau, epsx, dexdr, dexdgr, dexdtau )


! Written 12 November 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)
  real(REAL64), intent(in)           ::  tau                             !<  kinetic energy density (hartree/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdr                           !<  derivative with respect to rho of rho*epsx (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdgr                          !<  derivative with respect to grho of rho*epsx (hartree/bohr^2)
  real(REAL64), intent(out)          ::  dexdtau                         !<  derivative with respect to tau of rho*epsx (dimensionless)

! local variables

  real(REAL64)     ::  rs                                                !  Wigner-Seitz radius
  real(REAL64)     ::  exlda, d_exlda_dr                                 !  LDA exchange energy
  real(REAL64)     ::  xkf, d_xkf_dr                                     !  Fermi wave-vector and derivative
  real(REAL64)     ::  grloc                                             !  local value of grho
  real(REAL64)     ::  fx, d_fx_dgr, d_fx_dr, d_fx_ds, d_fx_dtau         !  enancement factor, and derivative

  real(REAL64)     ::  s, d_s_dr, d_s_dgr                                !  s parameter and derivatives

  real(REAL64)     ::  tausingle, d_tausingle_dr, d_tausingle_dgr        !  single orbital limit of tau
  real(REAL64)     ::  tauunif, d_tauunif_dr                             !  uniform limit of tau
  real(REAL64)     ::  alpha, d_alpha_dr, d_alpha_dgr, d_alpha_dtau      !  alpha parameter and derivatives

  real(REAL64)     ::  fxalpha, d_fxalpha_dalpha
  real(REAL64)     ::  d_fxalpha_dr, d_fxalpha_dgr, d_fxalpha_dtau
  real(REAL64)     ::  gxs, d_gxs_ds
  real(REAL64)     ::  h1x, d_h1x_ds

  real(REAL64)     ::  z, d_z_dalpha, d_z_ds                             !  Argument of Chebychev polynomial

! parameters

  real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  AKF = (9*PI/4)**(UM/(3*UM))

! TASK_x parameters

  real(REAL64), parameter  ::  AX = (4 / (9*PI))**(UM/3)

  real(REAL64), parameter  ::  H0X = 1.174_REAL64
  real(REAL64), parameter  ::  BX = 4.9479_REAL64

! TASK_x Chebyshev coeficients

  real(REAL64), parameter  ::  A0 = 0.938719_REAL64
  real(REAL64), parameter  ::  A1 =-0.076371_REAL64
  real(REAL64), parameter  ::  A2 =-0.0150899_REAL64

  real(REAL64), parameter  ::  B0 =-0.628591_REAL64
  real(REAL64), parameter  ::  B1 =-2.10315_REAL64
  real(REAL64), parameter  ::  B2 =-0.5_REAL64
  real(REAL64), parameter  ::  B3 = 0.103153_REAL64
  real(REAL64), parameter  ::  B4 = 0.128591_REAL64

! Avoids calculating irrelevant exponentials compared to unity

  real(REAL64), parameter  ::  ARGMAX = 50.0_REAL64


  rs = (3/(4*PI*rho))**(UM/3)

  grloc = grho
  if(grho < ZERO) grloc = ZERO

  exlda = -((3*UM) / (4*UM)) / (PI*AX*rs)
  d_exlda_dr = exlda / (3*rho)


  xkf = AKF / rs                                                         !  (3*PI*PI * rho)**(UM/3)
  d_xkf_dr = (UM/3) * xkf / rho                                          !  (UM/3)*(3*PI*PI * rho)**(UM/3) / rho

  s = grloc / (2 * xkf * rho)
  d_s_dr = -s/rho - (s/xkf)*d_xkf_dr
  d_s_dgr = UM / (2 * xkf * rho)

  tausingle = grloc*grloc / (8*rho)                                      !   xkf*xkf * rho * s*s /2
  d_tausingle_dr = -tausingle/rho
  d_tausingle_dgr = 2* grloc / (8*rho)
  tauunif = (3*UM / 10) * xkf*xkf * rho
  d_tauunif_dr = (3*UM / 10) * xkf*xkf + 2*(3*UM / 10) * xkf*d_xkf_dr * rho
  alpha = (tau - tausingle) / tauunif
  d_alpha_dr = -d_tausingle_dr/tauunif - alpha*d_tauunif_dr/tauunif
  d_alpha_dgr = -d_tausingle_dgr/tauunif
  d_alpha_dtau = UM/tauunif

  z = (alpha-UM) / (alpha+UM)
  d_z_dalpha = 2 / ((alpha+UM)*(alpha+UM))
  fxalpha = B0 + B1*z + B2*(2*Z*Z-UM) + B3*(4*z*z*z-3*z) + B4*(8*z*z*z*z-8*z*z+UM)
  d_fxalpha_dalpha = B1 + B2*4*z + B3*(12*z*z-3*UM) + B4*(32*z*z*z-16*z)
  d_fxalpha_dalpha = d_fxalpha_dalpha*d_z_dalpha
  d_fxalpha_dr = d_fxalpha_dalpha*d_alpha_dr
  d_fxalpha_dgr = d_fxalpha_dalpha*d_alpha_dgr
  d_fxalpha_dtau = d_fxalpha_dalpha*d_alpha_dtau

  if(s > (BX/ARGMAX)*(BX/ARGMAX)) then
    gxs = UM - exp(-BX/sqrt(s))
    d_gxs_ds = (gxs-UM) * BX / (2*sqrt(s)*s)
  else
    gxs = UM
    d_gxs_ds = ZERO
  endif

  z = (s*s-UM) / (s*s+UM)
  d_z_ds = 4*s / ((s*s+UM)*(s*s+UM))
  h1x = A0 + A1*Z + A2*(2*Z*Z-UM)
  d_h1x_ds = (A1 + A2*4*z)*d_z_ds

  fx = H0X*gxs + (UM - fxalpha)*(h1x - H0X)*gxs**10
  d_fx_ds = H0X*d_gxs_ds + (UM - fxalpha)*d_h1x_ds*gxs**10 + (UM - fxalpha)*(h1x - H0X)*10*d_gxs_ds*gxs**9

  d_fx_dr = d_fx_ds*d_s_dr - d_fxalpha_dr*(h1x - H0X)*gxs**10
  d_fx_dgr = d_fx_ds*d_s_dgr  - d_fxalpha_dgr*(h1x - H0X)*gxs**10
  d_fx_dtau = - d_fxalpha_dtau*(h1x - H0X)*gxs**10

  epsx = exlda * fx
  dexdr = epsx + rho*(d_exlda_dr*fx + exlda*d_fx_dr)
  dexdgr = rho * exlda * d_fx_dgr
  dexdtau = rho * exlda * d_fx_dtau

  return

end subroutine xc_mgga_x_task

!>  Furness, Kaplan, Ning, Perdew and Sun, J.Phys.Chem.Lett. 11, 8208 (2020)


subroutine xc_mgga_x_r2scan(rho, grho, tau, epsx, dexdr, dexdgr, dexdtau )


! Written 25 November 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)
  real(REAL64), intent(in)           ::  tau                             !<  kinetic energy density (hartree/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdr                           !<  derivative with respect to rho of rho*epsx (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdgr                          !<  derivative with respect to grho of rho*epsx (hartree/bohr^2)
  real(REAL64), intent(out)          ::  dexdtau                         !<  derivative with respect to tau of rho*epsx (dimensionless)

! local variables

  real(REAL64)     ::  rs                                                !  Wigner-Seitz radius
  real(REAL64)     ::  exlda, d_exlda_dr                                 !  LDA exchange energy
  real(REAL64)     ::  xkf, d_xkf_dr                                     !  Fermi wave-vector and derivative
  real(REAL64)     ::  grloc                                             !  local value of grho
  real(REAL64)     ::  fx, d_fx_dgr, d_fx_dr, d_fx_ds, d_fx_dtau         !  enancement factor, and derivative

  real(REAL64)     ::  s, d_s_dr, d_s_dgr                                !  s parameter and derivatives

  real(REAL64)     ::  tausingle, d_tausingle_dr, d_tausingle_dgr        !  single orbital limit of tau
  real(REAL64)     ::  tauunif, d_tauunif_dr                             !  uniform limit of tau
  real(REAL64)     ::  alpha, d_alpha_dr, d_alpha_dgr, d_alpha_dtau      !  alpha parameter and derivatives

  real(REAL64)     ::  fxalpha, d_fxalpha_dalpha
  real(REAL64)     ::  d_fxalpha_dr, d_fxalpha_dgr, d_fxalpha_dtau
  real(REAL64)     ::  gxs, d_gxs_ds
  real(REAL64)     ::  h1x, d_h1x_ds

  real(REAL64)     ::  arg, d_arg_ds
  real(REAL64)     ::  xp, d_xp_ds


! parameters

  real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  AKF = (9*PI/4)**(UM/(3*UM))

! r2SCAN_x parameters

  real(REAL64), parameter  ::  AX = (4 / (9*PI))**(UM/3)

  real(REAL64), parameter  ::  H0X = 1.174_REAL64
  real(REAL64), parameter  ::  BX = 4.9479_REAL64
  real(REAL64), parameter  ::  XMU = (10*UM) / (81*UM)
  real(REAL64), parameter  ::  XK1 = 0.065_REAL64
  real(REAL64), parameter  ::  DX = 1.24_REAL64
  real(REAL64), parameter  ::  C1X = 0.667_REAL64
  real(REAL64), parameter  ::  C2X = 0.8_REAL64
  real(REAL64), parameter  ::  C2 = -0.162742_REAL64
  real(REAL64), parameter  ::  DP2 = 0.361_REAL64

  real(REAL64), parameter  ::  ETA = 0.001_REAL64
  real(REAL64), parameter  ::  CETA = (20*UM) / (27*UM) + (ETA*5) / (3*UM)

! Polynomial coeficients

  real(REAL64), parameter  ::  CX0 = UM
  real(REAL64), parameter  ::  CX1 =-0.667_REAL64
  real(REAL64), parameter  ::  CX2 =-0.4445555_REAL64
  real(REAL64), parameter  ::  CX3 =-0.663086601049_REAL64
  real(REAL64), parameter  ::  CX4 = 1.45129704449_REAL64
  real(REAL64), parameter  ::  CX5 =-0.887998041597_REAL64
  real(REAL64), parameter  ::  CX6 = 0.234528941479_REAL64
  real(REAL64), parameter  ::  CX7 =-0.023185843322_REAL64

! Avoids calculating irrelevant exponentials compared to unity

  real(REAL64), parameter  ::  ARGMAX = 50.0_REAL64


  rs = (3/(4*PI*rho))**(UM/3)

  grloc = grho
  if(grho < ZERO) grloc = ZERO

  exlda = -((3*UM) / (4*UM)) / (PI*AX*rs)
  d_exlda_dr = exlda / (3*rho)


  xkf = AKF / rs                                                         !  (3*PI*PI * rho)**(UM/3)
  d_xkf_dr = (UM/3) * xkf / rho                                          !  (UM/3)*(3*PI*PI * rho)**(UM/3) / rho

  s = grloc / (2 * xkf * rho)
  d_s_dr = -s/rho - (s/xkf)*d_xkf_dr
  d_s_dgr = UM / (2 * xkf * rho)

  tausingle = grloc*grloc / (8*rho)                                      !   xkf*xkf * rho * s*s /2
  d_tausingle_dr = -tausingle/rho
  d_tausingle_dgr = 2* grloc / (8*rho)
  tauunif = (3*UM / 10) * xkf*xkf * rho
  d_tauunif_dr = (3*UM / 10) * xkf*xkf + 2*(3*UM / 10) * xkf*d_xkf_dr * rho
  alpha = (tau - tausingle) / (tauunif + ETA*tausingle)
  d_alpha_dr = -d_tausingle_dr * ( UM / (tauunif + ETA*tausingle)        &
               + ETA*alpha / (tauunif + ETA*tausingle) )                 &
               - alpha*d_tauunif_dr/ (tauunif + ETA*tausingle)
  d_alpha_dgr = -d_tausingle_dgr * ( UM / (tauunif + ETA*tausingle)      &
               + ETA*alpha / (tauunif + ETA*tausingle) )
  d_alpha_dtau = UM / (tauunif + ETA*tausingle)

  if(alpha < ZERO) then
    fxalpha = exp(-c1x*alpha / (UM-alpha))
    d_fxalpha_dalpha = -c1x*fxalpha / ((UM-alpha)*(UM-alpha))
  elseif(alpha > 2.5) then
    fxalpha = -DX*exp(c2x*alpha / (UM-alpha))
    d_fxalpha_dalpha = c2x*fxalpha / ((UM-alpha)*(UM-alpha))
  else
    fxalpha = CX0 + alpha*(CX1 + alpha*(CX2 + alpha*(CX3 + alpha*(CX4 +    &
              alpha*(CX5 + alpha*(CX6+alpha*CX7))))))
    d_fxalpha_dalpha = CX1 + alpha*(2*CX2 + alpha*(3*CX3 + alpha*(4*CX4 +  &
              alpha*(5*CX5 + alpha*(6*CX6+alpha*7*CX7)))))
  endif

  d_fxalpha_dr = d_fxalpha_dalpha*d_alpha_dr
  d_fxalpha_dgr = d_fxalpha_dalpha*d_alpha_dgr
  d_fxalpha_dtau = d_fxalpha_dalpha*d_alpha_dtau

  if(s > (BX/ARGMAX)*(BX/ARGMAX)) then
    gxs = UM - exp(-BX/sqrt(s))
    d_gxs_ds = (gxs-UM) * BX / (2*sqrt(s)*s)
  else
    gxs = UM
    d_gxs_ds = ZERO
  endif

  arg = -s*s*s*s / (DP2*DP2*DP2*DP2)
  d_arg_ds = -4*s*s*s / (DP2*DP2*DP2*DP2)
  xp = ( CETA*C2*exp(arg) + XMU) * s*s
  d_xp_ds = 2*( CETA*C2*exp(arg) + XMU)*s + CETA*C2*exp(arg) * s*s * d_arg_ds
  h1x = UM + XK1 - XK1 / (UM + xp / XK1)
  d_h1x_ds = ( UM / ((UM + xp / XK1)*(UM + xp / XK1)) ) * d_xp_ds

  fx = H0X*gxs + (UM - fxalpha)*(h1x - H0X)*gxs
  d_fx_ds = H0X*d_gxs_ds + (UM - fxalpha)*d_h1x_ds*gxs + (UM - fxalpha)*(h1x - H0X)*d_gxs_ds

  d_fx_dr = d_fx_ds*d_s_dr - d_fxalpha_dr*(h1x - H0X)*gxs

  d_fx_dgr = d_fx_ds*d_s_dgr  - d_fxalpha_dgr*(h1x - H0X)*gxs
  d_fx_dtau = - d_fxalpha_dtau*(h1x - H0X)*gxs

  epsx = exlda * fx
  dexdr = epsx + rho*(d_exlda_dr*fx + exlda*d_fx_dr)
  dexdgr = rho * exlda * d_fx_dgr
  dexdtau = rho * exlda * d_fx_dtau

  return

end subroutine xc_mgga_x_r2scan


subroutine xc_mgga_c_r2scan( rho, grho, tau, epsc, decdr, decdgr, decdtau )


! Written 20 November 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)
  real(REAL64), intent(in)           ::  tau                             !<  kinetic energy density (hartree/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdr                           !<  derivative with respect to rho of rho*epsc (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdgr                          !<  derivative with respect to grho of rho*epsc (hartree/bohr^2)
  real(REAL64), intent(out)          ::  decdtau                         !<  derivative with respect to tau of rho*epsc (dimensionless)

! local variables

  real(REAL64)     ::  rs                                                !  Wigner-Seitz radius
  real(REAL64)     ::  xkf, d_xkf_dr                                     !  Fermi wave-vector and derivative
  real(REAL64)     ::  grloc                                             !  local value of grho

  real(REAL64)     ::  s, d_s_dr, d_s_dgr                                !  s parameter and derivatives
  real(REAL64)     ::  d_rs_dr                                           !  d rs / d rho

  real(REAL64)     ::  tausingle, d_tausingle_dr, d_tausingle_dgr        !  single orbital limit of tau
  real(REAL64)     ::  tauunif, d_tauunif_dr                             !  uniform limit of tau
  real(REAL64)     ::  alpha, d_alpha_dr, d_alpha_dgr, d_alpha_dtau      !  alpha parameter and derivatives

  real(REAL64)     ::  eclda, vclda, d_eclda_dr
  real(REAL64)     ::  eclda0, d_eclda0_drs

  real(REAL64)     ::  tt, d_tt_drs, d_tt_ds
  real(REAL64)     ::  betrs, d_betrs_drs
  real(REAL64)     ::  w1, d_w1_drs

  real(REAL64)     ::  yfunc, d_yfunc_drs,  d_yfunc_ds
  real(REAL64)     ::  dely, d_dely_drs, d_dely_ds
  real(REAL64)     ::  gg, d_gg_drs, d_gg_ds

  real(REAL64)     ::  h1, d_h1_drs, d_h1_ds
  real(REAL64)     ::  w0, d_w0_drs
  real(REAL64)     ::  h0, d_h0_drs, d_h0_ds
  real(REAL64)     ::  ec0, d_ec0_dr, d_ec0_dgr
  real(REAL64)     ::  ec1, d_ec1_dr,d_ec1_dgr

  real(REAL64)     ::  fc, d_fc_dalpha, d_fc_dr, d_fc_dgr, d_fc_dtau

  real(REAL64)     ::  arg, d_arg_ds
  real(REAL64)     ::  d_eclda_drs, d2_eclda_drs2, d2_eclda0_drs2

  real(REAL64)     ::  b_pw, db_pwdrs, d2b_pwdrs2
  real(REAL64)     ::  c_pw, dc_pwdrs, d2c_pwdrs2

  real(REAL64)     ::  arg1, d_arg1_drs
  real(REAL64)     ::  arg2, d_arg2_drs
  real(REAL64)     ::  arg3, d_arg3_ds


! parameters

  real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  AKF = (9*PI/4)**(UM/(3*UM))

! r2scan_c parameters

  real(REAL64), parameter  ::  B1C = 0.0285764_REAL64
  real(REAL64), parameter  ::  B2C = 0.0889_REAL64
  real(REAL64), parameter  ::  B3C = 0.125541_REAL64
  real(REAL64), parameter  ::  CHI0 = 0.128025_REAL64  !  0.1280504?
  real(REAL64), parameter  ::  GAM = (UM-log(2*UM)) / (PI*PI)

  real(REAL64), parameter  ::  BETAMB = 0.066724550603149220_REAL64  ! 0.066725_REAL64
  real(REAL64), parameter  ::  C1C = 0.64_REAL64
  real(REAL64), parameter  ::  C2C = 1.5_REAL64
  real(REAL64), parameter  ::  DC = 0.7_REAL64

  real(REAL64), parameter  ::  DP2 = 0.361_REAL64
  real(REAL64), parameter  ::  ETA = 0.001_REAL64

  real(REAL64), parameter  ::  BETAB = 0.1_REAL64
  real(REAL64), parameter  ::  BETAC = 0.1778_REAL64

  real(REAL64), parameter  ::  AX = -(3*UM/4)*(3/PI)**(UM/3)
  real(REAL64), parameter  ::  AXC = (3*UM/(4*PI))**(UM/3)*AX
  real(REAL64), parameter  ::  CTP = (3*PI*PI/(16*UM))**(2*UM/3)

! Polynomial coeficients

  real(REAL64), parameter  ::  CC0 = UM
  real(REAL64), parameter  ::  CC1 =-0.64_REAL64
  real(REAL64), parameter  ::  CC2 =-0.4352_REAL64
  real(REAL64), parameter  ::  CC3 =-1.535685604549_REAL64
  real(REAL64), parameter  ::  CC4 = 3.061560252175_REAL64
  real(REAL64), parameter  ::  CC5 =-1.915710236206_REAL64
  real(REAL64), parameter  ::  CC6 = 0.516884468372_REAL64
  real(REAL64), parameter  ::  CC7 =-0.051848879792_REAL64

  real(REAL64), parameter  ::  DFC2 = CC1 + 2*CC2 + 3*CC3 + 4*CC4 + 5*CC5 + 6*CC6 + 7*CC7

! PW92 parameters

  real(REAL64), parameter  ::  BETA1_PW = 7.5957_REAL64
  real(REAL64), parameter  ::  BETA2_PW = 3.5876_REAL64
  real(REAL64), parameter  ::  BETA3_PW = 1.6382_REAL64
  real(REAL64), parameter  ::  BETA4_PW = 0.49294_REAL64
!  real(REAL64), parameter  ::  A_PW = 0.062182_REAL64                   !  original value
  real(REAL64), parameter  ::  A_PW = 0.0621814_REAL64                   !  modified value from libxc
  real(REAL64), parameter  ::  ALP_PW = 0.21370_REAL64


! Avoids calculating irrelevant exponentials compared to unity

  real(REAL64), parameter  ::  ARGMAX = 50.0_REAL64

! Avoids division by zero

  real(REAL64), parameter  ::  DBL_MIN = 10*tiny(UM), DBL_EPSILON = 4*epsilon(UM)


  rs = (3/(4*PI*rho))**(UM/3)
  d_rs_dr = - rs / (3*rho)

  grloc = grho
  if(grho < ZERO) grloc = ZERO

! Perdew Wang
!  call xc_lda_c_pw92(rho, eclda, vclda)
!  d_eclda_dr = (vclda-eclda) / rho

  b_pw = BETA1_pw*sqrt(rs) +   BETA2_pw*rs +  BETA3_pw*rs*sqrt(rs) + BETA4_pw*rs*rs
  db_pwdrs = BETA1_pw / (2*sqrt(rs)) + BETA2_pw + (3*BETA3_pw*sqrt(rs)) / 2    &
                     + 2*BETA4_PW*rs
  d2b_pwdrs2 = -BETA1_pw / (4*sqrt(rs)*rs) + 3*BETA3_pw / (4*sqrt(rs)) + 2*BETA4_PW

  c_pw = UM + UM / (A_PW*b_pw)
  dc_pwdrs = - (c_pw - UM)*db_pwdrs / b_pw
  d2c_pwdrs2 = -dc_pwdrs*db_pwdrs / b_pw - (c_pw - UM)*d2b_pwdrs2 / b_pw       &
      + (c_pw - UM)*db_pwdrs*db_pwdrs / (b_pw*b_pw)

  eclda = - A_PW*(UM + ALP_PW*rs)*log(c_pw)
  d_eclda_drs = -A_PW*ALP_PW*log(c_pw) - A_PW*(UM + ALP_PW*rs)*(UM /c_pw)*dc_pwdrs
  d2_eclda_drs2 = -2*A_PW*ALP_PW*(UM /c_pw)*dc_pwdrs +                    &
                   A_PW*(UM + ALP_PW*rs)*(UM / (c_pw*c_pw))*dc_pwdrs*dc_pwdrs -   &
                   A_PW*(UM + ALP_PW*rs)*(UM /c_pw)*d2c_pwdrs2
  d_eclda_dr = d_eclda_drs*d_rs_dr

  xkf = AKF / rs                                                         !  (3*PI*PI * rho)**(UM/3)
  d_xkf_dr = (UM/3) * xkf / rho                                          !  (UM/3)*(3*PI*PI * rho)**(UM/3) / rho

  s = grloc / (2 * xkf * rho)
  d_s_dr = -s/rho - (s/xkf)*d_xkf_dr
  d_s_dgr = UM / (2 * xkf * rho)

  tausingle = grloc*grloc / (8*rho)
  d_tausingle_dr = -tausingle/rho
  d_tausingle_dgr = 2* grloc / (8*rho)
  tauunif = (3*UM / 10) * xkf*xkf * rho
  d_tauunif_dr = (3*UM / 10) * xkf*xkf + 2*(3*UM / 10) * xkf*d_xkf_dr * rho
  alpha = (tau - tausingle) / (tauunif + ETA*tausingle)
  d_alpha_dr = -d_tausingle_dr * ( UM / (tauunif + ETA*tausingle)        &
               + ETA*alpha / (tauunif + ETA*tausingle) )                 &
               - alpha*d_tauunif_dr/ (tauunif + ETA*tausingle)
  d_alpha_dgr = -d_tausingle_dgr * ( UM / (tauunif + ETA*tausingle)      &
               + ETA*alpha / (tauunif + ETA*tausingle) )
  d_alpha_dtau = UM / (tauunif + ETA*tausingle)


  tt = CTP*s*s/rs
  d_tt_drs = -CTP*s*s/(rs*rs)
  d_tt_ds = 2*CTP*s/rs
  betrs = BETAMB*(UM + BETAB*rs) / (UM + BETAC*rs)
  d_betrs_drs = BETAMB*BETAB / (UM + BETAC*rs) -                         &
                BETAMB*BETAC*(UM + BETAB*rs) / ((UM + BETAC*rs)*(UM + BETAC*rs))

  w1 = exp(-eclda/GAM) - UM
  d_w1_drs = -exp(-eclda/GAM)*(d_eclda_dr / d_rs_dr) / GAM

  eclda0 = -B1C / (UM + B2C*sqrt(rs) + B3C*rs)
  d_eclda0_drs = eclda0*eclda0*(B3C + B2C / (2*sqrt(rs))) / B1C
  d2_eclda0_drs2 = 2*eclda0*eclda0*eclda0*(B3C + B2C / (2*sqrt(rs))) *    &
                      (B3C + B2C / (2*sqrt(rs))) / (B1C*B1C) -           &
                   eclda0*eclda0*(B2C / (4*rs*sqrt(rs))) / B1C

  yfunc = betrs*tt / (GAM*w1)
  d_yfunc_drs = d_betrs_drs*tt / (GAM*w1) - yfunc*d_w1_drs / w1 +        &
                 betrs*d_tt_drs / (GAM*w1)
  d_yfunc_ds = betrs*d_tt_ds / (GAM*w1)

  arg = -s*s*s*s / (DP2*DP2*DP2*DP2)
  d_arg_ds = -4*s*s*s / (DP2*DP2*DP2*DP2)

  arg1 = dfc2 / (27*GAM*w1)
  d_arg1_drs = -arg1*d_w1_drs / w1
  arg2 = 20*rs*(d_eclda0_drs - d_eclda_drs) - 45*ETA*(eclda0 - eclda)
  d_arg2_drs = 20*(d_eclda0_drs - d_eclda_drs) +                         &
               20*rs*(d2_eclda0_drs2 - d2_eclda_drs2) -                  &
               45*ETA*(d_eclda0_drs - d_eclda_drs)
  arg3 = s*s*exp(arg)
  d_arg3_ds = 2*s*exp(arg) + arg3*d_arg_ds
  dely = arg1*arg2*arg3
  d_dely_drs = d_arg1_drs*arg2*arg3 + arg1*d_arg2_drs*arg3
  d_dely_ds = arg1*arg2*d_arg3_ds

  gg = UM / sqrt(sqrt(UM+4*(yfunc - dely)))
  d_gg_drs = -gg * (d_yfunc_drs - d_dely_drs) / (UM+4*(yfunc - dely))
  d_gg_ds = -gg * (d_yfunc_ds - d_dely_ds) / (UM+4*(yfunc - dely))

  h1 = GAM*log( UM + w1*(UM - gg) )
  d_h1_drs = GAM*( d_w1_drs*(UM-gg) - w1*d_gg_drs ) / ( UM + w1*(UM-gg) )
  d_h1_ds = -GAM*w1*d_gg_ds / ( UM + w1*(UM-gg) )

  ec1 = eclda + h1
  d_ec1_dr = d_eclda_dr + d_h1_drs*d_rs_dr + d_h1_ds*d_s_dr
  d_ec1_dgr = d_h1_ds*d_s_dgr

  w0 = exp(-eclda0/B1C) - UM
  d_w0_drs = -exp(-eclda0/B1C)*d_eclda0_drs/B1C

  arg = UM + 4*CHI0*s*s
  d_arg_ds = 8*CHI0*s
  h0 = B1C * log( UM + w0*( UM - UM / (sqrt(sqrt(arg))) ) )
  d_h0_drs = B1C*d_w0_drs*( UM - UM / (sqrt(sqrt(arg))))         /             &
              ( UM + w0*( UM - UM / (sqrt(sqrt(arg))) ) )
  d_h0_ds = B1C*(w0/4)*(UM / (sqrt(sqrt(arg))))*(UM/arg)*d_arg_ds   /          &
              ( UM + w0*( UM - UM / (sqrt(sqrt(arg))) ) )

  ec0 = eclda0 + h0

  d_ec0_dr = d_eclda0_drs*d_rs_dr + d_h0_drs*d_rs_dr + d_h0_ds*d_s_dr
  d_ec0_dgr = d_h0_ds*d_s_dgr

  if(alpha < ZERO) then
    fc = exp(-c1c*alpha / (UM-alpha))
    d_fc_dalpha = -c1c*fc / ((UM-alpha)*(UM-alpha))
  elseif(alpha > 2.5) then
    fc = -DC*exp(c2c*alpha / (UM-alpha))
    d_fc_dalpha = c2c*fc / ((UM-alpha)*(UM-alpha))
  else
    fc = CC0 + alpha*(CC1 + alpha*(CC2 + alpha*(CC3 + alpha*(CC4 +    &
              alpha*(CC5 + alpha*(CC6+alpha*CC7))))))
    d_fc_dalpha = CC1 + alpha*(2*CC2 + alpha*(3*CC3 + alpha*(4*CC4 +  &
              alpha*(5*CC5 + alpha*(6*CC6+alpha*7*CC7)))))
  endif

  d_fc_dr = d_fc_dalpha*d_alpha_dr
  d_fc_dgr = d_fc_dalpha*d_alpha_dgr
  d_fc_dtau = d_fc_dalpha*d_alpha_dtau

  epsc = ec0 + (UM-fc)*(ec1-ec0)

  decdr = epsc + rho*( d_ec0_dr - d_fc_dr*(ec1-ec0)                 &
                    + (UM-fc)*(d_ec1_dr-d_ec0_dr) )
  decdgr = rho * ( d_ec0_dgr - d_fc_dgr*(ec1-ec0)                   &
                 + (UM-fc)*d_ec1_dgr - (UM-fc)*d_ec0_dgr )
  decdtau =-d_fc_dtau*(ec1-ec0)*rho

  return

end subroutine xc_mgga_c_r2scan
