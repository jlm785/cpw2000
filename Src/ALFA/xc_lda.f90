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
!>  in the local density approximation .
!>  Lengths in Bohr, energies in Hartrees.
!>  Adapted from a package by L.C. Balbas and J.M. Soler, Dec'96. Version 0.5.
!>
!>  \author       J. L. Martins, L.C. Balbas, J.M. Soler, J.M. Pacheco
!>  \version      5.12
!>  \date         23 February 1999. 13 November 2025.
!>  \copyright    GNU Public License v2

subroutine xc_lda(author, rho, epsx, epsc, vx, vc)

! Written 23 february 1999. jlm
! Modified by J.L. Martins and J.M. Pacheco.
! Separated in subroutines, 10 November 2025. JLM
! New parameters to reproduce libxc resuts. Bug in Wigner. 13 November 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author                          !<  type of xc wanted (ca=pz , pw92 , vwn, wi)
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vx                              !<  exchange potential (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vc                              !<  correlation potential (hartree/bohr^3)

! parameters

  real(REAL64), parameter  :: ZERO = 0.0_REAL64
  real(REAL64), parameter  :: EPS = 1.0E-24_REAL64


! initial stuff

  if(rho < EPS) then

    epsx = ZERO
    epsc = ZERO
    vx = ZERO
    vc = ZERO

  else

!   E x c h a n g e

    call xc_lda_x(rho, epsx, vx)

!   C o r r e l a t i o n

    if ( author == 'ca' .or. author == 'CA' .or.                         &
         author == 'pz' .or. author == 'PZ') then

!     P e r d e w    a n d     Z u n g e r

      call xc_lda_c_pz(rho, epsc, vc)


    elseif ( author == 'vwn' .or. author == 'VWN' ) then

!     V o s k o     W i l k     and     N u s a i r

      call xc_lda_c_vwn(rho, epsc, vc)


    elseif ( author == 'wi' .or. author == 'WI') then

!     W i g n e r

      call xc_lda_c_wigner(rho, epsc, vc)

    elseif ( author == 'pw92' .or. author == 'PW92' ) then

!     P e r d e w     &     W a n g     1 9 9 2

      call xc_lda_c_pw92(rho, epsc, vc)

    else

      write(6,'("  STOPPED  in xc_lda:   unknown correlation")')
      write(6,*) '     ',author

      stop

    endif

  endif

  return

end subroutine xc_lda


!>  Computes the local density exchange energy and potential

subroutine xc_lda_x(rho, epsx, vx)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vx                              !<  exchange potential (hartree/bohr^3)

! local variables

  real(REAL64)     :: rs                                                 !  Wigner-Seitz radius
  real(REAL64)     :: a0                                                 !  exchange constant

! constants

  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

  rs = (3/(4*PI*rho))**(UM/3)

  a0 = (4 / (9*PI))**(UM/3)
  vx = -UM / (PI*a0*rs)
  epsx = ( (3*UM) / 4) * vx

  return

end subroutine xc_lda_x


!>  Perdew-Zunger parameterization of Ceperley-Alder
!>  correlation. Ref: Perdew and Zunger, Phys. Rev. B 23 5075 (1981). jlm

subroutine xc_lda_c_pz(rho, epsc, vc)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vc                              !<  correlation potential (hartree/bohr^3)

! local variables

  real(REAL64)     :: rs                                                 !  Wigner-Seitz radius
  real(REAL64)     :: sqrs, rslog, rsdecdrs


! parameters

  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

! PZ parameters

  real(REAL64)             ::  be
  real(REAL64), parameter  ::  A1_PZ = 1.0529_REAL64
  real(REAL64), parameter  ::  A2_PZ = 0.3334_REAL64
  real(REAL64), parameter  ::  B_PZ = 0.1423_REAL64
  real(REAL64), parameter  ::  C1_PZ = 0.0311_REAL64
  real(REAL64), parameter  ::  C2_PZ = 0.002_REAL64
  real(REAL64), parameter  ::  C3_PZ = 0.048_REAL64
  real(REAL64), parameter  ::  C4_PZ = 0.0116_REAL64


  rs = (3/(4*PI*rho))**(UM/3)

  if (rs > UM) then
    sqrs = sqrt(rs)
    be = UM + A1_PZ*sqrs + A2_PZ*rs
    epsc = -B_PZ / be
    rsdecdrs = B_PZ * (A1_PZ*sqrs/2 + A2_PZ*rs) / (be*be)
    vc = epsc - rsdecdrs/3
  else
    rslog = log(rs)
    epsc = (C1_PZ + C2_PZ*rs)*rslog - C3_PZ - C4_PZ*rs
    rsdecdrs = (C1_PZ + C2_PZ*rs) + C2_PZ*rslog*rs - C4_PZ*rs
    vc = epsc - rsdecdrs/3
  endif

  return

end subroutine xc_lda_c_pz


!>  Vosko-Wilk-Nusair parameterization of Ceperley-Alder
!>  (this paper is for LDA and has the right formula, eq. b.11)
!>  Wilk and Vosko, J. Phys. C15 (82) 2139. jmp

subroutine xc_lda_c_vwn(rho, epsc, vc)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vc                              !<  correlation potential (hartree/bohr^3)

! local variables

  real(REAL64)     ::  rs                                                 !  Wigner-Seitz radius
  real(REAL64)     ::  sqrs, rsdecdrs

  real(REAL64)     ::  xd_vwn
  real(REAL64)     ::  arg1, arg2, arg3

! parameters

  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: UM = 1.0_REAL64

! VWN parameters

  real(REAL64), parameter  ::  X_VWN = 0.10498_REAL64                    !  Sign changed with respect to paper
  real(REAL64), parameter  ::  A_VWN = 0.0310907_REAL64
  real(REAL64), parameter  ::  B_VWN = 3.72744_REAL64
  real(REAL64), parameter  ::  C_VWN = 12.9352_REAL64

  real(REAL64), parameter  ::  QQ  = sqrt(4*C_VWN - B_VWN * B_VWN)
  real(REAL64), parameter  ::  C1 = 2*B_VWN / QQ
  real(REAL64), parameter  ::  C2 =-B_VWN * X_VWN / (X_VWN*X_VWN - B_VWN*X_VWN + C_VWN)
  real(REAL64), parameter  ::  C3 = 2*(B_VWN - 2*X_VWN) / QQ
  real(REAL64), parameter  ::  B1 = (B_VWN * X_VWN + C_VWN) / (X_VWN * C_VWN)
  real(REAL64), parameter  ::  B2 = (X_VWN + B_VWN) / (X_VWN * C_VWN)
  real(REAL64), parameter  ::  B3 = UM / (X_VWN * C_VWN)


  rs = (3/(4*PI*rho))**(UM/3)

  sqrs  = sqrt(rs)
  xd_vwn = rs + B_VWN*sqrs + C_VWN

  arg1 = rs / xd_vwn
  arg2 = QQ / (2*sqrs + B_VWN)
  arg3 = (sqrs + X_VWN) * (sqrs + X_VWN) / xd_vwn

  epsc = log(arg1) + (C1 - C2 * C3) * atan(arg2) - C2 * log(arg3)
  epsc = A_VWN * epsc

  rsdecdrs = A_VWN * (UM + B1*sqrs) / (UM + B1*sqrs + B2*rs + B3*rs*sqrs)
  vc = epsc - rsdecdrs / 3

  return

end subroutine xc_lda_c_vwn


!>  Wigner correlation energy
!>  Phys. Rev. 46, 1002 (1934)
!>  Trans. Faraday Soc. 34, 678 (1938)

subroutine xc_lda_c_wigner(rho, epsc, vc)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vc                              !<  correlation potential (hartree/bohr^3)

! local variables

  real(REAL64)     :: rs                                                 !  Wigner-Seitz radius
  real(REAL64)     :: rsdecdrs

! parameters

  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: UM = 1.0_REAL64

! Wigner parameters

!  real(REAL64), parameter  ::  A_WI = -0.4377645_REAL64                 !  My old value
!  real(REAL64), parameter  ::  A_WI = -0.29_REAL64                      !  Value of 1938 paper of Wigner
  real(REAL64), parameter  ::  A_WI = -0.44_REAL64                       !  Used in libxc
!  real(REAL64), parameter  ::  B_WI = 5.1_REAL64                        !  Value of 1938 paper of Wigner
  real(REAL64), parameter  ::  B_WI = 7.8_REAL64


  rs = (3/(4*PI*rho))**(UM/3)

  epsc = A_WI / (rs + B_WI)
  rsdecdrs = -rs*A_WI / ((rs + B_WI)*(rs + B_WI))
  vc = epsc - rsdecdrs / 3

  return

end subroutine xc_lda_c_wigner

!>  Implements the perdew-wang'92 local correlation (beyond rpa).
!>  J.P. Perdew and Y. Wang, PRB 45, 13244 (1992) using eq.(10). lcb, jms

subroutine xc_lda_c_pw92(rho, epsc, vc)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vc                              !<  correlation potential (hartree/bohr^3)

! local variables

  real(REAL64)     :: rs                                                 !  Wigner-Seitz radius
  real(REAL64)     :: sqrs
  real(REAL64)     :: rsdecdrs

! parameters

  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: UM = 1.0_REAL64

! PW92 parameters

  real(REAL64)             ::  b_pw, c_pw, dc_pwdrs, db_pwdrs
  real(REAL64), parameter  ::  BETA1_PW = 7.5957_REAL64
  real(REAL64), parameter  ::  BETA2_PW = 3.5876_REAL64
  real(REAL64), parameter  ::  BETA3_PW = 1.6382_REAL64
  real(REAL64), parameter  ::  BETA4_PW = 0.49294_REAL64
!  real(REAL64), parameter  ::  A_PW = 0.062182_REAL64                   !  original value
  real(REAL64), parameter  ::  A_PW = 0.0621814_REAL64                   !  modified value from libxc
  real(REAL64), parameter  ::  ALP_PW = 0.21370_REAL64


  rs = (3/(4*PI*rho))**(UM/3)

  sqrs = sqrt(rs)

  b_pw = BETA1_pw*sqrs +   BETA2_pw*rs +  BETA3_pw*rs*sqrs + BETA4_pw*rs*rs
  db_pwdrs = BETA1_pw / (2*sqrs) + BETA2_pw + (3*BETA3_pw*sqrs) / 2 + 2*BETA4_PW*rs
  c_pw = UM + UM / (A_PW*b_pw)
  dc_pwdrs = - (c_pw - UM)*db_pwdrs / b_pw

  epsc = - A_PW * (UM + ALP_PW*rs) *  log(c_pw)
  rsdecdrs = -A_PW*rs * (ALP_PW * log(c_pw) + (UM + ALP_PW*rs) * dc_pwdrs / c_pw )
  vc = epsc - rsdecdrs / 3

  return

end subroutine xc_lda_c_pw92
