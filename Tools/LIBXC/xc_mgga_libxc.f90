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
!>  in the Tran-Blaha meta-GGA
!>
!>  Version that interfaces with libxc
!>  needs module from libxcf03.f90
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.05
!>  \date         30 September 2022.
!>  \copyright    GNU Public License v2

subroutine xc_mgga(author, rho, grho, lap_rho, twotau,                   &
                        epsx, epsc, vx, vc, ctb09 )

!  Written 30 September 2022. JLM
!  Split from Carlos Loia Reis previous code.

  use xc_f03_lib_m

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author                          !<  type of correlation wanted.
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of charge density
  real(REAL64), intent(in)           ::  lap_rho                         !<  Laplacian of charge density
  real(REAL64), intent(in)           ::  twotau                          !<  Twice the kinetic energy density
  real(REAL64), intent(in)           ::  ctb09                           !<  TB constant

! output

  real(REAL64), intent(out)          ::  epsx                            !<  ZERO exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  epsc                            !<  ZERO correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vx                              !<  Tran-Blaha exchange potential (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vc                              !<  LDA correlation potential (hartree/bohr^3)

! libxc interface stuff

  TYPE(xc_f03_func_t)        ::  xc_func
!  TYPE(xc_f03_func_info_t)   ::  xc_info
  integer(8), parameter      ::  np = 1
  integer                    ::  func_id
  real(REAL64)               ::  rhotmp(1)
  real(REAL64)               ::  sigtmp(1)
  real(REAL64)               ::  tautmp(1)
  real(REAL64)               ::  laptmp(1)
!  real(REAL64)               ::  etmp(1)
  real(REAL64)               ::  vtmp(1)
  real(REAL64)               ::  dedgtmp(1)
  real(REAL64)               ::  vtautmp(1)
  real(REAL64)               ::  vlaptmp(1)
  real(REAL64)               :: set_ctb09(1)

! parameters

  real(REAL64), parameter :: ZERO = 0.0_REAL64


  rhotmp(1) = rho
  sigtmp(1) = grho*grho
  laptmp(1) = lap_rho
  tautmp(1) = twotau/2

  set_ctb09(1) = ctb09

  call xc_lda(author, rho, epsx, epsc, vx, vc )

  func_id = 208                             ! tb09
!  func_id = 206                         ! br89

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
  if(func_id == 208) call xc_f03_func_set_ext_params(xc_func, set_ctb09)
  call xc_f03_mgga_vxc(xc_func, np, rhotmp, sigtmp, laptmp, tautmp, vtmp, dedgtmp, vtautmp, vlaptmp )
  call xc_f03_func_end(xc_func)

  vx = vtmp(1)

  epsx = ZERO
  epsc = ZERO

  return

end subroutine xc_mgga

