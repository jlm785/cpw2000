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

!>     Computes the exchange correlation energy and potential
!>     in the generalized gradient approximation (GGA).
!>     non-spin-polarized
!>     Lengths in Bohr, energies in Hartrees.
!>
!>  Version that interfaces with libxc
!>  needs module from libxcf03.f90
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.05
!>  \date         29 September 2022.
!>  \copyright    GNU Public License v2

subroutine xc_gga( author, rho, grho,                                    &
                   epsx, epsc, dexdr, decdr, dexdgr, decdgr )


!  Written 29 September 2022. JLM
!  Split from Carlos Loia Reis previous code.
! Jorge Pacheco, Jose Luis Martins

  use xc_f03_lib_m

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author                          !<  type of correlation wanted.
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  ZERO exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  epsc                            !<  ZERO correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdr                           !<  derivative with respect to rho of rho*epsx (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdr                           !<  derivative with respect to rho of rho*epsc (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdgr                          !<  derivative with respect to grho of rho*epsx (hartree/bohr^2)
  real(REAL64), intent(out)          ::  decdgr                          !<  derivative with respect to grho of rho*epsc (hartree/bohr^2)

! libxc interface stuff

  TYPE(xc_f03_func_t)        ::  xc_func
!  TYPE(xc_f03_func_info_t)   ::  xc_info
  integer(8), parameter      ::  np = 1
  integer                    ::  func_id
  real(REAL64)               ::  rhotmp(1)
  real(REAL64)               ::  sigtmp(1)
  real(REAL64)               ::  etmp(1)
  real(REAL64)               ::  vtmp(1)
  real(REAL64)               ::  dedsigtmp(1)

  rhotmp(1) = rho
  sigtmp(1) = grho*grho


  if (author == 'pbe' .or. author == 'PBE') then

    func_id = 130

    call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_gga_exc_vxc(xc_func, np, rhotmp, sigtmp, etmp, vtmp, dedsigtmp)
    call xc_f03_func_end(xc_func)

    epsc = etmp(1)
    decdr = vtmp(1)
    decdgr = 2*dedsigtmp(1)*grho

  else

    write(6,'("   STOPPED in xc_gga:   unknown correlation")')
    write(6,*) '     ',author

    stop

  endif

! e x c h a n g e

  func_id = 101

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
  call xc_f03_gga_exc_vxc(xc_func, np, rhotmp, sigtmp, etmp, vtmp, dedsigtmp)
  call xc_f03_func_end(xc_func)

  epsx = etmp(1)
  dexdr = vtmp(1)
  dexdgr = 2*dedsigtmp(1)*grho

  return

end subroutine xc_gga
