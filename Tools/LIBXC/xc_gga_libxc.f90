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
!>  PBE uses the libxc constants (GGA_X_PBE and GGA_C_PBE),
!>  the native version uses those of the original paper.
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.13
!>  \date         29 September 2022. 10 November 2026.
!>  \copyright    GNU Public License v2

subroutine xc_gga( author, rho, grho,                                    &
                   epsx, epsc, dexdr, decdr, dexdgr, decdgr )


!  Written 29 September 2022. JLM
!  Split from Carlos Loia Reis previous code.
! Jorge Pacheco, Jose Luis Martins
!  chrsameinfo, rho < EPS and grho < 0 as in Src. 10 November 2026. JLM+claude

  use xc_f03_lib_m

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author                          !<  type of correlation wanted (PBE)
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdr                           !<  derivative with respect to rho of rho*epsx (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdr                           !<  derivative with respect to rho of rho*epsc (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdgr                          !<  derivative with respect to grho of rho*epsx (hartree/bohr^2)
  real(REAL64), intent(out)          ::  decdgr                          !<  derivative with respect to grho of rho*epsc (hartree/bohr^2)

! local variables

  real(REAL64)               ::  grloc                                   !  local value of grho

! libxc interface stuff

  type(xc_f03_func_t)        ::  xc_func
  integer(8), parameter      ::  NP = 1
  integer                    ::  idx, idc
  real(REAL64)               ::  rhotmp(1)
  real(REAL64)               ::  sigtmp(1)
  real(REAL64)               ::  etmp(1)
  real(REAL64)               ::  vtmp(1)
  real(REAL64)               ::  dedsigtmp(1)

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

  else

    if ( chrsameinfo(author, 'PBE' ) ) then

!     Perdew, Burke and Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)

      idx = 101                                                          !  GGA_X_PBE
      idc = 130                                                          !  GGA_C_PBE

    else

      write(6,'("   STOPPED in xc_gga (libxc):   unknown correlation")')
      write(6,*) '     ',author

      stop

    endif

    grloc = grho
    if(grho < ZERO) grloc = ZERO

    rhotmp(1) = rho
    sigtmp(1) = grloc*grloc

!   C o r r e l a t i o n

    call xc_f03_func_init(xc_func, idc, XC_UNPOLARIZED)
    call xc_f03_gga_exc_vxc(xc_func, NP, rhotmp, sigtmp, etmp, vtmp, dedsigtmp)
    call xc_f03_func_end(xc_func)

    epsc = etmp(1)
    decdr = vtmp(1)
    decdgr = 2*dedsigtmp(1)*grloc

!   E x c h a n g e

    call xc_f03_func_init(xc_func, idx, XC_UNPOLARIZED)
    call xc_f03_gga_exc_vxc(xc_func, NP, rhotmp, sigtmp, etmp, vtmp, dedsigtmp)
    call xc_f03_func_end(xc_func)

    epsx = etmp(1)
    dexdr = vtmp(1)
    dexdgr = 2*dedsigtmp(1)*grloc

  endif

  return

end subroutine xc_gga
