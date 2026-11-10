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
!>  in the meta generalized gradient approximation (mGGA).
!>  non-spin-polarized
!>  Lengths in Bohr, energies in Hartrees.
!>
!>  Version that interfaces with libxc
!>  needs module from libxcf03.f90
!>
!>  By default libxc enforces sigma <= 8 rho tau (Fermi hole curvature),
!>  which changes the results where tau < tau_W (alpha < 0).  This
!>  happens in numerical noise, and systematically for the deorbitalized
!>  functionals with the Thomas-Fermi tau (R2TF, LATF, TATF), where
!>  alpha = 1 - 5 s^2 / 3.  The native subroutines do not enforce it.
!>  With LFHC = .FALSE. the results are the same as the native ones,
!>  except where tau <= 0 (libxc uses its tau threshold).
!>
!>  \author       claude
!>  \version      5.13
!>  \date         10 November 2026.
!>  \copyright    GNU Public License v2

subroutine xc_mgga( author, rho, grho, tau,                              &
                   epsx, epsc, dexdr, decdr, dexdgr, decdgr,             &
                   dexdtau, decdtau )

! Written 10 November 2026, same interface as xc_mgga in Src. claude

  use xc_f03_lib_m

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author                          !<  type of xc wanted (LAK, TASK, R2SC)
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of the electron density (1/bohr^4)
  real(REAL64), intent(in)           ::  tau                             !<  kinetic energy density (hartree/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdr                           !<  derivative with respect to rho of rho*epsx (hartree/bohr^3)
  real(REAL64), intent(out)          ::  decdr                           !<  derivative with respect to rho of rho*epsc (hartree/bohr^3)
  real(REAL64), intent(out)          ::  dexdgr                          !<  derivative with respect to grho of rho*epsx (hartree/bohr^2)
  real(REAL64), intent(out)          ::  decdgr                          !<  derivative with respect to grho of rho*epsc (hartree/bohr^2)
  real(REAL64), intent(out)          ::  dexdtau                         !<  derivative with respect to tau of rho*epsx (dimensionless)
  real(REAL64), intent(out)          ::  decdtau                         !<  derivative with respect to tau of rho*epsc (dimensionless)

! local variables

  real(REAL64)               ::  grloc                                   !  local value of grho
  logical                    ::  lcorlda                                 !  correlation is LDA

! libxc interface stuff

  type(xc_f03_func_t)        ::  xc_func
  integer(8), parameter      ::  NP = 1
  integer                    ::  idx, idc
  real(REAL64)               ::  rhotmp(1)
  real(REAL64)               ::  sigtmp(1)
  real(REAL64)               ::  laptmp(1)
  real(REAL64)               ::  tautmp(1)
  real(REAL64)               ::  etmp(1)
  real(REAL64)               ::  vtmp(1)
  real(REAL64)               ::  dedsigtmp(1)
  real(REAL64)               ::  vlaptmp(1)
  real(REAL64)               ::  vtautmp(1)

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks

! parameters

  real(REAL64), parameter  :: ZERO = 0.0_REAL64
  real(REAL64), parameter  :: EPS = 1.0E-24_REAL64

  logical, parameter       :: LFHC = .FALSE.                             !  .TRUE. is the libxc default


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

    lcorlda = .FALSE.

    if (chrsameinfo(author, 'LAK' ) ) then

!     T. Lebeda, T. Aschebrock, and S. Kummel,
!     Phys. Rev. Lett. 133, 136402 (2024)

      idx = 342                                                          !  MGGA_X_LAK
      idc = 345                                                          !  MGGA_C_LAK

    elseif (chrsameinfo(author, 'TASK' ) ) then

!     T. Aschebrock, and S. Kummel,
!     Phys.Rev.Research. 1, 033082 (2019), with PW92 correlation

      idx = 707                                                          !  MGGA_X_TASK
      idc = 13                                                           !  LDA_C_PW_MOD
      lcorlda = .TRUE.

    elseif (chrsameinfo(author, 'R2SC' ) ) then

!     Furness, Kaplan, Ning, Perdew and Sun,
!     J.Phys.Chem.Lett. 11, 8208 (2020)

      idx = 497                                                          !  MGGA_X_R2SCAN
      idc = 498                                                          !  MGGA_C_R2SCAN

    else

      write(6,'("   STOPPED in xc_mgga (libxc):   unknown correlation")')
      write(6,*) '     ',author

      stop

    endif

    grloc = grho
    if(grho < ZERO) grloc = ZERO

    rhotmp(1) = rho
    sigtmp(1) = grloc*grloc
    laptmp(1) = ZERO
    tautmp(1) = tau

!   E x c h a n g e

    call xc_f03_func_init(xc_func, idx, XC_UNPOLARIZED)
    if(.not. LFHC) call xc_f03_func_set_fhc_enforcement(xc_func, .FALSE.)
    call xc_f03_mgga_exc_vxc(xc_func, NP, rhotmp, sigtmp, laptmp, tautmp,    &
                             etmp, vtmp, dedsigtmp, vlaptmp, vtautmp)
    call xc_f03_func_end(xc_func)

    epsx = etmp(1)
    dexdr = vtmp(1)
    dexdgr = 2*dedsigtmp(1)*grloc
    dexdtau = vtautmp(1)

!   C o r r e l a t i o n

    call xc_f03_func_init(xc_func, idc, XC_UNPOLARIZED)
    if(lcorlda) then
      call xc_f03_lda_exc_vxc(xc_func, NP, rhotmp, etmp, vtmp)
      dedsigtmp(1) = ZERO
      vtautmp(1) = ZERO
    else
      if(.not. LFHC) call xc_f03_func_set_fhc_enforcement(xc_func, .FALSE.)
      call xc_f03_mgga_exc_vxc(xc_func, NP, rhotmp, sigtmp, laptmp, tautmp,  &
                               etmp, vtmp, dedsigtmp, vlaptmp, vtautmp)
    endif
    call xc_f03_func_end(xc_func)

    epsc = etmp(1)
    decdr = vtmp(1)
    decdgr = 2*dedsigtmp(1)*grloc
    decdtau = vtautmp(1)

  endif

  return

end subroutine xc_mgga
