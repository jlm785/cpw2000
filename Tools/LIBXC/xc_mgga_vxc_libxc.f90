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

!>  Computes the exchange correlation potential for the
!>  meta-GGA that only have a potential (Tran-Blaha, Becke-Roussel).
!>  The energies are those of the LDA.
!>
!>  Version that interfaces with libxc
!>  needs module from libxcf03.f90
!>
!>  The Becke-Roussel hole potential is obtained from the
!>  Tran-Blaha potential with c = 2/3,  v_BR = (3/2) v_TB09(c=2/3).
!>
!>  By default libxc enforces sigma <= 8 rho tau (Fermi hole curvature),
!>  which changes the results where tau < tau_W (alpha < 0).  This
!>  happens in numerical noise, and systematically for the deorbitalized
!>  functionals with the Thomas-Fermi tau (R2TF, LATF, TATF), where
!>  alpha = 1 - 5 s^2 / 3.  The native subroutines do not enforce it.
!>  With LFHC = .FALSE. the results are the same as the native ones,
!>  except where tau <= 0 (libxc uses its tau threshold).
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.13
!>  \date         30 September 2022. 10 November 2026.
!>  \copyright    GNU Public License v2

subroutine xc_mgga_vxc(author_x, author_c, rho, grho, lap_rho, tau,      &
                     epsx, epsc, vx, vc, ctb09 )

!  Written 30 September 2022. JLM
!  Split from Carlos Loia Reis previous code.
!  Was xc_mgga, interface of xc_mgga_vxc in Src (tau instead of twotau),
!  both libxc parameters of TB09, BR. 10 November 2026. JLM+claude

  use xc_f03_lib_m

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author_x                        !<  type of exchange wanted (TB09 = TBL, BR)
  character(len = *), intent(in)     ::  author_c                        !<  type of correlation wanted (LDA)
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of charge density
  real(REAL64), intent(in)           ::  lap_rho                         !<  Laplacian of charge density
  real(REAL64), intent(in)           ::  tau                             !<  The kinetic energy density
  real(REAL64), intent(in)           ::  ctb09                           !<  TB constant

! output

  real(REAL64), intent(out)          ::  epsx                            !<  LDA exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  epsc                            !<  LDA correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vx                              !<  Tran-Blaha exchange potential (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vc                              !<  LDA correlation potential (hartree/bohr^3)

! local variables

  real(REAL64)               ::  cpar                                    !  parameter c of libxc TB09
  real(REAL64)               ::  fac                                     !  multiplying factor

! libxc interface stuff

  type(xc_f03_func_t)        ::  xc_func
  integer(8), parameter      ::  NP = 1
  real(REAL64)               ::  rhotmp(1)
  real(REAL64)               ::  sigtmp(1)
  real(REAL64)               ::  tautmp(1)
  real(REAL64)               ::  laptmp(1)
  real(REAL64)               ::  vtmp(1)
  real(REAL64)               ::  dedgtmp(1)
  real(REAL64)               ::  vtautmp(1)
  real(REAL64)               ::  vlaptmp(1)

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks

! parameters

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  :: EPS = 1.0E-24_REAL64

  logical, parameter       :: LFHC = .FALSE.                             !  .TRUE. is the libxc default


  epsx = ZERO
  epsc = ZERO

  if(rho < EPS) then

    vx = ZERO
    vc = ZERO

  else

    if ( chrsameinfo(author_x, 'TBL' ) .or. chrsameinfo(author_x, 'TB09' ) ) then

      cpar = ctb09
      fac = UM

    elseif ( chrsameinfo(author_x, 'BR' )  ) then

      cpar = (2*UM) / (3*UM)
      fac = (3*UM) / (2*UM)

    else

      write(6,'("  STOPPED  in xc_mgga_vxc (libxc):   unknown exchange-correlation")')
      write(6,*) '     ',author_x,'  ', author_c

      stop

    endif

    call xc_lda(author_c, rho, epsx, epsc, vx, vc )

    rhotmp(1) = rho
    sigtmp(1) = grho*grho
    laptmp(1) = lap_rho
    tautmp(1) = tau

    call xc_f03_func_init(xc_func, 208, XC_UNPOLARIZED)                  !  MGGA_X_TB09
    if(.not. LFHC) call xc_f03_func_set_fhc_enforcement(xc_func, .FALSE.)
    call xc_f03_func_set_ext_params_name(xc_func, 'c', cpar)
    call xc_f03_func_set_ext_params_name(xc_func, '_alpha', ZERO)
    call xc_f03_mgga_vxc(xc_func, NP, rhotmp, sigtmp, laptmp, tautmp,     &
                         vtmp, dedgtmp, vlaptmp, vtautmp )
    call xc_f03_func_end(xc_func)

    vx = fac*vtmp(1)

  endif

  return

end subroutine xc_mgga_vxc
