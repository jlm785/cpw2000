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
!>
!>  Version that interfaces with libxc
!>  needs module from libxcf03.f90
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.13
!>  \date         29 September 2022. 10 November 2026.
!>  \copyright    GNU Public License v2


subroutine xc_lda( author, rho, epsx, epsc, vx, vc )

!  Written 29 September 2022. JLM
!  Split from Carlos Loia Reis previous code.
!  chrsameinfo, rho < EPS as in Src, PW92 is LDA_C_PW_MOD (13). 10 November 2026. JLM+claude

  use xc_f03_lib_m

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author                          !<  type of correlation wanted (CA = PZ, VWN, WI, PW92)
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)

! output

  real(REAL64), intent(out)          ::  epsx                            !<  exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  epsc                            !<  correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vx                              !<  exchange potential (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vc                              !<  correlation potential (hartree/bohr^3)

! libxc interface stuff

  type(xc_f03_func_t)        ::  xc_func
  integer(8), parameter      ::  NP = 1
  integer                    ::  func_id
  real(REAL64)               ::  rhotmp(1)
  real(REAL64)               ::  etmp(1)
  real(REAL64)               ::  vtmp(1)

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks

! parameters

  real(REAL64), parameter  :: ZERO = 0.0_REAL64
  real(REAL64), parameter  :: EPS = 1.0E-24_REAL64


  if(rho < EPS) then

    epsx = ZERO
    epsc = ZERO
    vx = ZERO
    vc = ZERO

  else

    if ( chrsameinfo(author, 'CA' ) .or.  chrsameinfo(author, 'PZ' ) ) then

!     P e r d e w    a n d     Z u n g e r

      func_id = 9                                                        !  LDA_C_PZ

    elseif ( chrsameinfo(author, 'VWN' ) ) then

!     V o s k o     W i l k     and     N u s a i r

      func_id = 7                                                        !  LDA_C_VWN

    elseif ( chrsameinfo(author, 'WI' ) ) then

!     W i g n e r

      func_id = 2                                                        !  LDA_C_WIGNER

    elseif ( chrsameinfo(author, 'PW92' ) ) then

!     P e r d e w     &     W a n g     1 9 9 2,   same constants as xc_lda_c_pw92 in Src

      func_id = 13                                                       !  LDA_C_PW_MOD

    else

      write(6,'("  STOPPED  in xc_lda (libxc):   unknown correlation")')
      write(6,*) '     ',author

      stop

    endif

    rhotmp(1) = rho

!   E x c h a n g e

    call xc_f03_func_init(xc_func, 1, XC_UNPOLARIZED)                    !  LDA_X
    call xc_f03_lda_exc_vxc(xc_func, NP, rhotmp, etmp, vtmp)
    call xc_f03_func_end(xc_func)

    epsx = etmp(1)
    vx = vtmp(1)

!   C o r r e l a t i o n

    call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_lda_exc_vxc(xc_func, NP, rhotmp, etmp, vtmp)
    call xc_f03_func_end(xc_func)

    epsc = etmp(1)
    vc = vtmp(1)

  endif

  return

end subroutine xc_lda
