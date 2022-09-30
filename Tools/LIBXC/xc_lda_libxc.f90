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
!>  \version      5.05
!>  \date         29 September 2022.
!>  \copyright    GNU Public License v2


subroutine xc_lda( author, rho, epsx, epsc, vx, vc )

!  Written 29 September 2022. JLM
!  Split from Carlos Loia Reis previous code.

  use xc_f03_lib_m

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author                          !<  type of correlation wanted.
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)

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
  real(REAL64)               ::  etmp(1)
  real(REAL64)               ::  vtmp(1)

! E x c h a n g e

  func_id = 1
  rhotmp(1) = rho

   call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
     call xc_f03_lda_exc_vxc(xc_func, np, rhotmp, etmp, vtmp)
   call xc_f03_func_end(xc_func)

   epsx = etmp(1)
   vx = vtmp(1)

! C o r r e l a t i o n

  if ( author == 'ca' .or. author == 'CA' .or.                           &
       author == 'pz' .or. author == 'PZ') then

  func_id = 9

  elseif ( author == 'vwn' .or. author == 'VWN' ) then

!   V o s k o     W i l k     and     N u s a i r

  func_id = 7

  elseif ( author == 'wi' .or. author == 'WI') then

!   W i g n e r

  func_id = 2

  elseif ( author == 'pw92' .or. author == 'PW92' ) then

!   P e r d e w     &     W a n g     1 9 9 2

  func_id = 12

  else
    write(6,'("  STOPPED  in xc_lda (libxc):   unknown correlation")')
    write(6,*) '     ',author

    stop

  endif

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_lda_exc_vxc(xc_func, np, rhotmp, etmp, vtmp)
    epsc = etmp(1)
    vc = vtmp(1)
  call xc_f03_func_end(xc_func)

  return

end subroutine xc_lda
