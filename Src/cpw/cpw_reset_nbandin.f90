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

!>  checks if the requested
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         25 June 2024.
!>  \copyright    GNU Public License v2

subroutine cpw_reset_nbandin(io6, ntype, zv, nbandin,                    &
         mxdtyp)

! Written 25 June 2024. JLM

  implicit none

  integer, parameter        ::  REAL64 = selected_real_kind(12)

! input variables

  integer, intent(in)                ::  mxdtyp                          !<  maximum number of type of atoms

  integer, intent(in)                ::  io6                             !<  default output tape
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  real(REAL64), intent(in)           ::  zv(mxdtyp)                      !<  valence of atom with type i

! input and output

  integer, intent(inout)             ::  nbandin                         !<  target for number of bands

! local variables
  real(REAL64)          ::  zsum                                         !  sum of ionic charges

! parameters

  real(REAL64), parameter  :: ZERO = 0.0_REAL64

! counters, etc...

  integer  :: n


  zsum = ZERO
  do n = 1,ntype
    zsum = zsum + zv(n)
  enddo
  if(2*nbandin < nint(zsum))  then
    write(io6,*)
    write(io6,*)
    write(io6,*)
    write(io6,*)   "  WARNING       WARNING       WARNING       WARNING"
    write(io6,*)
    write(io6,*)   "  The value of nbandin, ", nbandin, " is too small"
    write(io6,*)
    nbandin = nint(1.1*zsum/2) + 5
    write(io6,*)   "  nbandin increased to ", nbandin
    write(io6,*)
    write(io6,*)
    write(io6,*)
  endif

  return

end subroutine cpw_reset_nbandin
