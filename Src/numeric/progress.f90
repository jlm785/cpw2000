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

subroutine progress(j,n)
!  a non-advancing status counter...
!  input n = number of steps
!  input j = current step 
!  note that achar(13) brings the cursor to the begining of line
  implicit none
  integer :: j, n
  write(6,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
  & " Percent Complete: ", (real(j)/real(n))*100.0, "%"
end subroutine progress
