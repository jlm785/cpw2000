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

!>  Number of states to use so that a degenerate level is not cut.
!>  It is the smallest number not smaller than neig_min that is
!>  followed by a gap, among the neig_max calculated states.  If there
!>  is none, it is the largest number smaller than neig_min followed
!>  by a gap.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         23 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_int_neig_level(neig_min, neig_max, ei, nuse)

! Written 23 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  neig_min                        !<  minimum number of states wanted
  integer, intent(in)                ::  neig_max                        !<  number of calculated states
  real(REAL64), intent(in)           ::  ei(neig_max)                    !<  eigenvalues in increasing order (hartree)

! output

  integer, intent(out)               ::  nuse                            !<  number of states to use

! parameters

  real(REAL64), parameter     ::  TOL = 1.0E-5_REAL64                    !  tolerance for degenerate energies (hartree)

! counters

  integer   ::  n


  nuse = 0

  do n = neig_min,neig_max-1
    if(ei(n+1) - ei(n) > TOL) then
      nuse = n

      exit

    endif
  enddo

  if(nuse == 0) then
    do n = min(neig_min,neig_max)-1,1,-1
      if(ei(n+1) - ei(n) > TOL) then
        nuse = n

        exit

      endif
    enddo
  endif

  if(nuse == 0) nuse = min(neig_min,neig_max)

  return

end subroutine ao_int_neig_level
