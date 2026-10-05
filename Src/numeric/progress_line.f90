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

!>  Progress of a loop j = 1, ..., n written in a single line,
!>  with at most LMAX characters.  It writes a label when j = 1,
!>  the value of j every nstep steps, and ends the line when j = n.
!>  nstep is chosen from n so the line is never longer than LMAX.
!>
!>  Unlike progress, it does not use a carriage return,
!>  so it is suitable for output redirected to a file.
!>  Nothing else should be written to the same unit inside the loop.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         4 October 2026.
!>  \copyright    GNU Public License v2

subroutine progress_line(j, n, label, io)

! Written 4 October 2026. claude

  implicit none

! input

  integer, intent(in)                   ::  j                            !<  current step
  integer, intent(in)                   ::  n                            !<  number of steps
  character(len=*), intent(in)          ::  label                        !<  label at the beginning of the line
  integer, intent(in)                   ::  io                           !<  tape number

! local variables

  integer                   ::  iw                                       !  width of each number
  integer                   ::  nfit                                     !  number of values that fit in the line
  integer                   ::  nstep                                    !  values are written every nstep steps
  character(len=20)         ::  fmt

! parameters

  integer, parameter        ::  LMAX = 120                               !  maximum length of the line


  if(n < 1 .or. j < 1 .or. j > n) return

! width of the numbers, including a blank

  iw = 2
  if(n > 9) iw = int(log10(real(n))) + 2

! the first and last steps are always written, so two values less fit

  nfit = (LMAX - 2 - len(label)) / iw - 2
  nfit = max(nfit, 1)
  nstep = (n + nfit - 1) / nfit

  if(j == 1) then
    write(io,'(2x,a)',advance='no') label
  endif

  if(mod(j,nstep) == 0 .or. j == 1 .or. j == n) then
    write(fmt,'("(i",i0,")")') iw
    write(io,fmt,advance='no') j
  endif

  if(j == n) then
    write(io,*)
  endif

  return

end subroutine progress_line
