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

!>  Calculates the values in kmscr, describing the representation
!>  of quantities (initially the screening potential) in real space.
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         24 September 2026.
!>  \copyright    GNU Public License v2

subroutine size_kmscr(kmax, flgdal, mrgdual, idshift, kmscr)

! Written 24 September 2026 in a cleanup of existing code. JLM
! Added mrgdual, size_fft uses kmscr(1:3) (dual). 24 September 2026. JLM
! Checks of idshift and mrgdual. 24 September 2026. JLM+claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  character(len=4), intent(in)       ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  integer, intent(in)                ::  mrgdual                         !<  margin added to kmax/2 in the dual approximation (2 is the default, 4 is used for folding)
  integer, intent(in)                ::  idshift                         !<  shift of the fft mesh, used /= 0 only in highly banked memory.

! output

  integer, intent(out)               ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

! local variables

  integer        ::  mxdfft                                              !  array dimension for fft transform
  integer        ::  mxdwrk                                              !  array dimension for fft transform workspace

  integer        ::  nsfft(3)
  integer        ::  idsh                                                !  idshift after the checks


! idshift enlarges the leading dimension of the mesh (id = n1 + idshift).
! The arrays that hold the mesh (vscr and other quantities in real space)
! are dimensioned with the value of mxdfft = (n1+1)*n2*n3 given by size_fft,
! which is only enough for idshift = 0 or 1.  Other values may cause
! problems elsewhere (array overflows), so they are reset to 0.

  idsh = idshift
  if(idshift /= 0 .and. idshift /= 1) then
    write(6,*)
    write(6,'("   WARNING in size_kmscr:  idshift = ",i6," is not 0 or 1, ",   &
         & "it was changed to 0")') idshift
    write(6,*)

    idsh = 0
  endif

! the margin is only a sanity check, it is not changed.

  if(mrgdual < 1 .or. mrgdual > 5) then
    write(6,*)
    write(6,'("   WARNING in size_kmscr:  mrgdual = ",i6," is outside the ",  &
         & "range 1 to 5")') mrgdual
    write(6,*)
  endif

  if(flgdal == 'DUAL') then
    kmscr(1) = kmax(1)/2 + mrgdual
    kmscr(2) = kmax(2)/2 + mrgdual
    kmscr(3) = kmax(3)/2 + mrgdual
  else
    kmscr(1) = kmax(1)
    kmscr(2) = kmax(2)
    kmscr(3) = kmax(3)
  endif

  call size_fft(kmscr, nsfft, mxdfft, mxdwrk)

  kmscr(4) = nsfft(1)
  kmscr(5) = nsfft(2)
  kmscr(6) = nsfft(3)
  kmscr(7) = nsfft(1) + idsh

  return

end subroutine size_kmscr
