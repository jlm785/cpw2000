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

!>  Calculates the auto-correlation function from a function projected
!>  in several intervals
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         8 March 2023.
!>  \copyright    GNU Public License v2

subroutine plot_z1D_local_corr(nn, ave, nmat, height, rleft, nrepeat, width)

! written 8 March 2023. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nn                              !<  number of points
  integer, intent(in)                ::  nmat                            !<  number of materials

  real(REAL64), intent(in)           ::  ave(nn)                         !<  layer average of function in the 3d direction of fft grid

  real(REAL64), intent(in)           ::  height                          !<  height of the cell
  real(REAL64), intent(in)           ::  rleft(nmat)                     !<  left boundary of material
  integer, intent(in)                ::  nrepeat(nmat)                   !<  number of repeating units for each material

! output

  real(REAL64), intent(out)          ::  width(nmat)                     !<  width for the double average

! allocatable arrays

  real(REAL64), allocatable   ::  prave(:)                               !  projected layer average
  real(REAL64), allocatable   ::  fac(:)                                 !  autocorrelation function

! local variables

  integer           ::  kstart                                           !  starting point to search for a maximum
  real(REAL64)      ::  xpeak                                            !  position of nearest maximum
  logical           ::  lfound                                           !  true if a peak was found away from extrema

  real(REAL64)      ::  rright
  integer           ::  nleft, nright, nnmat



! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer       ::  n, j


  allocate(prave(nn))
  allocate(fac(0:nn/2))

  do n = 1,nmat

    rright = rleft(mod(n,nmat) + 1)
    nleft = nint(rleft(n)*nn)
    nright = nint(rright*nn)

    prave = ZERO

    if(rright - rleft(n) > ZERO) then
      nnmat = nright-nleft
      kstart = nnmat / nrepeat(n)
      do j = 1,nnmat
        prave(j) = ave(j+nleft)
        WRITE(6,*) J,j+nleft
      enddo
    else
      nnmat = nright+nn-nleft
      kstart = nnmat / nrepeat(n)
      do j = 1,nnmat
        prave(j) = ave(mod(j+nleft-1,nn)+1)
        WRITE(6,*) J,mod(j+nleft-1,nn)+1
      enddo
    endif

    call plot_z1D_auto_corr(nnmat, prave, kstart, fac, xpeak, lfound)

    if(.not. lfound) then
      write(6,*)
      write(6,*) '  peak for material ',n,' may not be accurate'
      write(6,*)
    endif

    width(n) = xpeak*height/nn

  enddo


  deallocate(prave)
  deallocate(fac)

  return

end subroutine plot_z1D_local_corr
