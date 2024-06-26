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

subroutine plot_z1D_local_corr(nn, ave, nmat, height, rbottom, nrepeat, width)

!  written 8 March 2023. JLM
!  changed left/right to bottom/top. 5 September 2023

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nn                              !<  number of points
  integer, intent(in)                ::  nmat                            !<  number of materials

  real(REAL64), intent(in)           ::  ave(nn)                         !<  layer average of function in the 3d direction of fft grid

  real(REAL64), intent(in)           ::  height                          !<  height of the cell
  real(REAL64), intent(in)           ::  rbottom(nmat)                   !<  bottom boundary of material
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

  real(REAL64)      ::  rtop                                             !  top boundary of material
  integer           ::  nbottom, ntop, nnmat                             !  index of bottom of material, top of material and number of points



! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64
  real(REAL64), parameter  ::  EPS = 1.0E-6_REAL64

! counters

  integer       ::  n, j


  allocate(prave(nn))
  allocate(fac(0:nn/2))

  do n = 1,nmat

    rtop = rbottom(mod(n,nmat) + 1)
    nbottom = nint(rbottom(n)*nn)
    ntop = nint(rtop*nn)

    prave = ZERO

!   for nmat = 1 rtop=rbottom

    if(rtop - rbottom(n) > EPS) then
      nnmat = ntop-nbottom
      kstart = nnmat / nrepeat(n)
      do j = 1,nnmat
        prave(j) = ave(j+nbottom)
      enddo
    else
      nnmat = ntop+nn-nbottom
      kstart = nnmat / nrepeat(n)
      do j = 1,nnmat
        prave(j) = ave(mod(j+nbottom-1,nn)+1)
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
