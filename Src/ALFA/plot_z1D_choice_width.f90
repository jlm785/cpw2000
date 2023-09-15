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

!>  Chooses the widths for the plotting of double averages.
!>
!>  For the double average see  PRL 61, 734 (1988).
!>
!>  \author       Jose Luis Martins
!>  \version      5.07
!>  \date         13 September 2023.
!>  \copyright    GNU Public License v2


subroutine plot_z1d_choice_width(ioreplay, nwidth, width, lopt,          &
         ntype, natom, nameat,                                           &
         height, nmat, nrepeat, widthgeom, widthrho,                     &
         nn, xfirst, xpeak,                                              &
         mxdtyp, mxdwid)


! Extracted from plot_rho_v_average. 13 September 2023. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdwid                          !<  maximum number of widths

  integer, intent(in)                ::  nmat                            !<  number of materials

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i

  integer, intent(in)                ::  nrepeat(nmat)                   !<  number of repeat units for each material

  real(REAL64), intent(in)           ::  widthgeom(nmat)                 !<  width for the double average from material info
  real(REAL64), intent(in)           ::  widthrho(nmat)                  !<  width for the double average from locally projected charge density
  real(REAL64), intent(in)           ::  xpeak(mxdtyp)                   !<  position of first peak of auto-correlation for atom of type i

  real(REAL64), intent(in)           ::  height                          !<  height of the cell
  integer, intent(in)                ::  nn                              !<  number of plot points
  real(REAL64), intent(in)           ::  xfirst                          !<  first peak of electron autocorrelation function

! output

  integer, intent(out)               ::  nwidth                          !<  number of widths
  real(REAL64), intent(out)          ::  width(mxdwid)                   !<  width for the double average
  logical, intent(out)               ::  lopt                            !<  if true the width is optimized

! allocatable arrays

  integer, allocatable               ::  indx(:)                         !  index of species with most atoms
  real(REAL64), allocatable          ::  widthold(:)                     !  width for the double average in earlier code

! main variables

  integer                     ::  nab, ntotal, nwd
  integer                     ::  ichoice
  integer                     ::  nwidthold
  integer                     ::  nptot

! counters

  integer      :: i, j, k, nt



! identifies most abundant species and finds default widths
! old scheme kept for compatibility

  allocate(indx(mxdtyp))
  allocate(widthold(mxdtyp))

  call isort(ntype, natom, indx)

  nab = 0
  ntotal = 0
  do nt = 1,ntype
    ntotal = ntotal + natom(nt)
  enddo
  do nt = 1,ntype
    if(4*natom(nt) > ntotal) nab = nab + 1
  enddo

  nptot = 0
  do j = 1,nmat
    nptot = nptot + nrepeat(j)
  enddo

  if(nab == 0) then
    nwidthold = 1
    widthold(1) = height/nptot
  else
    nwidthold = nab
    do j=1,nab
      widthold(j) = xpeak(indx(ntype-j+1))*height/nn
    enddo
  endif

  write(6,*)
  write(6,*)
  write(6,*)
  write(6,*) ' You have to choose both the number and'
  write(6,*) ' width of the square well averages'
  write(6,*)
  write(6,'("  Width from total number of planes: ",g14.6)') height/nptot

  write(6,*)
  write(6,'("  Width from electron density autocorrelation: ",g14.6)')   &
         xfirst*height/nn
  write(6,*)

  do nt = 1,ntype
    write(6,'("  Width from ",a2," atomic autocorrelation: ",g14.6)')    &
         nameat(indx(nt)),xpeak(indx(nt))*height/nn
    write(6,*)
  enddo

  write(6,*)
  write(6,*)
  write(6,*)
  write(6,'("  Old algorithm suggest ",i2," widths with values: ",       &
        &     99f12.5)') nwidthold,(widthold(j),j=1,nwidthold)
  write(6,*)

  write(6,*)
  write(6,'("  Materials description (layer size) suggest ",i2,          &
        &     " widths with values: ",99f12.5)')                         &
               nmat,(widthgeom(j),j=1,nmat)
  write(6,*)

  write(6,*)
  write(6,'("  Density auto-correlation for each material suggest ",i2,  &
        &     " widths with values: ",99f12.5)')                         &
               nmat,(widthrho(j),j=1,nmat)
  write(6,*)


  write(6,*)
  write(6,*) "  Choose which values you want"
  write(6,*) "  1) Optimization by searching flatter results"
  write(6,*) "  2) Materials description"
  write(6,*) "  3) Density autocorrelation"
  write(6,*) "  4) Old algorithm"
  write(6,*) "  5) Enter other values of your choice"
  write(6,*)
  write(6,*) "  Enter your choice (1--5)"


  read(5,*) ichoice
  write(ioreplay,'(2x,i5,"   choice of widths")') ichoice

  if(ichoice < 1 .or. ichoice > 5) then
    write(6,*)
    write(6,*) "  Wrong value, enter again"
    read(5,*) ichoice
    write(ioreplay,'(2x,i5,"   choice of widths")') ichoice
    if(ichoice < 1 .or. ichoice > 4) then
      write(6,*)
      write(6,*) "  Wrong value, using optimization"
      write(6,*)
      ichoice = 1
    endif
  endif

  lopt = .FALSE.
  if(ichoice == 1) lopt = .TRUE.

  nwidth = nmat

  if(ichoice == 1 .or. ichoice == 2) then
    do i = 1,nwidth
      width(i) = widthgeom(i)
    enddo
  elseif(ichoice == 3) then
    do i = 1,nwidth
      width(i) = widthrho(i)
    enddo
  elseif(ichoice == 4) then
    do i = 1,nwidth
      width(i) = widthold(i)
    enddo
  else

!   user chosen values

    write(6,'(" Enter number of averages ")')
    write(6,'(" Must be positive and smaller or equal to ",i5)') mxdwid
    write(6,*)

    read(5,*)  nwd
    write(ioreplay,'(2x,i8,"   number of average widths")') nwd

    if(nwd < 1 .or. nwd > mxdwid) then

      write(6,*) "  Invalid answer enter number of averages again"

      read(5,*)  nwd
      write(ioreplay,'(2x,i8,"   number of average widths again")') nwd

      if(nwd < 1 .or. nwd > mxdwid) then

        write(6,*) "  STOP Invalid answer again."

        stop

      endif


    endif

!   only place where nwidth =/= nmat

    nwidth = nwd
    write(6,'(" Enter the ",i5," widths ")') nwidth
    read(5,*) (width(i),i=1,nwidth)
    write(ioreplay,'(2x,20f15.5)') (width(i),i=1,nwidth)

  endif

! After this point widths are set.

  do k=1,nwidth

    if(width(k) <= 0 .or. width(k) > height) then

      write(6,*) "  STOP Invalid widths  "
      write(6,'("  Width ",i5," is ",f14.6)') k, width(k)
      write(6,'("  Should be positive and smaller than ",f14.6)') height

      stop

    endif

  enddo



! deallocates the stuff

  deallocate(widthold)
  deallocate(indx)

  return

end subroutine plot_z1D_choice_width
