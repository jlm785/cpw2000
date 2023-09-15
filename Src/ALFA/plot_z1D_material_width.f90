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

!>  Asks the user to choose tha materials widths
!>  For new way of doing the double averages.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         2 April 2023.
!>  \copyright    GNU Public License v2


subroutine plot_z1D_material_width(ioreplay,                             &
       ntot, iptype, ipnatom, rat, adot,                                 &
       height, nptot, nrepeat, rbottom, widthgeom,                       &
       nmat, mxdtyp, mxdatm)

! extracted from plot_rho_v_average. 2 April 2023. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations

  integer, intent(in)                ::  nmat                            !<  number of materials
  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  ntot                            !<  total number of atoms
  integer, intent(in)                ::  iptype(mxdatm*mxdtyp)           !<  points to atom type
  integer, intent(in)                ::  ipnatom(mxdatm*mxdtyp)          !<  points to atom of a given type

! output


  real(REAL64), intent(out)          ::  height                          !<  height of the cell
  integer, intent(out)               ::  nptot                           !<  total number of repeat units (sum of nrepeat)

  integer, intent(out)               ::  nrepeat(nmat)                   !<  number of repeat units
  real(REAL64), intent(out)          ::  rbottom(nmat)                   !<  bottom atom (lowest z) of material (zero if nmat = 1)

  real(REAL64), intent(out)          ::  widthgeom(nmat)                 !<  width for the double average from material info

! other variables

  real(REAL64)        ::  rtop
  logical             ::  lrepeat

! allocatable array

  integer, allocatable               ::  kchoice(:)

! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  BOHR = 0.5291772109_REAL64

! counters

  integer      ::  i, j, k, k2


  allocate(kchoice(nmat))

  do j = 1,nmat

    write(6,*)
    write(6,'("  Enter the number of repeat units for material ",i3)') j
    write(6,*) '  (atomic or lattice planes, depends on material)'
    write(6,*)

    read(5,*) nrepeat(j)
    write(ioreplay,'(2x,i8,"   number of repeating units for material",i3)') nrepeat(j),j

    if(nrepeat(j) < 1 .or. nrepeat(j) > ntot) then

      write(6,*)
      write(6,*) '  Answer must be between 1 and ',ntot
      write(6,*) '  Enter the index again'
      write(6,*)

      read(5,*) nrepeat(j)
      write(ioreplay,'(2x,i8,"    repeating unit again")') nrepeat(j)

      if(nrepeat(j) < 1) then
        write(6,*)
        write(6,*) '  Wrong value exiting program'
        write(6,*)

        stop

      endif

    endif

!   for nmat = 1

    rbottom(j) = ZERO

    if(nmat > 1) then

      write(6,*)
      write(6,*) '  Use the above table to enter the index of the bottom'
      write(6,*) '  atom (lowest value of z, last column) in that material.'
      write(6,*) '  Take into account the periodicity!'
      write(6,*)

      read(5,*) k
      write(ioreplay,'(2x,i8,"    bottom atom")') k

      if(k < 1 .or. k > ntot) then

        write(6,*)
        write(6,*) '  Answer must be between 1 and ',ntot
        write(6,*) '  Enter the index again'
        write(6,*)

        read(5,*) k2
        write(ioreplay,'(2x,i8,"    bottom atom again")') k2

        if(k2 < 1 .or. k2 > ntot) then
          write(6,*)
          write(6,*) '  Wrong value exiting program'
          write(6,*)

          stop

        endif

        k = k2

      endif

      if(j > 1) then
        lrepeat = .FALSE.
        do i = 1,j-1
          if (k == kchoice(i)) lrepeat = .TRUE.
        enddo

        if(lrepeat) then

          write(6,*)
          write(6,*) '  Bottom atoms cannot be repeated'
          write(6,*) '  Enter the index again'
          write(6,*)

          read(5,*) k
          write(ioreplay,'(2x,i8,"    bottom atom repeated")') k

          if(k < 1 .or. k > ntot) then

            write(6,*)
            write(6,*) '  Answer must be between 1 and ',ntot
            write(6,*)
            write(6,*) '  Wrong value exiting program'
            write(6,*)

            stop

          endif

          do i = 1,j-1
            if (k == kchoice(i)) then

              write(6,*)
              write(6,*) '  Bottom atoms cannot be repeated'
              write(6,*)
              write(6,*) '  Wrong value exiting program'
              write(6,*)

              stop

            endif
          enddo

        endif

      endif

      kchoice(j) = k

!     rbottom should always be positive

      if(k > 1) then
        rbottom(j) = rat(3,ipnatom(k),iptype(k))                         &
                     - floor(rat(3,ipnatom(k),iptype(k))+0.000001)       &
                     + rat(3,ipnatom(k-1),iptype(k-1))                   &
                     - floor(rat(3,ipnatom(k-1),iptype(k-1))+0.000001)
        rbottom(j) = rbottom(j) / 2
      else
        rbottom(1) = rat(3,ipnatom(1),iptype(1))                         &
                     - floor(rat(3,ipnatom(1),iptype(1))+0.000001)       &
                     + rat(3,ipnatom(ntot),iptype(ntot))                 &
                     - (floor(rat(3,ipnatom(ntot),iptype(ntot))+0.000001)+1)
        rbottom(1) = rbottom(1) / 2
        if(rbottom(1) < ZERO) rbottom(1) =  rbottom(1) + UM
      endif

      write(6,*)
      write(6,'("  Bottom boundary of material ",i3," is at  ",f8.3)') j,rbottom(j)

    endif

  enddo

  nptot = 0
  do j = 1,nmat
    nptot = nptot + nrepeat(j)
  enddo

! rewrite code if you want to enable other direction!!!!
! As it is it averages in the plane of first two lattice vectors,
! that is, in the direction of the third reciprocal lattice vector.


  height = adot(3,3) - adot(2,3)*adot(3,2)/adot(2,2)                     &
                     - adot(1,3)*adot(3,1)/adot(1,1)                     &
           + 2*adot(1,2)*adot(2,3)*adot(3,1)/(adot(2,2)*adot(1,1))
  height = sqrt(height)

  write(6,*)
  write(6,'("  the height of the cell is: ",f12.4," Ang.")') height*BOHR
  write(6,*)

  if(nmat > 1) then

    do j = 1,nmat
      rtop = rbottom(mod(j,nmat) + 1)
      if(rtop - rbottom(j) > ZERO) then
        widthgeom(j) = (rtop- rbottom(j))*height / nrepeat(j)
      else
        widthgeom(j) = (rtop + UM - rbottom(j))*height / nrepeat(j)
      endif
    enddo

  else
    widthgeom(1) = height / nptot
  endif

  deallocate(kchoice)

  return

end subroutine plot_z1D_material_width
