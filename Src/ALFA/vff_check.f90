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

!>  checks if the existing bonds are allowed by the VFF parametrization
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.11 (1.7 of md)
!>  \date         November 2020.  5 June 2024.
!>  \copyright    GNU Public License v2

subroutine vff_check(iowrite, nbond, ibond, nangl, iangl,                &
                     dist, nameat, ityp, natotal, mxdtyp)

! Written November 2020. J.L.Martins, C.S.Loia
! Indentation, 5 June 2024. JLM

  implicit none

  integer, parameter        :: REAL64 = selected_real_kind(12)

! input variables

  integer, intent(in)       :: natotal                                   !<  maximum number of atoms
  integer, intent(in)       :: mxdtyp                                    !<  maximum number of type of atoms

  integer, intent(in)       :: iowrite                                   !<  tape number

  integer, intent(in)       :: nbond                                     !<  number of bonds (should be 2*natotal)
  integer, intent(in)       :: ibond(2,2*natotal)                        !<  atoms at each end of the bond
  integer, intent(in)       :: nangl                                     !<  number of angles (should be 6*natotal)
  integer, intent(in)       :: iangl(3,6*natotal)                        !<  atoms that define the

  real(REAL64), intent(in)  :: dist((mxdtyp*(mxdtyp+1))/2)               !<  equilibrium distance for bond.  It is < 0 if bond is not allowed
  integer, intent(in)       :: ityp(natotal)
  character(len=2),intent(in) :: nameat(mxdtyp)                          !<  chemical symbol of type of atom

! parameters

  real(REAL64), parameter  :: ZERO = 0.0_REAL64

! counters, etc...

  integer  :: n
  integer  :: ia, ja, ka, itia, itja, itka, itija, itika, itjka


! loop over bonds

  do n = 1,nbond

    ia = ibond(1,n)
    ja = ibond(2,n)
    itia = ityp(ia)
    itja = ityp(ja)
    if(itia > itja) then
      itija = (itia*(itia-1))/2 + itja
    else
      itija = (itja*(itja-1))/2 + itia
    endif

    if(dist(itija) < ZERO) then

      write(iowrite,*)
      write(iowrite,*) '  stopped in vff_check: bonds between '
      write(iowrite,'(3x,a2," and ",a2," are not parametrized")')        &
              nameat(itia), nameat(itja)

      stop

    endif

  enddo

! loop over angles

  do n = 1,nangl

    ia = iangl(1,n)
    ja = iangl(2,n)
    ka = iangl(3,n)
    itia = ityp(ia)
    itja = ityp(ja)
    itka = ityp(ka)
    if(itja > itka) then
      itjka = (itja*(itja-1))/2 + itka
    else
      itjka = (itka*(itka-1))/2 + itja
    endif
    if(itia > itja) then
      itija = (itia*(itia-1))/2 + itja
    else
      itija = (itja*(itja-1))/2 + itia
    endif
    if(itia > itka) then
      itika = (itia*(itia-1))/2 + itka
    else
      itika = (itka*(itka-1))/2 + itia
    endif

    if(dist(itija) < ZERO) then

      write(iowrite,*)
      write(iowrite,*) '  stopped in vff_check: bonds between '
      write(iowrite,'(2x,a2," and ",a2," are not parametrized")')        &
              nameat(itia), nameat(itja)

      stop

    endif

    if(dist(itika) < ZERO) then

      write(iowrite,*)
      write(iowrite,*) '  stopped in vff_check: bonds between '
      write(iowrite,'(2x,a2," and ",a2," are not parametrized")')        &
              nameat(itia), nameat(itka)

      stop

    endif

  enddo

  return

end subroutine vff_check
