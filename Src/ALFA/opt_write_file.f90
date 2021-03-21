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

!>writes the imaginary part of the dielectric function to a file.
!>It is commented for gnuplot or xmgrace!

  subroutine opt_write_file(filename, iotape, title, subtitle, func,     &
             nhist, ehist, dhist, vcell, ztot)

! Adapted from dos_print_file July 2020. CLR
! Modified, normalization, 1 October 2020. JLM
! Modified, float, 17 December 2020. JLM
! copyright Carlos Loia Reis,  INESC-MN.

! version 4.99 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  character(len=*) , intent(in)       ::  filename                       !<  input file
  integer, intent(in)                 ::  iotape                         !<  tape number
  character(len=50), intent(in)       ::  title                          !<  title for plots
  character(len=140), intent(in)      ::  subtitle                       !<  subtitle for plots
  character(len=50) , intent(in)      ::  func                           !<  name of function

  integer, intent(in)                 ::  nhist                          !<  number of points in histograms
  real(REAL64), intent(in)            ::  ehist(nhist)                   !<  energies of the histogram
  real(REAL64), intent(in)            ::  dhist(nhist)                   !<  density of states

  real(REAL64), intent(in)            ::  vcell                          !<  unit cell volume (in atomic units)
  real(REAL64), intent(in)            ::  ztot                           !<  total charge density (electrons/cell)

! local variables

  logical     ::  lfloat

! counters

  integer     ::  i

! constants

  real(REAL64), parameter :: HARTREE = 27.21138386_REAL64

  lfloat = .FALSE.

  do i=1,nhist
     if(dhist(i) > 10000.0) then
       lfloat = .TRUE.
       exit
     endif
  enddo
 
! printout to filename

  open(unit=iotape,file=filename,form='formatted')

  write(iotape,'("#   ",a50)') title
  write(iotape,'("#   ",a140)') subtitle
  write(iotape,'("#")')
  write(iotape,'("# ",f14.5,"  vcell (cell volume a.u. )")') vcell
  write(iotape,'("# ",f14.5,"  ztot (electrons/cell )")') ztot
  write(iotape,'("#")')
  write(iotape,'("#   Energy in eV",5x,a50)') func
  write(iotape,'("#")')

  if(lfloat) then
    do i=1,nhist
      write(iotape,'(2x,f15.8,2x,e15.8)') ehist(i)*HARTREE, dhist(i)
    enddo
  else
    do i=1,nhist
      write(iotape,'(2x,2(2x,f15.8))') ehist(i)*HARTREE, dhist(i)
    enddo
  endif

  close(unit=iotape)

  return
end subroutine opt_write_file
