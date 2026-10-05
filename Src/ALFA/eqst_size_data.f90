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

!>  Opens the data file and finds the array dimensions.
!>  The default file name is eqst.dat, if it does not exist
!>  asks for the file name. The file is left open and rewound.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         3 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_size_data(iodata, datafile, mxdstr, mxdnpt)

! Written from eqst_read_data. 3 October 2026. claude

  implicit none

! input

  integer, intent(in)                   ::  iodata                       !<  tape number of the data file

! output

  character(len=*), intent(out)         ::  datafile                     !<  name of the data file
  integer, intent(out)                  ::  mxdstr                       !<  number of structures
  integer, intent(out)                  ::  mxdnpt                       !<  maximum number of energy points (per structure)

! local variables

  character(len=100)        ::  fline
  logical                   ::  lexist
  integer                   ::  ioerr
  integer                   ::  nstr, npt

! counters

  integer                   ::  ns, mm


! finds the data file

  datafile = 'eqst.dat'

  inquire(file = trim(datafile), exist = lexist)

  if(.not. lexist) then
    write(6,*)
    write(6,*) '  The default data file eqst.dat was not found.'
    write(6,*) '  Enter the name of the data file'
    write(6,*)

    read(5,'(a)') datafile
    datafile = adjustl(datafile)

    inquire(file = trim(datafile), exist = lexist)

    if(.not. lexist) then
      write(6,*)
      write(6,*) '  STOPPED in eqst_size_data:  file ',trim(datafile),' not found'
      write(6,*)

      stop

    endif
  endif

  open(unit = iodata, file = trim(datafile), status = 'OLD', form = 'FORMATTED')

  write(6,*)
  write(6,*) '  Reading data from file ',trim(datafile)
  write(6,*)

! type of equation of state is checked later

  read(iodata,*,iostat=ioerr) fline
  if(ioerr /= 0) call eqst_size_data_stop('type of equation of state', 0)

  read(iodata,*,iostat=ioerr) nstr
  if(ioerr /= 0) call eqst_size_data_stop('number of structures', 0)
  if(nstr < 1) then
    write(6,*)
    write(6,*) '  STOPPED in eqst_size_data:  number of structures = ',nstr
    write(6,*)

    stop

  endif

  mxdstr = nstr
  mxdnpt = 0

  do ns = 1,nstr

    read(iodata,*,iostat=ioerr) npt
    if(ioerr /= 0) call eqst_size_data_stop('number of points', ns)
    if(npt < 4) then
      write(6,*)
      write(6,*) '  STOPPED in eqst_size_data:  structure ',ns,' has ',npt,' points.'
      write(6,*) '  At least 4 are needed for the fit.'
      write(6,*)

      stop

    endif

    if(npt > mxdnpt) mxdnpt = npt

!   volume factor and label, and the data lines

    read(iodata,'(a100)',iostat=ioerr) fline
    if(ioerr /= 0) call eqst_size_data_stop('volume factor and label', ns)

    do mm = 1,npt
      read(iodata,'(a100)',iostat=ioerr) fline
      if(ioerr /= 0) call eqst_size_data_stop('data point', ns)
    enddo

  enddo

! maximum pressure and number of points for the plots

  if(nstr > 1) then
    read(iodata,'(a100)',iostat=ioerr) fline
    if(ioerr /= 0) call eqst_size_data_stop('maximum pressure, number of points', 0)
  endif

  rewind(iodata)

  return

contains

!>  stops with an error message if the data file ends prematurely

  subroutine eqst_size_data_stop(what, nst)

    implicit none

    character(len=*), intent(in)        ::  what                         !<  what was being read
    integer, intent(in)                 ::  nst                          !<  structure being read, 0 if none

    write(6,*)
    write(6,*) '  STOPPED in eqst_size_data:  error reading ',what
    if(nst > 0) write(6,*) '  for structure ',nst
    write(6,*) '  in file ',trim(datafile)
    write(6,*)

    stop

  end subroutine eqst_size_data_stop

end subroutine eqst_size_data
