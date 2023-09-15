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

!>  Writes the restart file.
!>  For molecular dynamics trajectory should be the same.
!>  for l_bfgs minimization, algorithm is restarted.
!>
!>  \author       Jose Luis Martins, Benedito Costa Cabral
!>  \version      5.06
!>  \date         26 May 1999, 17 August 2004, 14 October 2022.
!>  \copyright    GNU Public License v2

subroutine move_restart_out(io7, filename, flgcal,                       &
    ntype, natom, nameat, atmass, rat, vat, adot, vadot,                 &
    istmd, tstep, beta, tempk, iseed, strext, press, celmas,             &
    rat1, frc1, adot1, frcel1,                                           &
    mxdatm, mxdtyp)

! Written 26 may 99. jlm
! Rewritten 19 november 2002. jlm
! Modified 14 october 2003. jlm
! Modified (formats) 17 august 2004. jlm
! Modified 6 January 2017, f90. JLM
! Modified, documentation, August 2019. JLM
! Modified, filename, 14 October 2022. JLM

  implicit none

  integer, parameter :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type

  character(len=6), intent(in)       ::  flgcal                          !<  type of calculation
  integer, intent(in)                ::  io7                             !<  tape number
  character(len=*), intent(in)       ::  filename                        !<  file name

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i

  real(REAL64), intent(in)           ::  atmass(mxdtyp)                  !<  atomic mass of atoms of type i

  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  vat(3,mxdatm,mxdtyp)            !<  d rat / d t  velocity in lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  real(REAL64), intent(in)           ::  vadot(3,3)                      !<  d adot / d t  rate of change of metric

  integer,intent(in)                 ::  istmd                           !<  md step. Equal to 1 in first step of molecular dynamics
  real(REAL64), intent(in)           ::  tstep                           !<  molecular dynamics time step (in a.u.)

  real(REAL64), intent(in)           ::  tempk                           !<  ionic temperature (in Kelvin)
  real(REAL64), intent(in)           ::  beta                            !<  friction coefficient/mass (in a.u.)

  integer, intent(in)                ::  iseed                           !<  seed for random number generator

  real(REAL64), intent(in)           ::  strext(3,3)                     !<  external applied stress
  real(REAL64), intent(in)           ::  press                           !<  external pressure
  real(real64), intent(in)           ::  celmas                          !<  fictitious cell mass

  real(REAL64), intent(in)           ::  rat1(3,mxdatm,mxdtyp)           !<  previous value of lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  frc1(3,mxdatm,mxdtyp)           !<  previous value of force on atom j of type i
  real(REAL64), intent(in)           ::  adot1(3,3)                      !<  previous value of adot
  real(REAL64), intent(in)           ::  frcel1(3,3)                     !<  previous cell "force" (covariant components)

! counters

  integer       ::  nt, i, j

  open(unit = io7, file = filename, status = 'unknown', form = 'formatted')

  write(io7,'(2x,a6)') flgcal

  write(io7,'(9(2x,e24.16))') ((adot(i,j),i=1,3),j=1,3)
  write(io7,'(2x,i10)') ntype

  do nt = 1,ntype
    write(io7,'(2x,i10,2x,a2,2x,e24.16)') natom(nt), nameat(nt), atmass(nt)

    do i = 1,natom(nt)
      write(io7,'(3(2x,e24.16))') (rat(j,i,nt),j=1,3)
    enddo
  enddo

  if(flgcal == 'MICRO ') then

    write(io7,'(2x,i10,2x,e24.16,"    iteration,etc")') istmd,tstep
    do nt=1,ntype
      do i=1,natom(nt)
        write(io7,'(9(2x,e24.16))') (vat(j,i,nt),j=1,3),                 &
                    (rat1(j,i,nt),j=1,3), (frc1(j,i,nt),j=1,3)
      enddo
    enddo

  elseif(flgcal == 'LANG  ') then

    write(io7,'(2x,i10,3(2x,e24.16),2x,i10,"    iteration,etc")')        &
                 istmd, tstep, beta, tempk, iseed
    do nt = 1,ntype
      do i = 1,natom(nt)
        write(io7,'(11(2x,e24.16))') (vat(j,i,nt),j=1,3),                &
                    (rat1(j,i,nt),j=1,3), (frc1(j,i,nt),j=1,3)
      enddo
    enddo

  elseif(flgcal == 'VCSLNG' .or. flgcal == 'EPILNG') then

    write(io7,'(2x,i10,3(2x,e24.16),2x,i10,"    iteration,etc")')        &
                 istmd, tstep, beta, tempk, iseed
    write(io7,'(11(2x,e24.16))') press,((strext(i,j),i=1,3),j=1,3), celmas
    do nt = 1,ntype
      do i = 1,natom(nt)
        write(io7,'(9(2x,e24.16))') (vat(j,i,nt),j=1,3),                 &
                    (rat1(j,i,nt),j=1,3), (frc1(j,i,nt),j=1,3)
      enddo
    enddo
    write(io7,'(9(2x,e24.16))') ((vadot(i,j),i=1,3),j=1,3)
    write(io7,'(9(2x,e24.16))') ((adot1(i,j),i=1,3),j=1,3)
    write(io7,'(9(2x,e24.16))') ((frcel1(i,j),i=1,3),j=1,3)

  elseif(flgcal == 'VCSMIC') then

    write(io7,'(2x,i10,2x,e24.16,"    iteration,etc")') istmd ,tstep
    write(io7,'(11(2x,e24.16))') press,((strext(i,j),i=1,3),j=1,3), celmas
    do nt = 1,ntype
      do i = 1,natom(nt)
        write(io7,'(9(2x,e24.16))') (vat(j,i,nt),j=1,3),                 &
                    (rat1(j,i,nt),j=1,3), (frc1(j,i,nt),j=1,3)
      enddo
    enddo
    write(io7,'(9(2x,e24.16))') ((vadot(i,j),i=1,3),j=1,3)
    write(io7,'(9(2x,e24.16))') ((adot1(i,j),i=1,3),j=1,3)
    write(io7,'(9(2x,e24.16))') ((frcel1(i,j),i=1,3),j=1,3)

  elseif(flgcal == 'VCSLBF' .or. flgcal == 'EPILBF') then

    write(io7,'(11(2x,e24.16))') press,((strext(i,j),i=1,3),j=1,3)

  endif

  close(unit = io7)

  return

end subroutine move_restart_out
