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

!>  Writes the atom.dat file for the pseudopotential generation code
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         20 February 2025
!>  \copyright    GNU Public License v2

subroutine pw2o_atom_dat(nameat,                                         &
    irel, icore, icorr, iray, psdtitle, zv)


! Written 20 February 2025.

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  character(len=3), intent(in)       ::  irel                            !<  type of calculation relativistic/spin
  character(len=4), intent(in)       ::  icore                           !<  type of partial core correction
  character(len=2), intent(in)       ::  icorr                           !<  type of correlation
  character(len=60), intent(in)      ::  iray                            !<  information about pseudopotential
  character(len=10), intent(in)      ::  psdtitle(20)                    !<  further information about pseudopotential

  character(len=2), intent(in)       ::  nameat                          !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  zv                              !<  valence of atom with type i


! local:

  integer               ::  io                                           !  tape number
  character(len = 20)   ::  filename
  character(len = 2)    ::  pcc
  character(len = 1)    ::  ispp

  real(REAL64)          ::  rmax                                         !  maximum mesh radius
  real(REAL64)          ::  aa                                           !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)          ::  bb                                           !  a = exp(-aa)/znuc, b = 1/bb

  integer                           ::  ncore                            !  canonical number of core orbitals
  integer                           ::  nval                             !  canonical number of interesting valence orbitals
  integer                           ::  jhard                            !  flag for accuracy/speed compromise
  integer                           ::  no(5), lo(5)                     !  configuration
  real(REAL64)                      ::  zo(5)                            !  occupation
  real(REAL64)                      ::  rc(5)                            !  core radii
  character(len=30)                 ::  status                           !  quality of

  character(len=10)     ::  tmptitle

  integer               ::  iz

  real(REAL64)          ::  zvsum
  integer               ::  ncan
  character(len=1)      ::  spd
  integer               ::  lcan
  real(REAL64)          ::  rcan

! constants

  real(REAL64), parameter    ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-6_REAL64

! counters

  integer       ::  n


   call p_tbl_charge(nameat,iz)


  io = 10


! convert the flags

  if(icore == 'nc  ') then
    pcc = 'pg'
  elseif(icore == 'pcec') then
    pcc = 'pe'
  else
    write(6,*)
    write(6,*) '   unrecognized flag for core correction: ',icore
    write(6,*)

    return

  endif

  if(irel == 'isp') then
    ispp = 's'
  elseif(irel == 'rel') then
    ispp = 'r'
  elseif(irel == 'nrl') then
    ispp = ' '
  else
    write(6,*)
    write(6,*) '   unrecognized flag for relativistic calculation: ',irel
    write(6,*)

    return

  endif

! default mesh, configuration, core radii

  call atom_atm_default_mesh(rmax, aa, bb)

  jhard = 0
  call atom_p_tbl_config(nameat, ncore, nval, no, lo, zo, jhard)

  call atom_p_tbl_psd_tm2(nameat, rc, status)

  zvsum = ZERO
  do n = 1,nval
    zvsum = zvsum + zo(n)
  enddo

  if(abs(zvsum - zv) > EPS) then
    write(6,*)
    write(6,'("  number of valence electrons is not canonical :",2f10.3)') zvsum, zv
    write(6,*) '  Not enough information to write atom.dat'
    write(6,*)

    return

  endif

  do n = 1,nval
    tmptitle = psdtitle(2*n-1)
    read(tmptitle(1:1),'(i1)') ncan

    if(ncan /= no(n)) then
      write(6,*)
      write(6,*) '  inconsistent principal quantum number :', ncan, no(n)
      write(6,*) '  Not enough information to write atom.dat'
      write(6,*)

      return

    endif

    read(tmptitle(2:2),'(a1)') spd

    if(spd == 's') then
      lcan = 0
    elseif(spd == 'p') then
      lcan = 1
    elseif(spd == 'd') then
      lcan = 2
    elseif(spd == 'f') then
      lcan = 3
    elseif(spd == 'g') then
      lcan = 4
    else
      write(6,*)
      write(6,*) '  inconsistent angular quantum character :', spd
      write(6,*) '  Not enough information to write atom.dat'
      write(6,*)

      return

    endif

    if(lcan /= lo(n)) then
      write(6,*)
      write(6,*) '  inconsistent angular quantum number :', lcan, lo(n)
      write(6,*) '  Not enough information to write atom.dat'
      write(6,*)

      return

    endif

    tmptitle = psdtitle(2*n)
    if(ispp == ' ' .or. ispp == 'r') then
      read(tmptitle(6:10),'(f5.2)') rcan
    elseif(ispp == 's') then
      read(tmptitle(7:10),'(f4.2)') rcan
    endif
    rc(n) = rcan

  enddo

! open file

  filename = adjustl(trim(nameat))//"_atom.dat"
  open(unit = io, file = trim(filename),status='UNKNOWN', form='FORMATTED')

  write(io,'(3x,a2,10x,a60)') pcc,iray
  write(io,'(8x,a3)') 'tm2'
  write(io,'(" n=",a2," c=",a2,a1)') nameat,icorr,ispp


  write(io,'(6f10.1)') UM*iz, ZERO, ZERO, rmax, aa, bb

  write(io,'(1x,i4,1x,i4)') ncore,nval

  do n = 1,nval
    write(io,'(1x,i4,1x,i4,1x,f8.2,2x,f8.2)') no(n), lo(n), zo(n), ZERO
  enddo

  write(io,'(1x,5(f7.2,2x))') (rc(n),n=1,nval)
  write(io,'(50x)')

  close(unit=io)

  return

end subroutine pw2o_atom_dat
