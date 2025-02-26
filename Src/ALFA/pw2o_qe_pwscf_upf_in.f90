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

!>  Writes the Xy_nc_upf.in pseudopotential file for Quantum Espresso.
!>  May not work if it includes spin-orbit.  Use the pseudopotential
!>  generation code instead.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         15 October 2018. 25 February 2025.
!>  \copyright    GNU Public License v2

subroutine pw2o_qe_pwscf_upf_in(nameat,                                  &
    irel, iray, psdtitle,                                                &
    nkb)


! Written 15 October 2018.
! Documentation, 12 february 2021. JLM
! Case llocal is not available, 13 November 2023. JLM
!

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  character(len=3), intent(in)       ::  irel                            !<  type of calculation relativistic/spin
  character(len=60), intent(in)      ::  iray                            !<  information about pseudopotential
  character(len=10), intent(in)      ::  psdtitle(20)                    !<  further information about pseudopotential

  character(len=2), intent(in)       ::  nameat                          !<  chemical symbol for the type i
  integer, intent(in)                ::  nkb(0:3,-1:1)                   !<   kb pseudo.  normalization for atom k, ang. mom. l

! local:

  integer               ::  io                                           !  tape number
  character(len = 20)   ::  filename
  real(REAL64)          ::  rc, occup
  character(len = 2)    ::  nl

  integer               ::  iz

  integer               ::  lmax,llocal

  character(len=80)     ::  config
  character(len=3)      ::  str
  character(len=80)     ::  fmt
  character(len=80)     ::  pseudofile

  character(len=10)     ::  tmptitle

  integer               ::  ioerror

! constants

  real(REAL64), parameter    ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer       ::  j, l


  call p_tbl_charge(nameat,iz)
  call pw2o_qe_default_conf(iz,config)

  io = 10

! open file

  filename = adjustl(trim(nameat))//"_nc_upf.in"
  open(unit = io, file = trim(filename),status='UNKNOWN', form='FORMATTED')

  write(io,'(a7)') " &input"

  write(io,'(9x,"title=''",a2,"'',")') nameat
  write(io,'(9x,"zed=",f12.4,",")') iz*UM
  if(irel == 'rel') then
    write(io,'(9x,"rel=2",",")')
  else
    write(io,'(9x,"rel=0",",")')
  endif
  write(io,'(9x,"beta=0.5,")')
  write(io,'(9x,"eminld=-4.0,")')
  write(io,'(9x,"emaxld=4.0,")')
  write(io,'(9x,"deld=0.012d0,")')
  write(io,'(9x,"nld=3,")')
  j = len(adjustl(trim(config)))
  if(j < 10) then
    write(str,'("a",i1)') j
  else
    write(str,'("a",i2)') j
  endif
  fmt = '(9x,"config=''",'//str//',"''",)'
  write(io,fmt) config
  write(io,'(9x,"iswitch=3,")')
  write(io,'(9x,"dft=''LDA'',")')

  write(io,'(" /")')

  write(io,'(a8)') " &inputp"

  write(io,'(9x,"pseudotype=1,")')

! hack for the case llocal is not available (old pseudos)

  read(iray(58:60),*,iostat=ioerror) llocal
  if(ioerror /= 0) then
    write(6,*)
    write(6,*) '   WARNING    WARNING    WARNING'
    write(6,*) '   unable to find the local pseudopotential'
    write(6,*) '   change lloc in ',filename
    write(6,*)
    llocal = -9
  endif

  if(llocal == -1) then
    write(6,*)
    write(6,*) '   WARNING    WARNING    WARNING'
    write(6,*) '   QE does not have the maximum of potentials flag'
    write(6,*) '   change lloc in ',filename
    write(6,*)
  endif



  lmax = 0
  do l = 0,3
  do j = -1,1
    if(nkb(l,j) /= 0) then
      if(l > lmax) lmax = l
    endif
  enddo
  enddo
  if(llocal > lmax) lmax = llocal
  write(io,'(9x,"lloc=",i2,",")') llocal
  write(io,'(9x,"tm=.true.,")')
  pseudofile = adjustl(trim(nameat))//'_LDA_ncpp.UPF'
  write(io,'(9x,"file_pseudopw=''",a15,"'',")') trim(pseudofile)

  write(io,'(" /")')

  if(irel == 'rel') then
    write(io,'(i3)') 2*lmax+1
    do l=0,lmax
      tmptitle = psdtitle(2*l+1)
      read(tmptitle(1:2),*) nl
      if(nl(2:2) == 's') nl(2:2) = 'S'
      if(nl(2:2) == 'p') nl(2:2) = 'P'
      if(nl(2:2) == 'd') nl(2:2) = 'D'
      if(nl(2:2) == 'f') nl(2:2) = 'F'
      if(nl(2:2) == 'g') nl(2:2) = 'G'
      read(tmptitle(4:9),*) occup
      if(occup < 0.001) then
        write(6,*)
        write(6,*) '   WARNING    WARNING    WARNING'
        write(6,*) '   QE does not like unbound orbitals.  You may have to'
        write(6,*) '   make changes in ',filename,' for channel ',nl(2:2)
        write(6,*) '   and atom ',nameat
        write(6,*)
      endif
      tmptitle = psdtitle(2*l+2)
      read(tmptitle(6:10),*) rc
      write(io,'(a2,2i3,4f6.2)') nl,l+1,l,occup,ZERO,rc,rc
      if(l /=0) then
        write(io,'(a2,2i3,4f6.2)') nl,l+1,l,occup,ZERO,rc,rc
      endif
    enddo
  else
    write(io,'(i3)') lmax+1
    do l=0,lmax
      tmptitle = psdtitle(2*l+1)
      read(tmptitle(1:2),*) nl
      if(nl(2:2) == 's') nl(2:2) = 'S'
      if(nl(2:2) == 'p') nl(2:2) = 'P'
      if(nl(2:2) == 'd') nl(2:2) = 'D'
      if(nl(2:2) == 'f') nl(2:2) = 'F'
      if(nl(2:2) == 'g') nl(2:2) = 'G'
      read(tmptitle(4:9),*) occup
      tmptitle = psdtitle(2*l+2)
      read(psdtitle(6:10),*) rc
      write(io,'(a2,2i3,4f6.2)') nl,l+1,l,occup,ZERO,rc,rc
    enddo
  endif

  close(unit=io)

  return

end subroutine pw2o_qe_pwscf_upf_in
