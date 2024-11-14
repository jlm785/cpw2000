!> writes the Xy_nc_upf.in pseudopotential file for Quantum Espresso

subroutine write_pwscf_upf_in(nameat,                                    &
    irel, iray, ititle,                                                  &
    nkb)

! subroutine write_pwscf_upf_in(nameat, zv,                                &
!     irel, icore, icorr, iray, ititle,                                    &
!     nqnl, delqnl, vkbraw, nkb, vloc, dcor, dval,                         &
!     mxdlqp)

! Written 15 October 2018.
! Documentation, 12 february 2021. JLM
! case llocal is not avaulable, 13 November 2023. JLM
! copyright  J.L.Martins, INESC-MN.

! version 5.11

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

!  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential

  character(len=3), intent(in)       ::  irel                            !<  type of calculation relativistic/spin
!  character(len=4), intent(in)       ::  icore                           !<  type of partial core correction
!  character(len=2), intent(in)       ::  icorr                           !<  type of correlation
  character(len=60), intent(in)      ::  iray                            !<  information about pseudopotential
  character(len=70), intent(in)      ::  ititle                          !<  further information about pseudopotential

  character(len=2), intent(in)       ::  nameat                          !<  chemical symbol for the type i
!  real(REAL64), intent(in)           ::  zv                              !<  valence of atom with type i

!  integer, intent(in)                ::  nqnl                            !<  number of points for pseudo interpolation for atom k
!  real(REAL64), intent(in)           ::  delqnl                          !<  step used in the pseudo interpolation for atom k
!  real(REAL64), intent(in)       :: vkbraw(-2:mxdlqp,0:3,-1:1)           !<  (1/q**l) * kb nonlocal pseudo. for atom k, ang. mom. l. (non normalized to vcell, hartree)
  integer, intent(in)                ::  nkb(0:3,-1:1)                   !<   kb pseudo.  normalization for atom k, ang. mom. l
!  real(REAL64), intent(in)           ::  vloc(-1:mxdlqp)                 !<  local pseudopotential for atom k (hartree)
!  real(REAL64), intent(in)           ::  dcor(-1:mxdlqp)                 !<  core charge density for atom k
!  real(REAL64), intent(in)           ::  dval(-1:mxdlqp)                 !<  valence charge density for atom k

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

  integer               ::  ioerror

! constants

  real(REAL64), parameter    ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer       ::  j, l


  call p_tbl_charge(nameat,iz)
  call default_conf(iz,config)

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
      read(ititle(l*20+1:l*20+2),*) nl
      if(nl(2:2) == 's') nl(2:2) = 'S'
      if(nl(2:2) == 'p') nl(2:2) = 'P'
      if(nl(2:2) == 'd') nl(2:2) = 'D'
      if(nl(2:2) == 'f') nl(2:2) = 'F'
      if(nl(2:2) == 'g') nl(2:2) = 'G'
      read(ititle(l*20+4:l*20+9),*) occup
      if(occup < 0.001) then
        write(6,*)
        write(6,*) '   WARNING    WARNING    WARNING'
        write(6,*) '   QE does not like unbound orbitals.  You may have to'
        write(6,*) '   make changes in ',filename,' for channel ',nl(2:2)
        write(6,*) '   and atom ',nameat
        write(6,*)
      endif
      read(ititle(l*20+16:l*20+20),*) rc
      write(io,'(a2,2i3,4f6.2)') nl,l+1,l,occup,ZERO,rc,rc
      if(l /=0) then
        write(io,'(a2,2i3,4f6.2)') nl,l+1,l,occup,ZERO,rc,rc
      endif
    enddo
  else
    write(io,'(i3)') lmax+1
    do l=0,lmax
      read(ititle(l*20+1:l*20+2),*) nl
      if(nl(2:2) == 's') nl(2:2) = 'S'
      if(nl(2:2) == 'p') nl(2:2) = 'P'
      if(nl(2:2) == 'd') nl(2:2) = 'D'
      if(nl(2:2) == 'f') nl(2:2) = 'F'
      if(nl(2:2) == 'g') nl(2:2) = 'G'
      read(ititle(l*20+4:l*20+9),*) occup
      read(ititle(l*20+16:l*20+20),*) rc
      write(io,'(a2,2i3,4f6.2)') nl,l+1,l,occup,ZERO,rc,rc
    enddo
  endif

  return
end subroutine write_pwscf_upf_in
