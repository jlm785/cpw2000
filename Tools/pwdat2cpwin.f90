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

!>  Converts the old PW.DAT file format to the cpw.in format.
!>  The result is written in cpw_from_pwdat.in.
!>  If a cpw.in file (old style, without the crystal structure)
!>  is present, the type of run, the Keating correction and
!>  the TBL choice of exchange and correlation are taken from it,
!>  the other options have to be copied by hand.
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         October 2018, 18 October 2026.
!>  \copyright    GNU Public License v2

program pwdat2cpwin

! Written October 2018. JLM
! New API of write_cpwout, options from an old cpw.in, size_mxdtyp_mxdatm
!   and read_data moved here from the library. 18 October 2026. JLM+claude

  use esdf

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! dimensions

  integer                            ::  mxdtyp                          !  array dimension of types of atoms
  integer                            ::  mxdatm                          !  array dimension of number of atoms of a given type

! crystal structure

  real(REAL64)                       ::  alatt                           !  lattice constant
  real(REAL64)                       ::  adot(3,3)                       !  metric in direct space in atomic units (Bohr radius)
  integer                            ::  ntype                           !  number of types of atoms
  integer, allocatable               ::  natom(:)                        !  number of atoms of type i
  real(REAL64), allocatable          ::  rat(:,:,:)                      !  lattice coordinates of atom j of type i
  real(REAL64), allocatable          ::  atmass(:)                       !  atomic mass of atoms of type i
  character(len=2), allocatable      ::  nameat(:)                       !  chemical symbol for the type i

! cutoff, bands and k-points

  real(REAL64)                       ::  emax                            !  kinetic energy cutoff of plane wave expansion (Hartree).
  integer                            ::  nbandin                         !  target for number of bands
  integer                            ::  nx, ny, nz                      !  size of the integration mesh in k-space (nx*ny*nz)
  real(REAL64)                       ::  sx, sy, sz                      !  offset of the integration mesh (usually 0.5)

! options

  character(len=6)                   ::  flgcal                          !  type of calculation
  logical                            ::  lkeat                           !  Keating correction
  character(len=4)                   ::  author                          !  exchange and correlation
  logical                            ::  ltbl                            !  TBL exchange and correlation

! other variables

  character(len=250)                 ::  meta_pwdat                      !  metadata from PW.DAT
  character(len=20)                  ::  fileout                         !  name of output file
  integer                            ::  ioerr                           !  input/output error
  logical                            ::  lcpwin                          !  an old cpw.in exists


  call size_mxdtyp_mxdatm('PW.DAT', 21, mxdtyp, mxdatm)

  open(unit = 5, file = 'PW.DAT', status = 'OLD', form = 'FORMATTED')

  read(5,'(a250)',iostat=ioerr) meta_pwdat
  if(ioerr == 0) then
    write(6,*) meta_pwdat(1:60)
    write(6,*) meta_pwdat(61:110)
    write(6,*) meta_pwdat(111:250)
  else
    write(6,*) '  pwdat2cpwin:   Unable to read titles'
  endif
  write(6,*)
  write(6,*)

  allocate(natom(mxdtyp))
  allocate(rat(3,mxdatm,mxdtyp))
  allocate(atmass(mxdtyp))
  allocate(nameat(mxdtyp))

! gets crystal data

  call read_data(adot, ntype, natom, nameat, rat, atmass, alatt,         &
      mxdtyp, mxdatm)

  read(5,'(1x,f9.4)') emax

  read(5,'(i5,15x,3i5,5x,3f5.2)') nbandin, nx, ny, nz, sx, sy, sz

  close(unit = 5)

! options from an old cpw.in

  flgcal = 'ONE   '
  lkeat = .FALSE.
  author = 'CA  '

  inquire(file = 'cpw.in', exist = lcpwin)

  if(lcpwin) then

    write(6,*)
    write(6,*) '  Type of run, Keating and TBL options from cpw.in'
    write(6,*)

    call esdf_init('cpw.in')

    flgcal = esdf_string('MD.TypeOfRun', flgcal)
    call chrcap(flgcal,6)
    lkeat = esdf_boolean('MD.UseKeatingCorrections', lkeat)
    author = esdf_string('XC.Authors', author)
    call chrcap(author,4)

    call esdf_close

  endif

  ltbl = (author == 'TBL ')

  fileout = 'cpw_from_pwdat.in'

  call write_cpwout(21, fileout, meta_pwdat, flgcal, 'pwdat2cpwin  ',    &
      adot, ntype, natom, nameat, rat, atmass, alatt,                    &
      emax, nbandin, nx, ny, nz, sx, sy, sz,                             &
      lkeat, ltbl,                                                       &
      mxdtyp, mxdatm)

  write(6,*)
  write(6,*) '  Wrote file ', trim(fileout)
  write(6,*)

  deallocate(natom)
  deallocate(rat)
  deallocate(atmass)
  deallocate(nameat)

end program pwdat2cpwin


!>  Finds the array sizes for the atomic structure, PW.DAT version.
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         2 May 2014, January 2020.
!>  \copyright    GNU Public License v2

subroutine size_mxdtyp_mxdatm(filename, iopw, mxdtyp, mxdatm)

! written 2 May 2014. JLM
! Modified, documentation, ioerr, January 2020. JLM
! Indentation. 28 September 2026. JLM+claude
! Moved from the library to Tools/pwdat2cpwin.f90. 18 October 2026. JLM+claude


  implicit none

! input


  character(len=*), intent(in)       ::  filename                        !<  file to be written
  integer, intent(in)                ::  iopw                            !<  tape number

! output

  integer, intent(out)               ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(out)               ::  mxdatm                          !<  array dimension of types of atoms

! local variables

  integer     ::  ntype, natom
  integer     ::  ioerr

! counters

  integer     ::  nt, j


  open(unit=iopw,file=filename,form='formatted',status='old',            &
        IOSTAT=ioerr)


  if(ioerr /= 0) then
    write(6,*)
    write(6,*)
    write(6,*) '  stopped in size_mxdtyp_mxdatm'
    write(6,*) '  Unable to open ',filename,' file'

    stop

  endif

! reads title

  read(iopw,*)

! read scale factor

  read(iopw,*)

! read the basis vectors

  read(iopw,*)
  read(iopw,*)
  read(iopw,*)

! read the number of atoms of each type and their names.

  read(iopw,*) ntype
  mxdtyp = ntype
  mxdatm = 0
  do nt=1,ntype
    read(iopw,*) natom
    mxdatm = max(mxdatm,natom)
    do j=1,natom
      read(iopw,*)
    enddo
  enddo

  close(unit=iopw)

  return

end subroutine size_mxdtyp_mxdatm


!>  Reads the crystal structure from a data file
!>  compatible with Berkeley/Sverre data files.
!>  Hence the complicated fixed format structure.
!>
!>  \author       Sverre Froyen, José Luís Martins, Carlos Loia Reis
!>  \version      5.13
!>  \date         15 October 1993, 8 December 2016.
!>  \copyright    GNU Public License v2

subroutine read_data(adot, ntype, natom, nameat, rat, atmass, alatt,     &
    mxdtyp, mxdatm)

! Adapted from Sverre Froyen plane wave program
! Written 15 October 1993. jlm
! Modified 15 January 1999. jlm
! Modified 2013/2014. CLR
! Modified,f90, near_rational, 30 May 2014. JLM
! Modified, near_rational bug trap, 21 June 2014. JLM
! Modified, ntry, printing, STRUCT, December 8, 2016. JLM
! Modified, documentation, August 2019. JLM
! Indentation. 28 September 2026. JLM+claude
! Moved from the library to Tools/pwdat2cpwin.f90. 18 October 2026. JLM+claude


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms

! output

  real(REAL64), intent(out)          ::  adot(3,3)                       !<  metric in direct space
  integer, intent(out)               ::  ntype                           !<  number of types of atoms
  integer, intent(out)               ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(out)      ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(out)          ::  atmass(mxdtyp)                  !<  atomic mass (in a.u.) of atom of type i

  real(REAL64), intent(out)          ::  alatt                           !<  lattice constant

! other variables

  real(REAL64)      ::  avec(3,3)

  character(len=1)  ::  iop(3,3),jop(3)
  real(REAL64)      ::  car(3)
  real(REAL64)      ::  sgn
  real(REAL64)      ::  tmp
  real(REAL64)      ::  difmax

  logical           ::  lf                          !  an approximation was found
  logical           ::  lsq                         !  it is the square of a rati onal
  integer           ::  nnum, nden                      !  x ~ nnum/ndem or x**2 ~ nnum/nden and sign of nnum is the same as x
  integer           ::  ntry                            !  tries denominators up to ntry


! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  :: AMU = 1822.888485_REAL64
  real(REAL64), parameter  :: EPS = 2.0E-6                               !  Criteria for rounding the values (depends on the f12.6 in input)
  real(REAL64), parameter  :: EPSMALL = 1.0E-14                          !  Criteria for worthwhile rounding

! counters

  integer    ::  i, j, ja, jmax, nt, ntt

! hardcoded value of ntry

  ntry = 36

! read scale factor

  read(5,'(f12.6)') alatt

  write(6,*)
  write(6,*)
  write(6,*)'  Crystal Structure:'
  write(6,*)
  write(6,*)
  write(6,'(4x,f14.6,"   reference lattice constant")') alatt
  write(6,*)

! read the basis vectors for the lattice (cartesian coordinates and
! atomic units divided by ascale). modify according to iop...
! multiply by alatt.

  do j=1,3
    read(5,'(3(2x,a1,f12.6))') (iop(i,j),avec(i,j),i=1,3)
  enddo

  do j=1,3
  do i=1,3
    sgn = UM
    if (avec(i,j) < ZERO) sgn = -sgn
    avec(i,j) = abs(avec(i,j))
    if (iop(i,j)  ==  'S') avec(i,j) = sqrt(avec(i,j))
    if (iop(i,j)  ==  'C') avec(i,j) = avec(i,j)**(UM/3)
    if (iop(i,j)  ==  'T') avec(i,j) = avec(i,j)/3
    if (iop(i,j)  ==  'H') avec(i,j) = avec(i,j)/sqrt(3*UM)
    avec(i,j) = sgn*avec(i,j)
  enddo
  enddo

!      checks if it is close to a rational number or the sqrt of a rational

  do i=1,3
  do j=1,3

    call near_rational(avec(j,i), lf, lsq, nnum, nden, ntry, 5.0*EPS)

    if(lf) then
      if(lsq) then
        tmp = sign(sqrt(UM*nnum/(UM*nden)),avec(j,i))
      else
        tmp = UM*nnum/(UM*nden)
      endif
      if(abs(tmp - avec(j,i)) < EPS .and.                                &
         abs(tmp - avec(j,i)) > EPSMALL) then
        write(6,*)
        write(6,'("   WARNING:   avec(",i3,",",i3,")")') j,i
        write(6,'("   changed from ",f20.12," to ",f20.12)')             &
                 avec(j,i),tmp
        write(6,*)
        avec(j,i) = tmp
      endif
    endif

  enddo
  enddo

  do j=1,3
  do i=1,3
    avec(i,j) = alatt*avec(i,j)
  enddo
  enddo

  write(6,*)
  write(6,*) '  Primitive Translation Vectors'
  write(6,'(24x,"in a.u.",28x,"in lattice units")')

  do j=1,3
    write(6,'("  a",i1,"=",3(1x,e13.6),5x,3(2x,f7.3))')                  &
                 j,(avec(i,j),i=1,3),(avec(i,j)/alatt,i=1,3)
  enddo

! read the number of atoms of each type and their names.
! read atomic positions in unit cell (in units of
! basis vectors). modify according to iop.

  read(5,'(i5)') ntype

  if(ntype > mxdtyp .or. ntype <= 0) then
    write(6,*)
    write(6,'("    STOPPED  in read_data:   ",                      &
        &        "  increase mxdtyp to ",i6)') ntype

    stop

  endif

  do nt=1,ntype

    read(5,'(i5,3x,a2)') natom(nt),nameat(nt)

    if(natom(nt) > mxdatm .or. natom(nt) <= 0) then
      write(6,*)
      write(6,'("    STOPPED in read_data:   ",                     &
          &        "  increase mxdatm to ",i6)') natom(nt)

      stop

    endif

    call p_tbl_mass(nameat(nt), atmass(nt))

    jmax = natom(nt)
    do ja=1,jmax
      read(5,'(3(2x,a1,f12.6))') (jop(i),rat(i,ja,nt),i=1,3)
      do i=1,3
        sgn = UM
        if (rat(i,ja,nt) < ZERO) sgn = -sgn
        rat(i,ja,nt) = abs(rat(i,ja,nt))
        if (jop(i)  ==  'S') rat(i,ja,nt) = sqrt(rat(i,ja,nt))
        if (jop(i)  ==  'T')  rat(i,ja,nt) = rat(i,ja,nt)/3
        rat(i,ja,nt) = sgn*rat(i,ja,nt)
      enddo
    enddo
  enddo

! checks if it is close to a rational number or the sqrt of a rational

  difmax = ZERO

  do nt=1,ntype
  do ja=1,natom(nt)
  do i=1,3

     call near_rational(rat(i,ja,nt), lf, lsq, nnum, nden, ntry, 5.0*EPS)

    if(lf) then
      if(lsq) then
        tmp = sign(sqrt(UM*nnum/(UM*nden)),rat(i,ja,nt))
      else
        tmp = UM*nnum/(UM*nden)
      endif
      if(abs(tmp - rat(i,ja,nt)) < EPS .and.                             &
         abs(tmp - rat(i,ja,nt)) > EPSMALL) then
!        write(6,*)
!        write(6,'("   WARNING:    rat(",i3,",",i3,",",i3,")")')     &
!     &                i,ja,nt
!        write(6,'("   changed from ",f20.12," to ",f20.12)')        &
!     &              rat(i,ja,nt),tmp
!        write(6,*)
        difmax = max(difmax,abs(tmp - rat(i,ja,nt)))
        rat(i,ja,nt) = tmp
      endif
    endif

  enddo
  enddo
  enddo

  if(difmax > sqrt(EPS*EPSMALL)) then
        write(6,*)
        write(6,'("   WARNING:    atomic positions changed in ",    &
            &           "read_data")')
        write(6,'("   Maximum correction ",f20.12)') difmax
        write(6,*)
  endif

! printout

  write(6,*)
  write(6,'("  No. Type     Position(lattice coord.)",10x,          &
      &          "Position(cartesian coord.)",15x,"Mass")')
  write(6,*)

  ntt = 0
  do nt=1,ntype
    jmax = natom(nt)
    do ja=1,jmax
      ntt = ntt + 1
      car(1) = avec(1,1)*rat(1,ja,nt) + avec(1,2)*rat(2,ja,nt) +         &
               avec(1,3)*rat(3,ja,nt)
      car(2) = avec(2,1)*rat(1,ja,nt) + avec(2,2)*rat(2,ja,nt) +         &
               avec(2,3)*rat(3,ja,nt)
      car(3) = avec(3,1)*rat(1,ja,nt) + avec(3,2)*rat(2,ja,nt) +         &
               avec(3,3)*rat(3,ja,nt)

      write(6,'(1x,i3,3x,a2,2x,3(2x,f7.3),3x,3(1x,e13.6),5x,f8.3)')      &
                ntt,nameat(nt),(rat(i,ja,nt),i=1,3),                     &
                (car(i),i=1,3),atmass(nt)

    enddo
  enddo

! convert mass

  do nt=1,ntype
    atmass(nt) = AMU*atmass(nt)
  enddo

! calculate metric

  adot(1,1) = avec(1,1)*avec(1,1) + avec(2,1)*avec(2,1) +                &
              avec(3,1)*avec(3,1)
  adot(2,2) = avec(1,2)*avec(1,2) + avec(2,2)*avec(2,2) +                &
              avec(3,2)*avec(3,2)
  adot(3,3) = avec(1,3)*avec(1,3) + avec(2,3)*avec(2,3) +                &
              avec(3,3)*avec(3,3)
  adot(1,2) = avec(1,1)*avec(1,2) + avec(2,1)*avec(2,2) +                &
              avec(3,1)*avec(3,2)
  adot(1,3) = avec(1,1)*avec(1,3) + avec(2,1)*avec(2,3) +                &
              avec(3,1)*avec(3,3)
  adot(2,3) = avec(1,2)*avec(1,3) + avec(2,2)*avec(2,3) +                &
              avec(3,2)*avec(3,3)
  adot(2,1) = adot(1,2)
  adot(3,1) = adot(1,3)
  adot(3,2) = adot(2,3)

!  lfile = .FALSE.
!  inquire(file='STRUCT.DAT',exist=lfile)
!  if (lfile) then
!    open(unit=222,file='STRUCT.DAT',form='unformatted')
!    write(6,*) ' In read_data str lfile=', lfile
!  endif
!  if(lfile) then
!    write(6,*) ' Found struct.dat file!'
!    write(6,*) ' Metric and atom positions will be overrriden.'
!    write(6,*)
!    do i=1,3
!      read(222) (adot(ja,i),ja=1,3)
!    enddo
!    read(222) ntype
!    do nt=1,ntype
!      read(222) natom(nt)
!      do i=1,natom(nt)
!        read(222)(rat(k,i,nt),k=1,3)
!      enddo
!    enddo
!    close(unit = 222)
!  endif


  return

end subroutine read_data


