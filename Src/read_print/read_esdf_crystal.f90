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

!>  reads the crystal geometry from a file
!>  opened by a previous edsf_init
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         2017-2021, 8 July 2025.
!>  \copyright    GNU Public License v2

subroutine read_esdf_crystal(ipr,                                        &
    adot, ntype, natom, nameat, rat,                                     &
    atmass, alatt, lgeom,                                                &
    mxdtyp, mxdatm)

! Written June 16, 2017. jlm
! Modified 6 June 2019, bug in atomic mass. JLM
! Modified, documentation, August 10 2019. JLM
! Modified, prints reciprocal lattice vectors, 16 April 2021. JLM
! Modified to allow "negative" volumes. 8 July 2025. JLM


  use esdf

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms

  integer, intent(in)                ::  ipr                             !<  print control: ipr < 3 no printing, except warnings.

! output


  real(REAL64), intent(out)          ::  adot(3,3)                       !<  metric in direct space
  integer, intent(out)               ::  ntype                           !<  number of types of atoms
  integer, intent(out)               ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(out)      ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(out)          ::  atmass(mxdtyp)                  !<  atomic mass (in a.u.) of atom of type i

  real(REAL64), intent(out)          ::  alatt                           !<  lattice constant
  logical, intent(out)               ::  lgeom                           !<  indicates if geometry was successfully read.

! local allocatable variables

  integer, allocatable       ::  indx(:)
  integer, allocatable       ::  iz(:)

! local variables

  integer                    ::  nlines
  real(REAL64)               ::  avec(3,3)
  real(REAL64)               ::  bvec(3,3), vcell
  real(REAL64)               ::  coor(3)

  logical           ::  lf                         !  an approximation was found
  logical           ::  lsq                        !  it is the square of a rati onal
  integer           ::  nnum, nden                 !  x ~ nnum/ndem or x**2 ~ nnum/nden and sign of nnum is the same as x
  integer           ::  ntry                       !  tries denominators up to ntry
  integer           ::  izz                        !  atomic number
  real(REAL64)      ::  tmp

  logical           ::  lperm                      !  for checking permutations
  integer           ::  ioerr

  integer           ::  nat_1, nat_2               !  total number of atoms

! parameters

  real(REAL64), parameter  :: AMU = 1822.888485_REAL64

  real(REAL64), parameter  ::  UM = 1.0_REAL64
  real(REAL64), parameter  ::  ZERO = 0.0_REAL64
  real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64

  real(REAL64), parameter  ::  EPS = 2.0E-6_REAL64                       !  Criteria for rounding the values (depends on the f12.6 in input)
  real(REAL64), parameter  ::  EPSMALL = 1.0E-14_REAL64                  !  Criteria for worthwhile rounding

! counters

  integer   ::  i, j, nt


! hardcoded value of ntry

  ntry = 36

! starts reading the geometry of the crystal

  lgeom = .TRUE.


! Lattice constant


  alatt = UM
  if(esdf_defined('LatticeConstant')) then
    alatt = esdf_physical('LatticeConstant',alatt,'bohr')
  else
    write(6,*)
    write(6,*) '   WARNING'
    write(6,*)
    write(6,*) '   Lattice constant assumed to be 1.000'
  endif

  if(alatt < 0.01) then
    write(6,*)
    write(6,*) '   WARNING'
    write(6,*)
    write(6,*) '   Lattice constant looks very small'
    write(6,'("    The value of alatt is: ",e12.3)') alatt
  endif

  if(alatt > 100.0) then
    write(6,*)
    write(6,*) '   WARNING'
    write(6,*)
    write(6,*)
    write(6,*) '   Lattice constant looks very large'
    write(6,'("    The value of alatt is: ",e12.3)') alatt
  endif


! lattice vectors

  if(esdf_block('LatticeVectors',nlines)) then
    if(nlines /= 3) then
      write(6,*)
      write(6,*) '  Stopped in read_esdf_crystal  '
      write(6,*) '  Wrong Number of Lines in LatticeVectors!'

      stop

    endif

    ioerr = 0
    do i=1,nlines
      read(block_data(i),*,iostat=ioerr) avec(1,i), avec(2,i), avec(3,i)

      if(ioerr /= 0) then
        write(6,*)
        write(6,*) '    Stopped in read_esdf_crystal'
        write(6,*) '    error reading LatticeVectors'

        stop

      endif

    enddo

!   checks if it is close to a rational number or the sqrt of a rational

    do i=1,3
    do j=1,3

      call near_rational(avec(j,i),lf,lsq,nnum,nden,ntry,5.0*EPS)

      if(lf) then
        if(lsq) then
          tmp = sign(sqrt(UM*nnum/(UM*nden)),avec(j,i))
        else
          tmp = UM*nnum/(UM*nden)
        endif
        if(abs(tmp - avec(j,i)) < EPS .and.                              &
           abs(tmp - avec(j,i)) > EPSMALL) then
          write(6,*)
          write(6,'("   WARNING:   avec(",i3,",",i3,")")') j,i
          write(6,'("   changed from ",f20.12," to ",f20.12)') avec(j,i),tmp
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
    write(6,*) '  Input Primitive Translation Vectors'
    write(6,'(24x,"in a.u.",28x,"in lattice units")')

    do j=1,3
      write(6,'("  a",i1,"=",3(1x,e13.6),5x,3(2x,f7.3))')                &
                 j,(avec(i,j),i=1,3),(avec(i,j)/alatt,i=1,3)
    enddo

!   calculate metric

    adot(1,1) = avec(1,1)*avec(1,1) + avec(2,1)*avec(2,1) + avec(3,1)*avec(3,1)
    adot(2,2) = avec(1,2)*avec(1,2) + avec(2,2)*avec(2,2) + avec(3,2)*avec(3,2)
    adot(3,3) = avec(1,3)*avec(1,3) + avec(2,3)*avec(2,3) + avec(3,3)*avec(3,3)
    adot(1,2) = avec(1,1)*avec(1,2) + avec(2,1)*avec(2,2) + avec(3,1)*avec(3,2)
    adot(1,3) = avec(1,1)*avec(1,3) + avec(2,1)*avec(2,3) + avec(3,1)*avec(3,3)
    adot(2,3) = avec(1,2)*avec(1,3) + avec(2,2)*avec(2,3) + avec(3,2)*avec(3,3)
    adot(2,1) = adot(1,2)
    adot(3,1) = adot(1,3)
    adot(3,2) = adot(2,3)

!   compute the reciprocal vectors and cell volume

    bvec(1,1) = avec(2,2)*avec(3,3) - avec(3,2)*avec(2,3)
    bvec(2,1) = avec(3,2)*avec(1,3) - avec(1,2)*avec(3,3)
    bvec(3,1) = avec(1,2)*avec(2,3) - avec(2,2)*avec(1,3)
    bvec(1,2) = avec(2,3)*avec(3,1) - avec(3,3)*avec(2,1)
    bvec(2,2) = avec(3,3)*avec(1,1) - avec(1,3)*avec(3,1)
    bvec(3,2) = avec(1,3)*avec(2,1) - avec(2,3)*avec(1,1)
    bvec(1,3) = avec(2,1)*avec(3,2) - avec(3,1)*avec(2,2)
    bvec(2,3) = avec(3,1)*avec(1,2) - avec(1,1)*avec(3,2)
    bvec(3,3) = avec(1,1)*avec(2,2) - avec(2,1)*avec(1,2)

!   cell volume

    vcell = bvec(1,1)*avec(1,1) + bvec(2,1)*avec(2,1) +                  &
            bvec(3,1)*avec(3,1)

    if(abs(vcell) < EPSMALL) then
      write(6,*)
      write(6,'("    STOPPED in read_esdf_crystal.    Cell volume ",     &
        & "squared= ",e12.4)') vcell

      stop

    endif

    if(vcell < ZERO) then
      write(6,*)
      write(6,'("    WARNING in read_esdf_crystal.    Negatively ",      &
        & "oriented lattice vectors ",e12.4)') vcell
      write(6,*)
    endif

    do j=1,3
      bvec(1,j) = 2*PI*bvec(1,j)/vcell
      bvec(2,j) = 2*PI*bvec(2,j)/vcell
      bvec(3,j) = 2*PI*bvec(3,j)/vcell
    enddo

    write(6,*)
    write(6,*) '  Input Reciprocal Lattice Vectors'
    write(6,'(24x,"in a.u.",28x,"in lattice units (2 pi / a)")')

    do j=1,3
      write(6,'("  b",i1,"=",3(1x,e13.6),5x,3(2x,f7.3))')                &
                 j,(bvec(i,j),i=1,3),(alatt*bvec(i,j)/(2*PI),i=1,3)
    enddo

    write(6,*)
    write(6,*) '  WARNING:  Inside the code other orientations of'
    write(6,*) '  lattice vectors may be used as they are arbitrary.'
    write(6,*) '  Internal orientation uses CONVENTIONS!'
    write(6,*)
    write(6,*) '  Always use lattice coordinates with this code...'
    write(6,*)

!   internal representation of lattice vectors

    call adot_to_avec_sym(adot,avec,bvec)

    write(6,*)
    write(6,*) '  Internal Primitive Translation Vectors'
    write(6,'(24x,"in a.u.",28x,"in lattice units")')

    do j=1,3
      write(6,'("  a",i1,"=",3(1x,e13.6),5x,3(2x,f7.3))')                &
                 j,(avec(i,j),i=1,3),(avec(i,j)/alatt,i=1,3)
    enddo

    write(6,*)
    write(6,*) '  Internal Reciprocal Lattice Vectors'
    write(6,'(24x,"in a.u.",28x,"in lattice units (2 pi / a)")')

    do j=1,3
      write(6,'("  b",i1,"=",3(1x,e13.6),5x,3(2x,f7.3))')                &
                 j,(bvec(i,j),i=1,3),(alatt*bvec(i,j)/(2*PI),i=1,3)
    enddo

  else
    lgeom = .FALSE.
  endif


! atomic species

  if(esdf_defined('NumberOfSpecies')) then

    ntype = esdf_integer('NumberOfSpecies',1)
    if(ntype > mxdtyp) then
      write(6,*)
      write(6,'("   Stopped in read_esdf_crystal:")')
      write(6,'("   The number of species ",i5," is larger than",        &
        & "  the allocated dimension",i5)') ntype,mxdtyp

      stop

    endif

    if(ntype < 1) then
      ntype = 0
      write(6,*)
      write(6,*) '   WARNING'
      write(6,*)
      write(6,*) '   Number of elements < 1.'
      write(6,*) '   Will be set by Chemical_Species_Label.'
    endif

    if(ntype > 120) then
      write(6,*)
      write(6,*) '   WARNING'
      write(6,*)
      write(6,*) '   Number of elements looks very high'
      write(6,'("    The value of ntype is: ",i8)') ntype
    endif

  else
    ntype = 0
  endif

  if(ntype == 0) then
    allocate(indx(100))
    allocate(iz(100))
  else
    allocate(indx(ntype))
    allocate(iz(ntype))
  endif

  if(esdf_block('Chemical_Species_Label',nlines)) then

    if(ntype == 0) then
      if(nlines > 100) then
        write(6,*)
        write(6,'("   Stopped in read_esdf_crystal:")')
        write(6,'("   Insuficient allocated space. Set ntype")')

        stop

      endif
      ntype = nlines
    endif

    if(nlines /= ntype) then
      write(6,*)
      write(6,*) '  Stopped in read_esdf_crystal.'
      write(6,*) '  Wrong Number of Linesin Chemical_Species_Label!'
      write(6,*) '  '

      stop

    endif

    ioerr = 0
    do i=1,nlines

      read(block_data(i),*,iostat=ioerr) indx(i),iz(i),nameat(indx(i))

      if(ioerr /= 0) then
        write(6,*)
        write(6,*) '    Stopped in read_esdf_crystal,'
        write(6,*) '    error reading Chemical_Species_Label.'

        stop

      endif

    enddo

!   checks if indx is a permutation

    call check_permutation(indx,ntype,lperm)

    if(.NOT. lperm) then
      write(6,*)
      write(6,*) '  Stopped in read_esdf_crystal.'
      write(6,*) '  Wrong species indexin Chemical_Species_Label!'

      stop

    endif

!   checks if the atomic numbers are correct

    do i = 1,ntype
      call p_tbl_charge(nameat(indx(i)),izz)
      if(izz /= iz(indx(i))) then
        write(6,*)
        write(6,*) '   WARNING   '
        write(6,*)
        write(6,'("   For atom type ",i4," element ",a2," using ",       &
          & "atomic number ",i3," instead of input value ",i3)')         &
          indx(i),nameat(indx(i)),izz,iz(indx(i))
      endif
    enddo

!   default atomic massas

    do nt = 1,ntype
      call p_tbl_mass(nameat(nt),atmass(nt))
    enddo

  else
    lgeom = .FALSE.
  endif

  deallocate(indx)
  deallocate(iz)

! atomic positions

  if(esdf_block('AtomicCoordinatesAndAtomicSpecies',nlines)) then

    do nt = 1,ntype
      natom(nt) = 0
    enddo

    ioerr = 0
    do i=1,nlines

      read(block_data(i),*,iostat=ioerr) (coor(j),j=1,3),nt

      if(ioerr /= 0) then
        write(6,*)
        write(6,*) '    Stopped in read_esdf_crystal'
        write(6,*) '    error reading AtomicCoordinates'
        write(6,*) '    in line ',i,' of that block'

        stop

      endif

      if(nt > mxdtyp) then
        write(6,*)
        write(6,*)   '   STOPPED in read_esdf_crystal'
        write(6,'("  There is an atom of type ",i5," but maximum",       &
           & " allowed value is ",i5)') nt,mxdtyp

        stop

      endif
      natom(nt) = natom(nt) + 1
      if(natom(nt) > mxdatm) then
        write(6,*)
        write(6,*)   '   STOPPED in read_esdf_crystal'
        write(6,'("  There are at least ",i5," atoms of ",a2,            &
          & " but mxdatm is ",i5)') natom(nt),nameat(nt),mxdatm

        stop

      endif

      rat(1,natom(nt),nt) = coor(1)
      rat(2,natom(nt),nt) = coor(2)
      rat(3,natom(nt),nt) = coor(3)

    enddo

  else
    lgeom = .FALSE.
  endif

! checks consistency

  if(esdf_defined('NumberOfAtoms')) then
    nat_1 = esdf_integer('NumberOfAtoms',1)
  endif

  nat_2 = 0
  do nt = 1,ntype
    if(natom(nt) == 0) then
      write(6,*)
      write(6,*) '    Stopped in read_esdf_crystal.'
      write(6,*) '    Could not find coordinates for atom type ', nameat(nt)
      write(6,*)

      stop

    endif
    nat_2 = nat_2 + natom(nt)
  enddo

  if(nat_1 /= nat_2) then
    write(6,*)
    write(6,*) '    Stopped in read_esdf_crystal.'
    write(6,*) '    Number of atoms found: ', nat_2
    write(6,*) '    Expected nymber of atoms: ',nat_1
    write(6,*)

    stop

  endif


! prints the results

  if(ipr > 2) then
    write(6,*)
    write(6,*)
    write(6,*)  '  The values set by read_esdf_crystal are:'
    write(6,*)
    write(6,'("   The value of alatt is: ",e12.3)') alatt
    write(6,'("   The value of ntype is: ",i5)') ntype
    write(6,'("   The value of natom, nameat, atmass")')
    write(6,'("   and respective atom coordinates (rat) are")')
    do nt = 1,ntype
      write(6,'(3x,i6,3x,a2,3x,f12.3)') natom(nt), nameat(nt), atmass(nt)
      do i = 1,natom(nt)
        write(6,'(3(3x,f12.3))') (rat(j,i,nt),j=1,3)
      enddo
    enddo
    write(6,*)
    write(6,*)
  endif

! convert mass

  do nt = 1,ntype
    atmass(nt) = AMU*atmass(nt)
  enddo

  return
end subroutine read_esdf_crystal
