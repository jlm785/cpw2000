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

!>  Performs the analysis of the crystal geometry
!>  from the cpw.in file
!>
!>  \author       Jose Luis Martins
!>  \version      5.04
!>  \date         22 April 2021. 23 January 2022.
!>  \copyright    GNU Public License v2

subroutine voronoi_sub(ioreplay)

! Added call to sym_space_group_name. 23 January 2022. JLM

  use esdf

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

! dimensions

  integer                            ::  mxdtyp                          !  array dimension of types of atoms
  integer                            ::  mxdatm                          !  array dimension of number of atoms of a given type

! crystal variables

  real(REAL64)                       ::  adot(3,3)                       !  metric in direct space in atomic units (Bohr radius)

  integer                            ::  ntype                           !  number of types of atoms
  integer, allocatable               ::  natom(:)                        !  number of atoms of type i
  real(REAL64), allocatable          ::  rat(:,:,:)                      !  lattice coordinates of atom j of type i
  real(REAL64), allocatable          ::  atmass(:)                       !  atomic mass of atoms of type i
  character(len=2), allocatable      ::  nameat(:)                       !  chemical symbol for the type i

! metric variables

  real(REAL64)          ::  adotnig(3,3), adotsym(3,3)
  integer               ::  ibravais, mtotal(3,3)
  real(REAL64)          ::  avec(3,3)
  real(REAL64)          ::  aconv(3,3), avecnig(3,3)

! space group

  integer                            ::  ntrans                          !<  number of symmetry operations in the factor group
  integer                            ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64)                       ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! neighbour variables

  integer, allocatable               ::  nneighb(:,:)                    !<  number of neighbors for each atom
  integer, allocatable               ::  neighbtype(:,:,:)               !<  type of nth neighbors for each atom
  real(REAL64), allocatable          ::  rneighb(:,:,:,:)                !<  relative position of the n-th neighbor
  real(REAL64), allocatable          ::  wneighb(:,:,:)                  !<  strength of the n-th neighbor bond

! other variables

  integer                 ::  ipr, iprmet

  character(len=60)       ::  filename                                   !
  integer                 ::  io

  character(len=1)        ::  yesno

  logical                 ::  lmetric, lspace

! unused variables

  integer        ::  mxdgvein, mxdnstin, mxdlqp, mxdlao

! parameter

  real(REAL64), parameter ::  TOL = 1.0E-6_REAL64
  integer, parameter      ::  MXDNB = 40

! counters

  integer        ::  i, j, k

  filename = 'PW_RHO_V.DAT'
  io = 21

  call pw_rho_v_in_size(filename, io,                                    &
     mxdtyp, mxdatm, mxdgvein, mxdnstin, mxdlqp, mxdlao)

  allocate(natom(mxdtyp))
  allocate(rat(3,mxdatm,mxdtyp))
  allocate(atmass(mxdtyp))
  allocate(nameat(mxdtyp))

! copyed/adapted from pw_rho_v_in_crystal_calc

  open(unit=io,file=trim(filename),status='old',form='UNFORMATTED')

  read(io) ntype
  if(ntype > mxdtyp) then
    write(6,*)
    write(6,'("  STOPPED in voronoi_sub")')
    write(6,'("  ntype = ",i8," is greater than mxdtyp =",i8)')          &
         ntype,mxdtyp

    stop

  endif

  read(io) (natom(i),i=1,ntype)
  do i=1,ntype
    if(natom(i) > mxdatm) then
      write(6,*)
      write(6,'("  STOPPED in voronoi_sub")')
      write(6,'("  for atom type ",i4," natom = ",i8," is ",             &
        & "greater than mxdatm =",i8)') i,natom(i),mxdatm

      stop

    endif
  enddo

  read(io)                                                               !  bdate
  read(io)                                                               !  author
  read(io)                                                               !  emax
  read(io)                                                               !  meta

  read(io) ((adot(i,j),i=1,3),j=1,3)

  read(io) (nameat(i),i=1,ntype)

  do i=1,ntype
    read(io) ((rat(j,k,i),j=1,3),k=1,natom(i))
  enddo

  write(6,*)
  write(6,'("  Analysis of the crystal structure ")')
  write(6,*)


  write(6,*)
  write(6,*) '  Do you want information on the crystal lattice?'
  read(5,*) yesno
  write(ioreplay,*) yesno,'   crystal lattice'

  if(yesno == 'y' .or. yesno == 'Y') then

    iprmet = 1
    call metric_ident(adot, adotnig, adotsym, ibravais, mtotal,          &
                  avec, aconv, avecnig, TOL, iprmet)

    lmetric = .TRUE.

  else
    lmetric = .FALSE.
  endif


  write(6,*)
  write(6,*) '  Do you want a list of symmetry operations'
  read(5,*) yesno
  write(ioreplay,*) yesno,'   symmetry operations'

  if(yesno == 'y' .or. yesno == 'Y') then

    ipr = 1
    call sym_identify(1, ipr, TOL,                                       &
          ntrans, mtrx, tnp,                                             &
          adot, ntype, natom, rat,                                       &
          mxdtyp, mxdatm)

    lspace = .TRUE.

  else
    lspace = .FALSE.
  endif


  write(6,*)
  write(6,*) '  Do you want a list of possible space group names'
  read(5,*) yesno
  write(ioreplay,*) yesno,'   space group names'

  if(yesno == 'y' .or. yesno == 'Y') then

    if(.not. lmetric) then
      iprmet = 0
      call metric_ident(adot, adotnig, adotsym, ibravais, mtotal,        &
                  avec, aconv, avecnig, TOL, iprmet)
    endif
    if(.not. lspace) then
      ipr = 0
      call sym_identify(1, ipr, TOL,                                     &
          ntrans, mtrx, tnp,                                             &
          adot, ntype, natom, rat,                                       &
          mxdtyp, mxdatm)
    endif

    call sym_space_group_name(ibravais, ntrans, mtrx, tnp)

  endif


  write(6,*)
  write(6,*) '  Do you want a list of atomic neighbors'
  read(5,*) yesno
  write(ioreplay,*) yesno,'   atomic neighbors'

  if(yesno == 'y' .or. yesno == 'Y') then

    allocate(nneighb(mxdatm,mxdtyp))
    allocate(neighbtype(MXDNB,mxdatm,mxdtyp))
    allocate(rneighb(3,MXDNB,mxdatm,mxdtyp))
    allocate(wneighb(MXDNB,mxdatm,mxdtyp))

    call voronoi_neighb(adot, ntype, natom, nameat, rat,                 &
                     nneighb, neighbtype, rneighb, wneighb,              &
                     mxdtyp, mxdatm, MXDNB)

    call voronoi_neighb_print(adot, ntype, natom, nameat, rat,           &
                     nneighb, neighbtype, rneighb, wneighb,              &
                     mxdtyp, mxdatm, MXDNB)

    deallocate(nneighb)
    deallocate(neighbtype)
    deallocate(rneighb)
    deallocate(wneighb)

  endif


  write(6,*)
  write(6,*) '  Do you want a CIF file with crystal structure?'
  read(5,*) yesno
  write(ioreplay,*) yesno,'   atomic neighbors'

  if(yesno == 'y' .or. yesno == 'Y') then

   write(6,*)
   write(6,*) '  Preparing a CIF file from crystal data.'
   write(6,*) '  Uploading it to FINDSYM web site'
   write(6,*) '  (https://stokes.byu.edu/iso/findsym.php)'
   write(6,*) '  Will identify the space group of the crystal'
   write(6,*) '  and you can download a CIF file with that information.'
   write(6,*) '  The downloaded file can be converted to input for other'
   write(6,*) '  electronic structure codes with cif2cell.'
   write(6,*) '  (https://pypi.org/project/cif2cell/)'
   write(6,*)

   call write_cif(adot, ntype, natom, nameat, rat,                       &
                       mxdtyp, mxdatm)

  endif

  close(unit=io)

  deallocate(natom)
  deallocate(rat)
  deallocate(atmass)
  deallocate(nameat)

  return
end subroutine voronoi_sub
