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
!>  \date         22 April 2021, 23 January 2022.
!>  \copyright    GNU Public License v2

program cpwin_geom

! Added call to sym_space_group_name. 23 January 2022. JLM

  use esdf

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! dimensions

  integer                            ::  mxdtyp                          !<  array dimension of types of atoms
  integer                            ::  mxdatm                          !<  array dimension of number of atoms of a given type

! crystal variables

  real(REAL64)                       ::  alatt                           !<  lattice constant

  real(REAL64)                       ::  adot(3,3)                       !<  metric in direct space in atomic units (Bohr radius)

  integer                            ::  ntype                           !<  number of types of atoms
  integer, allocatable               ::  natom(:)                        !<  number of atoms of type i
  real(REAL64), allocatable          ::  rat(:,:,:)                      !<  lattice coordinates of atom j of type i
  real(REAL64), allocatable          ::  atmass(:)                       !<  atomic mass of atoms of type i
  character(len=2), allocatable      ::  nameat(:)                       !<  chemical symbol for the type i

  character(len=250)                 ::  meta_pwdat                      !<  metadata

! metric variables

  real(REAL64)                       ::  adotnig(3,3), adotsym(3,3)
  integer                            ::  ibravais, mtotal(3,3)
  real(REAL64)                       ::  avec(3,3)
  real(REAL64)                       ::  aconv(3,3), avecnig(3,3)

! space group

  integer                            ::  ntrans                          !<  number of symmetry operations in the factor group
  integer                            ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64)                       ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

  integer                            ::  code_group                      !  point group


! neighbour variables

  integer, allocatable               ::  nneighb(:,:)                    !<  number of neighbors for each atom
  integer, allocatable               ::  neighbtype(:,:,:)               !<  type of nth neighbors for each atom
  real(REAL64), allocatable          ::  rneighb(:,:,:,:)                !<  relative position of the n-th neighbor
  real(REAL64), allocatable          ::  wneighb(:,:,:)                  !<  strength of the n-th neighbor bond

! other variables

  integer               ::  ipr
  character(len=6)      ::  fcpwin                                       !<  filename with input data (new style)
  logical               ::  lgeom, lgeom2

! parameter

  real(REAL64), parameter ::  TOL = 1.0E-6_REAL64
  integer, parameter      ::  MXDNB = 40

  fcpwin = 'cpw.in'
  ipr = 3

  call size_mxdtyp_mxdatm_esdf(ipr, fcpwin, lgeom,                       &
              mxdtyp, mxdatm)

  if(lgeom) then

    allocate(natom(mxdtyp))
    allocate(rat(3,mxdatm,mxdtyp))
    allocate(atmass(mxdtyp))
    allocate(nameat(mxdtyp))

    call esdf_init(fcpwin)

    call read_esdf_crystal(ipr,                                          &
        adot, ntype, natom, nameat, rat,                                 &
        atmass, alatt, lgeom2,                                           &
        mxdtyp, mxdatm)

    call esdf_close

    if(lgeom2) then

      ipr = 1
      call metric_ident(adot, adotnig, adotsym, ibravais, mtotal,        &
                      avec, aconv, avecnig, TOL, ipr)

      call sym_identify(1, ipr, TOL,                                     &
          ntrans, mtrx, tnp,                                             &
          adot, ntype, natom, rat,                                       &
          mxdtyp, mxdatm)

      call sym_point_group_name(adot, 1, code_group, ntrans, mtrx)

      call sym_space_group_name(ibravais, code_group, ntrans, mtrx, tnp)

      allocate(nneighb(mxdatm,mxdtyp))
      allocate(neighbtype(MXDNB,mxdatm,mxdtyp))
      allocate(rneighb(3,MXDNB,mxdatm,mxdtyp))
      allocate(wneighb(MXDNB,mxdatm,mxdtyp))

      call voronoi_neighb(adot, ntype, natom, nameat, rat,               &
                       nneighb, neighbtype, rneighb, wneighb,            &
                       mxdtyp, mxdatm, MXDNB)

      call voronoi_neighb_print(adot, ntype, natom, nameat, rat,         &
                       nneighb, neighbtype, rneighb, wneighb,            &
                       mxdtyp, mxdatm, MXDNB)

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

      call write_cif(adot, ntype, natom, nameat, rat,                    &
                       mxdtyp, mxdatm)

      deallocate(nneighb)
      deallocate(neighbtype)
      deallocate(rneighb)

    else
      write(6,*)
      write(6,*) '  could not get enough information from ',fcpwin
      write(6,*) '  to figure out the geometry'
      write(6,*)
    endif

    deallocate(natom)
    deallocate(rat)
    deallocate(atmass)
    deallocate(nameat)

  else
    write(6,*)
    write(6,*) '  could not get enough information from ',fcpwin
    write(6,*) '  to find mxdtyp, mxdatm '
    write(6,*)
  endif

  stop
end program cpwin_geom
