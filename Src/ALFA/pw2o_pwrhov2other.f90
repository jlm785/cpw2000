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

!>  Reads the  PW_RHO_V.DAT file from cpw and
!>  writes the input file for pwSCF (quantum espresso) and abinit
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         12 October 2018. 25 February 2025
!>  \copyright    GNU Public License v2

subroutine pw2o_pwrhov2other(ioreplay)

! Written October 12, 2018 from in_rho_v.f90. JLM
! Modernized, February 12, 2021. JLM
! Removed unused variables, 14 November 2024. JLM
! Converted to a subroutine called by cpw_post_process. 25 February 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input


  integer, intent(in)                     ::  ioreplay                   !<  tape number for reproducing calculations

! atomic structure variables

  real(REAL64)                       ::  adot(3,3)                       !  metric in direct space
  integer                            ::  ntype                           !  number of types of atoms
  integer,allocatable                ::  natom(:)                        !  number of atoms of type i
  character(len=2),allocatable       ::  nameat(:)                       !  chemical symbol for the type i
  real(REAL64),allocatable           ::  rat(:,:,:)                      !  k-th component (in lattice coordinates) of the position of the n-th atom of type i

! information about the calculation

  character(len=4)                   ::  author                          !  type of xc wanted (CA=PZ , PW92 , PBE)
  character(len=6)                   ::  flgscf                          !  type of self consistent field and diagonalizatioN
  character(len=4)                   ::  flgdalin                        !  whether the dual approximation is used
  real(REAL64)                       ::  emaxin                          !  kinetic energy cutoff of plane wave expansion (hartree).
  real(REAL64)                       ::  teleck                          !  electronic temperature (in kelvin)

  character(len=60)                  ::  pwline                          !  identifier of the calculation.  May contain miscellaneous information!
  character(len=50)                  ::  title                           !  title for plots
  character(len=140)                 ::  subtitle                        !  title for plots
  character(len=250)                 ::  meta_cpw2000                    !  metadata from cpw2000

  integer                            ::  nbandin                         !  target for number of bands
  integer                            ::  nx,ny,nz                        !  divisions of Brillouin zone for integration (Monkhorst-Pack)
  real(REAL64)                       ::  sx,sy,sz                        !  shift of points in division of Brillouin zone for integration (Monkhorst-Pack)
  real(REAL64)                       ::  alatt                           !  lattice constant
  real(REAL64)                       ::  efermi                          !  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

! Pseudopotential variables

  integer, allocatable               ::  nqnl(:)                         !  number of points for pseudo interpolation for atom k
  real(REAL64), allocatable          ::  delqnl(:)                       !  step used in the pseudo interpolation for atom k
  real(REAL64), allocatable          ::  vkb(:,:,:,:)                    !  kb nonlocal pseudo. for atom k, ang. mom. l. (not normalized to vcell, hartree)
  integer, allocatable               ::  nkb(:,:,:)                      !   kb pseudo.  normalization for atom k, ang. mom. l
  real(REAL64), allocatable          ::  vloc  (:,:)                     !  local pseudopotential for atom k (hartree)
  real(REAL64), allocatable          ::  dcor(:,:)                       !  core charge density for atom k
  real(REAL64), allocatable          ::  dval (:,:)                      !  valence charge density for atom k
  integer, allocatable               ::  norbat(:)                       !  number of atomic orbitals for atom k
  integer, allocatable               ::  nqwf(:)                         !  number of points for wavefunction interpolation for atom k
  real(REAL64), allocatable          ::  delqwf(:)                       !  step used in the wavefunction interpolation for atom k
  integer, allocatable               ::  lorb(:,:)                       !  angular momentum of orbital n of atom k
  real(REAL64), allocatable          ::  wvfao(:,:,:)                    !  wavefunction for atom k, ang. mom. l
  logical                            ::  latorb                          !  indicates if all atoms have information about atomic orbitals
  real(REAL64), allocatable          ::  zv(:)                           !  valence of atom with type i
  real(REAL64)                       ::  ztot                            !  total charge density (electrons/cell)

  character(len=3), allocatable      ::  irel(:)                         !  type of calculation relativistic/spin
  character(len=4), allocatable      ::  icore(:)                        !  type of partial core correction
  character(len=2), allocatable      ::  icorr(:)                        !  type of correlation
  character(len=60), allocatable     ::  iray(:)                         !  information about pseudopotential
  character(len=10), allocatable     ::  psdtitle(:,:)                   !  further information about pseudopotential

! input effective potential and charge density

  integer                            ::  ngin                            !  input value of ng
  integer                            ::  kmaxin(3)                       !  max value of kgvin(i,n)
  integer, allocatable               ::  kgvin(:,:)                      !  input values of kgv
  complex(REAL64), allocatable       ::  phasein(:)                      !  phase factor of G-vector n
  real(REAL64), allocatable          ::  conjin(:)                       !  is -1 if one must take the complex conjugate of x*phase
  integer                            ::  nsin                            !  input value of ns
  integer, allocatable               ::  mstarin(:)                      !  input values of kgv
  complex(REAL64), allocatable       ::  veffin(:)                       !  input effective potential (local+Hartree+Xc) for the prototype g-vector
  complex(REAL64), allocatable       ::  denin(:)                        !  input charge density for the prototype g-vector
  complex(REAL64), allocatable       ::  denbondin(:)                    !  input bonding charge density for the prototype g-vector

! space group information


  integer                            ::  ntrans                          !  number of symmetry operations in the factor group
  integer                            ::  mtrx(3,3,48)                    !  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation
  real(REAL64)                       ::  tnp(3,48)                       !  2*pi* i-th component (in lattice coordinates) of the fractional translation vector of the k-th symmetry operation

! dimensions

  integer                            ::  mxdtyp                          !  array dimension of types of atoms
  integer                            ::  mxdatm                          !  array dimension of number of atoms of a given type

  integer                            ::  mxdgvein                        !  array dimension for input g-space vectors
  integer                            ::  mxdnstin                        !  array dimension for input g-space stars

  integer                            ::  mxdlqp                          !  array dimension for local potential
  integer                            ::  mxdlao                          !  array dimension of orbital per atom type

! other variables

  integer                 ::  ipr
  real(REAL64)            ::  ealraw
  integer                 ::  iotape
  real(REAL64)            ::  emax

  character(len=60)       ::  filename

  logical                 ::  lso                                        !  presence of relativistic pseudopotentials
  character(len=1)        ::  yesno

  character(len=30)       ::  fileband

! counters

  integer       ::  nt


! open file and reads data

  filename = 'PW_RHO_V.DAT'
  iotape = 21

  fileband = 'BAND_LINES.DAT'

  call pw_rho_v_in_size(filename, iotape,                                &
         mxdtyp, mxdatm, mxdgvein, mxdnstin, mxdlqp, mxdlao)

  allocate(natom(mxdtyp))
  allocate(nameat(mxdtyp))
  allocate(rat(3,mxdatm,mxdtyp))

  allocate(kgvin(3,mxdgvein))
  allocate(phasein(mxdgvein))
  allocate(conjin(mxdgvein))
  allocate(mstarin(mxdnstin))
  allocate(denin(mxdnstin))
  allocate(denbondin(mxdnstin))
  allocate(veffin(mxdnstin))


  allocate(nqnl(mxdtyp))
  allocate(delqnl(mxdtyp))
  allocate(vkb(-2:mxdlqp,0:3,-1:1,mxdtyp))
  allocate(nkb(0:3,-1:1,mxdtyp))
  allocate(vloc(-1:mxdlqp,mxdtyp))
  allocate(dcor(-1:mxdlqp,mxdtyp))
  allocate(dval(-1:mxdlqp,mxdtyp))
  allocate(norbat(mxdtyp))
  allocate(nqwf(mxdtyp))
  allocate(delqwf(mxdtyp))
  allocate(lorb(mxdlao,mxdtyp))
  allocate(wvfao(-2:mxdlqp,mxdlao,mxdtyp))
  allocate(zv(mxdtyp))

  allocate(iray(mxdtyp))
  allocate(psdtitle(20,mxdtyp))
  allocate(irel(mxdtyp))
  allocate(icore(mxdtyp))
  allocate(icorr(mxdtyp))


  ipr = 1
  icorr = '  '

  call pw_rho_v_in(filename, iotape, ipr,                                &
         pwline, title, subtitle, meta_cpw2000,                          &
         author, flgscf, flgdalin,                                       &
         emaxin, teleck, nx,ny,nz, sx,sy,sz, nbandin, alatt, efermi,     &
         adot, ntype, natom, nameat, rat,                                &
         ntrans, mtrx, tnp,                                              &
         ngin, kmaxin, kgvin, phasein, conjin, nsin, mstarin,            &
         veffin, denin, denbondin,                                       &
         irel, icore, icorr, iray, psdtitle,                             &
         ealraw, zv, ztot,                                               &
         nqnl, delqnl, vkb, nkb, vloc, dcor, dval,                       &
         norbat, nqwf, delqwf, wvfao, lorb, latorb,                      &
         mxdtyp, mxdatm, mxdgvein, mxdnstin, mxdlqp, mxdlao)

  emax = emaxin

  write(6,*)
  write(6,*) '   Do you want an input file for abinit? (y/n)'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'     input file for abinit'

  if(yesno == 'Y' .or. yesno == 'y') then

    call pw2o_abinit_in(adot, ntype, natom, nameat, rat, alatt,          &
        emax, nx,ny,nz,                                                  &
        mxdtyp, mxdatm)

    write(6,*)
    write(6,*) '   Finished writing input file for abinit'
    write(6,*) '   Check if it is correct before using'
    write(6,*)

  endif

  write(6,*)
  write(6,*) '   Do you want an input file for quantum espresso? (y/n)'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'     input file for QE'

  if(yesno == 'Y' .or. yesno == 'y') then

    call pw2o_qe_pwscf_in(meta_cpw2000,                                  &
        adot, ntype, natom, nameat, rat, alatt,                          &
        emax, nbandin, nx,ny,nz, sx,sy,sz,                               &
        mxdtyp, mxdatm)

    write(6,*)
    write(6,*) '   Finished writing scf input file for Quantum Espresso'
    write(6,*) '   Check if it is correct before using'
    write(6,*)

    lso = .FALSE.
    do nt = 1,ntype
      if(irel(nt) == 'rel') lso = .TRUE.
    enddo

    call pw2o_qe_bands_in(meta_cpw2000, lso, fileband,                   &
        adot, ntype, natom, nameat, rat, alatt,                          &
        emax, nbandin,                                                   &
        mxdtyp, mxdatm)

    write(6,*)
    write(6,*) '   Finished writing input file for Quantum Espresso'
    write(6,*) '   band structure.  Check if it is correct before using'
    write(6,*)
    do nt = 1,ntype

      call pw2o_qe_pwscf_upf_in(nameat(nt),                              &
          irel(nt), iray(nt), psdtitle(:,nt),                            &
          nkb(0:3,-1:1,nt))

      write(6,*)
      write(6,*) '   Finished writing pseudopotential generation file'
      write(6,*) '   that sometimes can be used with Quantum Espresso'
      write(6,*) '   for atom ', nameat(nt)

    enddo

  endif

  write(6,*)
  write(6,*)
  write(6,*)
  write(6,*) '   Do you want an input file for pseudopotential code? (y/n)'
  write(6,*)
  write(6,*) '   You can use it (with eventual modifications)'
  write(6,*) '   to generate an equivalent pseudopotential file'
  write(6,*) '   for Quantum Espresso'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'     input file for pseudopotential generation'

  if(yesno == 'Y' .or. yesno == 'y') then

    do nt = 1,ntype

      call pw2o_atom_dat(nameat(nt),                                     &
        irel(nt), icore(nt), icorr(nt), iray(nt), psdtitle(:,nt), zv(nt) )

      write(6,*)
      write(6,*) '   Finished writing pseudopotential generation file'
      write(6,*) '   that can be used with atom pseudopotential code'
      write(6,*) '   to generate a Quantum Espresso'
      write(6,*) '   pseudopotential for ', nameat(nt)

    enddo

  endif

  write(6,*)


  deallocate(natom)
  deallocate(nameat)
  deallocate(rat)

  deallocate(kgvin)
  deallocate(phasein)
  deallocate(conjin)
  deallocate(mstarin)
  deallocate(denin)
  deallocate(denbondin)
  deallocate(veffin)


  deallocate(nqnl)
  deallocate(delqnl)
  deallocate(vkb)
  deallocate(nkb)
  deallocate(vloc)
  deallocate(dcor)
  deallocate(dval)
  deallocate(norbat)
  deallocate(nqwf)
  deallocate(delqwf)
  deallocate(lorb)
  deallocate(wvfao)
  deallocate(zv)

  deallocate(iray)
  deallocate(psdtitle)
  deallocate(irel)
  deallocate(icore)
  deallocate(icorr)

  return

end subroutine pw2o_pwrhov2other

