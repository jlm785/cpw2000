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

!>  This subroutines reads the file "filename",  (default PW_RHO_V.DAT),
!>  data from a self-consistent calculation
!>
!>  \author       Jose Luis Martins
!>  \version      5.02
!>  \date         May 16, 2014, 31 December 2020. JLM
!>  \copyright    GNU Public License v2

subroutine pw_rho_v_in(filename, io, ipr,                                &
         pwline, title, subtitle, meta_cpw2000,                          &
         author, flgscf, flgdal, emax, teleck,                           &
         nx,ny,nz, sx,sy,sz, nband,                                      &
         ntrans, mtrx, tnp,                                              &
         alatt, adot, ntype, natom, nameat, rat,                         &
         ng, kmax, kgv, phase, conj, ns, mstar,                          &
         veff, den, denbond,                                             &
         irel, icore, icorr, iray, ititle,                               &
         ealraw, zv, ztot,                                               &
         nqnl, delqnl, vkbraw, nkb, vloc, dcor, dval,                    &
         norbat, nqwf, delqwf, wvfao, lorb, latorb,                      &
         mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdlao)

! Written January 12, 2014. JLM
! Modified (split) May 27, 2014. JLM
! Modified pwline,title,subtitle, 1 August 2014. JLM
! Modified, mxdlao, December 1, 2015.  JLM
! Modified, meta_cpw2000, author,nx,etc, January 10, 2017. JLM
! Modified October 15 2018 to pass generation information for pwSCF.

! API CHANGED FROM EARLIER VERSIONS!!!!!!!!!!!!!!!!!!

! Modified, documentation, February 4 2020. JLM

  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  character(len=*), intent(in)       ::  filename                        !<  name of file
  integer, intent(in)                ::  io                              !<  number of tape to which the pseudo is added.

  integer, intent(in)                ::  ipr                             !<  printing level

! output

  character(len=60), intent(out)     ::  pwline                          !<  identification of the calculation
  character(len=50), intent(out)     ::  title                           !<  title for plots
  character(len=140), intent(out)    ::  subtitle                        !<  subtitle for plots
  character(len=250), intent(out)    ::  meta_cpw2000                    !<  metadata from cpw2000

  character(len=3), intent(out)      ::  author                          !<  type of xc wanted (CA=PZ , PW92 , PBE)
  character(len=6), intent(out)      ::  flgscf                          !<  type of self consistent field and diagonalization
  character(len=4), intent(out)      ::  flgdal                          !<  whether the dual approximation is used
  real(REAL64), intent(out)          ::  emax                            !<  kinetic energy cutoff of plane wave expansion (Hartree).
  real(REAL64), intent(out)          ::  teleck                          !<  electronic temperature (in Kelvin)

  integer, intent(out)               ::  nband                           !<  target for number of bands
  integer, intent(out)               ::  nx,ny,nz                        !<  divisions of Brillouin zone for integration (Monkhorst-Pack)
  real(REAL64), intent(out)          ::  sx,sy,sz                        !<  shift of points in division of Brillouin zone for integration (Monkhorst-Pack)
  real(REAL64), intent(out)          ::  alatt                           !<  lattice constant

  integer, intent(out)               ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(out)               ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(out)          ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

  real(REAL64), intent(out)          ::  adot(3,3)                       !<  metric in direct space
  integer, intent(out)               ::  ntype                           !<  number of types of atoms
  integer, intent(out)               ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(out)      ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(out)               ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(out)               ::  kmax(3)                         !<  max value of kgv(i,n)
  integer, intent(out)               ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(out)       ::  phase(mxdgve)                   !<   phase factor of G-vector n
  real(REAL64), intent(out)          ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(out)               ::  ns                              !<  number os stars with length less than gmax
  integer, intent(out)               ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

  complex(REAL64), intent(out)       ::  veff(mxdnst)                    !<  effective potential (local+Hartree+Xc) for the prototype g-vector in star j
  complex(REAL64), intent(out)       ::  den(mxdnst)                     !<  valence charge density for the prototype g-vector in star j
  complex(REAL64), intent(out)       ::  denbond(mxdnst)                 !<  bonding charge density for the prototype g-vector in star j

  character(len=3), intent(out)      ::  irel(mxdtyp)                    !<  type of calculation relativistic/spin
  character(len=4), intent(out)      ::  icore(mxdtyp)                   !<  type of partial core correction
  character(len=2), intent(out)      ::  icorr(mxdtyp)                   !<  type of correlation
  character(len=60), intent(out)     ::  iray(mxdtyp)                    !<  information about pseudopotential
  character(len=70), intent(out)     ::  ititle(mxdtyp)                  !<  further information about pseudopotential

  real(REAL64), intent(out)          ::  ealraw                          !<  G=0 contrib. to the total energy. (non norm. to vcell,hartree)
  integer, intent(out)               ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(out)          ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(out)       ::  vkbraw(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * kb nonlocal pseudo. for atom k, ang. mom. l. (non normalized to vcell, hartree)
  integer, intent(out)               ::  nkb(0:3,-1:1,mxdtyp)            !<   kb pseudo.  normalization for atom k, ang. mom. l
  real(REAL64), intent(out)          ::  vloc(-1:mxdlqp,mxdtyp)          !<  local pseudopotential for atom k (hartree)
  real(REAL64), intent(out)          ::  dcor(-1:mxdlqp,mxdtyp)          !<  core charge density for atom k
  real(REAL64), intent(out)          ::  dval(-1:mxdlqp,mxdtyp)          !<  valence charge density for atom k
  integer, intent(out)               ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(out)               ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(out)          ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  integer, intent(out)               ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k
  real(REAL64), intent(out)          ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  (1/q**l) * wavefunction for atom k, ang. mom. l (unnormalized to vcell)
  logical, intent(out)               ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  real(REAL64), intent(out)          ::  zv(mxdtyp)                      !<  valence of atom with type i
  real(REAL64), intent(out)          ::  ztot                            !<  total charge density (electrons/cell)

! counters

  integer    ::  i, j, k




  open(unit=io,file=trim(filename),status='old',form='UNFORMATTED')

! reads the first part of the file up to the geometry

  call pw_rho_v_in_crystal_calc(io,                                      &
         pwline, title, subtitle, meta_cpw2000,                          &
         author, flgscf, flgdal, emax, teleck,                           &
         nx,ny,nz, sx,sy,sz, nband,                                      &
         ng ,ns,                                                         &
         ntrans, mtrx, tnp,                                              &
         alatt, adot, ntype, natom, nameat, rat,                         &
         mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp)

! reads the self-consistent charge and potential

  read(io) (kmax(j),j=1,3)
  read(io) ((kgv(j,k),j=1,3),k=1,ng)
  read(io) (phase(i),conj(i),i=1,ng)
  read(io) (mstar(i),i=1,ns)

  read(io) (den(i),i=1,ns)
  read(io) (denbond(i),i=1,ns)
  read(io) (veff(i),i=1,ns)

! reads the pseudopotentials

  call pw_rho_v_in_pseudo(io, ipr, ealraw, author,                       &
         irel, icore, icorr, iray, ititle,                               &
         nqnl, delqnl, vkbraw, nkb, vloc, dcor, dval,                    &
         norbat, nqwf, delqwf, wvfao, lorb, latorb,                      &
         ntype, natom, nameat, zv, ztot,                                 &
         mxdtyp, mxdlqp, mxdlao)

  close(unit = io)

  return

end subroutine pw_rho_v_in
