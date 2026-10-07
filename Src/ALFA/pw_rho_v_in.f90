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
!>  data from a self-consistent calculation.
!>  Files written before October 2026 are recognized by the shorter
!>  first record (no mxdset).  They have one atomic basis set, no core
!>  kinetic energy density and no vtau (set to zero, with a warning
!>  for the generalized Kohn-Sham meta-GGA).
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         May 16, 2014, 7 October 2026.
!>  \copyright    GNU Public License v2

subroutine pw_rho_v_in(filename, io, ipr,                                &
         pwline, title, subtitle, meta_cpw2000,                          &
         author, flgscf, flgdal,                                         &
         emax, teleck, nx,ny,nz, sx,sy,sz, nband, alatt, efermi,         &
         adot, ntype, natom, nameat, rat,                                &
         ntrans, mtrx, tnp,                                              &
         ng, kmax, kgv, phase, conj, ns, mstar,                          &
         veff, vtau, den, denbond,                                       &
         irel, icore, icorr, iray, psdtitle,                             &
         ealraw, zv, ztot,                                               &
         nqnl, delqnl, vkbraw, nkb, vloc, dcor, dval, tauc_q,            &
         n_bsets, norbat, nqwf, delqwf, wvfao, lorb, latorb,             &
         mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdlao, mxdset)

! Written January 12, 2014. JLM
! Modified (split) May 27, 2014. JLM
! Modified pwline,title,subtitle, 1 August 2014. JLM
! Modified, mxdlao, December 1, 2015.  JLM
! Modified, meta_cpw2000, author,nx,etc, January 10, 2017. JLM
! Modified October 15 2018 to pass generation information for pwSCF.
! Modified, documentation, February 4 2020. JLM
! Modified, efermi order of ntrans, mtrx, tnp, 29 November 2021. JLM
! Modified, size of author, 13 January 2024.
! Modified, ititle -> psdtitle. 20 February 2025. JLM
! Documentation, one argument per declaration. 28 September 2026. JLM+claude
! vtau, several atomic basis sets, core tau, old files recognized. 7 October 2026. JLM+claude


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type
  integer, intent(in)                ::  mxdset                          !<  array dimension for number of atomic basis sets

  character(len=*), intent(in)       ::  filename                        !<  name of file
  integer, intent(in)                ::  io                              !<  number of tape to which the pseudo is added.

  integer, intent(in)                ::  ipr                             !<  printing level

! output

  character(len=60), intent(out)     ::  pwline                          !<  identification of the calculation
  character(len=50), intent(out)     ::  title                           !<  title for plots
  character(len=140), intent(out)    ::  subtitle                        !<  subtitle for plots
  character(len=250), intent(out)    ::  meta_cpw2000                    !<  metadata from cpw2000

  character(len=4), intent(out)      ::  author                          !<  type of xc wanted (CA=PZ , PW92 , PBE)
  character(len=6), intent(out)      ::  flgscf                          !<  type of self consistent field and diagonalization
  character(len=4), intent(out)      ::  flgdal                          !<  whether the dual approximation is used
  real(REAL64), intent(out)          ::  emax                            !<  kinetic energy cutoff of plane wave expansion (Hartree).
  real(REAL64), intent(out)          ::  teleck                          !<  electronic temperature (in Kelvin)

  integer, intent(out)               ::  nband                           !<  target for number of bands
  integer, intent(out)               ::  nx                              !<  divisions of Brillouin zone for integration (Monkhorst-Pack), direction 1
  integer, intent(out)               ::  ny                              !<  divisions of Brillouin zone for integration (Monkhorst-Pack), direction 2
  integer, intent(out)               ::  nz                              !<  divisions of Brillouin zone for integration (Monkhorst-Pack), direction 3
  real(REAL64), intent(out)          ::  sx                              !<  shift of points in division of Brillouin zone for integration (Monkhorst-Pack), direction 1
  real(REAL64), intent(out)          ::  sy                              !<  shift of points in division of Brillouin zone for integration (Monkhorst-Pack), direction 2
  real(REAL64), intent(out)          ::  sz                              !<  shift of points in division of Brillouin zone for integration (Monkhorst-Pack), direction 3
  real(REAL64), intent(out)          ::  alatt                           !<  lattice constant
  real(REAL64), intent(out)          ::  efermi                          !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

  real(REAL64), intent(out)          ::  adot(3,3)                       !<  metric in direct space
  integer, intent(out)               ::  ntype                           !<  number of types of atoms
  integer, intent(out)               ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(out)      ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(out)               ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(out)               ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(out)          ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

  integer, intent(out)               ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(out)               ::  kmax(3)                         !<  max value of kgv(i,n)
  integer, intent(out)               ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(out)       ::  phase(mxdgve)                   !<   phase factor of G-vector n
  real(REAL64), intent(out)          ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(out)               ::  ns                              !<  number os stars with length less than gmax
  integer, intent(out)               ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

  complex(REAL64), intent(out)       ::  veff(mxdnst)                    !<  effective potential (local+Hartree+Xc) for the prototype g-vector in star j
  complex(REAL64), intent(out)       ::  vtau(mxdnst)                    !<  d (rho eps_xc) / d tau (generalized Kohn-Sham meta-GGA, zero otherwise) for the prototype g-vector in star j
  complex(REAL64), intent(out)       ::  den(mxdnst)                     !<  valence charge density for the prototype g-vector in star j
  complex(REAL64), intent(out)       ::  denbond(mxdnst)                 !<  bonding charge density for the prototype g-vector in star j

  character(len=3), intent(out)      ::  irel(mxdtyp)                    !<  type of calculation relativistic/spin
  character(len=4), intent(out)      ::  icore(mxdtyp)                   !<  type of partial core correction
  character(len=2), intent(out)      ::  icorr(mxdtyp)                   !<  type of correlation
  character(len=60), intent(out)     ::  iray(mxdtyp)                    !<  information about pseudopotential
  character(len=10), intent(out)     ::  psdtitle(20,mxdtyp)             !<  further information about pseudopotential

  real(REAL64), intent(out)          ::  ealraw                          !<  G=0 contrib. to the total energy. (non norm. to vcell,hartree)
  integer, intent(out)               ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(out)          ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(out)       ::  vkbraw(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * kb nonlocal pseudo. for atom k, ang. mom. l. (non normalized to vcell, hartree)
  integer, intent(out)               ::  nkb(0:3,-1:1,mxdtyp)            !<   kb pseudo.  normalization for atom k, ang. mom. l
  real(REAL64), intent(out)          ::  vloc(-1:mxdlqp,mxdtyp)          !<  local pseudopotential for atom k (hartree)
  real(REAL64), intent(out)          ::  dcor(-1:mxdlqp,mxdtyp)          !<  core charge density for atom k
  real(REAL64), intent(out)          ::  dval(-1:mxdlqp,mxdtyp)          !<  valence charge density for atom k
  real(REAL64), intent(out)          ::  tauc_q(-1:mxdlqp,mxdtyp)        !<  partial core kinetic energy density  (hartree/bohr^3)
  integer, intent(out)               ::  n_bsets(mxdtyp)                 !<  number of basis sets for each atom k
  integer, intent(out)               ::  norbat(mxdset,mxdtyp)           !<  number of atomic orbitals for basis set nb and atom k
  integer, intent(out)               ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(out)          ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  integer, intent(out)               ::  lorb(mxdlao,mxdset,mxdtyp)      !<  angular momentum of orbital n of basis nb of atom k
  real(REAL64), intent(out)          ::  wvfao(-2:mxdlqp,mxdlao,mxdset,mxdtyp)  !<  (1/q**l) * wavefunction for atom k, basis nb, ang. mom. l (unnormalized to vcell)
  logical, intent(out)               ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  real(REAL64), intent(out)          ::  zv(mxdtyp)                      !<  valence of atom with type i
  real(REAL64), intent(out)          ::  ztot                            !<  total charge density (electrons/cell)

! local variables

  logical             ::  lnew                                           !  file written after October 2026
  integer             ::  i1, i2, i3, i4, i5, i6, i7
  integer             ::  ioerr
  character(len=4)    ::  xcbase                                         !  meta-GGA functional used in xc_mgga
  character(len=4)    ::  tausrc                                         !  source of tau, 'PSI ' for generalized Kohn-Sham

! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer    ::  i, j, k




  open(unit=io,file=trim(filename),status='old',form='UNFORMATTED')

! newer files have mxdset in the first record

  read(io,iostat=ioerr) i1, i2, i3, i4, i5, i6, i7
  lnew = ioerr == 0
  rewind(io)

! reads the first part of the file up to the geometry

  call pw_rho_v_in_crystal_calc(io,                                      &
         pwline, title, subtitle, meta_cpw2000,                          &
         author, flgscf, flgdal, emax, teleck,                           &
         nx,ny,nz, sx,sy,sz, nband, alatt, efermi,                       &
         ng ,ns,                                                         &
         ntrans, mtrx, tnp,                                              &
         adot, ntype, natom, nameat, rat,                                &
         mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp)

! reads the self-consistent charge and potential

  read(io) (kmax(j),j=1,3)
  read(io) ((kgv(j,k),j=1,3),k=1,ng)
  read(io) (phase(i),conj(i),i=1,ng)
  read(io) (mstar(i),i=1,ns)

  read(io) (den(i),i=1,ns)
  read(io) (denbond(i),i=1,ns)
  read(io) (veff(i),i=1,ns)

! vtau for the generalized Kohn-Sham meta-GGA

  call xc_author_tau(author, xcbase, tausrc)

  if(tausrc == 'PSI ' .and. lnew) then
    read(io) (vtau(i),i=1,ns)
  else
    do i = 1,ns
      vtau(i) = C_ZERO
    enddo
    if(tausrc == 'PSI ') then
      write(6,*)
      write(6,'("   WARNING in pw_rho_v_in:  file ",a," is in the old ", &
         &      "format, vtau of the ",a4," meta-GGA set to zero")')     &
             trim(filename), author
      write(6,*)
    endif
  endif

! reads the pseudopotentials

  call pw_rho_v_in_pseudo(io, ipr, ealraw, author, lnew,                 &
         irel, icore, icorr, iray, psdtitle,                             &
         nqnl, delqnl, vkbraw, nkb, vloc, dcor, dval, tauc_q,            &
         n_bsets, norbat, nqwf, delqwf, wvfao, lorb, latorb,             &
         ntype, natom, nameat, zv, ztot,                                 &
         mxdtyp, mxdlqp, mxdlao, mxdset)

  close(unit = io)

  return

end subroutine pw_rho_v_in
