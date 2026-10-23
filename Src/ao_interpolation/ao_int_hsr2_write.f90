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

!>  Writes the header of the file hsr2.dat with the modified tight-binding
!>  (MTB) data.  The header has all the information needed to use the
!>  file without PW_RHO_V.DAT:  titles, crystal structure, symmetry,
!>  the chosen atomic basis set with its radial functions (only up to
!>  the wave-vector needed by the plane-wave cutoff) and the setup of
!>  the MTB calculation.  The hamiltonian and overlap follow, written by
!>  NonOrthoInterpWriteToFile.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         24 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_int_hsr2_write(io, meta_cpw2000, title, subtitle,          &
    adot, ztot, emax, ntype, natom, nameat, rat,                         &
    ntrans, mtrx, tnp,                                                   &
    n_bsets, jbset, norbat, lorb, nqwf, delqwf, wvfao,                   &
    lso, loptical, imtb, mp_grid, nband,                                 &
    mxdtyp, mxdatm, mxdlqp, mxdlao, mxdset)

! Written 24 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type
  integer, intent(in)                ::  mxdset                          !<  array dimension of number of atomic basis sets

  integer, intent(in)                ::  io                              !<  tape number of hsr2.dat (open, unformatted)

  character(len=250), intent(in)     ::  meta_cpw2000                    !<  metadata from cpw2000 (superlattice information for unfolding)
  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  emax                            !<  kinetic energy cutoff of plane wave expansion (hartree)
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

  integer, intent(in)                ::  n_bsets(mxdtyp)                 !<  number of atomic basis sets for atom type nt
  integer, intent(in)                ::  jbset(mxdtyp)                   !<  basis set chosen for atom type nt
  integer, intent(in)                ::  norbat(mxdset,mxdtyp)           !<  number of atomic orbitals in basis set ib of atom type nt
  integer, intent(in)                ::  lorb(mxdlao,mxdset,mxdtyp)      !<  angular momentum of orbital n in basis set ib of atom type nt
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom type nt
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom type nt
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdset,mxdtyp)  !<  (1/q**l) * wavefunction for atom type nt, basis set ib, orbital n

  integer, intent(in)                ::  lso                             !<  1 with spin-orbit, 0 without
  integer, intent(in)                ::  loptical                        !<  1 with optical matrix elements, 0 without
  integer, intent(in)                ::  imtb                            !<  MTB method, 1 projection, 2 matching
  integer, intent(in)                ::  mp_grid(3)                      !<  grid of k-points of the MTB calculation
  integer, intent(in)                ::  nband                           !<  number of basis functions (spinors with spin-orbit)

! local variables

  integer              ::  ib(mxdtyp)                                    !  basis set used, as in atomic_orbital_c16
  integer              ::  nq(mxdtyp)                                    !  number of points of the radial functions written
  integer              ::  natmax, norbmax, nqmax

! parameters

  character(len=24), parameter  ::  HSR2VER = 'cpw2000 hsr2 version 1  '
  integer, parameter            ::  NMARGIN = 10                         !  points beyond the plane-wave cutoff

! counters

  integer   ::  nt, j, n


  natmax = 1
  norbmax = 1
  nqmax = 1
  do nt = 1,ntype
    ib(nt) = max(1,min(jbset(nt),n_bsets(nt)))
    nq(nt) = min(nqwf(nt), int(sqrt(2*emax)/delqwf(nt)) + NMARGIN)
    natmax = max(natmax,natom(nt))
    norbmax = max(norbmax,norbat(ib(nt),nt))
    nqmax = max(nqmax,nq(nt))
  enddo

  write(io) HSR2VER
  write(io) ntype, natmax, norbmax, nqmax

  write(io) meta_cpw2000, title, subtitle

  write(io) adot, ztot, emax
  write(io) (natom(nt),nt=1,ntype), (nameat(nt),nt=1,ntype)
  write(io) (((rat(j,n,nt),j=1,3),n=1,natom(nt)),nt=1,ntype)

  write(io) ntrans
  write(io) mtrx(:,:,1:ntrans), tnp(:,1:ntrans)

  write(io) (norbat(ib(nt),nt),nt=1,ntype), (nq(nt),nt=1,ntype),         &
            (delqwf(nt),nt=1,ntype)
  write(io) ((lorb(n,ib(nt),nt),n=1,norbat(ib(nt),nt)),nt=1,ntype)
  do nt = 1,ntype
    do n = 1,norbat(ib(nt),nt)
      write(io) (wvfao(j,n,ib(nt),nt),j=-2,nq(nt)-1)
    enddo
  enddo

  write(io) lso, loptical, imtb, mp_grid, nband

  return

end subroutine ao_int_hsr2_write
