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

!>  Reads the header of hsr2.dat, after ao_int_hsr2_size.
!>  The atomic basis set is the one chosen when the file was written,
!>  and it is returned as the only basis set (mxdset = 1).
!>
!>  \author       claude
!>  \version      5.13
!>  \date         24 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_int_hsr2_read(io, meta_cpw2000, title, subtitle,           &
    adot, ztot, emax, ntype, natom, nameat, rat,                         &
    ntrans, mtrx, tnp,                                                   &
    norbat, lorb, nqwf, delqwf, wvfao,                                   &
    lso, loptical, imtb, mp_grid, nband,                                 &
    mxdtyp, mxdatm, mxdlqp, mxdlao)

! Written 24 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdlqp                          !<  array dimension of the radial functions
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  integer, intent(in)                ::  io                              !<  tape number of hsr2.dat (open, unformatted, after ao_int_hsr2_size)
  integer, intent(in)                ::  ntype                           !<  number of types of atoms

! output

  character(len=250), intent(out)    ::  meta_cpw2000                    !<  metadata from cpw2000 (superlattice information for unfolding)
  character(len=50), intent(out)     ::  title                           !<  title for plots
  character(len=140), intent(out)    ::  subtitle                        !<  subtitle for plots

  real(REAL64), intent(out)          ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(out)          ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(out)          ::  emax                            !<  kinetic energy cutoff of plane wave expansion (hartree)
  integer, intent(out)               ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(out)      ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i

  integer, intent(out)               ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(out)               ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(out)          ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

  integer, intent(out)               ::  norbat(1,mxdtyp)                !<  number of atomic orbitals of atom type nt
  integer, intent(out)               ::  lorb(mxdlao,1,mxdtyp)           !<  angular momentum of orbital n of atom type nt
  integer, intent(out)               ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom type nt
  real(REAL64), intent(out)          ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom type nt
  real(REAL64), intent(out)          ::  wvfao(-2:mxdlqp,mxdlao,1,mxdtyp)   !<  (1/q**l) * wavefunction for atom type nt, orbital n

  integer, intent(out)               ::  lso                             !<  1 with spin-orbit, 0 without
  integer, intent(out)               ::  loptical                        !<  1 with optical matrix elements, 0 without
  integer, intent(out)               ::  imtb                            !<  MTB method, 1 projection, 2 matching
  integer, intent(out)               ::  mp_grid(3)                      !<  grid of k-points of the MTB calculation
  integer, intent(out)               ::  nband                           !<  number of basis functions (spinors with spin-orbit)

! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64

! counters

  integer   ::  nt, j, n


  read(io) meta_cpw2000, title, subtitle

  read(io) adot, ztot, emax
  read(io) (natom(nt),nt=1,ntype), (nameat(nt),nt=1,ntype)

  rat(:,:,:) = ZERO
  read(io) (((rat(j,n,nt),j=1,3),n=1,natom(nt)),nt=1,ntype)

  read(io) ntrans
  mtrx(:,:,:) = 0
  tnp(:,:) = ZERO
  read(io) mtrx(:,:,1:ntrans), tnp(:,1:ntrans)

  read(io) (norbat(1,nt),nt=1,ntype), (nqwf(nt),nt=1,ntype),             &
           (delqwf(nt),nt=1,ntype)
  lorb(:,:,:) = 0
  read(io) ((lorb(n,1,nt),n=1,norbat(1,nt)),nt=1,ntype)
  wvfao(:,:,:,:) = ZERO
  do nt = 1,ntype
    do n = 1,norbat(1,nt)
      read(io) (wvfao(j,n,1,nt),j=-2,nqwf(nt)-1)
    enddo
  enddo

  read(io) lso, loptical, imtb, mp_grid, nband

  return

end subroutine ao_int_hsr2_read
