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

!>  Writes a file with the self consistent
!>  potential and charge density for later processing
!>  by cpw_post_process
!>
!>  \author       Jose Luis Martins
!>  \version      5.03
!>  \date         22 April 2021, 21November 2021.
!>  \copyright    GNU Public License v2

subroutine pw_rho_v_out(filename, io, author, tblaha, flgscf, flgdal,    &
         meta_pwdat, meta_cpw2000,                                       &
         emax, teleck, nx,ny,nz, sx,sy,sz, nband, alatt, efermi,         &
         adot, ntype, natom, nameat, rat,                                &
         ntrans, mtrx, tnp,                                              &
         ng, kmax, kgv, phase, conj, ns, mstar,                          &
         veff, den, dens,                                                &
         mxdtyp, mxdatm, mxdgve, mxdnst)

! Written December 18-22, 2013. jlm
! Modified January 9, 2014. jlm
! Modified, pwline bug, July 10, 2014.  JLM
! Modified, pseudopotential output, April 16, 2004. JLM
! Increased size pwline. August 1 2014. JLM
! Modified c16, 24 October 2015. JLM
! Modified pwline length, 14 December 2016. JLM
! pwline replaced by meta_pwdat, meta_cpw2000 in input, January 5, 2017. JLM
! Modified, documentation, August 2019. JLM
! Modified, mxdlao, ntrans, 13 September 2021. JLM
! Modified, efermi, 29 November 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  integer, intent(in)                ::  io                              !<  number of tape to which the pseudo is added.

  character(len=250), intent(in)     ::  meta_pwdat                      !<  metadata from PW.DAT
  character(len=250), intent(in)     ::  meta_cpw2000                    !<  metadata from cpw2000

  character(len=*), intent(in)       ::  filename                        !<  input file
  character(len=4), intent(in)       ::  author                          !<  type of xc wanted (CA=PZ , PW92 , PBE)
  real(REAL64), intent(in)           ::  tblaha                          !<  Tran-Blaha constant
  character(len=6), intent(in)       ::  flgscf                          !<  type of self consistent field and diagonalizatioN
  character(len=4), intent(in)       ::  flgdal                          !<  whether the dual approximation is used
  real(REAL64), intent(in)           ::  emax                            !<  kinetic energy cutoff of plane wave expansion (hartree).
  real(REAL64), intent(in)           ::  teleck                          !<  electronic temperature (in kelvin)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

  integer, intent(in)                ::  nband                           !<  target for number of bands
  integer, intent(in)                ::  nx,ny,nz                        !<  divisions of Brillouin zone for integration (Monkhorst-Pack)
  real(REAL64), intent(in)           ::  sx,sy,sz                        !<  shift of points in division of Brillouin zone for integration (Monkhorst-Pack)
  real(REAL64), intent(in)           ::  alatt                           !<  lattice constant
  real(REAL64), intent(in)           ::  efermi                          !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kmax(3)                         !<  max value of kgv(i,n)
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

  complex(REAL64), intent(in)        ::  veff(mxdnst)                    !<  ionic potential (hartree) for the prototype g-vector in star j
  complex(REAL64), intent(in)        ::  den(mxdnst)                     !<  valence charge density for the prototype g-vector in star j
  complex(REAL64), intent(in)        ::  dens(mxdnst)                    !<  spherical atomic valence charge density for the prototype g-vector in star j

! local variables

  character(len=9 )   ::  bdate
  character(len=8)    ::  btime
  integer             ::  mxdl, mxdlao

! counters

  integer    ::  i, j, k


  call size_mxdlqp_lao(ntype, nameat,                                    &
         mxdtyp, mxdl, mxdlao)

  open(unit=io, file=filename, status='UNKNOWN', form='UNFORMATTED')

! quantities relevant for array size (+ntrans)

  write(io) ntype, ng, ns, mxdl, ntrans, mxdlao
  write(io) (natom(i),i=1,ntype)

! date and calculation identification and parameters

  call zedate(bdate)
  call zetime(btime)
  write(io) bdate,btime

  write(io) author, flgscf, flgdal
  write(io) emax, teleck, nx,ny,nz, sx,sy,sz, nband, alatt, efermi

  write(io) meta_pwdat, meta_cpw2000

! crystal structure (other than quantities relevant for array size)

  write(io) ((adot(i,j),j=1,3),i=1,3),                                   &
            ( ((mtrx(i,j,k),j=1,3),i=1,3), k=1,ntrans ),                 &
            ( (tnp(i,k),i=1,3), k=1,ntrans )
  write(io) (nameat(i),i=1,ntype)

  do i=1,ntype
    write(io) ((rat(j,k,i),j=1,3),k=1,natom(i))
  enddo

! G- space structure

  write(io) (kmax(j),j=1,3)
  write(io) ((kgv(j,k),j=1,3),k=1,ng)
  write(io) (phase(i),conj(i),i=1,ng)
  write(io) (mstar(i),i=1,ns)

! charge density, bonding charge density, effective (Hartree+XC) potential

  write(io) (den(i),i=1,ns)
  write(io) (den(i)-dens(i),i=1,ns)
  write(io) (veff(i),i=1,ns)

  call read_write_pseudo(io, ntype, nameat,                              &
      mxdtyp)

  close(unit = io)

  return
end subroutine pw_rho_v_out
