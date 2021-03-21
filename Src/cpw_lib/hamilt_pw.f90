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

!>     calculates the representation of the hamiltonian
!>     on a plane-wave basis for one k-point
!>     consolidates 3 subroutines tha are often called together.

       subroutine hamilt_pw(emax, rkpt, lkpg, veffr1, nanl,             &
     & mtxd, isort, qmod, ekpg, hdiag,                                  &
     & ng,kgv,                                                          &
     & nqnl,delqnl,vkb,nkb,                                             &
     & ntype,natom,rat,adot,                                            &
     & anlga,xnlkb,                                                     &
     & mxdtyp,mxdatm,mxdlqp,mxddim,mxdanl,mxdgve)

!      Adapted 2 June 2020 from h_kb_dia.f90. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.98


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type
       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdlqp                     !<  array dimension for local potential
       integer, intent(in)                ::  mxdanl                     !<  array dimension of number of projectors
       integer, intent(in)                ::  mxdgve                     !<  array dimension of G-space vectors
       
       real(REAL64), intent(in)           ::  emax                       !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
       real(REAL64), intent(in)           ::  rkpt(3)                    !<  k-point reciprocal lattice coordinates 
       logical, intent(in)                ::  lkpg                       !<  If true use the previous G-vectors (same mtxd and isort)
       real(REAL64), intent(in)           ::  veffr1                     !<  Average value (veff(1)) of the effective potential.
       integer, intent(in)                ::  nanl                       !<  number of projectors without spin 

       integer, intent(in)                ::  ng                         !<  size of g-space
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  G-vectors in reciprocal lattice coordinates 

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  lattice coordinates of atom j of type i
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in real space

       integer, intent(in)                ::  nqnl(mxdtyp)               !<  number of points for pseudo interpolation for atom k
       real(REAL64), intent(in)           ::  delqnl(mxdtyp)             !<  step used in the pseudo interpolation for atom k
       real(REAL64), intent(in)    ::    vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. NOT normalized to vcell, hartree
       integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)       !<  KB pseudo.  normalization for atom k, ang. mom. l

!      input and output

       integer, intent(inout)             ::  mtxd                       !<  wavefunction dimension
       integer, intent(inout)             ::  isort(mxddim)              !<  G-vector corresponding to coefficient i of wavefunction 

!      output

       real(REAL64), intent(out)          ::  qmod(mxddim)               !<  length of k+g-vector of row/column i
       real(REAL64), intent(out)          ::  ekpg(mxddim)               !<  kinetic energy (hartree) of k+g-vector of row/column i

       real(REAL64), intent(out)          ::  hdiag(mxddim)              !<  hamiltonian diagonal

       complex(REAL64), intent(out)       ::  anlga(mxddim,mxdanl)       !<  KB projectors without spin-orbit
       real(REAL64), intent(out)          ::  xnlkb(mxdanl)              !<  KB normalization without spin-orbit


       call hamilt_struct(emax,rkpt,mtxd,isort,qmod,ekpg,lkpg,           &
     & ng,kgv,adot,                                                      &
     & mxdgve,mxddim)

       call hamilt_diag_kb(mtxd,hdiag,qmod,ekpg,veffr1,                  &
     & nqnl,delqnl,vkb,nkb,                                              &
     & ntype,natom,adot,                                                 &
     & mxdtyp,mxdlqp,mxddim)

       call proj_nl_kb_c16(rkpt,mtxd,isort,nanl,                         &
     & ng,kgv,                                                           &
     & nqnl,delqnl,vkb,nkb,                                              &
     & ntype,natom,rat,adot,                                             &
     & anlga,xnlkb,                                                      &
     & mxdtyp,mxdatm,mxdlqp,mxddim,mxdanl,mxdgve)


       return

       end subroutine hamilt_pw
