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

!>  Calculates the matrices for the k.p method and for
!>  The calculation of the oscillator strength
!>
!>  \author       José Luís Martins
!>  \version      5.09
!>  \date         April 14, 2014. 30 November 2023.
!>  \copyright    GNU Public License v2

subroutine kdotp_matrix(mtxd, neig, psi, ei, rkpt, isort, nder,          &
      h0, dh0drk, d2h0drk2,                                              &
      ng, kgv,                                                           &
      ntype, natom, rat, adot,                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)


! Written April 14, 2014, from previous code. JLM
! Modified, documentation, nder, 20 February 2020. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  wavevector
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalues (Hartree) for rk0

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction

  integer, intent(in)                ::  nder                            !<  order of derivatives to be calculated. kdotp nder =2, oscillator nder = 1.

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

! output

  complex(REAL64), intent(out)       ::  h0(mxdbnd,mxdbnd)               !<  <Psi|H|Psi> without spin-orbit
  complex(REAL64), intent(out)       ::  dh0drk(mxdbnd,mxdbnd,3)         !<  d <Psi|H|Psi> d k
  complex(REAL64), intent(out)       ::  d2h0drk2(mxdbnd,mxdbnd,3,3)     !<  d^2 <Psi|H|Psi> d k^2

! local allocatable arrays

  complex(REAL64), allocatable       ::  pmat(:,:,:)                     !  <Psi|p_j|Psi>

  complex(REAL64), allocatable       ::  anlga(:,:)                      !  KB projectors without spin-orbit
  real(REAL64), allocatable          ::  xnlkb(:)                        !  KB normalization without spin-orbit
  complex(REAL64), allocatable       ::  danlgadrk(:,:,:)                !  d anlga / d rkpt
  complex(REAL64), allocatable       ::  d2anlgadrk2(:,:,:,:)            !  d 2 anlga / d 2 rkpt

  complex(REAL64), allocatable       ::  vnl0(:,:)                       !  <Psi|V_NL|Psi>
  complex(REAL64), allocatable       ::  dvnl0drk(:,:,:)                 !  d <Psi|V_NL|Psi> d k
  complex(REAL64), allocatable       ::  d2vnl0drk2(:,:,:,:)             !  d^2 <Psi|V_NL|Psi> d k^2

! local variables

  integer           ::  mxdanl         !  array dimension of number of projectors
  integer           ::  nanl, nanlso, nanlspin   !  number of KB projectors
  real(REAL64)      ::  vcell, bdot(3,3)


! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64
  complex(REAL64), parameter  :: C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer    ::  i, j, m, n


  if(nder < 0 .or. nder > 2) then
    write(6,*)
    write(6,*) '  Stopped in kdotp_matrix:  nder = ',nder
    write(6,*)

    stop

  elseif(nder == 0) then
    write(6,*)
    write(6,*) '  WARNING in kdotp_matrix. nder = 0',                    &
                   ' does not make sense'
    write(6,*)
  endif

  call adot_to_bdot(adot, vcell, bdot)


! momentum matrix elements

  allocate(pmat(3,mxdbnd,mxdbnd))

  call psi_p_psi(mtxd, neig, psi, pmat, rkpt, isort, ng, kgv, .FALSE.,   &
           mxddim, mxdbnd, mxdgve)

! non-local pseudopotential matrix elements

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso,  nanlspin,       &
      mxdtyp)

  mxdanl = nanl

  allocate(anlga(mxddim,mxdanl))
  allocate(xnlkb(mxdanl))

  allocate(vnl0(mxdbnd,mxdbnd))

  allocate(danlgadrk(mxddim,mxdanl,3))
  allocate(dvnl0drk(mxdbnd,mxdbnd,3))

  if(nder == 2) then
    allocate(d2anlgadrk2(mxddim,mxdanl,3,3))
    allocate(d2vnl0drk2(mxdbnd,mxdbnd,3,3))
  else
    allocate(d2anlgadrk2(1,1,3,3))
    allocate(d2vnl0drk2(1,2,3,3))
  endif


  call proj_nl_kb_der_c16(rkpt, mtxd, isort, nanl, nder,                 &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlga, xnlkb, danlgadrk, d2anlgadrk2,                              &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)


  call psi_vnl_psi_der(mtxd, neig, nanl, psi, nder,                      &
      vnl0, dvnl0drk, d2vnl0drk2,                                        &
      anlga, xnlkb, danlgadrk, d2anlgadrk2,                              &
      mxddim, mxdbnd, mxdanl)

  do i=1,neig
  do j=1,neig
    h0(j,i) = C_ZERO
  enddo
  enddo

! eigenvalues

  do j=1,neig
    h0(j,j) =  h0(j,j) + cmplx(ei(j),ZERO,REAL64)
  enddo

! kinetic and non-local energy

  do i=1,neig
  do j=1,neig
    do m=1,3
    dh0drk(j,i,m) = dvnl0drk(j,i,m)
      do n=1,3
        dh0drk(j,i,m)  =  dh0drk(j,i,m) + bdot(m,n)*pmat(n,j,i)
      enddo
    enddo
  enddo
  enddo

  if(nder == 2) then

    do i=1,neig
    do j=1,neig
      do n=1,3
      do m=1,3
        d2h0drk2(j,i,m,n) = d2vnl0drk2(j,i,m,n)
      enddo
      enddo
    enddo
    enddo

    do j=1,neig
      do n=1,3
      do m=1,3
        d2h0drk2(j,j,m,n) = d2h0drk2(j,j,m,n) + bdot(m,n)
      enddo
      enddo
    enddo

  endif


  deallocate(anlga)
  deallocate(xnlkb)
  deallocate(vnl0)

  deallocate(danlgadrk)
  deallocate(pmat)

  deallocate(dvnl0drk)

  deallocate(d2vnl0drk2)
  deallocate(d2anlgadrk2)

  return

end subroutine kdotp_matrix

