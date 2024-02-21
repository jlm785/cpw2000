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

!>  Calculates the force and stress
!>
!>  \author       Jose Luis Martins
!>  \version      5.10
!>  \date         20 October 93, 21 February 2024.
!>  \copyright    GNU Public License v2

  subroutine force_stress_kb(force, stress, energy,                      &
      forcew, strew, strxc, ealpha, flgpsd,                              &
      ntype, natom, nameat, rat, adot,                                   &
      ntrans, mtrx, tnp,                                                 &
      ng, kgv, phase, conj, ns, mstar, ek,                               &
      nqnl, delqnl, vkb, nkb,                                            &
      vion, vhar, vxc, den,                                              &
      mtxd_allk, isort_allk, psi_allk, occ_allk,                         &
      vql, dnc, dvql, ddc,                                               &
      nrk, nband, rk, wgk,                                               &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve, mxdnst, mxdnrk)

! Version 4.0. 20 october 93. jlm
! Version 4.40 24/07/2002. jlm
! Modified to symmetrize strxc (gga) 28/9/2004, version 4.42 jlm
! Modified f90, January 2017.  JLM
! Modified, documentation, January 2020. JLM
! Modified, indentation, remove print, 21 February 2024. JLM


  implicit none

  integer, parameter :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdnrk                          !<  array dimension of k-points

  real(REAL64), intent(in)           ::  energy                          !<  Total electronic energy (Hartree)
  real(REAL64), intent(in)           ::  forcew(3,mxdatm,mxdtyp)         !<  d enerew / d rat,  Ewald contribution to force (contravariant components)
  real(REAL64), intent(in)           ::  strew(3,3)                      !<  d enerew / d adot,  Ewald contribution to stress tensor (contravariant components)
  real(REAL64), intent(in)           ::  ealpha                          !<  G=0 contribution to the total energy (Hartree)
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  real part of the phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of G-vectors in the j-th star
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  complex(REAL64), intent(in)        ::  vion(mxdnst)                    !<  ionic potential for the prototype G-vector in star j
  complex(REAL64), intent(in)        ::  den(mxdnst)                     !<  total charge density for the prototype G-vector
  complex(REAL64), intent(in)        ::  vhar(mxdnst)                    !<  Hartree potential for the prototype G-vector
  complex(REAL64), intent(in)        ::  vxc(mxdnst)                     !<  exchange+correlation potential for the prototype G-vector

  integer, intent(in)                ::  mtxd_allk(mxdnrk)               !<  dimension of the hamiltonian for k-point n
  integer, intent(in)                ::  isort_allk(mxddim,mxdnrk)       !<  G-vector associated with k+G vector i of hamiltonian for k-point n
  complex(REAL64), intent(in)        ::  psi_allk(mxddim,mxdbnd,mxdnrk)  !<  eigenvectors for all k-points
  real(REAL64), intent(in)           ::  occ_allk(mxdnrk*mxdbnd)         !<  fractional ocupation of level j, for all the k-points

  real(REAL64), intent(in)           ::  vql(mxdtyp,mxdnst)              !<  local pseudopotential for atom type i and prototype g-vector in star j       real*8 floc(3,mxdatm,mxdtyp)
  real(REAL64), intent(in)           ::  dnc(mxdtyp,mxdnst)              !<  core charge for atom type i and prototype g-vector in star j
  complex(REAL64), intent(in)        ::  dvql(mxdnst)                    !<  derivative of the local pseudopotential for the prototype g-vector in star j
  complex(REAL64), intent(in)        ::  ddc(mxdnst)                     !<  derivative of the core charge for the prototype g-vector in star j

  integer, intent(in)                ::  nrk                             !<  number of k-points for integration in the irreducible wedge of the brillouin zone
  integer, intent(in)                ::  nband(mxdnrk)                   !<  number of bands for each k-points
  real(REAL64), intent(in)           ::  rk(3,mxdnrk)                    !<  component in lattice coordinates of the k-point in the mesh
  real(REAL64), intent(in)           ::  wgk(mxdnrk)                     !<  weight in the integration of k-point

! input and output

  real(REAL64), intent(inout)        ::  strxc(3,3)                      !<  contribution of xc to the stress tensor (contravariant,Hartree).  It is only symmetrized here.

! output

  real(REAL64), intent(out)          ::  force(3,mxdatm,mxdtyp)          !<  d energy / d rat,  force on the n-th atom of type i (contravarian components, hartree/bohr)
  real(REAL64), intent(out)          ::  stress(3,3)                     !<  d energy / d adot,  stress tensor (contravariant)

! local variables

  real(REAL64)         ::  vcell, bdot(3,3), adotm1(3,3)
  real(REAL64)         ::  rkpt(3)
  real(REAL64)         ::  sunsym(3,3), ssym(3,3)
  real(REAL64)         ::  strhl(3,3), strnlkb(3,3), strkin(3,3)
  integer              ::  iel, neig, ipr, mtxd

! allocatable variables

  real(REAL64), allocatable          ::  occp(:)                         !  ocupation*weight*spin deg. of eigenvector j

  real(REAL64), allocatable          ::  floc(:,:,:)
  real(REAL64), allocatable          ::  funsym(:,:,:)
  real(REAL64), allocatable          ::  fsym(:,:,:)
  real(REAL64), allocatable          ::  fnlkb(:,:,:)

! constants

  real(REAL64), parameter   :: ZERO = 0.0_REAL64
  real(REAL64), parameter   :: PI = 3.14159265358979323846_REAL64

! counters

  integer       ::   i, j, k, l, irk



  allocate(occp(mxdbnd))

  allocate(floc(3,mxdatm,mxdtyp))
  allocate(funsym(3,mxdatm,mxdtyp))
  allocate(fsym(3,mxdatm,mxdtyp))
  allocate(fnlkb(3,mxdatm,mxdtyp))


  if(flgpsd /= 'PSEUKB') then
    write(6,'("   STOPPED in force_stress_kb:   wrong pseudopotencial ",a6)') flgpsd

    stop

  endif

  call adot_to_bdot(adot,vcell,bdot)

  do i = 1,3
  do j = 1,3
    adotm1(i,j) = bdot(i,j) / (4*PI*PI)
  enddo
  enddo

! local contribution to force (covariant)

  call for_str_local_force(floc,                                         &
      ng, kgv, phase, conj, ns, mstar,                                   &
      vxc, den,                                                          &
      ntype, natom, rat,                                                 &
      vql, dnc,                                                          &
      mxdgve, mxdnst, mxdtyp, mxdatm)

! local contribution to stress (covariant)

  call for_str_local_stress(ealpha, strhl,                               &
      ng, kgv, ns, mstar, ek,                                            &
      vion, vhar, vxc, den,                                              &
      adot,                                                              &
      dvql, ddc,                                                         &
      mxdgve, mxdnst)

! strxc in GGA is not strictly symmetric.
! It is not elegant but each contribution is kept separate

  do i = 1,3
  do j = 1,3
    sunsym(i,j) = ZERO
    do k=1,3
    do l=1,3
      sunsym(i,j) = sunsym(i,j) + adot(i,k)*strxc(k,l)*adot(l,j)
    enddo
    enddo
  enddo
  enddo

  do i = 1,ntype
    do j = 1,natom(i)
    do k = 1,3
      funsym(k,j,i) = ZERO
    enddo
    enddo
  enddo

! symetrizes

  call for_str_sym(funsym, sunsym, fsym, ssym,                           &
      ntype, natom, rat,                                                 &
      ntrans, mtrx, tnp,                                                 &
      mxdtyp ,mxdatm)

  do i = 1,3
  do j = 1,3
    strxc(i,j) = ZERO
    do k = 1,3
    do l = 1,3
      strxc(i,j) = strxc(i,j) + adotm1(i,k)*ssym(k,l)*adotm1(l,j)
    enddo
    enddo
  enddo
  enddo

  do i = 1,3
  do j = 1,3
    sunsym(j,i) = ZERO
  enddo
  enddo

  do i = 1,ntype
    do j = 1,natom(i)
    do k = 1,3
      funsym(k,j,i) = ZERO
    enddo
    enddo
  enddo

  iel = 0
  do irk = 1,nrk

!   loop over k-points

    rkpt(1) = rk(1,irk)
    rkpt(2) = rk(2,irk)
    rkpt(3) = rk(3,irk)

    neig = nband(irk)
    mtxd = mtxd_allk(irk)

    do j = 1,neig
      iel = iel + 1
      occp(j) = 2*wgk(irk)*occ_allk(iel)
    enddo
!
    call for_str_kinetic_stress(strkin,                                  &
        mtxd, rkpt, neig, occp,                                          &
        isort_allk(:,irk) ,psi_allk(:,:,irk),                            &
        kgv,                                                             &
        mxdgve, mxddim, mxdbnd)

    do i = 1,3
    do j = 1,3
      sunsym(j,i) = sunsym(j,i) + strkin(j,i)
    enddo
    enddo

    call for_str_nl_kb(fnlkb, strnlkb,                                   &
        mtxd, rkpt, neig, occp,                                          &
        isort_allk(:,irk), psi_allk(:,:,irk),                            &
        kgv,                                                             &
        nqnl, delqnl, vkb, nkb,                                          &
        ntype, natom, rat, adot,                                         &
        mxdtyp, mxdatm, mxdgve, mxdlqp, mxddim, mxdbnd)

    do i = 1,3
    do j = 1,3
      sunsym(j,i) = sunsym(j,i) + strnlkb(j,i)
    enddo
    enddo

    do i = 1,ntype
      do j = 1,natom(i)
      do k = 1,3
        funsym(k,j,i) = funsym(k,j,i) + fnlkb(k,j,i)
      enddo
      enddo
    enddo

  enddo

  call for_str_sym(funsym, sunsym, fsym, ssym,                           &
      ntype, natom, rat,                                                 &
      ntrans, mtrx, tnp,                                                 &
      mxdtyp, mxdatm)

! sums up contributions

  do i = 1,3
  do j = 1,3
    stress(i,j) = strew(i,j) + strxc(i,j)
    do k=1,3
    do l=1,3
      stress(i,j) = stress(i,j) + adotm1(i,k) *                          &
                       (strhl(k,l) + ssym(k,l)) * adotm1(l,j)
    enddo
    enddo
  enddo
  enddo

  do i = 1,ntype
    do j = 1,natom(i)
    do k = 1,3
      force(k,j,i) = forcew(k,j,i)
      do l = 1,3
        force(k,j,i) = force(k,j,i) + adotm1(k,l)                        &
                        * (fsym(l,j,i) + floc(l,j,i))
      enddo
    enddo
    enddo
  enddo

  deallocate(occp)

  deallocate(floc)
  deallocate(funsym)
  deallocate(fsym)
  deallocate(fnlkb)

  return

end subroutine force_stress_kb

