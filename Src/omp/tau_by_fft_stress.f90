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

!>  Computes, symmetrizes, and adds the "kinetic energy density"
!>  from the eigenvectors at a given k-point of the irreducible wedge.
!>  The sum is done explicitly over all the points of the full
!>  integration mesh that are equivalent to that k-point (kmap from int_pnt),
!>  with the wave-functions obtained by symmetry (psi_rot_inv_shift),
!>  so that the derivatives with respect to the metric can be calculated.
!>  If requested, d tau / d bdot is calculated on the real space mesh of the
!>  potential (kmscr), without folding to G-space.
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.13
!>  \date         14 June 2026.
!>  \copyright    GNU Public License v2

subroutine tau_by_fft_stress(tauk, ldtau, kmscr, dtau_dbdot,             &
    nx, ny, nz, sx, sy, sz,                                              &
    mtxd, neig, occ_x_wgk, isort, psi,                                   &
    rkpt, irk, adot,                                                     &
    ntrans, mtrx, tnp,                                                   &
    nrk, rk, wgk, kmap,                                                  &
    ng, kgv, phase, conj, ns, inds, kmax, mstar,                         &
    mxddim, mxdbnd, mxdnrk, mxdgve, mxdnst, mxdscr)

! Adapted from early version of tau_by_fft. 14 June 2026.
! Debugged: reference point, ipoint for rotated basis, rotated basis
! built from the reference basis instead of hamilt_struct,
! weights from occ_x_wgk. 29 September 2026. JLM+claude
! d tau / d bdot on the real space mesh of the potential. 29 September 2026. JLM+claude


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdnrk                          !<  array dimension of k-points

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdscr                          !<  array dimension of dtau_dbdot (vscr)

  integer, intent(in)                ::  nx                              !<  size of the integration mesh in k-space (nx*ny*nz), direction 1
  integer, intent(in)                ::  ny                              !<  size of the integration mesh in k-space (nx*ny*nz), direction 2
  integer, intent(in)                ::  nz                              !<  size of the integration mesh in k-space (nx*ny*nz), direction 3
  real(REAL64), intent(in)           ::  sx                              !<  offset of the integration mesh (0.5 for Monkhorst-Pack, 0.0 for DOS), direction 1
  real(REAL64), intent(in)           ::  sy                              !<  offset of the integration mesh (0.5 for Monkhorst-Pack, 0.0 for DOS), direction 2
  real(REAL64), intent(in)           ::  sz                              !<  offset of the integration mesh (0.5 for Monkhorst-Pack, 0.0 for DOS), direction 3

  logical, intent(in)                ::  ldtau                           !<  also calculates d tauk / d bdot for use in stress
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  neig                            !<  number of eigenvectors
  real(REAL64), intent(in)           ::  occ_x_wgk(mxdbnd)               !<  occupation*k-weight*spin deg. of eigenvector j
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  |psi>

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  irk                             !<  k-point index
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

  integer, intent(in)                ::  nrk                             !<  number of k-points for integration in the irreducible wedge of the brillouin zone
  real(REAL64), intent(in)           ::  rk(3,mxdnrk)                    !<  component in lattice coordinates of the k-point in the mesh
  real(REAL64), intent(in)           ::  wgk(mxdnrk)                     !<  weight in the integration of k-point
  integer, intent(in)                ::  kmap(3,nx,ny,nz)                !<  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

! output

  complex(REAL64), intent(out)       ::  tauk(mxdnst)                    !<  symmetrized kinetic energy density (in stars of G)
  real(REAL64), intent(out)          ::  dtau_dbdot(3,3,mxdscr)          !<  d tau / d bdot on the real space mesh of the potential (hartree/bohr^3), bdot includes the 2 pi factors

! local allocatable arrays

  complex(REAL64), allocatable       ::  tauu(:)                         !  unsymmetrized density (in G)

  complex(REAL64), allocatable       ::  taumsh(:)
  real(REAL64), allocatable          ::  wrkfft(:)
  complex(REAL64), allocatable       ::  chd(:,:)
  integer, allocatable               ::  ipoint(:)

  complex(REAL64), allocatable       ::  chds(:,:)                       !  i (k+G) psi on the real space mesh of the potential
  real(REAL64), allocatable          ::  wrks(:)                         !  work array for the fft on that mesh
  integer, allocatable               ::  ipoints(:)                      !  gather-scatter index for that mesh

  complex(REAL64), allocatable       ::  psi_tr(:,:)                     !  transformed psi
  integer, allocatable               ::  isort_tr(:)                     !  isort for transformed k-point
  integer, allocatable               ::  igofkg(:,:,:)                   !  index of G-vector with components (k1,k2,k3)

! local variables

  integer         ::  mxdfft, mxdwrk
  integer         ::  jmin, jmax
  integer         ::  nsfft(3)
  integer         ::  n1, n2, n3, id
  integer         ::  nn1, nn2, nn3
  integer         ::  kd1, kd2, kd3
  integer         ::  k1, k2, k3
  integer         ::  ntot, it
  complex(REAL64) ::  xp

  real(REAL64)    ::  qk(3)
  real(REAL64)    ::  vcell, bdot(3,3)

  real(REAL64)    ::  rkpt_tr(3)                                         !  k-point obtained by symmetry (transformed)
  integer         ::  jsym                                               !  symmetry operation that brings rkpt to rkpt_tr
  integer         ::  kmap2                                              !  1 if there is an additional inversion (time reversal)
  integer         ::  mtxd_tr                                            !  mtxd for transformed point
  integer         ::  kgshift(3)
  integer         ::  kgv_rot(3)

  real(REAL64)    ::  umxyz                                              !  weight of a point of the full mesh
  real(REAL64)    ::  fac
  real(REAL64)    ::  dx, dy, dz
  real(REAL64)    ::  rk_loc(3)                                          !  k-point corresponding to ix, iy, iz
  real(REAL64)    ::  diff
  integer         ::  npoint                                             !  number of mesh points equivalent to irk

  integer         ::  mxdffs, mxdwrs                                     !  dimensions for the fft on the mesh of the potential
  integer         ::  ns1, ns2, ns3, ids                                 !  mesh of the potential, packing (ids,ns2,ns3)
  integer         ::  nns1, nns2, nns3
  integer         ::  ks1, ks2, ks3
  integer         ::  ind
  real(REAL64)    ::  facv

! counters

  integer         ::  i, j, n, m
  integer         ::  ix, iy, iz
  integer         ::  i1, i2, i3

! parameters

  real(REAL64), parameter :: SMALL = 1.0E-12_REAL64
  real(REAL64), parameter :: EPS = 1.0E-8_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_I = cmplx(ZERO,UM,REAL64)


  call adot_to_bdot(adot, vcell, bdot)

  do i = 1,ns
    tauk(i) = C_ZERO
  enddo

  if(ldtau) then
    if(kmscr(7)*kmscr(5)*kmscr(6) > mxdscr) then
      write(6,*)
      write(6,*) '   STOPPED in tau_by_fft_stress:  mesh of the potential ',     &
                 kmscr(7)*kmscr(5)*kmscr(6), ' larger than mxdscr = ', mxdscr

      stop

    endif
    dtau_dbdot(:,:,:) = ZERO
  endif

  if (neig < 1) return

! find min and max band with nonzero occupation

  jmin = 0
  jmax = 0
  do i = 1,neig
    if (abs(occ_x_wgk(i)) > SMALL .and. jmin == 0) jmin = i
    if (abs(occ_x_wgk(i)) > SMALL .and. jmin /= 0) jmax = i
  enddo

  if (jmin == 0) return

! paranoid check of the reference k-point

  if(abs( (rkpt(1)-rk(1,irk))*(rkpt(1)-rk(1,irk)) +                      &
          (rkpt(2)-rk(2,irk))*(rkpt(2)-rk(2,irk)) +                      &
          (rkpt(3)-rk(3,irk))*(rkpt(3)-rk(3,irk)) ) > EPS) then
    write(6,*)
    write(6,*) '   STOPPED in tau_by_fft_stress:  irk = ', irk
    write(6,*) '   rk(:,irk) = ', rk(1,irk), rk(2,irk), rk(3,irk)
    write(6,*) '   rkpt =      ', rkpt(1), rkpt(2), rkpt(3)

    stop

  endif

  allocate(psi_tr(mxddim,mxdbnd))
  allocate(isort_tr(mxddim))

! index of G-vectors, to rotate the basis of the reference k-point

  allocate(igofkg(-kmax(1):kmax(1),-kmax(2):kmax(2),-kmax(3):kmax(3)))

  igofkg(:,:,:) = 0
  do i = 1,ng
    igofkg(kgv(1,i),kgv(2,i),kgv(3,i)) = i
  enddo

! find n for fast fourier transform
! ni is the number of points used in direction i.

  call size_fft(kmax, nsfft, mxdfft, mxdwrk)

  allocate(chd(mxdfft,3))
  allocate(wrkfft(mxdwrk))

  n1 = nsfft(1)
  n2 = nsfft(2)
  n3 = nsfft(3)
  id = n1
  ntot = id * n2 * n3
  nn1 = (n1-1) / 2
  nn2 = (n2-1) / 2
  nn3 = (n3-1) / 2

  allocate(ipoint(mxddim))

  allocate(taumsh(mxdfft))

  taumsh(:) = C_ZERO

! mesh of the potential for d tau / d bdot

  if(ldtau) then
    call size_fft(kmscr, nsfft, mxdffs, mxdwrs)
    ns1 = kmscr(4)
    ns2 = kmscr(5)
    ns3 = kmscr(6)
    ids = kmscr(7)
    nns1 = (ns1-1) / 2
    nns2 = (ns2-1) / 2
    nns3 = (ns3-1) / 2
    mxdffs = max(mxdffs, ids*ns2*ns3)
    allocate(chds(mxdffs,3))
    allocate(wrks(mxdwrs))
    allocate(ipoints(mxddim))
  endif

! each point of the full mesh has weight umxyz = 1/(nx*ny*nz).
! occ_x_wgk / wgk(irk) is occupation*spin degeneracy.

  umxyz = UM / (nx*ny*nz)

  dx = UM / nx
  dy = UM / ny
  dz = UM / nz

! start loop over equivalent k-points of the full mesh

  npoint = 0

  do ix = 1,nx
  do iy = 1,ny
  do iz = 1,nz

    if(kmap(1,ix,iy,iz) /= irk) cycle

    npoint = npoint + 1

!   same expression as in int_pnt

    rk_loc(1) = ((ix-1) + sx) * dx
    rk_loc(2) = ((iy-1) + sy) * dy
    rk_loc(3) = ((iz-1) + sz) * dz

    if(abs( (rk_loc(1)-rk(1,irk))*(rk_loc(1)-rk(1,irk)) +                &
            (rk_loc(2)-rk(2,irk))*(rk_loc(2)-rk(2,irk)) +                &
            (rk_loc(3)-rk(3,irk))*(rk_loc(3)-rk(3,irk)) ) < EPS) then

!     it is the reference point

      rkpt_tr(:) = rkpt(:)
      mtxd_tr = mtxd
      isort_tr(1:mtxd) = isort(1:mtxd)
      psi_tr(1:mtxd,1:neig) = psi(1:mtxd,1:neig)

    else

      jsym = kmap(3,ix,iy,iz)
      kmap2 = kmap(2,ix,iy,iz)

      if(jsym < 1 .or. jsym > ntrans) then
        write(6,*)
        write(6,*) '   STOPPED in tau_by_fft_stress:  wrong symmetry operation ', jsym
        write(6,*) '   ix, iy, iz = ', ix, iy, iz

        stop

      endif

!     transformed point, k_tr = mtrx k   (or  -mtrx k  with time reversal)

      do j = 1,3
        rkpt_tr(j) = mtrx(j,1,jsym)*rkpt(1) + mtrx(j,2,jsym)*rkpt(2) +   &
                     mtrx(j,3,jsym)*rkpt(3)
      enddo
      if(kmap2 == 1) rkpt_tr(:) = -rkpt_tr(:)

!     paranoid check, rkpt_tr and rk_loc should differ by a reciprocal lattice vector

      do j = 1,3
        diff = rkpt_tr(j) - rk_loc(j)
        if(abs(diff - nint(diff)) > EPS) then
          write(6,*)
          write(6,*) '   STOPPED in tau_by_fft_stress:  inconsistent kmap'
          write(6,*) '   ix, iy, iz = ', ix, iy, iz, '  jsym, kmap2 = ', jsym, kmap2
          write(6,*) '   rkpt_tr =    ', rkpt_tr(1), rkpt_tr(2), rkpt_tr(3)
          write(6,*) '   rk_loc =     ', rk_loc(1), rk_loc(2), rk_loc(3)

          stop

        endif
      enddo

!     the basis of the transformed point is the rotated basis of the reference point,
!     so the set of k+G is exactly the same (independent of cutoff rounding or fixed k+G)

      mtxd_tr = mtxd
      do i = 1,mtxd
        do j = 1,3
          kgv_rot(j) = mtrx(j,1,jsym)*kgv(1,isort(i)) +                  &
                       mtrx(j,2,jsym)*kgv(2,isort(i)) +                  &
                       mtrx(j,3,jsym)*kgv(3,isort(i))
        enddo
        if(kmap2 == 1) kgv_rot(:) = -kgv_rot(:)
        if(abs(kgv_rot(1)) > kmax(1) .or. abs(kgv_rot(2)) > kmax(2) .or.  &
           abs(kgv_rot(3)) > kmax(3)) then
          isort_tr(i) = 0
        else
          isort_tr(i) = igofkg(kgv_rot(1),kgv_rot(2),kgv_rot(3))
        endif
        if(isort_tr(i) == 0) then
          write(6,*)
          write(6,*) '   STOPPED in tau_by_fft_stress:  rotated G-vector not found'
          write(6,*) '   G = ', kgv(1,isort(i)), kgv(2,isort(i)), kgv(3,isort(i))
          write(6,*) '   rotated = ', kgv_rot(1), kgv_rot(2), kgv_rot(3)

          stop

        endif
      enddo

      kgshift(:) = 0

      call psi_rot_inv_shift(mtrx(:,:,jsym), tnp(:,jsym),                &
          kmap2, kgshift, neig,                                          &
          rkpt, mtxd, isort, psi,                                        &
          rkpt_tr, mtxd_tr, isort_tr, psi_tr,                            &
          ng, kgv,                                                       &
          mxdgve, mxddim, mxdbnd)

    endif

!   fills array ipoint (gather-scatter index) for the transformed basis

    kd1 = 0
    kd2 = 0
    kd3 = 0
    do i = 1,mtxd_tr
      it = isort_tr(i)
      k1 = kgv(1,it)
      if (iabs(k1) > kd1) kd1 = iabs(k1)
      if (k1 < 0) k1 = n1 + k1
      k2 = kgv(2,it)
      if (iabs(k2) > kd2) kd2 = iabs(k2)
      if (k2 < 0) k2 = n2 + k2
      k3 = kgv(3,it)
      if (iabs(k3) > kd3) kd3 = iabs(k3)
      if (k3 < 0) k3 = n3 + k3
      ipoint(i) = (k3*n2 + k2)*id + k1 + 1
    enddo
    if (kd1 > nn1 .or. kd2 > nn2 .or. kd3 > nn3) then
      write(6,*)
      write(6,'("     STOPPED in tau_by_fft_stress.  size of matrix:    ", &
         &      i7," fft mesh ",3i5)') mtxd_tr, n1, n2, n3

      stop

    endif

!   same for the mesh of the potential.  The wave-function should fit
!   (dual approximation), products are then exact at the mesh points.

    if(ldtau) then
      ks1 = 0
      ks2 = 0
      ks3 = 0
      do i = 1,mtxd_tr
        it = isort_tr(i)
        k1 = kgv(1,it)
        if (iabs(k1) > ks1) ks1 = iabs(k1)
        if (k1 < 0) k1 = ns1 + k1
        k2 = kgv(2,it)
        if (iabs(k2) > ks2) ks2 = iabs(k2)
        if (k2 < 0) k2 = ns2 + k2
        k3 = kgv(3,it)
        if (iabs(k3) > ks3) ks3 = iabs(k3)
        if (k3 < 0) k3 = ns3 + k3
        ipoints(i) = (k3*ns2 + k2)*ids + k1 + 1
      enddo
      if (ks1 > nns1 .or. ks2 > nns2 .or. ks3 > nns3) then
        write(6,*)
        write(6,'("     STOPPED in tau_by_fft_stress.  wave-function does not ", &
           &      "fit in the mesh of the potential ",3i5)') ns1, ns2, ns3

        stop

      endif
    endif

!   start loop over eigenvectors

    do j = jmin,jmax

!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        chd(i,1) = C_ZERO
        chd(i,2) = C_ZERO
        chd(i,3) = C_ZERO
      enddo
!$omp end parallel do

!     components of i (k+G) psi in reciprocal lattice coordinates

!$omp parallel do default(shared) private(i,qk)
      do i = 1,mtxd_tr
        qk(1) = rkpt_tr(1) + real(kgv(1,isort_tr(i)),REAL64)
        qk(2) = rkpt_tr(2) + real(kgv(2,isort_tr(i)),REAL64)
        qk(3) = rkpt_tr(3) + real(kgv(3,isort_tr(i)),REAL64)

        chd(ipoint(i),1) = C_I*qk(1)*psi_tr(i,j)
        chd(ipoint(i),2) = C_I*qk(2)*psi_tr(i,j)
        chd(ipoint(i),3) = C_I*qk(3)*psi_tr(i,j)
      enddo
!$omp end parallel do

!     Fourier transform to real space

      call cfft_wf_c16(chd(:,1), id,n1,n2,n3, kd1,kd2,kd3, -1, wrkfft,mxdwrk)
      call cfft_wf_c16(chd(:,2), id,n1,n2,n3, kd1,kd2,kd3, -1, wrkfft,mxdwrk)
      call cfft_wf_c16(chd(:,3), id,n1,n2,n3, kd1,kd2,kd3, -1, wrkfft,mxdwrk)

!     |grad psi|^2 using the reciprocal space metric (includes 2 pi factors).
!     occ_x_wgk(j)/wgk(irk) is occupation * spin degeneracy, the 1/2 is 1/2 p^2 / m

      fac = occ_x_wgk(j) * umxyz / (2*wgk(irk))

!$omp parallel do default(shared) private(i,xp)
      do i = 1,ntot
        xp = conjg(chd(i,1))*(bdot(1,1)*chd(i,1) + bdot(1,2)*chd(i,2) + bdot(1,3)*chd(i,3)) +   &
             conjg(chd(i,2))*(bdot(2,1)*chd(i,1) + bdot(2,2)*chd(i,2) + bdot(2,3)*chd(i,3)) +   &
             conjg(chd(i,3))*(bdot(3,1)*chd(i,1) + bdot(3,2)*chd(i,2) + bdot(3,3)*chd(i,3))
        taumsh(i) = taumsh(i) + fac*real(xp,REAL64)
      enddo
!$omp end parallel do

!     d tau / d bdot(n,m) = fac Re( conjg(d_n psi) d_m psi ) / vcell  on the mesh of the potential.
!     tau is linear in bdot, so sum_nm bdot(n,m) d tau / d bdot(n,m) = tau.

      if(ldtau) then

        facv = fac / vcell

        do n = 1,3

!$omp parallel do default(shared) private(i)
          do i = 1,ids*ns2*ns3
            chds(i,n) = C_ZERO
          enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i,qk)
          do i = 1,mtxd_tr
            qk(n) = rkpt_tr(n) + real(kgv(n,isort_tr(i)),REAL64)
            chds(ipoints(i),n) = C_I*qk(n)*psi_tr(i,j)
          enddo
!$omp end parallel do

          call cfft_wf_c16(chds(:,n), ids,ns1,ns2,ns3, ks1,ks2,ks3, -1,  &
              wrks, mxdwrs)

        enddo

!$omp parallel do default(shared) private(i1,i2,i3,ind,n,m)
        do i3 = 1,ns3
        do i2 = 1,ns2
        do i1 = 1,ns1
          ind = i1 + (i2-1)*ids + (i3-1)*ids*ns2
          do m = 1,3
          do n = 1,3
            dtau_dbdot(n,m,ind) = dtau_dbdot(n,m,ind) +                  &
                facv*real(conjg(chds(ind,n))*chds(ind,m),REAL64)
          enddo
          enddo
        enddo
        enddo
        enddo
!$omp end parallel do

      endif

    enddo                  !  end loop over eigenvectors

  enddo
  enddo
  enddo                    !  end loop over equivalent k-points

! paranoid check, the number of equivalent points should agree with the weight

  if(abs(npoint*umxyz - wgk(irk)) > EPS) then
    write(6,*)
    write(6,*) '   STOPPED in tau_by_fft_stress:  number of equivalent points'
    write(6,*) '   irk = ', irk, '  npoint = ', npoint, '  wgk*nx*ny*nz = ', wgk(irk)*nx*ny*nz

    stop

  endif

! Fourier transform to momentum space

  call rfft_c16(taumsh, id,n1,n2,n3, 1, wrkfft,mxdwrk)

  allocate(tauu(mxdgve))

  call gvec_mesh_fold(tauu, taumsh, id,n1,n2,n3,                         &
      ng, kgv,                                                           &
      mxdgve, mxdfft)

  deallocate(ipoint)
  deallocate(chd)
  deallocate(wrkfft)
  deallocate(taumsh)

  deallocate(psi_tr)
  deallocate(isort_tr)
  deallocate(igofkg)

  if(ldtau) then
    deallocate(chds)
    deallocate(wrks)
    deallocate(ipoints)
  endif

  call gvec_star_of_g_fold(tauk, tauu, .FALSE.,                          &
       ng, phase, conj, ns, inds, mstar,                                 &
       mxdgve, mxdnst)

  deallocate(tauu)

  return

end subroutine tau_by_fft_stress
