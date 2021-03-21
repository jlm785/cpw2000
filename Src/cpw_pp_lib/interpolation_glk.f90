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

!>  Interpolates the bands between nrk_int points.
!>  k=0 is the point to be interpolated, k = 1,2,...,nrk_int
!>  are the points used in interpolation


subroutine interpolation_glk(nrk_int, emax, neig, xsvd, csvd,            &
    xw_all, rkpt_all, mtxd_all, neig_all, isort_all, psi_all,            &
    ei, psi, hpsi, mtxd, isort, qmod, ekpg,                              &
    ng, kgv,                                                             &
    ntype, natom, rat, adot,                                             &
    nqnl, delqnl, vkb, nkb,                                              &
    vscr, kmscr,                                                         &
    mxdtyp, mxdatm, mxddim, mxdlqp, mxdbnd, mxdgve, mxdscr)

! Written 23 August 2020. JLM
! Based on interpolation_svd and out_band_dos
! Modified, qmod-->ekpg in hk_psi. 13 February 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

! version 4.989

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr


  integer, intent(in)                ::  nrk_int                         !<  number of points used in the integration
  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  integer, intent(in)                ::  neig                            !<  number of eigenvectors required
  real(REAL64), intent(in)           ::  xsvd                            !<  singular values < xsvd are discarded
  real(REAL64), intent(in)           ::  csvd                            !<  use maximum of csvd*neig states after decomposition

  real(REAL64), intent(in)           ::  xw_all(nrk_int)                 !<  relative coordinates of the k-point to be interpolated
  real(REAL64), intent(in)           ::  rkpt_all(3,nrk_int)             !<  coordinates of k-points used for the interpolation

  integer, intent(in)                ::  mtxd_all(nrk_int)               !<  wavefunction dimension (basis size) for k-point
  integer, intent(in)                ::  neig_all(nrk_int)               !<  number of wavefunctions for k-point
  integer, intent(in)                ::  isort_all(mxddim,nrk_int)       !<  g-vector associated with row/column i of hamiltonian

  complex(REAL64), intent(in)        ::  psi_all(mxddim,mxdbnd,nrk_int)  !<  wavevector of points used in the interpolation




  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. NOT normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

! output

  real(REAL64)                       ::  ei(mxdbnd)
  complex(REAL64), intent(out)       ::  psi(mxddim,mxdbnd)              !<  wavevector
  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)             !<  |hpsi> =  V_NL |psi>
  integer, intent(out)               ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(out)               ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  real(REAL64), intent(out)          ::  qmod(mxddim)                    !<  length of k+g-vector of row/column i
  real(REAL64), intent(out)          ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i


! allocatable arrays

  complex(REAL64), allocatable       ::  anlga(:,:)                      !  KB projectors without spin-orbit
  real(REAL64), allocatable          ::  xnlkb(:)                        !  KB normalization without spin-orbit

  complex(REAL64), allocatable       ::  psi_svd(:,:)
  complex(REAL64), allocatable       ::  hpsi_svd(:,:)
  real(REAL64), allocatable          ::  eig_svd(:)

  complex(REAL64), allocatable       ::  vec(:,:)                        !  overlap matrix S
  complex(REAL64), allocatable       ::  hred(:,:)                       !  reduced hamiltonian

  real(REAL64), allocatable          ::  singval(:)                      !  singular values

! local variables

  real(REAL64)           ::  rkpt(3)                                     !  k-point to be interpolated

  real(REAL64)           ::  xsum

  integer                ::  nanl                                        !  number of Kleinman-Bylander projectorswithout spin-orbit.
  integer                ::  mxdanl
  integer                ::  nanlso                                      !  number of Kleinman-Bylander projectors with spin-orbit.

  integer                ::  neig_tot
  integer                ::  i1, i2

  integer                ::  nglk, nglkmax                               !  number of states used from SVD output

  integer                ::  info

  logical                ::  lnewanl

! constants

  real(REAL64), parameter     :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  :: C_UM = cmplx(UM,ZERO,REAL64)
  real(REAL64), parameter     :: EPS = 0.000001_REAL64                   !  criteria for error

! counters

  integer    ::  i, j, n, k


! some checks

  if(nrk_int < 2) then
    write(6,'("   stopped in interpolation_glk: nrk = ",i7)') nrk_int

    stop

  endif

  xsum = ZERO
  do j = 1,nrk_int
    xsum = xsum + xw_all(j)
  enddo

  if(abs(xsum - UM) > EPS) then
    write(6,'("   stopped in interpolation_glk:")')
    write(6,'("   sum of weights is not 1")')

    stop

  endif

  do j=1,3
    rkpt(j) = ZERO
    do n = 1,nrk_int
      rkpt(j) = rkpt(j) + xw_all(n)*rkpt_all(j,n)
    enddo
  enddo

  call hamilt_struct(emax, rkpt, mtxd, isort, qmod, ekpg, .FALSE.,       &
      ng, kgv, adot,                                                     &
      mxdgve, mxddim)

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, mxdtyp)

  mxdanl = nanl

  allocate(xnlkb(mxdanl))
  allocate(anlga(mxddim,mxdanl))

  call proj_nl_kb_c16(rkpt, mtxd, isort, nanl,                           &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlga, xnlkb,                                                      &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

! constructs psi_svd

  allocate(psi_svd(mxddim,mxdbnd*nrk_int))
  allocate(hpsi_svd(mxddim,mxdbnd*nrk_int))
  allocate(eig_svd(mxdbnd*nrk_int))
  allocate(singval(mxdbnd*nrk_int))

  neig_tot = 0
  do k = 1,nrk_int

    neig_tot = neig_tot + neig_all(k)

    if(neig_tot > mxdbnd*nrk_int) then
      write(6,'("   Stopped in interpolation_glk,  ",                    &
              "neig_tot > mxdbnd*nrk_int: ",2i5)')                       &
              neig_tot, mxdbnd*nrk_int

      stop

    endif
  enddo

  neig_tot = 0
  do k = 1,nrk_int

    i1 = neig_tot + 1
    i2 = neig_tot + neig_all(k)

    call psi_convert(neig_all(k), mtxd_all(k), isort_all(:,k),           &
        psi_all(:,:,k), mtxd , isort , psi_svd(:,i1:i2),                 &
        mxddim, mxdbnd)

    do j = i1,i2
    do i = 1,mtxd
      psi_svd(i,j) = xw_all(k) * psi_svd(i,j)
    enddo
    enddo

    neig_tot = neig_tot + neig_all(k)

  enddo

! singular value decomposition

  call svd_c16(psi_svd, singval, mtxd, neig_tot, info,                   &
      mxddim, mxdbnd*nrk_int)

  if(info /= 0) then
    write(6,'("   Stopped in interpolation_glk:  ",                      &
               "svd_c16 info = ",i5)') info

    stop

  endif


! finds how many states should be used

  nglkmax = nint(csvd*(neig))
  nglkmax = max(nglkmax,neig)
  nglkmax = min(nglkmax,2*neig)

  nglk = neig
  if(nglkmax > neig) then
    do j = neig+1,nglkmax
      if(singval(j) < xsvd) exit
      nglk = nglk + 1
    enddo
  endif

! diagonalizes in SVD sub-space

  lnewanl = .TRUE.

  call hk_psi_c16(mtxd, nglk, psi_svd, hpsi_svd, lnewanl,                &
      ng,  kgv,                                                          &
      ekpg, isort, vscr, kmscr,                                          &
      anlga, xnlkb, nanl,                                                &
      mxddim, mxdbnd*nrk_int, mxdanl, mxdgve, mxdscr)

  deallocate(xnlkb)
  deallocate(anlga)

  allocate(vec(nglk,nglk))
  allocate(hred(nglk,nglk))

  call zgemm('c','n', nglk, nglk, mtxd, C_UM, psi_svd, mxddim,           &
      hpsi_svd, mxddim, C_ZERO ,hred, nglk)

  call diag_c16(nglk, hred, eig_svd, vec, nglk, info)

  if(info /= 0) then
    write(6,'("   Stopped in interpolation_glk:  ",                      &
               "diag_c16 info = ",i5)') info

    stop

  endif

  do n = 1,neig
    ei(n) = eig_svd(n)
  enddo

  call zgemm('n','n', mtxd, neig, nglk, C_UM, psi_svd, mxddim,           &
      vec, nglk, C_ZERO, psi, mxddim)
  call zgemm('n','n', mtxd, neig, nglk, C_UM, hpsi_svd, mxddim,          &
      vec, nglk, C_ZERO, hpsi, mxddim)

  deallocate(vec)
  deallocate(hred)

  deallocate(psi_svd)
  deallocate(hpsi_svd)
  deallocate(singval)

  return
end subroutine interpolation_glk
