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

!>  Solves the Sternheimer equation with a RCI (reverse communication interface)
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.09
!>  \date         13 January 2023, 15 december 2023.
!>  \copyright    GNU Public License v2


subroutine berry_stern_solve_spin(mtxd, neig, psi_sp, ei, dhdkpsi_sp,    &
    dpsi_sp, tol,                                                        &
    nlevel, levdeg, leveigs,                                             &
    isort, ekpg,                                                         &
    vscr_sp, kmscr, nsp,                                                 &
    ng, kgv,                                                             &
    nanlsp, anlsp, xnlkbsp,                                              &
    mxddim, mxdbnd, mxdgve, mxdscr, mxdasp, mxdlev, mxddeg, mxdnsp)

! adapted from psi_vnl_psi_der, psi_p_psi and CLR phonon hk_psi_nl_lr_c16
! spin version 15 December 2023. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension for screening potential
  integer, intent(in)                ::  mxdasp                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdlev                          !<  array dimension for number of levels
  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels
  integer, intent(in)                ::  mxdnsp                          !<  array dimension for number of spin components (1,2,4)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (not counting spin)
  integer, intent(in)                ::  neig                            !<  number of bands (including spin)

  real(REAL64), intent(in)           ::  tol                             !<  tolerance solution

  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates

  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for
  integer, intent(in)                ::  nsp                             !<  number of spin components ox xc-potential (1,2,4)
  real(REAL64), intent(in)           ::  vscr_sp(mxdscr,mxdnsp)          !<  screened potential in the fft real space mesh with spin components

  integer, intent(in)                ::  nanlsp                          !<  number of projectors with spin
  complex(REAL64), intent(in)        ::  anlsp(2*mxddim,mxdasp)          !<  KB projectors with spin-orbit
  real(REAL64), intent(in)           ::  xnlkbsp(mxdasp)                 !<  KB normalization with spin-orbit

  integer, intent(in)                ::  nlevel                          !<  number of energy levels
  integer, intent(in)                ::  levdeg(mxdlev)                  !<  degeneragy of level
  integer, intent(in)                ::  leveigs(mxdlev,mxddeg)          !<  points to degenerate level

  complex(REAL64), intent(in)        ::  psi_sp(2*mxddim, mxdbnd)        !<  |psi_sp> (in principle eigen-functions)
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalue no. i. (hartree)

  complex(REAL64), intent(in)        ::  dhdkpsi_sp(2*mxddim,mxdbnd,3)   !<  (d H /d k) |Psi>

! output

  complex(REAL64), intent(out)       ::  dpsi_sp(2*mxddim, mxdbnd, 3)    !<  d |psi_sp> / d k


! local allocatable variables

  REAL(REAL64), allocatable          ::  x(:), b(:)
  REAL(REAL64), allocatable          ::  tmp(:,:)
  complex(REAL64), allocatable       ::  ac(:,:), bc(:,:), acguess(:,:)

  integer, allocatable               ::  ipar(:)                         !  should not be changed during iterations
  REAL(REAL64), allocatable          ::  dpar(:)                         !  should not be changed during iterations

! local variables

  complex(REAL64)      ::  zz

  real(REAL64)      ::  xx, xxb

  integer           ::  job
  integer           ::  niter, nitmax
  logical           ::  lfound

  real(REAL64)      ::  xpre

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer    ::  j, k, n, m
  integer    ::  nl, nk, ml, mk

! external

  complex(REAL64) ,external   :: zdotc
  real(REAL64) ,external   :: dnrm2


! solves the equation

  nitmax = max(150, 2*2*mtxd)

  allocate(b(2*2*mtxd), x(2*2*mtxd))
  allocate(tmp(2*2*mtxd,4))
  allocate(ac(2*mtxd,1), bc(2*mtxd,1), acguess(2*mtxd,1))

! for compatibility with intel MKL

  allocate(ipar(128),dpar(128))

  do nl = 1,nlevel
  do nk = 1,levdeg(nl)
    n = leveigs(nl,nk)

    do j = 1,3

!     initial guess from perturbation theory
!     doesn't help much...

      ac(:,1) = C_ZERO

      do ml = 1,nlevel
        if(ml /= nl) then
          do mk = 1,levdeg(ml)
            m = leveigs(ml,mk)
            zz = zdotc(2*mtxd, psi_sp(:,m), 1, dhdkpsi_sp(:,n,j), 1)
            zz = zz / (ei(n) - ei(m))
            call zaxpy(2*mtxd, zz, psi_sp(:,m), 1, ac(:,1), 1)
          enddo
        endif
      enddo

      acguess(:,1) = ac(:,1)

      ac(:,1) = C_ZERO

      do k = 1,2*mtxd
        x(2*k-1) = real(ac(k,1),REAL64)
        x(2*k  ) = dimag(ac(k,1))
      enddo

      bc(:,1) = dhdkpsi_sp(1:2*mtxd,n,j)

      call berry_project_one('O', psi_sp, bc(:,1), 2*mtxd, neig,         &
            2*mxddim, mxdbnd)

      do k = 1,2*mtxd
        b(2*k-1) = real(bc(k,1),REAL64)
        b(2*k) = dimag(bc(k,1))
      enddo
      xxb = dnrm2(2*2*mtxd, b, 1)

      lfound = .FALSE.

      call dcg_init (2*2*mtxd, x, b, job, ipar, dpar, tmp )

      do niter = 1,nitmax

        call dcg (2*2*mtxd, x, b, job, ipar, dpar, tmp )

        if(job == 1) then

!         matrix-vector multiply

          do k = 1,2*mtxd
            ac(k,1) = cmplx(tmp(2*k-1, 1 ),tmp(2*k, 1 ),REAL64)
          enddo

          call berry_project_one('O', psi_sp, ac(:,1), 2*mtxd, neig,     &
            2*mxddim, mxdbnd)

          call hk_psi_spin_c16(mtxd, 1, ac, bc, .TRUE.,                  &
             ng, kgv,                                                    &
             ekpg, isort, vscr_sp, kmscr, nsp,                           &
             anlsp, xnlkbsp, nanlsp,                                     &
             mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

          call zaxpy(2*mtxd, cmplx(-ei(n),ZERO,REAL64), ac(:,1), 1, bc(:,1), 1)
          call zscal(2*mtxd, -C_UM, bc(:,1), 1)

!         |bc> = (E_n - H) |ac>. now orthogonalize within same energy

          call berry_project_one('O', psi_sp, bc(:,1), 2*mtxd, neig,     &
            2*mxddim, mxdbnd)

          do k = 1,2*mtxd
            tmp(2*k-1, 2 ) = real(bc(k,1),REAL64)
            tmp(2*k, 2 )   = dimag(bc(k,1))
          enddo

        elseif(job == 2) then

!         check convergence

!          xx = dnrm2(2*mtxd, tmp(:, 3 ), 1)

          xx = sqrt(dpar(5))

          if ( xxb == ZERO ) then
            if ( xx <= tol ) then
              lfound = .TRUE.
!dbg              write(6,'(i5,5x,"  error = ",e14.4)') niter, xx

              exit

            end if
          else
            if ( xx <= tol * xxb ) then
              lfound = .TRUE.
!dbg              write(6,'(i5,5x,"  error = ",e14.4)') niter, xx

              exit

            end if
          end if

        elseif(job == 3) then

!         preconditioning

          do k = 1,2*mtxd

            xpre = ekpg((k+1)/2)

            if(xpre > UM) then
              xpre = UM / xpre
            else
              xpre = UM
            endif

            ac(k,1) = 1.3*cmplx(-xpre*tmp(2*k-1, 3 ),-xpre*tmp(2*k, 3 ),REAL64)
!            ac(k,1) = -cmplx(tmp(2*k-1, 3 ),tmp(2*k, 3 ),REAL64)

          enddo

          call berry_project_one('O', psi_sp, ac(:,1), 2*mtxd, neig,     &
            2*mxddim, mxdbnd)


          do k = 1,2*mtxd
            tmp(2*k-1, 4 ) = real(ac(k,1),REAL64)
            tmp(2*k, 4 )   = dimag(ac(k,1))
          enddo

         endif

      enddo

      if(.not. lfound) then
        write(6,*) '   maximum number of iterations exceeded'

        write(6,*) ' n = ',n,' j = ',j
        write(6,'("  xx = ",e14.3)') xx

        stop

      endif

      do k = 1,2*mtxd
        dpsi_sp(k,n,j) =  cmplx(x(2*k-1), x(2*k),REAL64)
      enddo

      call berry_project_one('O', psi_sp, dpsi_sp(:,n,j), 2*mtxd, neig,  &
            2*mxddim, mxdbnd)

      dpsi_sp(1:2*mtxd,n,j) = dpsi_sp(1:2*mtxd,n,j) + acguess(1:2*mtxd,1)

    enddo

  enddo
  enddo

  deallocate(b, x)
  deallocate(tmp)
  deallocate(ac, bc, acguess)

  deallocate(ipar,dpar)

  return

end subroutine berry_stern_solve_spin

