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
!>  States are treated all at once.
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.06
!>  \date         13 January 2023.
!>  \copyright    GNU Public License v2


subroutine berry_stern_solve_all(mtxd, neig, psi, ei, dhdkpsi, dpsi, tol,  &
    isort, ekpg,                                                         &
    vscr, kmscr,                                                         &
    ng, kgv,                                                             &
    nanl, anlga, xnlkb,                                                  &
    mxddim, mxdbnd, mxdgve, mxdscr, mxdanl)

! adapted from psi_vnl_psi_der, psi_p_psi and CLR phonon hk_psi_nl_lr_c16


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: LMAX = 3                                !   hard coded max ang. mom.

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension for screening potential
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  wavefunction dimension

  real(REAL64), intent(in)           ::  tol                             !<  tolerance solution

  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates

  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for
  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh

  integer, intent(in)                ::  nanl                            !<  half of number of projectors without spin
  complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)            !<  KB projectors without spin-orbit
  real(REAL64), intent(in)           ::  xnlkb(mxdanl)                   !<  KB normalization without spin-orbit

  complex(REAL64), intent(in)        ::  psi(mxddim, mxdbnd)             !<  |psi> (in principle eigen-functions)
  real(REAL64), intent(in)           ::  ei(mxddim)                      !<  eigenvalue no. i. (hartree)

  complex(REAL64), intent(in)        ::  dhdkpsi(mxddim,mxdbnd,3)        !<  (d H /d k) |Psi>

! output

  complex(REAL64), intent(out)       ::  dpsi(mxddim, mxdbnd, 3)         !<  d |psi> / d k


! local allocatable variables

  REAL(REAL64), allocatable          ::  x(:), b(:)
  REAL(REAL64), allocatable          ::  tmp(:,:)
  complex(REAL64), allocatable       ::  ac(:,:), bc(:,:)

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

! external

  complex(REAL64) ,external   :: zdotc
  real(REAL64) ,external   :: dnrm2


! solves the equation

  nitmax = max(150, 2*mtxd)

  allocate(b(2*mtxd*neig), x(2*mtxd*neig))
  allocate(tmp(2*mtxd*neig,4))
  allocate(ac(mtxd,1), bc(mtxd,1))

! for compatibility with intel MKL

  allocate(ipar(128),dpar(128))

! loop over directions

  do j = 1,3


    do n = 1,neig

!     initial guess from perturbation theory
!     doesn't help much...

      ac(:,1) = C_ZERO

      do m = 1,neig
        if(n /= m) then
          if(abs(ei(n)-ei(m)) < 0.0001) then
            write(6,*) '   stopped in berry_stern_solve_all.'
            write(6,*) '   not safe for degenerate states'

            stop

          endif
          zz = zdotc(mtxd, psi(:,m), 1, dhdkpsi(:,n,j), 1)
          zz = zz / (ei(n) - ei(m))
          call zaxpy(mtxd, zz, psi(:,m), 1, ac(:,1), 1)
        endif
      enddo

      do k = 1,mtxd
        x(2*k-1) = real(ac(k,1),REAL64)
        x(2*k  )   = dimag(ac(k,1))
      enddo

      bc(:,1) = dhdkpsi(1:mtxd,n,j)

      zz = zdotc(mtxd, psi(:,n), 1, bc(:,1), 1)
      call zaxpy(mtxd, -zz, psi(:,n), 1, bc(:,1), 1)

      do k = 1,mtxd
        b(2*k-1) = real(bc(k,1),REAL64)
        b(2*k) = dimag(bc(k,1))
      enddo
      xxb = dnrm2(2*mtxd, b, 1)

      lfound = .FALSE.

      call dcg_init (2*mtxd, x, b, job, ipar, dpar, tmp )

      do niter = 1,nitmax

        call dcg (2*mtxd, x, b, job, ipar, dpar, tmp )

        if(job == 1) then

!         matrix-vector multiply

          do k = 1,mtxd
            ac(k,1) = cmplx(tmp(2*k-1, 1 ),tmp(2*k, 1 ),REAL64)
          enddo

          zz = zdotc(mtxd, psi(:,n), 1, ac(:,1), 1)
          call zaxpy(mtxd, -zz, psi(:,n), 1, ac(:,1), 1)

          call hk_psi_c16(mtxd, 1, ac, bc, .TRUE.,                       &
             ng, kgv,                                                    &
             ekpg, isort, vscr, kmscr,                                   &
             anlga, xnlkb, nanl,                                         &
             mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

          call zaxpy(mtxd, cmplx(-ei(n),ZERO,REAL64), ac(:,1), 1, bc(:,1), 1)
          call zscal(mtxd, -C_UM, bc(:,1), 1)

          zz = zdotc(mtxd, psi(:,n), 1, bc(:,1), 1)
          call zaxpy(mtxd, -zz, psi(:,n), 1, bc(:,1), 1)

          do k = 1,mtxd
            tmp(2*k-1, 2 ) = real(bc(k,1),REAL64)
            tmp(2*k, 2 )   = dimag(bc(k,1))
          enddo

        elseif(job == 2) then

!         check convergence

!          xx = dnrm2(2*mtxd, tmp(:, 3 ), 1)

          xx = sqrt(dpar(5))

!          if ( xxb == ZERO ) then
            if ( xx <= tol ) then
              lfound = .TRUE.
              write(6,'(i5,5x,"  error = ",e14.4)') niter, xx

              exit

            end if
!           else
!             if ( xx <= tol * xxb ) then
!               lfound = .TRUE.
! !dbg              write(6,'(i5,5x,"  error = ",e14.4)') niter, xx
!
!               exit
!
!             end if
!           end if

        elseif(job == 3) then

!         preconditioning

          do k = 1,mtxd

            xpre = ekpg(k)

            if(xpre > UM) then
              xpre = UM / xpre
            else
              xpre = UM
            endif

            ac(k,1) = 1.3*cmplx(-xpre*tmp(2*k-1, 3 ),-xpre*tmp(2*k, 3 ),REAL64)

          enddo

          zz = zdotc(mtxd, psi(:,n), 1, ac(:,1), 1)
          call zaxpy(mtxd, -zz, psi(:,n), 1, ac(:,1), 1)


          do k = 1,mtxd
            tmp(2*k-1, 4 ) = real(ac(k,1),REAL64)
            tmp(2*k, 4 )   = dimag(ac(k,1))
          enddo

         endif

      enddo

      if(.not. lfound) then
        write(6,*) '   stopped in berry_stern_solve_all.'
        write(6,*) '   maximum number of iterations exceeded'

        stop

      endif

      do k = 1,mtxd
        dpsi(k,n,j) =  cmplx(x(2*k-1), x(2*k),REAL64)
      enddo

    enddo

  enddo

! end of loop over directions

  deallocate(b, x)
  deallocate(tmp)
  deallocate(ac, bc)

  deallocate(ipar,dpar)

  return

end subroutine berry_stern_solve_all

