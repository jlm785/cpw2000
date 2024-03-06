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

!>  Calculates the inverse effective mass d^2 E_i / d k^2 in a given direction.
!>  Can be used for other rank 2 (not counting eigenstate index) Berry quantities
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         4 November 2023. 5 March 2024.
!>  \copyright    GNU Public License v2

subroutine berry_effective_mass(xk, adot, psidhdkpsi, tmass,             &
    d2eidxk2,                                                            &
    nlevel, levdeg, leveigs,                                             &
    mxdbnd, mxdlev, mxddeg)

! Written 4 November 2023. JLM
! Modified, first order perturbation must be diagonalized.  5 March 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdlev                          !<  array dimension for number of levels
  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels

  real(REAL64), intent(in)           ::  xk(3)                           !<  k-direction in reciprocal lattice coordinates

  integer, intent(in)                ::  nlevel                          !<  number of energy levels
  integer, intent(in)                ::  levdeg(mxdlev)                  !<  degeneragy of level
  integer, intent(in)                ::  leveigs(mxdlev,mxddeg)          !<  points to degenerate level

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  complex(REAL64), intent(in)     ::  psidhdkpsi(mxddeg,mxddeg,3,mxdlev) !<  <psi_n| d H / d k |psi_m> for each energy level (lattice coordinates)
  complex(REAL64), intent(in)        ::  tmass(mxddeg,mxddeg,3,3,mxdlev) !<  Tensor associated with effective mass d^2 E_i /d k_1 d k_2 (lattice coordinates)

! output

  real(REAL64), intent(out)          ::  d2eidxk2(mxdbnd)                !<  d^2 E_i /d xk^2  in xk direction

! allocatable arrays

  complex(REAL64), allocatable       ::  h1(:,:)                         !  first order perturbation
  complex(REAL64), allocatable       ::  h1eigvec(:,:)                   !  eigenvectors of first order
  real(REAL64), allocatable          ::  h1eig(:)                        !  eigenvalues of first order

  complex(REAL64), allocatable       ::  h2(:,:)                         !  second order perturbation
  complex(REAL64), allocatable       ::  h2sub(:,:)                      !  second order perturbation on sublevel
  complex(REAL64), allocatable       ::  h2eigvec(:,:)                   !  eigenvectors
  real(REAL64), allocatable          ::  h2eig(:)                        !  eigenvalues

  integer, allocatable               ::  levdegh1(:)                     !  degeneracy of sublevel
  integer, allocatable               ::  leveigh1(:,:)                   !  points to degenerate level

  complex(REAL64), allocatable       ::  tmp(:,:)                        !  temporary matrix

! local variables

  real(REAL64)      ::  vcell, bdot(3,3)
  real(REAL64)      ::  yk(3), xknorm
  integer           ::  info

  integer           ::  nlh1                         !  number of sublevels
  integer           ::  nh1, nnh1, maxdeg

! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
  real(REAL64), parameter     ::  EPS = 1.0E-12_REAL64
  real(REAL64), parameter     ::  TOL = 1.0E-5_REAL64

! counters

  integer    ::  i, j, n
  integer    ::  nl, nk, mk
  integer    ::  jl


! normalizes xk

  call adot_to_bdot(adot, vcell, bdot)

  xknorm = ZERO
  do i = 1,3
  do j = 1,3
    xknorm = xknorm + xk(i)*bdot(i,j)*xk(j)
  enddo
  enddo

  if(xknorm < EPS) then
    write(6,*)
    write(6,*) "   STOPPED in berry_efective_mass, null direction  ",xk(1),xk(2),xk(3)
    write(6,*)

    stop

  endif

  do i = 1,3
    yk(i) = xk(i) / sqrt(xknorm)
  enddo

! mxddeg is small no point in optimizing

  allocate(h1(mxddeg,mxddeg))
  allocate(h1eigvec(mxddeg,mxddeg))
  allocate(h1eig(mxddeg))

  allocate(h2(mxddeg,mxddeg))
  allocate(h2sub(mxddeg,mxddeg))
  allocate(h2eigvec(mxddeg,mxddeg))
  allocate(h2eig(mxddeg))

  allocate(levdegh1(mxddeg))
  allocate(leveigh1(mxddeg,mxddeg))

  allocate(tmp(mxddeg,mxddeg))

! loop over energy levels

  do nl = 1,nlevel

    if(levdeg(nl) == 1) then

!     no degeneracy

      n = leveigs(nl,1)
      d2eidxk2(n) = ZERO
      do i = 1,3
      do j = 1,3
        d2eidxk2(n) = d2eidxk2(n) + yk(i)*real(tmass(1,1,i,j,nl),REAL64)*yk(j)
      enddo
      enddo

    else

!     Degenerate case.  First diagonalize the first order perturbation

      do nk = 1,levdeg(nl)
      do mk = 1,levdeg(nl)
        h1(nk,mk) = C_ZERO
        do j = 1,3
          h1(nk,mk) = h1(nk,mk) + psidhdkpsi(nk,mk,j,nl)*yk(j)
        enddo
      enddo
      enddo

      call diag_c16(levdeg(nl), h1, h1eig, h1eigvec, mxddeg, info)
      if(info /= 0) then

        write(6,*) "   STOPPED in berry_effective_mass first order  info = ",info

        stop

      endif

      nh1 = levdeg(nl)
      nnh1 = nh1

      call berry_degeneracy(.FALSE., nh1, nnh1, h1eig, TOL,              &
          nlh1, maxdeg, levdegh1, leveigh1,                              &
          mxddeg, mxddeg, mxddeg)

!     calculates second order matrix

      do nk = 1,nh1
      do mk = 1,nh1
        h2(nk,mk) = C_ZERO
        do i = 1,3
        do j = 1,3
          h2(nk,mk) = h2(nk,mk) + yk(i)*tmass(nk,mk,i,j,nl)*yk(j)
        enddo
        enddo
      enddo
      enddo

!     rotates second order matrix according to the eigenvectors of the first order

      call zgemm('n','n', nh1, nh1, nh1, C_UM, h2, mxddeg,               &
                 h1eigvec, mxddeg, C_ZERO, tmp, mxddeg)
      call zgemm('c','n', nh1, nh1, nh1, C_UM, h1eigvec, mxddeg,         &
                 tmp, mxddeg, C_ZERO, h2, mxddeg)

!     iterates over sublevels

      do jl = 1,nlh1

        if(levdegh1(jl) == 1) then

          j = leveigh1(jl,1)
          n = leveigs(nl,j)

          d2eidxk2(n) = h2(j,j)

        else

          do i = 1,levdegh1(jl)
          do j = 1,levdegh1(jl)
            h2sub(i,j) = h2(leveigh1(jl,i),leveigh1(jl,j))
          enddo
          enddo

          call diag_c16(levdegh1(jl), h2sub, h2eig, h2eigvec, mxddeg, info)
          if(info /= 0) then

            write(6,*) "   STOPPED in berry_effective_mass second order  info = ",info

            stop

          endif

          do nk = 1,levdegh1(jl)
            j = leveigh1(jl,nk)
            n = leveigs(nl,j)
            d2eidxk2(n) = h2eig(nk)
          enddo

        endif

      enddo

!     end of loop over sub-levels

    endif

  enddo

! end of loop over levels

  deallocate(h1)
  deallocate(h1eigvec)
  deallocate(h1eig)

  deallocate(h2)
  deallocate(h2sub)
  deallocate(h2eigvec)
  deallocate(h2eig)

  deallocate(levdegh1)
  deallocate(leveigh1)

  deallocate(tmp)

  return

end subroutine berry_effective_mass
