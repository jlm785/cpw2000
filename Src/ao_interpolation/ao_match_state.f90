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

!>  imatch points to the old state that is most similar to the new states
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         29 October 2024
!>  \copyright    GNU Public License v2

subroutine ao_match_state(imatch, xover, mtxd,                           &
      neig_old, psi_old, ei_old, neig_new, psi_new, ei_new,              &
      mxddim, mxdbnd_old, mxdbnd_new)


! Adapted from out_band_match_state, 29 October 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd_old                      !<  array dimension for number of bands
  integer, intent(in)                ::  mxdbnd_new                      !<  array dimension for number of bands

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian (basis set size)

  integer, intent(in)                ::  neig_old                        !<  number of old eigenvectors
  complex(REAL64), intent(in)        ::  psi_old(mxddim, mxdbnd_old)     !<  component j of old eigenvectors
  real(REAL64), intent(in)           ::  ei_old(mxdbnd_old)              !<  eigenvalue of old eigenvectors

  integer, intent(in)                ::  neig_new                        !<  number of new eigenvectors
  complex(REAL64), intent(in)        ::  psi_new(mxddim, mxdbnd_new)     !<  component j of new eigenvectors
  real(REAL64), intent(in)           ::  ei_new(mxdbnd_new)              !<  eigenvalue of new eigenvectors

! output

  integer, intent(out)               ::  imatch(mxdbnd_new)              !<  points to the old state that is most similar to the new state
  real(REAL64), intent(out)          ::  xover(mxdbnd_new)               !<  overlap of new state with old state

! local allocatable arrays

  complex(REAL64), allocatable       ::  prod(:,:)                       !  old and new eigenvector products
  real(REAL64), allocatable          ::  prodsq(:,:)                     !  old and new eigenvector products, modulo square

  integer, allocatable               ::  levdeg(:)                       !  degeneracy of energy level
  integer, allocatable               ::  leveigs(:,:)                    !  states belonging to level

  integer, allocatable               ::  indlev(:)                       !  index within a level of new state

! local variables

  real(REAL64)      ::  xprod

  integer           ::  neig
  integer           ::  mxdlev                                           !  array dimension for number of levels
  integer           ::  mxddeg                                           !  array dimension for number of levels

  integer           ::  nlevel                                           !  number of energy levels (not counting degeneracies)
  integer           ::  maxdeg                                           !  maximum number of degeneracies

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
  real(REAL64), parameter     ::  EPS = 1.0E-8_REAL64
  real(REAL64), parameter     ::  TOL = 1.0E-5_REAL64
  real(REAL64), parameter     ::  EV = 27.21138505_REAL64

! counters

  integer   :: i, j, nl, nc

! The structure of the basis is the same, but the number of bands is different.
! This is the difference between this subroutine and out_band_match_state.


! calculate the products

  allocate(prod(neig_old,neig_new))
  allocate(prodsq(neig_old,neig_new))

  call zgemm('C', 'N', neig_old, neig_new, mtxd, C_UM, psi_old, mxddim,  &
          psi_new, mxddim, C_ZERO, prod, neig_old)

  do i = 1,neig_old
  do j = 1,neig_new
    prodsq(i,j) = prod(i,j)*conjg(prod(i,j))
  enddo
  enddo

! paranoid check things are consistent

  do j = 1,neig_new
    xprod = ZERO
    do i = 1,neig_old
      xprod = xprod + prodsq(i,j)
    enddo

    if(xprod > UM+EPS) then

      write(6,*)
      write(6,*)
      write(6,'("   STOPPED in ao_match_state: ")')
      write(6,'("   projection is greater than one! ")')
      write(6,*)

      exit

    endif

  enddo

! takes into account the degeneracy of old states.


  allocate(levdeg(1))
  allocate(leveigs(1,1))

  neig = neig_old

  call berry_degeneracy(.TRUE., neig_old, neig, ei_old, TOL,             &
         nlevel, maxdeg, levdeg, leveigs,                                &
         mxdbnd_old, 1, 1)

  mxdlev = nlevel
  mxddeg = maxdeg

  deallocate(levdeg)
  deallocate(leveigs)

  allocate(levdeg(mxdlev))
  allocate(leveigs(mxdlev,mxddeg))

  allocate(indlev(mxddeg))

! fills the information

  call berry_degeneracy(.FALSE., neig_old, neig, ei_old, TOL,            &
         nlevel, maxdeg, levdeg, leveigs,                                &
         mxdbnd_old, mxdlev, mxddeg)

! Finds projections greater than 1/2

  do i = 1,neig_new
    imatch(i) = 0
    do nl = 1,nlevel
      xprod = prodsq(leveigs(nl,1),i)
      if(levdeg(nl) > 1) then
        do j = 2,levdeg(nl)
          xprod = xprod + prodsq(leveigs(nl,j),i)
        enddo
      endif

      xover(i) = xprod
      if(xprod > UM/2) then
        imatch(i) = leveigs(nl,1)

        exit

      endif

    enddo

  enddo

! distribute within a degenerate level (probably not necessary...)

  do nl = 1,nlevel
    if(levdeg(nl) > 1) then

      nc = 0
      do i = 1,neig_new
        if(imatch(i) == leveigs(nl,1)) then
          nc = nc + 1
!         paranoid check
          if(nc > levdeg(nl)) then
            write(6,*)
            write(6,*) "  WARNING   inconsistent degeneracies in ao_match_state"
            write(6,*) "  check atomic orbital state ",i
            write(6,*)
          else
            indlev(nc) = i
          endif
        endif
      enddo

      if(nc > 1) then
        do j = 2,nc
          imatch(indlev(j)) = leveigs(nl,j)
        enddo
      endif

    endif
  enddo

! warns if there is a big shift in eigenvalues

  do i = 1,neig_new
    if(imatch(i) > 0) then
      if(abs(ei_new(i) - ei_old(imatch(i))) > UM/EV) then
            write(6,*)
            write(6,*) "  WARNING   correction larger than 1 eV in ao_match_state"
            write(6,*) "  check atomic orbital state ",i," and pw state ",imatch(i)
            write(6,'(2(5x,f12.5))') ei_new(i)*EV, ei_old(imatch(i))*EV
            write(6,*)
      endif
    endif
  enddo

  deallocate(prod)
  deallocate(prodsq)

  deallocate(levdeg)
  deallocate(leveigs)

  deallocate(indlev)

  return

end subroutine ao_match_state
