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
!>  \date         18 October 2013. 10 June 2020
!>  \copyright    GNU Public License v2

subroutine out_band_match_state(imatch, neig,                            &
      mtxdold, psiold, isold, mtxdnew, psi, isnew,                       &
      mxddim, mxdbnd, mxdold)

! version 4.53. 18 October 2013. jlm
! modified, bug with icmplx, 7 November 2018.
! modified fast version, documentation, 28 May 2019. JLM
! Dimension of psiold, 1 August 2019. JLM
! After implementing diagonalization in real case
! revert to old code! Or not!
! Removed icmplx, reintroduce for rpw, 10 June 2020. JLM
! Indentation, 29 October 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdold                          !<  psiold mxddim

  integer, intent(in)                ::  neig                            !<  number of eigenvectors
  integer, intent(in)                ::  mtxdold                         !<  number of coefficients in old vectors
  complex(REAL64), intent(in)        ::  psiold(mxdold, mxdbnd)          !<  component j of old eigenvectors
  integer, intent(in)                ::  isold(mxddim)                   !<  g-vector associated with row/column i of old hamiltonian
  integer, intent(in)                ::  mtxdnew                         !<  number of coefficients in new vectors
  complex(REAL64), intent(in)        ::  psi(mxddim, mxdbnd)             !<  component j of eigenvector i
  integer, intent(in)                ::  isnew(mxddim)                   !<  g-vector associated with row/column i of new hamiltonian

! output

  integer, intent(out)               ::  imatch(mxdbnd)                  !<  points to the old state that is most similar to the new state

! local allocatable arrays

  integer, allocatable               ::  ibold(:),ibnew(:)               !  reverse mapping of isold and isnew
  complex(REAL64), allocatable       ::  prod(:,:)                       !  old and new eigenvector products
  real(REAL64), allocatable          ::  prodsq(:)                       !  old and new eigenvector products, modulo square
  integer, allocatable               ::  ind(:)                          !  index array
  integer, allocatable               ::  irow(:),icol(:)                 !  indicates which states match

  complex(REAL64), allocatable       ::  psicomm(:,:)                    !  psi with alterated G vector ordering
  complex(REAL64), allocatable       ::  psioldcomm(:,:)                 !  psiold with alterated G vector ordering


! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer   :: i, j, ij, nmax, nc

! The structure of the basis may not be the same, so it is necessary
! to make them compatible.

  nmax = 0
  do i = mtxdold,1,-1
    if(isold(i) > nmax) nmax = isold(i)
  enddo

  do i = min(mtxdold,mtxdnew),1,-1
    if(isnew(i) > nmax) nmax = isnew(i)
  enddo

  allocate(ibold(nmax),ibnew(nmax))

  do i=1,nmax
    ibold(i) = 0
    ibnew(i) = 0
  enddo
  do i=1,mtxdold
    ibold(isold(i)) = i
  enddo
  do i=1,min(mtxdold,mtxdnew)
    ibnew(isnew(i)) = i
  enddo

  do i=1,nmax
    j = i
    if( ibold(i) == 0  .or. ibnew(i) == 0 ) exit
  enddo
  nmax = j - 1

! calculate the products

  allocate(prod(neig,neig))
  allocate(prodsq(neig*neig))

  allocate(psicomm(nmax,neig))
  allocate(psioldcomm(nmax,neig))

  do i = 1,neig
  do j = 1,nmax
    psicomm(j,i) = psi(ibnew(j),i)
  enddo
  enddo

  do i = 1,neig
  do j = 1,nmax
    psioldcomm(j,i) = psiold(ibold(j),i)
  enddo
  enddo

  call zgemm('C', 'N', neig, neig, nmax, C_UM, psioldcomm, nmax,         &
          psicomm, nmax, C_ZERO, prod, neig)

  do i = 1,neig
  do j = 1,neig
    ij = (i-1)*neig + j
    prodsq(ij) = prod(i,j)*conjg(prod(i,j))
  enddo
  enddo

  deallocate(psicomm)
  deallocate(psioldcomm)

! sorts the products

  allocate(ind(neig*neig))
  allocate(irow(neig),icol(neig))

  call sort(neig*neig,prodsq,ind)

  do i = 1,neig
    irow(i) = 0
    icol(i) = 0
  enddo

  nc = 0
  do ij = neig*neig,1,-1
    i = (ind(ij)-1)/neig + 1
    j = mod(ind(ij)-1,neig) + 1
    if(irow(i) == 0 .and. icol(j) == 0) then
      nc = nc + 1
      irow(i) = j
      icol(j) = i
    endif
    if(nc == neig) exit
  enddo

  do i = 1,neig
    imatch(icol(i)) = i
  enddo

! paranoid check that imatch is a permutation

  do i = 1,neig
    icol(i) = 0
  enddo
  do i = 1,neig
    icol(imatch(i)) = icol(imatch(i)) + 1
  enddo
  do i = 1,neig
    if(icol(i) /=1) then
      write(6,*)
      write(6,*)
      write(6,'("   WARNING   WARNING in out_band_match_state: ")')
      write(6,'("   imatch is not a permutation ")')
      write(6,*)
      write(6,*)
      do j = 1,neig
        imatch(j) = j
      enddo

      exit

    endif
  enddo

  deallocate(ibold,ibnew)
  deallocate(prod,prodsq)
  deallocate(ind)
  deallocate(irow,icol)

  return

end subroutine out_band_match_state
