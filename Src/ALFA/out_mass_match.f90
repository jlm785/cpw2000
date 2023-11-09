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

!>  Matches the left bands to the right side bands across a degeneracy
!>  Used for interpolation of multivalued functions
!>
!>  \author       Jose Luis Martins
!>  \version      5.08
!>  \date         7 November 2023.
!>  \copyright    GNU Public License v2

subroutine out_mass_match(neig, npt, ei_l, ipl,                          &
    mxdbnd)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  npt                             !<  2*npt+1 is the total number of interpolation points

  integer, intent(in)                ::  neig                            !<  wavefunction dimension

  real(REAL64), intent(in)           ::  ei_l(mxdbnd,-npt:npt)           !<  eigenvalue no. i. in the line (hartree).  ei_l(j,0) is the reference point

! output

  integer, intent(out)               ::  ipl(mxdbnd)                     !<  pointer to band at left that matches the one on right

! allocatable arrays with larger scope


  integer, allocatable               ::  levdeg(:)                       !  degeneracy of energy level
  integer, allocatable               ::  leveigs(:,:)                    !  states belonging to level

  real(REAL64), allocatable          ::  xin(:), yin(:)

  integer, allocatable               ::  ipk(:), iperm(:)        !  most similar levels.  Permutation.

! local variables


  integer           ::  nlevel                                           !  number of energy levels (not counting degeneracies)
  integer           ::  maxdeg                                           !  maximum number of degeneracies
  integer           ::  mxdlev                                           !  array dimension for number of levels
  integer           ::  mxddeg                                           !  array dimension for number of levels

  real(REAL64)      ::  y(0:2), dy(0:2)

  real(REAL64)      ::  dx
  logical           ::  lperm

  integer           ::  neigtmp

! constants

  real(REAL64), parameter     ::  TOL = 1.0E-8_REAL64

! counters

  integer    ::  j, n, m
  integer    ::  nl, nk, mk

! external

  complex(REAL64) ,external   :: zdotc


  do n = 1,neig
    ipl(n) = n
  enddo

  neigtmp = neig

  if(npt > 0) then

!   checks for degeneracies

    allocate(levdeg(1))
    allocate(leveigs(1,1))

    call berry_degeneracy(.TRUE., neig, neigtmp, ei_l(:,0), TOL,           &
           nlevel, maxdeg, levdeg, leveigs,                                &
           mxdbnd, 1, 1)

    mxdlev = nlevel
    mxddeg = maxdeg

    deallocate(levdeg)
    deallocate(leveigs)

    allocate(levdeg(mxdlev))
    allocate(leveigs(mxdlev,mxddeg))

    call berry_degeneracy(.FALSE., neig, neigtmp, ei_l(:,0), TOL,          &
           nlevel, maxdeg, levdeg, leveigs,                                &
           mxdbnd, mxdlev, mxddeg)

!   tries to match left and right sides

    allocate(xin(-npt:npt))
    allocate(yin(-npt:npt))

    allocate(ipk(mxddeg))
    allocate(iperm(mxddeg))

    do j = -npt,0
      xin(j) = (2*npt+j)
    enddo

    do n = 1,neig
      ipl(n) = n
    enddo

    do nl = 1,nlevel
      if(levdeg(nl) > 1) then
        do nk = 1,levdeg(nl)

          ipk(nk) = nk
          n = leveigs(nl,nk)
          do j = -npt,0
            yin(j) = ei_l(n,j+npt) - ei_l(n,0)
          enddo

          call poly_interp(y, dy, xin, yin, npt, 0)

          dx = abs(y(0) - (ei_l(n,-npt) - ei_l(n,0)))
          do mk = 1,levdeg(nl)
            if(mk /= nk) then
              m = leveigs(nl,mk)
              if(abs(y(0) - (ei_l(m,-npt)- ei_l(n,0))) + TOL < dx) then
                dx = abs(y(0) - ei_l(m,-npt))
                ipk(nk) = mk
               endif
            endif
          enddo

        enddo

!       checks it is a permutation

        lperm = .TRUE.
        do nk = 1,levdeg(nl)
          iperm(nk) = 0
        enddo
        do nk = 1,levdeg(nl)
          iperm(ipk(nk)) = iperm(ipk(nk)) + 1
        enddo
        do nk = 1,levdeg(nl)
          if(iperm(nk) /= 1) then
            lperm = .FALSE.
            exit
          endif
        enddo

        if(lperm) then
          do nk = 1,levdeg(nl)
            n = leveigs(nl,nk)
            m = leveigs(nl,ipk(nk))
            ipl(n) = m
          enddo
        endif

      endif
    enddo

    deallocate(levdeg)
    deallocate(leveigs)

    deallocate(xin,yin)
    deallocate(ipk,iperm)

  endif

  return

end subroutine out_mass_match





