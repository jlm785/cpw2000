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

!>  Calculates the degeneracy structure of eigenvalues
!>
!>  \author       Jose Luis Martins
!>  \version      5.08
!>  \date         9 October 2023.
!>  \copyright    GNU Public License v2

subroutine berry_degeneracy(ldims, neigin, neig, ei, tol,                &
     nlevel, maxdeg, levdeg, leveigs,                                    &
     mxdbnd, mxdlev, mxddeg)


! resets neig according to degeneracy.  9 November 2023. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdlev                          !<  array dimension for number of levels
  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels

  logical, intent(in)                ::  ldims                           !<  if true only calculates the needed dimensions

  integer, intent(in)                ::  neigin                          !<  number of eigenvalues available
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalue no. i. (hartree)
  real(REAL64), intent(in)           ::  tol                             !<  tolerance for identical eigenvalues

! output

  integer, intent(out)               ::  nlevel                          !<  number of energy levels
  integer, intent(out)               ::  maxdeg                          !<  maximum degeneragy present
  integer, intent(out)               ::  levdeg(mxdlev)                  !<  degeneragy of level
  integer, intent(out)               ::  leveigs(mxdlev,mxddeg)          !<  points to degenerate level

! input and output

  integer, intent(inout)             ::  neig                            !<  on input the number of desired states.  On output a larger number that takes into account the degeneracies

! counters

  integer       ::  n, ndegcount


  if(neigin > mxdbnd) then
    write(6,*)
    write(6,*) '   STOPPED in berry degeneracy'
    write(6,'("   neigin = ",i5,"  mxdbnd = ",i5)') neigin, mxdbnd
    write(6,*)

    stop

  endif

  if(ldims) then

!   finds dimensions.  only returns nlevel and maxdeg

    nlevel = 1
    maxdeg = 1
    ndegcount = 1

    if(neigin > 1) then
      do n = 2,neigin

        if(ei(n-1) - ei(n) > TOL) then
          write(6,*)
          write(6,*) '  STOPPED in berry degeneracy'
          write(6,*) '  Eigenvalues are not in increasing order'
          write(6,*)

          stop

        endif

        if(ei(n) - ei(n-1) > TOL) then
          nlevel = nlevel + 1
          ndegcount = 1
        else
          ndegcount = ndegcount + 1
          if(ndegcount > maxdeg) maxdeg = ndegcount
        endif

      enddo
    endif

  else

!   fills the information

    if(neig > neigin) then
      write(6,*)
      write(6,*) '   STOPPED in berry degeneracy'
      write(6,'("   neigin = ",i5,"  neig = ",i5)') neigin, neig
      write(6,*)

      stop

    endif

    levdeg(1) = 1
    leveigs(1,1) = 1
    nlevel = 1
    ndegcount = 1

    if(neigin > 1) then
      do n = 2,neigin

        if(ei(n) - ei(n-1) > TOL) then
          nlevel = nlevel + 1
          ndegcount = 1
        else
          ndegcount = ndegcount + 1
        endif

        if(nlevel > mxdlev) then
          write(6,*)
          write(6,*) '   STOPPED in berry degeneracy'
          write(6,'("   nlevel = ",i5,"  mxdlev = ",i5)') nlevel, mxdlev
          write(6,*)

          stop

        endif

        if(ndegcount > mxddeg) then
          write(6,*)
          write(6,*) '   STOPPED in berry degeneracy'
          write(6,'("   ndegcount = ",i5,"  mxddeg = ",i5)') ndegcount, mxddeg
          write(6,*)

          stop

        endif

        levdeg(nlevel) = ndegcount
        leveigs(nlevel,ndegcount) = n

      enddo
    endif

!   recalculates neig and nlevel

    ndegcount = 0
    do n = 1,nlevel
      ndegcount = ndegcount + levdeg(n)
      if(ndegcount >= neig) then
        neig = ndegcount
        nlevel = n

        exit

      endif
    enddo

  endif

  return

end subroutine berry_degeneracy
