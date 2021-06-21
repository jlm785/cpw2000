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

!>  tests the space group of the crystal

subroutine sym_test(ipr, tol, istatus,                                   &
  ntrans, mtrx, tnp,                                                     &
  ntype, natom, rat, adot,                                               &
  mxdtyp,mxdatm)

! Written 18 March 2004. JLM
! Modified, f90, 7 June 2014. JLM
! Modified, documentation, December 2019. JLM
! Modified, stop decision. 29 December 2020. JLM
! copyright INESC-MN/Jose Luis Martins

! version 4.99


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms

  integer, intent(in)                ::  ipr                             !<  print flag
  real(REAL64), intent(in)           ::  tol                             !<  tolerance for symmetry recognition

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

! output

  integer, intent(out)               ::  istatus                         !<  istatus = 0, successful; 1 not closed; 2 no inverse; 3 inconsistent with atomic positions

! local variables

  integer              ::  irotdir(3,3,48)
  real(REAL64)         ::  vcell,bdot(3,3)
  integer              ::  irot(3,3)
  real(REAL64)         ::  frac(3),tmp(3,3)

  integer              ::  icount,ires(48,48),jsucc
  real(REAL64)         ::  rat2(3),sdif,dif
 
! parameters

  real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counters

  integer i, j, k, n, n1, n2, n3, j2


! initial calculations

  istatus = 0

  call adot_to_bdot(adot,vcell,bdot)

  do i=1,3
  do j=1,3
    bdot(i,j) = bdot(i,j)/ (2*PI*2*PI)
  enddo
  enddo

  do n=1,ntrans
    do i=1,3
    do j=1,3
      tmp(i,j) = mtrx(i,1,n)*adot(1,j) +                                 &
                 mtrx(i,2,n)*adot(2,j) +                                 &
                 mtrx(i,3,n)*adot(3,j)
    enddo
    enddo

    do i=1,3
    do j=1,3
      irotdir(i,j,n) = nint(bdot(1,i)*tmp(1,j) +                         &
                            bdot(2,i)*tmp(2,j) +                         &
                            bdot(3,i)*tmp(3,j))
    enddo
    enddo
  enddo

! constructs multiplication table

  do n1=1,ntrans
  do n2=1,ntrans

    do i=1,3
    do j=1,3
      irot(i,j) = irotdir(i,1,n1)*irotdir(1,j,n2) +                      &
                  irotdir(i,2,n1)*irotdir(2,j,n2) +                      &
                  irotdir(i,3,n1)*irotdir(3,j,n2)
    enddo
    enddo

    do i=1,3
      frac(i) = irotdir(i,1,n1)*tnp(1,n2) +                              &
                irotdir(i,2,n1)*tnp(2,n2) +                              &
                irotdir(i,3,n1)*tnp(3,n2) + tnp(i,n1)
      frac(i) = frac(i) / (2*PI)
    enddo

    ires(n1,n2) = 0
    do n3=1,ntrans
      icount = 0
      do i=1,3
      do j=1,3
        if(irot(i,j) == irotdir(i,j,n3)) then
          icount = icount + 1
        endif
      enddo
      enddo
      sdif = ZERO
      do i=1,3
        dif = (tnp(i,n3) / (2*PI) - frac(i))
        dif = dif - nint(dif)
        sdif = sdif + abs(dif)
      enddo
      if(icount == 9 .and. sdif < 10.0*tol) then
        ires(n1,n2) = n3
      endif
    enddo
  enddo
  enddo

  if(ipr > 0) then
    write(6,*)
    write(6,'(12x,"multiplication table for factor group")')
    write(6,*)
    do n1=1,ntrans
      write(6,'(48i4)') (ires(n1,n2),n2=1,ntrans)
    enddo
  endif

! checks if it is well defined

  do n1=1,ntrans
  do n2=1,ntrans
    if(ires(n1,n2) == 0) then
      write(6,'("    WARNING in sym_test, the product of ",              &
        & "operations ",2i5," is not in the list")') n1,n2
      istatus = 1

      return

    endif
  enddo
  enddo

! checks there is an inverse

  if(ntrans > 1) then
    do n1=1,ntrans
    do n2=1,ntrans-1
    do n3=n2+1,ntrans
      if(ires(n1,n2) == ires(n1,n3)) then
        write(6,'("    WARNING in sym_test,  the product of ",           &
          & "operations ",2i5,"and ",2i5," is the same")') n1,n2,n1,n3
        istatus = 2

        return

      endif
      if(ires(n2,n1) == ires(n3,n1)) then
        write(6,'("    WARNING in sym_test,  the product of ",           &
          & "operations ",2i5,"and ",2i5," is the same")') n2,n1,n3,n1
        istatus = 2

        return

      endif
    enddo
    enddo
    enddo
  endif

! unity and associativity verified by construction


! checks the transformed atoms are in the right place

  do n1=1,ntrans

    do i=1,ntype
    do j=1,natom(i)

      do k=1,3
        rat2(k) = irotdir(k,1,n1)*rat(1,j,i) +                           &
                  irotdir(k,2,n1)*rat(2,j,i) +                           &
                  irotdir(k,3,n1)*rat(3,j,i) + tnp(k,n1) / (2*PI)
      enddo

      jsucc = 0
      do j2=1,natom(i)
        sdif = ZERO
        do k=1,3
          dif = rat2(k) - rat(k,j2,i)
          dif = dif - nint(dif)
          sdif = sdif + abs(dif)
        enddo
        if(sdif < tol) then
          jsucc = 1

          exit

        endif
      enddo

      if(jsucc == 0) then
        write(6,'("    WARNING in sym_test for operation ",i5,           &
          & "the",i4,"th atom of type",i5," has no match")') n1,j,i
        istatus = 3

        return

      endif
    enddo
    enddo

  enddo

  return
end subroutine sym_test
