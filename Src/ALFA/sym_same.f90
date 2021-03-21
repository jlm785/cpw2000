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

!> checks if the symmetry has changed.  Orphaned subroutine...

subroutine sym_same(lnewsym, isym, ipr, tol,                           &
                      ntrans, mtrx, tnp,                                 &
                      adot, ntype, natom, rat,                           &
                      mxdtyp, mxdatm)



! Written August 2019. jlm
! Modified, 11 November 2020. Does not stop. JLM
! Modified, does not modify adot,rat. 29 December 2020. JLM
! copyright INESC-MN/Jose Luis Martins

! version 4.99


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms

  integer, intent(in)                ::  isym                            !<  if 1 calculate symmetry, if 0 only identity is used
  integer, intent(in)                ::  ipr                             !<  print flag
  real(REAL64), intent(in)           ::  tol                             !<  tolerance for symmetry recognition

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! output

  logical , intent(out)              ::  lnewsym                         !<  indicates if the symmetry changed

! local variables

  real(REAL64)                       ::  adotnew(3,3)                    !  local metric in direct space
  integer                            ::  ntransnew                       !  new number of symmetry operations in the factor group
  real(REAL64), allocatable          ::  ratnew(:,:,:)                   !  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer                            ::  mtrxnew(3,3,48)                 !  new rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation
  real(REAL64)                       ::  tnpnew(3,48)                    !  new 2*pi* i-th component (in lattice coordinates) of the fractional translation vector of the k-th symmetry operation


  real(REAL64)                 ::  dif
 
! parameters

  real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter :: EPS = 1.0E-8_REAL64

! counters

  integer i, j, n


  allocate(ratnew(3,mxdatm,mxdtyp))
  
  adotnew(:,:) = adot(:,:)
  ratnew(:,:,:) = rat(:,:,:)

  call sym_identify(isym, ipr, tol,                                      &
                    ntransnew, mtrxnew, tnpnew,                          &
                    adotnew, ntype, natom, ratnew,                       &
                    mxdtyp, mxdatm)

  deallocate(ratnew)

  lnewsym = .FALSE.

  if(ntransnew /= ntrans) then

    if(ipr > 0) then
      write(6,*)
      write(6,'("    sym_same:  ntransnew /= ntrans",2i6)') ntrans, ntransnew
      write(6,*)
    endif

    lnewsym = .TRUE.

  else

    do n = 1,ntrans
      do i = 1,3
        dif = abs(tnpnew(i,n) - tnp(i,n)) / (2*PI)
        dif = dif - int(dif)
        if(dif > EPS) then
          if(ipr > 0) then
            write(6,*)
            write(6,'("   sym_same: tnp /= tnpnew",2i6,3f12.4)')           &
       &                 i,n,tnp(i,n),tnpnew(i,n),dif
            write(6,*)
          endif

          lnewsym = .TRUE.

        endif
        do j = 1,3
          if(mtrxnew(j,i,n) /= mtrx(j,i,n)) then
            if(ipr > 0) then
              write(6,*)
              write(6,'("   sym_same: mtrx /= mtrxnew",3i6,2f12.4)')     &
       &                    j,i,n,mtrx(j,i,n),mtrxnew(j,i,n)
              write(6,*)
            endif

            lnewsym = .TRUE.

          endif
        enddo
      enddo
    enddo

  endif

  return

end subroutine sym_same
