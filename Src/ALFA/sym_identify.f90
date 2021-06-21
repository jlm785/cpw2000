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

!>  Identifies the space group of the crystal.
!>  imposes higher symmetry in case of slight symmetry breaking

!>  It has a parameter, tol, that decides if two real numbers
!>  are within a noise tolerance.
!>  If the structure is a perfect superlattice, it may give (but shouldn't) buggy results!
!>  It imposes the identified symmetry on the crystal structure (see the use of tol).

subroutine sym_identify(isym, ipr, tol,                                  &
  ntrans, mtrx, tnp,                                                     &
  adot, ntype, natom, rat,                                               &
  mxdtyp, mxdatm)

! Written november 2003. jlm
! Modified april 2004 (superlattice bug). JLM
! Modified, f90, June 8 2014. JLM
! Modified December 8, 2016. sym_rat. JLM
! Modified, documentation, December 2019. JLM
! Modified, sym_test, 29 December 2020. JLM
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

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i

! input and output

  real(REAL64), intent(inout)        ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(inout)        ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

! output

  integer, intent(out)               ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(out)               ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(out)          ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! allocatable arrays

  integer, allocatable                ::  nattmp(:)
  real(REAL64), allocatable           ::  rattmp(:,:,:)
  integer, allocatable                ::  jsup(:)                        !  indicates the nsup atoms not related by superlattice translations (atoms with minimal number)

! local variables

  real(REAL64)       ::  avec(3,3), bvec(3,3)
  real(REAL64)       ::  vcell, adotnig(3,3), adotsym(3,3)
  integer            ::  iprmet
  integer            ::  ibravais, mtotal(3,3)
  integer            ::  nrot, isup, iprint
  real(REAL64)       ::  rot(3,3,48), temp(3,3)
  integer            ::  irotdir(3,3,48), irotrec(3,3,48)
  real(REAL64)       ::  aconv(3,3), avecnig(3,3)
  real(REAL64)       ::  dif, difmax

  integer            ::  istatus                                         !  istatus = 0, successful; 1 not closed; 2 no inverse; 3 inconsistent with atomic positions

  integer            ::  nsup                                            !  number of atoms not related by superlattice translations (atoms with minimal number)
 
! parameters

  real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counters

  integer i, j, n, nt


  if(isym == 0) then

!   does not use symmetry

    ntrans = 1
    tnp(1,1) = ZERO
    tnp(2,1) = ZERO
    tnp(3,1) = ZERO
    mtrx(1,1,1) = 1
    mtrx(1,2,1) = 0
    mtrx(1,3,1) = 0
    mtrx(2,1,1) = 0
    mtrx(2,2,1) = 1
    mtrx(2,3,1) = 0
    mtrx(3,1,1) = 0
    mtrx(3,2,1) = 0
    mtrx(3,3,1) = 1

  else

!   finds first the bravais lattice and "conventional" lattice vectors

    iprmet = 0
    if(ipr > 1) iprmet = 1

    call metric_ident(adot, adotnig, adotsym, ibravais, mtotal,          &
                      avec, aconv, avecnig, tol, iprmet)


!   replaces by symmetrized adot

    do i=1,3
    do j=1,3
      adot(j,i) = adotsym(j,i)
    enddo
    enddo

!   compute the reciprocal vectors and cell volume

    bvec(1,1) = avec(2,2)*avec(3,3) - avec(3,2)*avec(2,3)
    bvec(2,1) = avec(3,2)*avec(1,3) - avec(1,2)*avec(3,3)
    bvec(3,1) = avec(1,2)*avec(2,3) - avec(2,2)*avec(1,3)
    bvec(1,2) = avec(2,3)*avec(3,1) - avec(3,3)*avec(2,1)
    bvec(2,2) = avec(3,3)*avec(1,1) - avec(1,3)*avec(3,1)
    bvec(3,2) = avec(1,3)*avec(2,1) - avec(2,3)*avec(1,1)
    bvec(1,3) = avec(2,1)*avec(3,2) - avec(3,1)*avec(2,2)
    bvec(2,3) = avec(3,1)*avec(1,2) - avec(1,1)*avec(3,2)
    bvec(3,3) = avec(1,1)*avec(2,2) - avec(2,1)*avec(1,2)

!   cell volume

    vcell = bvec(1,1)*avec(1,1) + bvec(2,1)*avec(2,1) + bvec(3,1)*avec(3,1)

!   reciprocal lattice vectors without the 2*pi !!!!!!!!!!!!!!!!!!!!!!!!!

    do j=1,3
      bvec(1,j) = bvec(1,j)/vcell
      bvec(2,j) = bvec(2,j)/vcell
      bvec(3,j) = bvec(3,j)/vcell
    enddo

!   generates the rotation matrices in cartesian coordinates
!   of the point group of the lattice

    call sym_rot(ibravais, nrot, rot)

!   transforms into lattice (direct and reciprocal coordinates)

    do n=1,nrot
      do i=1,3
      do j=1,3
        temp(i,j) = rot(i,1,n)*avec(1,j) +                               &
                    rot(i,2,n)*avec(2,j) +                               &
                    rot(i,3,n)*avec(3,j)
      enddo
      enddo

      do i=1,3
      do j=1,3
        irotdir(i,j,n) = nint(bvec(1,i)*temp(1,j) +                      &
                              bvec(2,i)*temp(2,j) +                      &
                              bvec(3,i)*temp(3,j))
      enddo
      enddo

      do i=1,3
      do j=1,3
        temp(i,j) = avec(1,i)*rot(1,j,n) +                               &
                    avec(2,i)*rot(2,j,n) +                               &
                    avec(3,i)*rot(3,j,n)
      enddo
      enddo

      do i=1,3
      do j=1,3
        irotrec(i,j,n) = nint(temp(i,1)*bvec(1,j) +                      &
                              temp(i,2)*bvec(2,j) +                      &
                              temp(i,3)*bvec(3,j))
      enddo
      enddo

    enddo

!   checks if there is a fractional translaton
!   aplly unelegant fix for this case

    iprint = 1
    allocate(jsup(mxdatm))

    call sym_sup_lat(iprint, tol, isup, jsup, nsup,                      &
       adot, ntype, natom, rat,                                          &
       mxdtyp, mxdatm)

    if(isup > 0) then

      allocate(nattmp(mxdtyp+nsup))
      allocate(rattmp(3,mxdatm,mxdtyp+nsup))

      do nt = 1,ntype
        nattmp(nt) = natom(nt)
        do j = 1,natom(nt)
          rattmp(1,j,nt) = rat(1,j,nt)
          rattmp(2,j,nt) = rat(2,j,nt)
          rattmp(3,j,nt) = rat(3,j,nt)
        enddo
      enddo
      nattmp(ntype+1) = nsup
      do n = 1,nsup
        rattmp(1,n,ntype+1) = rat(1,jsup(n),isup)
        rattmp(2,n,ntype+1) = rat(2,jsup(n),isup)
        rattmp(3,n,ntype+1) = rat(3,jsup(n),isup)
      enddo

      write(6,'("  WARNING         ",i5," atoms of type ",i3,            &
        & "are considered different to fix superlattice bug")')          &
                 nsup,isup

      call sym_space_group(ntrans, mtrx, tnp, tol,                       &
         nrot, irotdir, irotrec,                                         &
         ntype+1, nattmp, rattmp,                                        &
         mxdtyp+1, mxdatm)

      deallocate(nattmp)
      deallocate(rattmp)

    else

      call sym_space_group(ntrans, mtrx, tnp, tol,                       &
         nrot, irotdir, irotrec,                                         &
         ntype, natom, rat,                                              &
         mxdtyp, mxdatm)

    endif

    deallocate(jsup)

  endif

  if(ipr > 0) then
    write(6,*)
    write(6,*)
    write(6,'(4x,"rotation matrices and fractional",                     &
        & " translations in lattice coordinates")')
    write(6,*)
    do n=1,ntrans
      write(6,'(i5,3(3x,3i3),4x,3f17.10,"    symmetry op.")')            &
         n,((mtrx(i,j,n),j=1,3),i=1,3),(tnp(j,n)/(2*PI),j=1,3)
    enddo
  endif

! tests it is a group

  iprmet = ipr-1
  call sym_test(iprmet, tol, istatus,                                    &
     ntrans, mtrx, tnp,                                                  &
     ntype, natom, rat, adot,                                            &
     mxdtyp, mxdatm)

  if(istatus /= 0) then
    write(6,*)
    write(6,*) '   STOPPED in sym_identify. Failed symmetry test.'

    stop

  endif

! imposes symmetry

  allocate(rattmp(3,mxdatm,mxdtyp))

  do nt = 1,ntype
    do j = 1,natom(nt)
      rattmp(1,j,nt) = rat(1,j,nt)
      rattmp(2,j,nt) = rat(2,j,nt)
      rattmp(3,j,nt) = rat(3,j,nt)
    enddo
  enddo

  call sym_rat(ntype, natom, rattmp, ntrans, mtrx, tnp,                  &
     mxdtyp, mxdatm)

  difmax = ZERO

  do nt = 1,ntype
    do j = 1,natom(nt)
      do i=1,3
        dif = abs(rattmp(i,j,nt)-rat(i,j,nt))
        difmax = max(dif,difmax)
      enddo
    enddo
  enddo

  if(ipr > 1) then
    write(6,*)
    write(6,*)  '  The atomic positions were symmetrized'
    write(6,'("   The maximum value of the correction was",f18.8)') difmax
    write(6,*)
  endif

  if(difmax > 0.001) then
    write(6,*)
    write(6,*)  '  Stopped in sym_identify:   difmax = ',difmax
    write(6,*)

    stop

  endif

  do nt = 1,ntype
    do j = 1,natom(nt)
      rat(1,j,nt) = rattmp(1,j,nt)
      rat(2,j,nt) = rattmp(2,j,nt)
      rat(3,j,nt) = rattmp(3,j,nt)
    enddo
  enddo

  deallocate(rattmp)

  return

end subroutine sym_identify

