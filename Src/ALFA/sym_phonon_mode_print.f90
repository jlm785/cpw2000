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

!>  Given a q-vector and a crystal structure finds a basis
!>  for phonon modes according to the irreducible representations
!>  of the little group.  Based on Quantum Espresso.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         6 December 2022.
!>  \copyright    GNU Public License v2

subroutine sym_phonon_mode_print(qvec, vpol, nrep, irep, io,             &
       adot, ntype, natom, nameat, rat,                                  &
       ntrans_q, mtrx_q, tnp_q,                                          &
       mxdtyp, mxdatm, mxdnat)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdnat                          !<  array dimension of total number of atoms

  integer, intent(in)                ::  io                              !<  default output

  real(REAL64), intent(in)           ::  qvec(3)                         !<  wave-vector in lattice coordinates

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  ntrans_q                        !<  number of symmetry operations in the little group of q
  integer, intent(in)                ::  mtrx_q(3,3,48)                  !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the little group of q
  real(REAL64), intent(in)           ::  tnp_q(3,48)                     !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the little group of q

  complex(REAL64), intent(in)        ::  vpol(3,mxdatm,mxdtyp,mxdnat)    !<  proper modes in lattice coordinates.

  integer, intent(in)                ::  nrep                            !<  number of representations present
  integer, intent(in)                ::  irep(mxdnat)                    !<  representation of each mode

! local allocatable arrays

  real(REAL64), allocatable          ::  ratcar(:,:,:)                   !  cartesian coordinates
  complex(REAL64), allocatable       ::  vpolcar(:,:,:)                  !  cartesian coordinates

! local variables


  real(REAL64)          ::  avec(3,3), bvec(3,3)
  real(REAL64)          ::  qcart(3)
  integer               ::  jmode, code_group

! constants

  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer             ::  nt, n, i, j, k, nr, jr



  call adot_to_avec_sym(adot,avec,bvec)

  do i = 1,3
    qcart(i) = zero
    do j = 1,3
      qcart(i) = qcart(i) + bvec(i,j)*qvec(j)
    enddo
  enddo

  allocate(ratcar(3,mxdatm,mxdtyp))
  allocate(vpolcar(3,mxdatm,mxdtyp))

  do nt = 1,ntype
    do j = 1,natom(nt)
      do i = 1,3
        ratcar(i,j,nt) = ZERO
        do k = 1,3
          ratcar(i,j,nt) = ratcar(i,j,nt) + avec(i,k)*rat(k,j,nt)
        enddo
      enddo
    enddo
  enddo

! vector

  write(io,*)
  write(io,*)
  write(io,'("   Basis for phonon modes according to the irreducible representations")')
  write(io,*)

  write(io,'("   The coordinates of the q-vector are")')
  write(io,*)
  write(io,'(20x,"Lattice",40x,"Cartesian")')
  write(io,*)
  write(io,'(3x,3f12.4,12x,3f12.4)') (qvec(j),j=1,3), (qcart(j),j=1,3)
  write(io,*)

! symmetry

  write(io,*)
  write(io,*)
  write(io,'(4x,"rotation matrices and fractional",                      &
      & " translations in lattice coordinates")')
  write(io,'(4x,"of the little group of q-vector.")')
  write(io,*)
  do n = 1,ntrans_q
    write(io,'(i5,3(3x,3i3),4x,3f17.10,"    symmetry op.")')             &
       n, ((mtrx_q(i,j,n),j=1,3),i=1,3), (tnp_q(j,n)/(2*PI),j=1,3)
  enddo

  call sym_point_group_name(adot, 1, code_group, ntrans_q, mtrx_q)


! representations

  do nr = 1,nrep
    write(io,*)
    write(io,*)
    write(io,'("   representation",i3)') nr
    write(io,*)
    jmode = 0
    do jr = 1,mxdnat
      if(irep(jr) == nr) then

        jmode = jmode + 1
        write(io,*)
        write(io,'("   mode",i3,5x,"atom coord (Lat)",7x,                &
            &   "Displacement (Lat)",6x,"atom coord (Cart)",9x,          &
            &   "Displacement (Cart)")') jmode
        write(io,*)
        do nt = 1,ntype
          do j = 1,natom(nt)
            do i = 1,3
              vpolcar(i,j,nt) = C_ZERO
              do k = 1,3
                vpolcar(i,j,nt) = vpolcar(i,j,nt) + avec(i,k)*vpol(k,j,nt,jr)
              enddo
            enddo
          enddo
        enddo

        do nt = 1,ntype
        do j = 1,natom(nt)
          write(io,'(3x,a2,1x,i3,2x,3f7.3,3x,3f7.3,3x,3f7.3,6x,3f7.3)')  &
                nameat(nt), j,                                           &
                (rat(k,j,nt),k=1,3), (real(vpol(k,j,nt,jr)),k=1,3),      &
                (ratcar(k,j,nt),k=1,3), (real(vpolcar(k,j,nt)),k=1,3)
          write(io,'(35x,3f7.3,30x,3f7.3)') (aimag(vpol(k,j,nt,jr)),k=1,3),     &
                                            (aimag(vpolcar(k,j,nt)),k=1,3)
        enddo
        enddo

      endif
    enddo
  enddo

  deallocate(ratcar)
  deallocate(vpolcar)


  return

end subroutine sym_phonon_mode_print

