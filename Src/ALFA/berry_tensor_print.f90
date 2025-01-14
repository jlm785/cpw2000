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

!>  Auxiliary subroutine to print the topological tensors
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         20 April 2024, 13 January 2025.
!>  \copyright    GNU Public License v2

subroutine berry_tensor_print(levdegnl, adot, t_lat, csym, cident, nlev, &
       mxddeg)

! Extracted and modified from out_qgeom_print 20 April 2024. JLM
! dimensions t_lat. 23 November 2024. JLM
! improve printing for better greping, 13 January 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels

  integer, intent(in)                ::  levdegnl                        !<  degeneracy of level
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  t_lat(mxddeg,mxddeg,3,3)        !<  Tensor in lattice coordinates
  character(len=1), intent(in)       ::  csym                            !<  indicates wether the tensor is symetric (S) or antisymetric (A)
  character(len=10), intent(in)      ::  cident                          !<  identifier for easier greping
  integer, intent(in)                ::  nlev                            !<  level number (used only for greping)

! local allocatable arrays

  real(REAL64), allocatable          ::  t_car(:,:,:,:)                  !   Tensor in cartesian coordinates
  real(REAL64), allocatable          ::  trace_r(:,:)                    !   Trace on coordinates of tensor
  real(REAL64), allocatable          ::  ev(:)
  real(REAL64), allocatable          ::  vec(:,:)


! local variables

  real(REAL64)       ::  trace_n(3,3)
  real(REAL64)       ::  pvec(3)
  logical            ::  lcheck
  integer            ::  info
  real(REAL64)       ::  det, xsum
  real(REAL64)       ::  ev3(3), vec3(3,3)

  character(len=50)  ::  cformat


! constants

  real(REAL64), parameter     ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer    ::  j, n, k, m



  if(levdegnl < 1) then
    write(6,*)
    write(6,*) '   error in berry_tensor_print'
    write(6,*) '   degenearcy is not a positive number: ',levdegnl
    write(6,*)

    stop

  endif

  allocate(t_car(mxddeg,mxddeg,3,3))
  allocate(trace_r(mxddeg,mxddeg))
  allocate(ev(mxddeg))
  allocate(vec(mxddeg,mxddeg))

  call berry_tensor_lat2car(adot, t_lat, t_car, levdegnl,                &
         mxddeg)

  write(6,*)
  write(6,*) "   primitive (reciprocal) lattice coordinates"
  write(6,*)

  if(levdegnl < 10) then
    write(cformat,'("(",i1,"(3f12.3,4x))")') levdegnl
  else
    cformat = '(8(3f12.3,4x))'
  endif

  do n = 1,levdegnl
    do j = 1,3
      write(6,cformat) ((t_lat(n,m,j,k), k = 1,3), m = 1,levdegnl)
    enddo
    write(6,*)
  enddo

  write(6,*)
  write(6,*) "   Cartesian lattice coordinates"
  write(6,*)

  if(levdegnl < 10) then
    write(cformat,'("(",i1,"(3f12.3,4x),20x,a10,3i5)")') levdegnl
  else
    cformat = '(8(3f12.3,4x))'
  endif

  do n = 1,levdegnl
    do j = 1,3
      write(6,cformat) ((t_car(n,m,j,k), k = 1,3), m = 1,levdegnl),      &
           cident, n, j, nlev
    enddo
    write(6,*)
  enddo
  write(6,*)

  call  berry_tensor_trace(t_car, trace_n, trace_r, levdegnl, csym, lcheck, mxddeg)

  if(levdegnl == 1) then
    write(6,*)
    write(6,'("   Trace over coordinates",f12.3)') trace_r(1,1)
    write(6,*)
  else
    write(6,*)
    write(6,*) "   Trace over coordinates"
    write(6,*)
    do n = 1,levdegnl
      write(6,'(20(f12.3,4x))') (trace_r(n,m), m = 1,levdegnl)
    enddo
  endif
  write(6,*)
  if(levdegnl > 1) then
    xsum = ZERO
    do n = 1,levdegnl
      xsum= xsum + trace_r(n,n)
    enddo
    write(6,'("    Trace over trace of coordinates",f12.3)') xsum
    write(6,*)
  endif


  if(csym == 'S' .or. csym == 's') then

    if(lcheck .and. levdegnl > 1) then

      call diag_r8(levdegnl, trace_r, ev, vec, mxddeg, info)

      if(info /= 0) then
        write(6,*)
        write(6,*) "   Error in diagonalization in berry_tensor_print"
        write(6,*)
      else
        write(6,*)
        write(6,'("   Eigenvalues: ", 30f12.3)') (ev(n),n = 1,levdegnl)
        det = UM
        do n = 1,levdegnl
          det = det*ev(n)
        enddo
        write(6,'(3x,e14.3,"   Determinant ",a10,i5)') det, cident, nlev
        write(6,*)
      endif

    endif

  endif

  if(levdegnl > 1) then
    write(6,*)
    write(6,*) "   Trace over energy levels"
    write(6,*)
    do n = 1,3
      write(6,'(5x,3f12.3)') (trace_n(n,m), m = 1,3)
    enddo
    write(6,*)
  endif

  if(csym == 'A' .or. csym == 'a') then

    if(lcheck) then

      pvec(1) = trace_n(2,3)
      pvec(2) = trace_n(3,1)
      pvec(3) = trace_n(1,2)

      if(levdegnl == 1) then
        write(6,*)
        write(6,*) "   Associated pseudo-vector"
        write(6,*)
      else
        write(6,*)
        write(6,*) "   Pseudo-vector associated with the trace."
        write(6,*)
      endif
      write(6, '(5x,3f12.3,20x,a10,"vector",i5)') (pvec(k),k=1,3), cident, nlev
      xsum = pvec(1)*pvec(1) + pvec(2)*pvec(2) + pvec(3)*pvec(3)
      xsum = sqrt(xsum)
      write(6,*)
      write(6,'("    Module: ", f12.3)') xsum
      write(6,*)

    endif

  endif

  if(csym == 'S' .or. csym == 's') then

    if(lcheck) then
      xsum = ZERO
      do n = 1,3
        xsum= xsum + trace_n(n,n)
      enddo
      write(6,'("   Trace over trace of levels",f12.3)') xsum
      write(6,*)

      call diag_r8(3, trace_n, ev3, vec3, 3, info)

      if(info /= 0) then
        write(6,*)
        write(6,*) "   Error in diagonalization in berry_tensor_print 2"
        write(6,*)
      else
        write(6,*)
        write(6,'("   Eigenvalues: ", 30f12.3)') (ev3(n),n = 1,3)
        det = UM
        do n = 1,3
          det = det*ev3(n)
        enddo
        write(6,'(3x,e14.3,"   Determinant ",a10,i5)') det, cident, nlev
        write(6,*)
      endif

    endif

  endif


  deallocate(t_car)
  deallocate(trace_r)
  deallocate(ev)
  deallocate(vec)

  return

end subroutine berry_tensor_print
