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


!>  This subroutine calculates the bands on an uniform grid
!>  and the oscillator strengths for later processing.
!>  It uses fourier interpolations to compute matrix elements
!>  Adapted from out_dos
!>
!>  \author       Carlos Loia Reis
!>  \version      5.11
!>  \date         November 2018, 7 October 2024.
!>  \copyright    GNU Public License v2

subroutine ao_interpolation_out_ie(noiData,ztot,adot,ntrans,mtrx)


! Written November 2018. CLR
! Modified, indentation, documentation, 7 October 2024. JLM


  use NonOrthoInterp

  implicit none
  type(noiData_t)                    :: noiData

  integer, parameter                 :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group


! local allocatable arrays

  integer, allocatable               ::  kmap2(:,:,:,:)                  !  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
  integer, allocatable               ::  nband2(:)                       !  number of bands for each k-points
  integer, allocatable               ::  indk2(:,:)                      !  index of the six k-points neighbouring k-point i
  real(REAL64),allocatable           ::  rk2(:,:)                        !  component in lattice coordinates of the k-point in the mesh
  real(REAL64),allocatable           ::  w2(:)                           !  weight in the integration of k-point


  real(REAL64), allocatable          ::  efmei(:)                        ! E_f - E_i
  real(REAL64), allocatable          ::  efmei_allk(:,:)                 ! E_f - E_i for all k
  integer, allocatable               ::  indx(:)                         ! index for energy difference sorting

! oscillator strength stuff

  complex(REAL64), allocatable       ::  dh0drk(:,:,:)                   !  <i|d H / d k|f>
  complex(REAL64), allocatable       ::  ematif(:)                       ! |<i|d H / d k|f>|**2
  complex(REAL64), allocatable       ::  ematif_allk(:,:)                ! |<i|d H / d k|f>|**2

! interpolation stuff
  integer nband

  real(REAL64), allocatable          ::  ev(:), ev_wrk(:)

  complex(REAL64) , allocatable      ::  Hk(:,:), Sk(:,:)
  complex(REAL64), allocatable       ::  S12(:,:), S12_inv(:,:)
  complex(REAL64) , allocatable      ::  vec(:,:), U(:,:)
  complex(REAL64) , allocatable      ::  wrk(:,:)

  integer                            ::  mxdpnt

  real(REAL64)                       ::  rkpt(3)                         !  j-th component in lattice coordinates of the k-point

  integer                            ::  irk
  integer                            ::  ipr,nrk2
  logical                            ::  lfile
  integer                            ::  nbandi2,nx2,ny2,nz2
  real(REAL64)                       ::  sx2,sy2,sz2
  real(REAL64)                       ::  vcell, bdot(3,3)

  integer                            ::  ncond,nval

  integer                            ::  ioerr                           !  error in opening/reading/writing files

  integer  ::  info

  integer                            ::  io, io_so                       !  tape numbers

! parameters

  real(REAL64)   , parameter         :: HARTREE = 27.21138386_REAL64
  real(REAL64)   , parameter         :: PI = 3.14159265358979323846_REAL64

  real(REAL64)   , parameter         :: ZERO = 0.0_REAL64
  real(REAL64)   , parameter         :: UM = 1.0_REAL64
  complex(REAL64), parameter         :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter         :: C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer                            ::  i, j, k, n, jmax, ij, m, idir


  ipr = 2

  lfile = .FALSE.
  open(unit=11,file='DOS_MESH.DAT',status='old',iostat=ioerr, form = 'formatted')

  if(ioerr == 0) lfile = .TRUE.

  if(lfile) then
    read(11,*) nbandi2, nx2,ny2,nz2, sx2,sy2,sz2

    close(unit=11)

!    mxdbnd = nbandi2

  else
    nbandi2 = nint(ztot) + 4
    nx2 = 8
    ny2 = 8
    nz2 = 8
    sx2 = 0.0
    sy2 = 0.0
    sz2 = 0.0
  endif


  nband = noiData %nband

  mxdpnt = nx2*ny2*nz2

  allocate(kmap2(3,nx2,ny2,nz2))
  allocate(nband2(mxdpnt))
  allocate(indk2(6,mxdpnt))
  allocate(rk2(3,mxdpnt))
  allocate(w2(mxdpnt))


  call int_pnt(nbandi2, nx2, ny2, nz2, sx2, sy2, sz2, ipr,               &
      adot,                                                              &
      ntrans, mtrx,                                                      &
      nrk2, rk2, w2, nband2, indk2, kmap2,                               &
      mxdpnt, nband)


  if(noiData %lso ==0) then
   nval = nint(0.5*ztot)
  else
   nval = nint(ztot)
  endif


  ncond = nband - nval

  if(noiData %lso == 0) then
   open(unit=21, file="opt_ie.dat", form="unformatted")
  else
   open(unit=21, file="opt_ie_so.dat", form="unformatted")
  endif

   write(21) noiData %lso
   write(21) nx2,ny2,nz2
   write(21) ntrans
   write(21) mtrx
   write(21) nrk2
   write(21) rk2
   write(21) w2
   write(21) indk2
   write(21) kmap2

   write(21) ztot
   write(21) nband
   write(21) nval
   write(21) ncond
   write(21) adot



!dbg  write(6,*) 'NVAL,NCOND',NVAL,NCOND

! writes in a file for later processing

  call adot_to_bdot(adot,vcell,bdot)

  if(noiData %lso ==0) then

    io = 12
    open(unit=io,file='PW_IE_AA.DAT',form='formatted')
    write(io,'(2x,i6,2x,f20.10,2x,f20.10,2x,i1)') nrk2,ztot,vcell,2
    do irk=1,nrk2,20
      jmax = irk+19
      if(jmax .gt. nrk2) jmax = nrk2
      write(io,'(2x,20i6)') (nval*ncond,j=irk,jmax)
    enddo
    do irk=1,nrk2,10
      jmax = irk+9
      if(jmax > nrk2) jmax = nrk2
      write(io,'(2x,10f14.8)') (w2(j),j=irk,jmax)
    enddo
    do irk=1,nrk2
      write(io,'(2x,6i6)') (indk2(j,irk),j=1,6)
    enddo

  else

    io_so = 13
    open(unit=io_so,file='PW_IE_SO_AA.DAT',form='formatted')
    write(io_so,'(2x,i6,2x,f20.10,2x,f20.10,2x,i1)') nrk2, ztot, vcell,1
    do irk=1,nrk2,20
      jmax = irk+19
      if(jmax > nrk2) jmax = nrk2
      write(io_so,'(2x,20i6)') (nval*ncond,j=irk,jmax)
    enddo
    do irk=1,nrk2,10
      jmax = irk+9
      if(jmax > nrk2) jmax = nrk2
      write(io_so,'(2x,10f14.8)') (w2(j),j=irk,jmax)
    enddo
    do irk=1,nrk2
      write(io_so,'(2x,6i6)') (indk2(j,irk),j=1,6)
    enddo

  endif


  allocate(Hk(nband,nband), Sk(nband,nband), vec(nband,nband), U(nband,nband))
  allocate(wrk(nband,nband), S12(nband,nband), S12_inv(nband,nband))

  allocate(ev(nband))
  allocate(ev_wrk(nband))

  allocate(efmei(nval*ncond))
  allocate(efmei_allk(nval*ncond,nrk2))

  allocate(ematif(nval*ncond))
  allocate(ematif_allk(nval*ncond,nrk2))

  allocate(indx(4*nval*ncond))

! loop over k-points

  allocate( dh0drk(nband,nband,3)      )        ! d <Psi|H|Psi> d k

  do irk=1,nrk2

    rkpt(1) = rk2(1,irk)
    rkpt(2) = rk2(2,irk)
    rkpt(3) = rk2(3,irk)

!   Interpolate H and S to the current k-point

    call fi_hamiltonian_get_hk(noiData%fiData,rkpt,Hk, nband)
    call fi_hamiltonian_get_sk(noiData%fiData,rkpt,Sk, nband)

!   Compute S^(-1/2) and S^(1/2) matrices

    call GetS12(Sk,S12,S12_inv,wrk,ev_wrk,nband)

! Solve the eigenvalue problem and get U matrix

!    call ZMul('N','N',Hk,S12,wrk,nband)
    call zgemm('N','N', nband,nband,nband, C_UM,Hk,nband, S12,nband,   &
               C_ZERO,wrk,nband)
!    call ZMul('N','N',S12,wrk,Hk,nband)
    call zgemm('N','N', nband,nband,nband, C_UM,S12,nband, wrk,nband,  &
               C_ZERO,Hk,nband)


    call diag_c16(nband,Hk,ev,U,nband,info)


    if(info /= 0) then

      write(6,*)
      write(6,*) '    STOPPED in ao_interpolation_out_ie'
      write(6,*) '    Could not diagonalize hk, info = ',info
      write(6,*)

      stop

    endif


!   Get the eigenvectors vec
!    call ZMul('N','N',S12,U,vec,nband)
    call zgemm('N','N', nband,nband,nband, C_UM,S12,nband, U,nband,    &
               C_ZERO,vec,nband)

    ! Interpolate dh0drk to the current k-point
    call fi_hamiltonian_get_Op_k(noiData%fiData,rkpt,noiData%fiData%dh0drk_GridR,dh0drk, nband)

    ! bring dhdrk from the AO representation to H representation
    do idir=1, 3
!     call ZMul('C','N',vec,dh0drk(:,:,idir),wrk,nband)
     call zgemm('C','N', nband,nband,nband, C_UM,vec,nband, dh0drk(:,:,idir),nband, &
                C_ZERO,wrk,nband)
!     call ZMul('N','N',wrk,vec,dh0drk(:,:,idir),nband)
     call zgemm('N','N', nband,nband,nband, C_UM,wrk,nband, vec,nband,              &
                C_ZERO,dh0drk(:,:,idir),nband)
    enddo

    write(*,'("dh0drk",i5, 3f8.5, i5)') irk,rkpt, nband

    write(21) ev
    write(21) dh0drk


    do i = 1,nval
    do j = 1,ncond
      ij = (i-1)*ncond + j
      efmei(ij)= ev(nval+j) - ev(i)
      ematif(ij) = C_ZERO

!---- !!!!
      do n = 1,3
      do m = 1,3
      ematif(ij) = ematif(ij) + dh0drk(nval+j,i,n)*     &
                           adot(n,m)*dh0drk(i,nval+j,m)
      enddo
      enddo
      ematif(ij) = ematif(ij) / (4*PI*PI)
!---- !!!!


    enddo
    enddo

    call sort(nval*ncond,efmei,indx)

    do ij = 1,nval*ncond
      efmei_allk(ij,irk)  = efmei(indx(ij))
      ematif_allk(ij,irk) = ematif(indx(ij))
    enddo


  enddo

  do irk = 1,nrk2
    do i = 1,nval*ncond,8
      jmax = i+7
      if (jmax > nval*ncond) jmax = nval*ncond

      if(noiData %lso ==0) then
         write(io,'(15x,8f12.6)') (efmei_allk(j,irk),j=i,jmax)
         write(io,'(15x,8f12.6)') (real(ematif_allk(j,irk)),j=i,jmax)
      else
         write(io_so,'(15x,8f12.6)') (efmei_allk(j,irk),j=i,jmax)
         write(io_so,'(15x,8f12.6)') (real(ematif_allk(j,irk)),j=i,jmax)
      endif

    enddo
  enddo

  if(noiData %lso ==0) then
     write(io,'(3i6)') nx2,ny2,nz2
     do k=1,nz2
     do j=1,ny2
     do i=1,nx2
       write(io,'(3i8)') (kmap2(n,i,j,k),n=1,3)
     enddo
     enddo
     enddo
     close(unit=io)
  else
     write(io_so,'(3i6)') nx2,ny2,nz2
     do k=1,nz2
     do j=1,ny2
     do i=1,nx2
       write(io_so,'(3i8)') (kmap2(n,i,j,k),n=1,3)
     enddo
     enddo
     enddo
     close(unit=io_so)
  endif


  close(unit=21)

  deallocate(kmap2)
  deallocate(nband2)
  deallocate(indk2)
  deallocate(rk2)
  deallocate(w2)

  deallocate(Hk, Sk, vec, U)
  deallocate(wrk, S12, S12_inv)
  deallocate(efmei,efmei_allk,ematif,ematif_allk)

  deallocate(ev)
  deallocate(ev_wrk)
  deallocate(dh0drk)


  return

end subroutine ao_interpolation_out_ie
