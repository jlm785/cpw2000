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

!>     This subroutine calculates the bands on an uniform grid
!>     for later processing using the atomic orbital interpolation

      subroutine ao_interpolation_out_dos(noiData, ztot, adot, ntrans, mtrx)                                                

!      Adapted to use ao interpolatiion package July, 2014. clr
!      modified 25 November 2015. JLM
!      Modified documentation May 2020. JLM
!      copyright  Jose Luis Martins, Carlos Loia Reis/INESC-MN

!      version 4.98

       use NonOrthoInterp

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       type(noiData_t)                     ::  noiData                    !<  see NonOrthoInterp

       real(REAL64), intent(in)            ::  ztot                       !<  total charge density (electrons/cell)

       real(REAL64), intent(in)            ::  adot(3,3)                  !<  metric in direct space
       integer, intent(in)                 ::  ntrans                     !<  number of symmetry operations in the factor group
       integer, intent(in)                 ::  mtrx(3,3,48)               !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

!      local allocatable arrays

       integer, allocatable                ::  kmap2(:,:,:,:)             !  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
       integer, allocatable                ::  nband2(:)                  !  number of bands for each k-points
       integer, allocatable                ::  indk2(:,:)                 !  index of the six k-points neighbouring k-point i
       real(REAL64),allocatable            ::  rk2(:,:)                   !  component in lattice coordinates of the k-point in the mesh
       real(REAL64),allocatable            ::  w2(:)                      !  weight in the integration of k-point

!      local variables

       integer           ::  mxdbnd         !  array dimension for the number of bands
       integer           ::  mxdpnt         !  dimensions for dos k-points

       integer           ::  neig           !  number of eigenvectors required (maybe modified on output)
       real(REAL64)      ::  rkpt(3)        !  j-th component in lattice coordinates of the k-point

       integer           ::  irk
       integer           ::  ipr,nrk2
       logical           ::  lfile
       integer           ::  nbandi2,nx2,ny2,nz2
       real(REAL64)      ::  sx2,sy2,sz2
       real(REAL64)      ::  vcell, bdot(3,3)

       integer           ::  io, io_so     !  tape numbers 
       character(len=128)  ::  syscom

!      counters
       integer    ::  i, j, k, n, jmax

       real(REAL64), allocatable :: ev_interp(:)


!      calculates local potential in fft mesh
       
       ipr = 2

       lfile = .false.
       open(unit=11,file='DOS_MESH.DAT',status='old',err=10, form = 'formatted')
       lfile = .true.
 10    continue

       if(lfile) then
         read(11,*) nbandi2,nx2,ny2,nz2,sx2,sy2,sz2

         close(unit=11)
       
         mxdbnd = nbandi2

       else
       
         mxdbnd = nint(ztot) + 4

         nbandi2 = mxdbnd
         nx2 = 8
         ny2 = 8
         nz2 = 8
         sx2 = 0.0
         sy2 = 0.0
         sz2 = 0.0
       endif

       mxdpnt = nx2*ny2*nz2

       allocate(kmap2(3,nx2,ny2,nz2))
       allocate(nband2(mxdpnt))
       allocate(indk2(6,mxdpnt))
       allocate(rk2(3,mxdpnt))
       allocate(w2(mxdpnt))
   
       call int_pnt(nbandi2,nx2,ny2,nz2,sx2,sy2,sz2,ipr,                 &
     & adot,                                                             &
     & ntrans,mtrx,                                                      &
     & nrk2,rk2,w2,nband2,indk2,kmap2,                                   &
     & mxdpnt,mxdbnd)

       neig = nbandi2

!      writes in a file for later processing

       io = 12
       open(unit=io,file='PW_DOS.DAT',form='formatted')

       call adot_to_bdot(adot,vcell,bdot)

       write(io,'(2x,i6,2x,f20.10,2x,f20.10,2x,i1)') nrk2,ztot,vcell,2
       do irk=1,nrk2,20
         jmax = irk+19
         if(jmax .gt. nrk2) jmax = nrk2
         write(io,'(2x,20i6)') (nband2(j),j=irk,jmax)
       enddo
       do irk=1,nrk2,10
         jmax = irk+9
         if(jmax > nrk2) jmax = nrk2
         write(io,'(2x,10f14.8)') (w2(j),j=irk,jmax)
       enddo
       do irk=1,nrk2
         write(io,'(2x,6i6)') (indk2(j,irk),j=1,6)
       enddo

!      repeats for spin-orbit

       io_so = 13
       open(unit=io_so,file='PW_DOS_SO.DAT',form='formatted')

       call adot_to_bdot(adot,vcell,bdot)

       write(io_so,'(2x,i6,2x,f20.10,2x,f20.10,2x,i1)') nrk2, ztot,      &
     &              vcell,1
       do irk=1,nrk2,20
         jmax = irk+19
         if(jmax > nrk2) jmax = nrk2
         write(io_so,'(2x,20i6)') (2*nband2(j),j=irk,jmax)
       enddo
       do irk=1,nrk2,10
         jmax = irk+9
         if(jmax > nrk2) jmax = nrk2
         write(io_so,'(2x,10f14.8)') (w2(j),j=irk,jmax)
       enddo
       do irk=1,nrk2
         write(io_so,'(2x,6i6)') (indk2(j,irk),j=1,6)
       enddo
              
       allocate(ev_interp(noiData%nband))

!      loop over k-points

       do irk=1,nrk2

         rkpt(1) = rk2(1,irk)
         rkpt(2) = rk2(2,irk)
         rkpt(3) = rk2(3,irk)
         
         call NonOrthoInterpRun(noiData,rkpt,ev_interp)

         if(noiData%lso==1) then
          do i=1,2*neig,8
            jmax = i+7
            if (jmax > 2*neig) jmax = 2*neig
            write(io_so,'(15x,8f12.6)') (ev_interp(j),j=i,jmax)
          enddo
         else
          do i=1,neig,8
            jmax = i+7
            if (jmax > neig) jmax = neig
            write(io,'(15x,8f12.6)') (ev_interp(j),j=i,jmax)
          enddo
         endif
         
         call progress(irk,nrk2)
         call flush(6)
         
       enddo

       write(io,'(3i6)') nx2,ny2,nz2
       do k=1,nz2
       do j=1,ny2
       do i=1,nx2
         write(io,'(3i8)') (kmap2(n,i,j,k),n=1,3)
       enddo
       enddo
       enddo
       
       close(unit=io)

       write(io_so,'(3i6)') nx2,ny2,nz2
       do k=1,nz2
       do j=1,ny2
       do i=1,nx2
         write(io_so,'(3i8)') (kmap2(n,i,j,k),n=1,3)
       enddo
       enddo
       enddo
       
       close(unit=io_so)
       
       if(noiData%lso==1) then
!          call system("rm PW_DOS.DAT");
         syscom = "rm PW_DOS.DAT"
       else
!          call system("rm PW_DOS_SO.DAT");      
         syscom = "rm PW_DOS_SO.DAT"
       endif
       call execute_command_line(trim(syscom))

!      Cleanup
       
       deallocate(kmap2)
       deallocate(nband2)
       deallocate(indk2)
       deallocate(rk2)
       deallocate(w2)
       
       deallocate(ev_interp)
       
       return
       end subroutine ao_interpolation_out_dos
