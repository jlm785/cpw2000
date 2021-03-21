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

!>     calculates the number of integration k-points for Monkhorst-Pack mesh

       subroutine size_mxdnrk(nx,ny,nz,sx,sy,sz,adot,ntrans,mtrx,        &
     & mxdnrk)

!      written May 5, 2014. JLM
!      adapted from intpnt
!      modified to be faster O(N^2) -> O(NlogN). 31 July 2014. JLM
!      Modified, documentation, January 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.94

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  nx, ny, nz                 !<  size of the integration mesh in k-space (nx*ny*nz)
       real(REAL64), intent(in)           ::  sx, sy, sz                 !<  offset of the integration mesh (usually 0.5)

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in real space
       integer, intent(in)                ::  ntrans                     !<  number of symmetry operations in the factor group
       integer, intent(in)                ::  mtrx(3,3,48)               !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

!      output

       integer, intent(out)               ::  mxdnrk                     !<  array dimension for k-points
       
!      local allocatable arrays

       real(REAL64), allocatable          ::  xk(:,:)
       real(REAL64), allocatable          ::  rk(:,:)
       integer, allocatable               ::  kmap(:)
       real(REAL64), allocatable          ::  dkg(:)
       integer, allocatable               ::  indx(:), indxinv(:)
      
!      local variables

       real(REAL64)      ::  dx, dy, dz
       real(REAL64)      ::  rktran(3)
       integer           ::  nrk, kdel
       real(REAL64)      ::  diff
       logical           ::  lfound, lfoundinv
       real(REAL64)      ::  bdot(3,3), vcell
       real(REAL64)      ::  xt(3),xtmod

       integer           ::  mxdpnt

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter :: EPS = 0.000001_REAL64

!      counters

       integer   ::   i, j, k, l, n
       integer   ::  i1, i2, i3, im1, im2, im3
       integer   ::  kmax, kmin, kmc


       call adot_to_bdot(adot,vcell,bdot)

       mxdpnt = nx*ny*nz

       allocate(xk(3,mxdpnt))
       allocate(kmap(mxdpnt))

       allocate(rk(3,mxdpnt))
       allocate(dkg(mxdpnt))

!      set up uniform array of k points
!      kmap is used to mark reducible k points and also to
!      map reducible to irreducible k points

       dx = UM / nx
       dy = UM / ny
       dz = UM / nz
       n = 0
       do i = 1,nx
       do j = 1,ny
       do k = 1,nz
         n = n+1
         xk(1,n) = (i-1 + sx) * dx
         xk(2,n) = (j-1 + sy) * dy
         xk(3,n) = (k-1 + sz) * dz
         kmap(n) = n
       enddo
       enddo
       enddo

!      finds distance to closest G-point

       im1 = nint((abs(adot(1,2))+abs(adot(1,3)))/adot(1,1)) + 1
       im2 = nint((abs(adot(2,1))+abs(adot(2,3)))/adot(2,2)) + 1
       im3 = nint((abs(adot(3,1))+abs(adot(3,2)))/adot(3,3)) + 1
       do n = 1,mxdpnt
         dkg(n) = ZERO
         do i=1,3
         do j=1,3
           dkg(n) = dkg(n) + xk(i,n)*bdot(i,j)*xk(j,n)
         enddo
         enddo
         do i1 = -im1,im1
         do i2 = -im2,im2
         do i3 = -im3,im3
           xt(1) = xk(1,n) - i1*UM
           xt(2) = xk(2,n) - i2*UM
           xt(3) = xk(3,n) - i3*UM
           xtmod = ZERO
           do i=1,3
           do j=1,3
             xtmod = xtmod + xt(i)*bdot(i,j)*xt(j)
           enddo
           enddo
           if(xtmod < dkg(n)) dkg(n) = xtmod
         enddo
         enddo
         enddo
       enddo
       
!      creates a list of increasing lengths

       allocate(indx(mxdpnt))
       allocate(indxinv(mxdpnt))

       call sort(mxdpnt,dkg,indx)
       
       do i = 1,mxdpnt
         indxinv(indx(i)) = i
       enddo
       

!      reduce to irreducible zone

       nrk = 0
       n = mxdpnt

       do i = 1,n
         if (kmap(i) == i) then

!          new irreducible point
!          mark with negative kmap

           nrk = nrk+1
           rk(1,nrk) = xk(1,i)
           rk(2,nrk) = xk(2,i)
           rk(3,nrk) = xk(3,i)
           kmap(i) = -nrk
           if (i /= n) then

!            first identifies those which have the same shortest distance to a G-vector

             kmin = indxinv(i)
             do k = indxinv(i),1,-1

               if(abs(dkg(i) - dkg(indx(k))) > EPS) exit
               
               kmin = k
             enddo

             kmax = indxinv(i)
             do k = indxinv(i),n

               if(abs(dkg(i) - dkg(indx(k))) > EPS) exit
               
               kmax = k
             enddo

!            operate on irreducible rk with the symmetry operations

             do j=1,ntrans
               rktran(1) = mtrx(1,1,j)*rk(1,nrk) + mtrx(1,2,j)*rk(2,nrk) &
     &                   + mtrx(1,3,j)*rk(3,nrk)
               rktran(2) = mtrx(2,1,j)*rk(1,nrk) + mtrx(2,2,j)*rk(2,nrk) &
     &                   + mtrx(2,3,j)*rk(3,nrk)
               rktran(3) = mtrx(3,1,j)*rk(1,nrk) + mtrx(3,2,j)*rk(2,nrk) &
     &                   + mtrx(3,3,j)*rk(3,nrk)

!              translate to interval 0-1.

               do k=1,3
                 kdel = int(rktran(k))
                 if (rktran(k) < 0) kdel = kdel-1
                 rktran(k) = rktran(k) - kdel
               enddo

!              remove (mark) k points related to irreducible rk by symmetry

               do kmc = kmin,kmax
               
                 k = indx(kmc)

                 if (kmap(k) == k) then

!                  both the transformed rk ...

                   lfound = .TRUE.
                   do l=1,3
                     diff = rktran(l) - xk(l,k)
                     diff = abs(diff - nint(diff))
                     if (diff > eps) then
                       lfound = .FALSE.

                       exit

                     endif
                   enddo
                   
                   if(lfound) then
                     kmap(k) = nrk
                   else

!                    ... and its inverse (time reversal)

                     lfoundinv = .TRUE.
                     do l=1,3
                       diff = rktran(l) + xk(l,k)
                       diff = abs(diff - nint(diff))
                       if (diff > eps) then
                         lfoundinv = .FALSE.

                         exit

                       endif
                     enddo
                     if(lfoundinv) then
                       kmap(k) = nrk
                     endif
                   endif

                 endif                                                   !  kmap(k) == k

               enddo                                                     !  do k = i+1,n

             enddo                                                       !  do j = 1,ntrans

           endif
         endif
       enddo                                                             !  do i = 1,n

       mxdnrk = nrk

       deallocate(xk)
       deallocate(kmap)
       deallocate(rk)
       deallocate(dkg)
       deallocate(indx)
       deallocate(indxinv)

       return
       end subroutine size_mxdnrk
