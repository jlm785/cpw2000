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

module IrredBZM
      integer, parameter,private          :: REAL64 = selected_real_kind(12)  

  type IrredBZ_t
       integer, allocatable  ::  kmap(:,:,:,:)         !  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
       integer, allocatable  ::  wrk_nband(:)              !  number of bands for each k-points
       integer, allocatable  ::  indk(:,:)             !  index of the six k-points neighbouring k-point i
       real(REAL64),allocatable   ::  rk(:,:)          !  component in lattice coordinates of the k-point in the mesh
       real(REAL64),allocatable   ::  w(:)             !  weight in the integration of k-point
       integer nx,ny,nz
       integer nrk  
       
  end type
  
  contains
  
    subroutine UnpackEv(this, ev_irred, ev_full, nband, nkpt)
       implicit none
       type(IrredBZ_t) :: this
       integer nband
       integer nkpt
       real(REAL64) ev_irred(nband,nkpt), ev_full(nband,nkpt)
       integer i,j,k,iband,irk,idk
       
       irk=1
       do i=1,this % nz
       do j=1,this % ny
       do k=1,this % nx
        idk = iabs(this % kmap(1,i,j,k))
!        write(*,*) i,j,k, idk
         do iband=1,nband
           ev_full(iband,irk) = ev_irred(iband,idk)
         enddo
         irk=irk+1
       enddo
       enddo
       enddo
       
       
    end subroutine UnpackEv
  
    subroutine IrredBZInit(this,adot,ntrans,mtrx,nx,ny,nz,sx,sy,sz)
       implicit none
       type(IrredBZ_t) :: this
       integer mxdpnt
       real(REAL64)         ::  adot(3,3)
       integer              ::  ntrans                     !  number of symmetry operations in the factor group
       integer              ::  mtrx(3,3,48)               !  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
       integer nx,ny,nz
       real(REAL64)         ::  sx, sy, sz

       integer, parameter:: mxdbnd = 1
       integer, parameter:: nbandi=1
              
       mxdpnt = nx*ny*nz
       
       this % nx = nx
       this % ny = ny
       this % nz = nz
       
       allocate(this % kmap(3,nx,ny,nz))
       allocate(this % wrk_nband(mxdpnt))
       allocate(this % indk(6,mxdpnt))
       allocate(this % rk(3,mxdpnt))
       allocate(this % w(mxdpnt))
   
       call int_pnt(nbandi,nx,ny,nz,sx,sy,sz,2,                          &
     & adot,                                                             &
     & ntrans,mtrx,                                                      &
     & this % nrk,this % rk,this % w,this % wrk_nband,this % indk,this % kmap, &
     & mxdpnt,mxdbnd)
    
    
    end subroutine IrredBZInit

end module
