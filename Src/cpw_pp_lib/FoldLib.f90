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

! FoldLib. A set of subroutines used in folding/unfolding of k-points. 
! ver 2.0. CLR, July 2014.
! Style modifications, July 15, 2014. JLM
! copyright  Jose Luis Martins, Carlos Loia Reis /INESC-MN.

subroutine Fold_GetPkn(Pkn, iMinv, idet, irk, kgv, isort, xvec, nkpt, nband, mtxd, ng, mxdgve, mxddim, mxdbnd)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                     !  array dimension of plane-waves
  integer, intent(in)                ::  mxdgve                     !  array dimension of G-space vectors
  integer, intent(in)                ::  mxdbnd                     !  array dimension for the number of bands

  integer, intent(in)                ::  irk

  integer, intent(in)                ::  mtxd                       !  wavefunction dimension
  integer, intent(in)                ::  isort(mxddim)              !  G-vector corresponding to coefficient i of wavefunction 
  integer, intent(in)                ::  ng                         !  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)              !  G-vectors in reciprocal lattice coordinates 

  integer, intent(in)                ::  nkpt,nband
  complex(REAL64), intent(in)        ::  xvec(mxddim,mxdbnd)

  integer, intent(in)                ::  iMinv(3,3)
  integer, intent(in)                ::  idet

! output

  real(REAL64), intent(out)          ::  Pkn(nkpt,nband) 

! local variables

  integer  ::  q(3), p(3), r(3)

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counters
  
  integer  :: i,j


  Pkn(irk,:) = ZERO

! q are the reciprocal lattice vectors components in the supercell: 
! G_i = q(1)B(1) + q(2)B(2) + q(3)B(3) with q(j) integers.
! p are the corresponding reciprocal lattice vectors components in the
! reference cell. g_i = p(1)b(1) + p(2)b(2) + p(3)b(3).
! When the p(j) are integers then the g_i are reciprocal lalttice vectors of
! the reference cell. The components are given by: p() = {[iMinv]/idet} q().
! [iMinv] q() = p() idet.
! if [iMinv] q() / idet has a remainder zero then p is an integer.

  do i=1,mtxd
    q(1) = kgv(1,isort(I))
    q(2) = kgv(2,isort(I))
    q(3) = kgv(3,isort(I))
    call Fold_MultV3I(iMinv,q,p)
    r(1) = mod(abs(p(1)),abs(idet))    
    r(2) = mod(abs(p(2)),abs(idet))    
    r(3) = mod(abs(p(3)),abs(idet))    
    if(r(1)==0 .and. r(2)==0 .and. r(3)==0 ) then
      do j=1,nband
        Pkn(irk,j)= Pkn(irk,j) + real(xvec(i,j)*conjg(xvec(i,j)),REAL64)
      enddo         
    endif
  enddo

  return
end subroutine Fold_GetPkn


subroutine Fold_GetPknSO(Pkn, iMinv, idet, irk, kgv, isort, xvec, nkpt, nband, mtxd, ng, mxdgve, mxddim, mxdbnd)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                     !  array dimension of plane-waves
  integer, intent(in)                ::  mxdgve                     !  array dimension of G-space vectors
  integer, intent(in)                ::  mxdbnd                     !  array dimension for the number of bands

  integer, intent(in)                ::  irk

  integer, intent(in)                ::  mtxd                       !  wavefunction dimension
  integer, intent(in)                ::  isort(mxddim)              !  G-vector corresponding to coefficient i of wavefunction 
  integer, intent(in)                ::  ng                         !  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)              !  G-vectors in reciprocal lattice coordinates 

  integer, intent(in)                ::  nkpt,nband
  complex(REAL64), intent(in)        ::  xvec(mxddim,mxdbnd)

  integer, intent(in)                ::  iMinv(3,3)
  integer, intent(in)                ::  idet

! output

  real(REAL64), intent(out)          ::  Pkn(nkpt,nband) 

! local variables

  integer  ::  q(3), p(3), r(3)

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counters
  
  integer  :: i,j


  Pkn(irk,:) = ZERO

  
  do i=1,mtxd
    q(1) = kgv(1,isort(I))
    q(2) = kgv(2,isort(I))
    q(3) = kgv(3,isort(I))
    call Fold_MultV3I(iMinv,q,p)
    r(1) = mod(abs(p(1)),abs(idet))    
    r(2) = mod(abs(p(2)),abs(idet))    
    r(3) = mod(abs(p(3)),abs(idet))    
    if(r(1)==0 .and. r(2)==0 .and. r(3)==0 ) then
        do j=1,nband
          Pkn(irk,j)= Pkn(irk,j) + real( xvec(2*i-1,j)*conjg(xvec(2*i-1,j))       &
                                       + xvec(2*i,j)*conjg(xvec(2*i,j)) ,REAL64)
        enddo         
    endif
  enddo

  return
end subroutine Fold_GetPknSO

subroutine Fold_Get_adot_pc(pwline, avec, adot_pc)

   implicit none

   integer, parameter          :: REAL64 = selected_real_kind(12)
   
!  Input

    character(len=60),intent(in)      ::  pwline
    real(REAL64), intent(in)  :: avec(3,3)

!  Output
    real(REAL64), intent(out)  :: adot_pc(3,3)

!  Local

   integer ::  idetpw
   integer :: iMpw(3,3),IMpwinv(3,3)

!  JLM begin

   character(len=6)  ::  fccSL
   integer :: ioerr

!  JLM end
   
   real(REAL64)   ::  avec_pc(3,3)
   real(REAL64)   ::  avec_sc(3,3)

!  constants

   real(REAL64), parameter :: ZERO = 0.0_REAL64

!  counters
   
   integer  ::  i, j, k


   avec_sc(:,:) = avec(:,:)

!  JLM begin
        
   read(pwline,'(3(3i4,2x),2x,a6)',IOSTAT=ioerr) ((iMpw(i,j),i=1,3),j=1,3),fccSL
   
   if(fccSL /= 'fcc SL' .or. ioerr /= 0) then
     write(6,*)
     write(6,*) '  STOPPED in FoldLib:   unable to find supercell information'
     write(6,*)

!  JLM end
     
     stop

   endif

   call Fold_Inv33I(iMpw,iMpwinv, idetpw)
      
   do i=1,3
   do j=1,3
     avec_pc(j,i) = ZERO
     do k=1,3
       avec_pc(j,i) = avec_pc(j,i) + avec_sc(j,k)*iMpwinv(k,i)
     enddo
     avec_pc(j,i) = avec_pc(j,i)/idetpw
   enddo
   enddo
   
       adot_pc(1,1) = avec_pc(1,1)*avec_pc(1,1) + avec_pc(2,1)*avec_pc(2,1) +           &
     &             avec_pc(3,1)*avec_pc(3,1)
       adot_pc(2,2) = avec_pc(1,2)*avec_pc(1,2) + avec_pc(2,2)*avec_pc(2,2) +           &
     &             avec_pc(3,2)*avec_pc(3,2)
       adot_pc(3,3) = avec_pc(1,3)*avec_pc(1,3) + avec_pc(2,3)*avec_pc(2,3) +           &
     &             avec_pc(3,3)*avec_pc(3,3)
       adot_pc(1,2) = avec_pc(1,1)*avec_pc(1,2) + avec_pc(2,1)*avec_pc(2,2) +           &
     &             avec_pc(3,1)*avec_pc(3,2)
       adot_pc(1,3) = avec_pc(1,1)*avec_pc(1,3) + avec_pc(2,1)*avec_pc(2,3) +           &
     &             avec_pc(3,1)*avec_pc(3,3)
       adot_pc(2,3) = avec_pc(1,2)*avec_pc(1,3) + avec_pc(2,2)*avec_pc(2,3) +           &
     &             avec_pc(3,2)*avec_pc(3,3)
       adot_pc(2,1) = adot_pc(1,2)
       adot_pc(3,1) = adot_pc(1,3)
       adot_pc(3,2) = adot_pc(2,3)

   return
end subroutine Fold_Get_adot_pc


subroutine Fold_Prep(pwline, avec,rk_in, iMinv, idet, rk_out, nkpt)

   implicit none

   integer, parameter          :: REAL64 = selected_real_kind(12)
   
!  Input 

    character(len=60),intent(in)      ::  pwline
    real(REAL64), intent(in)          :: avec(3,3)
    integer, intent(in)               :: nkpt
    real(REAL64), intent(in)          :: rk_in(3,nkpt)
    
!   Output

    integer, intent(out)         :: iMinv(3,3)
    integer, intent(out)         :: idet
    real(REAL64), intent(out)    :: rk_out(3,nkpt)
    
!   Local    
   
   integer        :: iMpw(3,3)

   real(REAL64)   :: avec_sc(3,3)
   real(REAL64)   :: avec_sc_inv(3,3)

   real(REAL64)   :: avec_pc(3,3)
   real(REAL64)   :: avec_pc_inv(3,3)
   
   real(REAL64)   :: M(3,3)
   real(REAL64)   :: Minv(3,3)
   
   real(REAL64)   :: bvec_pc(3,3)
   real(REAL64)   :: bvec_sc(3,3)
   
   integer        :: iM(3,3)
      

!  JLM start

   integer        :: iMpwinv(3,3),idetpw
   character(len=6)   ::  fccSL
   integer :: ioerr
   
!  JLM end

!      counters

   integer i, j, k, irk

!  constants

   real(REAL64), parameter :: ZERO = 0.0_REAL64
   real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64


   avec_sc(:,:) = avec(:,:)

!  JLM start
  
   read(pwline,'(3(3i4,2x),2x,a6)',IOSTAT=ioerr) ((iMpw(i,j),i=1,3),j=1,3),fccSL  

   if(fccSL /= 'fcc SL' .or. ioerr /= 0)  then
     write(6,*)
     write(6,*) '  STOPPED in FoldLib:   unable to find supercell information'
     write(6,*)
     
     stop

   endif


   call Fold_Inv33I(iMpw,iMpwinv, idetpw)
      
   do i=1,3
   do j=1,3
     avec_pc(j,i) = ZERO
     do k=1,3
       avec_pc(j,i) = avec_pc(j,i) + avec_sc(j,k)*iMpwinv(k,i)
     enddo
     avec_pc(j,i) = avec_pc(j,i)/idetpw
   enddo
   enddo
   
!  JLM end

   call Fold_Tr33(avec_pc)

   call Fold_Tr33(avec_sc)

   call Fold_Inv33(avec_pc, avec_pc_inv)

   bvec_pc(:,:) = 2*PI*avec_pc_inv(:,:)
   
   call Fold_Tr33(bvec_pc)
   
   call Fold_Inv33(avec_sc, avec_sc_inv)
   
   bvec_sc(:,:) = 2*PI*avec_sc_inv(:,:)
   
   call Fold_Tr33(bvec_sc)

   call Fold_Mult33T2(avec_pc,bvec_pc,M)
   
   call Fold_Mult33(avec_sc,bvec_sc,M)

   call Fold_Mult33(avec_sc,avec_pc_inv,M)
   
   call Fold_Inv33(M, Minv)
   
   call Fold_MtoIM(M,iM)

   call Fold_Inv33I(iM, iMinv, idet)
   
   do irk=1,nkpt
      call Fold_MultV3(M,rk_in(:,irk),rk_out(:,irk))
   enddo

   return
end subroutine Fold_Prep



subroutine Fold_MultV3(M,v1,v2)

   implicit none

   integer, parameter          :: REAL64 = selected_real_kind(12)

!  input

   REAL(REAL64), intent(in)    ::  M(3,3),v1(3)
   
!  output

   REAL(REAL64), intent(out)   ::  v2(3)

!  local

   REAL(REAL64)  ::  x

!  constants

   real(REAL64), parameter :: ZERO = 0.0_REAL64

!  counters
   integer  ::  i, k

   do i=1,3
      x = ZERO
      do k=1,3
         x = x + M(i,k)*v1(k)
      enddo
      v2(i) = x
   enddo
   
   return
end subroutine Fold_MultV3



subroutine Fold_MultV3I(M,v1,v2)

   implicit none

!  input

   integer, intent(in)    ::  M(3,3),v1(3)
   
!  output

   integer, intent(out)   ::  v2(3)

!  local

   integer  :: x

!  counters

   integer  ::  i, k

   do i=1,3
      x = 0
      do k=1,3
         x = x + M(i,k)*v1(k)
      enddo
      v2(i) = x
   enddo
   
   return
end subroutine Fold_MultV3I



subroutine Fold_MultV32(M,v1,v2)

   implicit none

   integer, parameter          :: REAL64 = selected_real_kind(12)

!  input

   REAL(REAL64), intent(in)   ::  M(3,3),v1(3)

!  output

   REAL(REAL64), intent(out)  ::  v2(3)

!  local

   REAL(REAL64)  ::  x

!  constants

   real(REAL64), parameter :: ZERO = 0.0_REAL64

!  counters

   integer  ::  i, k

   do i=1,3
      x = ZERO
      do k=1,3
         x = x + M(k,i)*v1(k)
      enddo
      v2(i) = x
   enddo
   
   return
end subroutine Fold_MultV32



subroutine Fold_Mult33(M1,M2,M)

   implicit none

   integer, parameter          :: REAL64 = selected_real_kind(12)

!  input

   REAL(REAL64), intent(in)  ::  M1(3,3),M2(3,3)

!  output

   REAL(REAL64), intent(out) ::  M(3,3)

!  local

   REAL(REAL64)  ::  x

!  constants

   real(REAL64), parameter :: ZERO = 0.0_REAL64

!  counters

   integer  ::  i, j, k

   do i=1,3
      do j=1,3
         x = ZERO
         do k=1,3
            x = x + M1(i,k)*M2(k,j)
         enddo
         M(i,j) = x
      enddo
   enddo
   
   return
end subroutine Fold_Mult33



subroutine Fold_Mult33T(M1,M2,M)

   implicit none

   integer, parameter          :: REAL64 = selected_real_kind(12)

!  input

   REAL(REAL64), intent(in)  ::  M1(3,3),M2(3,3)

!  output

   REAL(REAL64), intent(out)  ::  M(3,3)

!  local

   REAL(REAL64)  ::  x

!  constants

   real(REAL64), parameter :: ZERO = 0.0_REAL64

!  counters

   integer  ::  i, j, k

   do i=1,3
      do j=1,3
         x = ZERO
         do k=1,3
            x = x + M1(k,i)*M2(k,j)
         enddo
         M(i,j) = x
      enddo
   enddo
   
   return
end subroutine Fold_Mult33T



subroutine Fold_Mult33T2(M1,M2,M)

   implicit none

   integer, parameter          :: REAL64 = selected_real_kind(12)

!  input

   REAL(REAL64), intent(in)  ::  M1(3,3),M2(3,3)

!  output
 
   REAL(REAL64), intent(out)  ::  M(3,3)

!  local

  REAL(REAL64)  ::  x

!  constants

   real(REAL64), parameter :: ZERO = 0.0_REAL64

!  counters

   integer  ::  i, j, k
   
   do i=1,3
      do j=1,3
         x = ZERO
         do k=1,3
            x = x + M1(i,k)*M2(j,k)
         enddo
         M(i,j) = x
      enddo
   enddo
   
   return
end subroutine Fold_Mult33T2



subroutine Fold_Tr33(M)

   implicit none

   integer, parameter          :: REAL64 = selected_real_kind(12)

!  input and output

   REAL(REAL64), intent(inout)  ::  M(3,3)

!  local

   REAL(REAL64)  ::   Mtr(3,3)

!  counters

   integer ::  i, j
   
   do i=1,3
    do j=1,3
      Mtr(i,j) = M(j,i)
    enddo
   enddo
   do i=1,3
    do j=1,3
      M(i,j) = Mtr(i,j)
    enddo
   enddo
   
   return
end subroutine Fold_Tr33



subroutine Fold_MtoIM(M,iM)

  implicit none

   integer, parameter          :: REAL64 = selected_real_kind(12)

!  input

 real(REAL64), intent(in)  ::   M(3,3)
  
!  output
  
  integer, intent(out)     ::  iM(3,3)

!  counters

  integer  ::  i, j

   do i=1,3
    do j=1,3
      iM(i,j) = NINT(M(i,j))
    enddo
   enddo
   
   return
end subroutine Fold_MtoIM



subroutine Fold_Inv33I(iM,iMinv2, idet)

   implicit none

!  input

   integer, intent(in)      :: iM(3,3)
  
!  output
  
   integer, intent(out)     :: iMinv2(3,3)  
   integer, intent(out)     :: idet

!  local

   integer   :: iMinv(3,3)

!  counters

   integer   :: i, j

   iMinv(1,1) = iM(2,2)*iM(3,3) - iM(3,2)*iM(2,3)
   iMinv(2,1) = iM(3,2)*iM(1,3) - iM(1,2)*iM(3,3)
   iMinv(3,1) = iM(1,2)*iM(2,3) - iM(2,2)*iM(1,3)
   iMinv(1,2) = iM(2,3)*iM(3,1) - iM(3,3)*iM(2,1)
   iMinv(2,2) = iM(3,3)*iM(1,1) - iM(1,3)*iM(3,1)
   iMinv(3,2) = iM(1,3)*iM(2,1) - iM(2,3)*iM(1,1)
   iMinv(1,3) = iM(2,1)*iM(3,2) - iM(3,1)*iM(2,2)
   iMinv(2,3) = iM(3,1)*iM(1,2) - iM(1,1)*iM(3,2)
   iMinv(3,3) = iM(1,1)*iM(2,2) - iM(2,1)*iM(1,2)

   idet = iMinv(1,1)*iM(1,1) + iMinv(2,1)*iM(2,1) + iMinv(3,1)*iM(3,1)
     
   do i=1,3
    do j=1,3
      iMinv2(i,j) = iMinv(j,i)
    enddo
   enddo
   
   return
end subroutine Fold_Inv33I



subroutine Fold_Inv33(M,Minv2)

   implicit none

   integer, parameter          :: REAL64 = selected_real_kind(12)

!  input

   REAL(REAL64)  ::  M(3,3)
  
!  output

   REAL(REAL64)  ::  Minv2(3,3)

!  local

   REAL(REAL64)  ::  det
   REAL(REAL64)  ::  Minv(3,3)

!  counters

   integer  ::  i,j
   
   Minv(1,1) = M(2,2)*M(3,3) - M(3,2)*M(2,3)
   Minv(2,1) = M(3,2)*M(1,3) - M(1,2)*M(3,3)
   Minv(3,1) = M(1,2)*M(2,3) - M(2,2)*M(1,3)
   Minv(1,2) = M(2,3)*M(3,1) - M(3,3)*M(2,1)
   Minv(2,2) = M(3,3)*M(1,1) - M(1,3)*M(3,1)
   Minv(3,2) = M(1,3)*M(2,1) - M(2,3)*M(1,1)
   Minv(1,3) = M(2,1)*M(3,2) - M(3,1)*M(2,2)
   Minv(2,3) = M(3,1)*M(1,2) - M(1,1)*M(3,2)
   Minv(3,3) = M(1,1)*M(2,2) - M(2,1)*M(1,2)

   det = Minv(1,1)*M(1,1) + Minv(2,1)*M(2,1) + Minv(3,1)*M(3,1)
     
   do i=1,3
    do j=1,3
      Minv2(i,j) = Minv(j,i)/det
    enddo
   enddo

   return
end subroutine Fold_Inv33
