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

subroutine ao_interpolation_write_hr(adot)
  use NonOrthoInterp

  implicit none
  integer, parameter           :: REAL64 = selected_real_kind(12)
  real(REAL64), parameter      :: bohr_to_ang= 0.52917720859_REAL64  
  real(REAL64), intent(in)     :: adot(3,3)
  real(REAL64)                 :: avec(3,3), bvec(3,3)
  
  complex(REAL64), allocatable :: H(:,:),S(:,:), Hw(:,:)
    
  real(REAL64) :: rkpt(3)
  integer nband, ikpt, i, j, k, n1 ,n2, n3, nws_1,nws_2,nws_3,jmax
  
  type (noiData_t) :: mtb 
  type(fiData_t)  :: ortho   

  character(len=1) :: yesno

  call NonOrthoInterpReadFromFile(mtb)
  
  nband = mtb%nband
  
  allocate (H(nband,nband) )
  allocate (S(nband,nband) )
  allocate (Hw(nband,nband) )

  write(*,*) 'current interpolation grid is ', mtb%fiData%mp_grid(:)

  write(6,*) '  do you want to change it (y/n) ?'
  read(5,*) yesno
  if(yesno == 'y' .or. yesno == 'Y') then
    write(*,*) 'enter new grid'
    read(5,*) n1,n2,n3  
  else
    n1 = mtb%fiData%mp_grid(1)
    n2 = mtb%fiData%mp_grid(2)
    n3 = mtb%fiData%mp_grid(3)
    write(*,*) 'using ' , n1,n2,n3
  endif
        
  write(*,*) 'current ws search size is ', mtb%fiData%ws_search_size(:)
  write(6,*) '  do you want to change it (y/n) ?'
  read(5,*) yesno
  if(yesno == 'y' .or. yesno == 'Y') then
    write(*,*) 'enter new ws_search_size'
    read(5,*) nws_1,nws_2,nws_3  
  else
    nws_1 = mtb%fiData%ws_search_size(1)
    nws_2 = mtb%fiData%ws_search_size(2)
    nws_3 = mtb%fiData%ws_search_size(3)
    write(*,*) 'using ' ,nws_1,nws_2,nws_3
  endif
        
  call fi_hamiltonian_setup(ortho,nband,n1,n2,n3,nws_1,nws_2,nws_3,mtb%fiData%real_metric)
  
  write(*,*) "orthogonalizing"

  ikpt = 1
  do i=1,ortho%mp_grid(1)
  do j=1,ortho%mp_grid(2)
  do k=1,ortho%mp_grid(3)
  
    rkpt(1) = ortho%kpt_latt(1,ikpt)
    rkpt(2) = ortho%kpt_latt(2,ikpt)
    rkpt(3) = ortho%kpt_latt(3,ikpt)
    
    call fi_hamiltonian_get_hk(mtb%fiData,rkpt,H, nband)
    call fi_hamiltonian_get_sk(mtb%fiData,rkpt,S, nband)
    
    call OrthoH(H,S,Hw,nband)
  
    ortho%Ham_k(:,:,ikpt) = Hw
    ortho%S_k(:,:,ikpt) = S
            
!    write(*,'("kpt", i5,3f8.5)') ikpt, rkpt

    call progress(ikpt,ortho%mp_grid(1)*ortho%mp_grid(2)*ortho%mp_grid(3))
    
    ikpt = ikpt +1
    
  
  enddo
  enddo
  enddo
    
  call fi_hamiltonian_get_hr(ortho)
                
!  write mtb_hr.dat file comptible with Wannier90   

  write(*,*) 
  write(*,*) "writing mtb_hr.dat and mtb_tb.dat"
  open(unit = 256 , file = "mtb_hr.dat", form="formatted")
  open(unit = 257 , file = "mtb_tb.dat", form="formatted")
  
  write(256,*) "generated by mtb"
  write(256,*) nband
  write(256,*) ortho%nrpts

  write(257,*) "generated by mtb"
  call adot_to_avec_sym(adot,avec,bvec)

  write(257,'(3e22.12)') (avec(k,1)*bohr_to_ang,k=1,3)
  write(257,'(3e22.12)') (avec(k,2)*bohr_to_ang,k=1,3)
  write(257,'(3e22.12)') (avec(k,3)*bohr_to_ang,k=1,3)

  write(257,*) nband
  write(257,*) ortho%nrpts
  
  do i=1, ortho%nrpts, 15      
    jmax = i + 14
    if(jmax > ortho%nrpts) jmax = ortho%nrpts        
    write(256,'(15i5)') (ortho%ndegen(j), j=i, jmax)
    write(257,'(15i5)') (ortho%ndegen(j), j=i, jmax)
  enddo
  

  write(257,*)
  do i=1, ortho%nrpts
        write(257,'(3i5)') ortho%irvec(:,i) 
    do j=1, nband
      do k=1, nband
        write(256,'(3i5, 2i5, 2e22.12)') ortho%irvec(:,i), j, k, 27.21138505_REAL64*ortho%Ham_r(j,k,i)      
        write(257,'(2i5, 2e22.12)')                        j, k, 27.21138505_REAL64*ortho%Ham_r(j,k,i)      
      enddo
    enddo
    write(257,*)

    call progress(i,ortho%nrpts)
    call flush(6)
  enddo     
  
  
 close(256)
 close(257)
 
 deallocate(H,S,Hw)
      
 write(*,*) "done"
  
  
  
end subroutine