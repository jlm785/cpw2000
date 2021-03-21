       subroutine compact_expand(array,na,nb,lda,cop)

!      compacts or expands a 2D array

!      written 29 October 2015. JLM
!      copyright INESC-MN/Jose Luis Martins

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  na                         !  first array data size
       integer, intent(in)                ::  nb                         !  second array data size
       integer, intent(in)                ::  lda                        !  leading dimension of array
       character(len=1), intent(in)       ::  cop                        !  toggles compact/expand

!      input and output

       complex(REAL64), intent(inout)     ::  array(lda*nb)

!      counters

       integer    ::  i, j, ic, ie

       if(nb > 1) then
         if(cop == 'C' .or. cop == 'c') then

           do i = 2,nb
             ie = (i-1)*lda
             ic = (i-1)*na
             do j = 1,na
               array(ic + j) = array(ie + j)
             enddo
           enddo

         elseif(cop == 'E' .or. cop == 'e') then

           do i = nb,2,-1
             ie = (i-1)*lda
             ic = (i-1)*na
             do j = na,1,-1
               array(ie + j) = array(ic + j)
             enddo
           enddo

         endif
       endif

       return
       end subroutine compact_expand
       
       


       
