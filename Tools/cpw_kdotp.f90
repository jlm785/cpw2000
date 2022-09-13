       program cpw_kdotp


       implicit none

       integer              :: ioreplay = 35                             !<  hard coded tape unit for replay
       character(len = 15)  :: fname = 'replay_post.dat'                 !<  hard coded filename for replay


       open(unit = ioreplay, file = fname,  status = 'UNKNOWN',          &
     &                   form = 'FORMATTED')

       if(ioreplay < 7) stop   'ioreplay < 7'

       write(6,*)
       write(6,*) '  Welcome to cpw_analysis'
       write(6,*)

!      loop over tasks


           call kdotp_sub(ioreplay)

       close(unit = ioreplay)

       stop
       end program cpw_kdotp
