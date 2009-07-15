      program conv_dat

! Program to convert big<->little endian data in files

       integer in_dat, out_dat
       character*80 infile,outfile


!- 1 - Get the user's input

       indx = iargc( )
       if(indx.ne.2)then
        print*,' Usage: conv_dat infile outfile'
        stop
       endif

       call getarg(1,infile)
       call getarg(2,outfile)


!- 2 - Open files

       open (10,file= infile,status='old',form='unformatted'
!WORD     &          ,access='direct',recl=1)
     &          ,access='direct',recl=4)
       open (11,file=outfile,status='unknown',form='unformatted'
!WORD     &          ,access='direct',recl=1)
     &          ,access='direct',recl=4)

!- 3 - Practically infinite loop:  reads, converts and writes words.

        infnt = int(1.e9)
        do i=1,infnt

           read(10,rec=i,err=1001) in_dat

           out_dat = int_conv(in_dat)

           write(11,rec=i) out_dat

        enddo

!- 4 - finish the program

        print*, 'The conversion did not complete.'
        print*, 'The file is too big.'

 1001 continue

      end

!--------------------------------------------------

      function int_conv(in_dat)

       integer i1,i2
       character ch1(4), ch2(4)
       equivalence(i1,ch1)
       equivalence(i2,ch2)

        i1=in_dat

          do k=1,4
           ch2(k) = ch1(5-k)
          enddo

        int_conv = i2

      end 
