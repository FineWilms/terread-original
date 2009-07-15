      program t

      character*1 ns,ew
      character*12 file
      integer slon,slat

      do i = 1,9
       slon=-140+(i-1)*40
       ew="E"
       if ( slon.gt.0 ) ew="W"

       do j = 1,3
        slat=90-(j-1)*50 
        ns="N"
        if ( slon.lt.0 ) ns="S"
        write(file,'(a1,i3.3,a1,i3.3,".DEM")') ew,abs(slon),ns,abs(slat)
        write(6,*) file
       enddo ! j

      enddo ! i

      do i = 1,6
       slon=-120+(i-1)*60
       ew="E"
       if ( slon.ge.0 ) ew="W"
       write(file,'(a1,i3.3,"S60.DEM")') ew,abs(slon)
       write(6,*) file
      enddo ! i

      stop
      end
