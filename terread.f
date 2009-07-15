      program terread

      parameter (rtd=180.0/3.14159265 )
      parameter ( nter = 37 )

      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'map.h'   ! em
      include 'xyzinfo.h'   ! rlat,rlong
      include 'indices.h'
      include 'parm.h'   ! to pass rlong0,rlat0,schmidt  to setxyz
      logical olam

      include 'rwork.h' ! rmsk,inum,inumx,zss,almsk,tmax,tmin,tsd,rlatd,rlond,grid,id1km,rlonx,rlonn,rlatx,rlatn

      logical debug, do1km, okx, oky, do250
      character*60 fileout
      character*9 formout

      integer nx,ny
      integer i,j
      character*1 ns,ew
      character*11 file
      character*8 ausfile(nter)

      real lons,lats
      real lone,late

      data debug/.true./
      data do1km/.true./
      data do250/.true./
      data idia/24/,jdia/72/
      data id/2/,jd/2/

      namelist / topnml / ds, du, tanl, rnml, stl1, stl2, debug
     &  ,luout, fileout, olam, wbd, sbd, dlon, dlat
     &  ,idia, jdia ,rlong0, rlat0, schmidt
     &  ,do1km, do250, id, jd

      open ( unit=5,file='top.nml',status='unknown' )
      read ( 5,topnml, end=5 )
 5    write ( 6,topnml )

      if ( do250 ) then
        do n = 1,37
          read(88,*) ausfile(n)
          write(6,*) ausfile(n)
        enddo ! n = 1,37
      endif

      call setxyz

      gridx=-1.
      rlatx=-999.
      rlonx=-999.
      gridn=+999.
      rlatn=+999.
      rlonn=+999.
      !imax=-999
      !jmax=-999
      !imin=+999
      !jmin=+999
      do j=1,jl
       do i=1,il
         n=i+(j-1)*il
         grid(i,j)=(ds/em(i,j))/1.e3 ! km
         rlond(i,j)=rlong(n)*rtd
         if(rlond(i,j).gt.180.)rlond(i,j)=rlond(i,j)-360.
         rlatd(i,j)=rlat(n)*rtd
         gridx=max(gridx,grid(i,j))
         gridn=min(gridn,grid(i,j))
         if ( grid(i,j).lt.20. ) then    ! use 1km when grid spc < 20 km
           rlatx=max(rlatx,rlatd(i,j))
           rlonx=max(rlonx,rlond(i,j))
           rlatn=min(rlatn,rlatd(i,j))
           rlonn=min(rlonn,rlond(i,j))
           !jmax=max(jmax,j)
           !imax=max(imax,i)
           !jmin=min(jmin,j)
           !imin=min(imin,i)
           !write(6,*)i,j,n,grid(i,j),rlond(i,j),rlatd(i,j)
         endif ! grid < 20 km
         if ( j.eq.49 )write(6,*)i,j,rlond(i,j),rlatd(i,j),grid(i,j)
         if ( (il+1.lt.j .and. j.lt.2*il) .and. i.eq.25 )
     &         write(6,*)i,j,rlond(i,j),rlatd(i,j),grid(i,j)
         if(rlond(i,j).gt.180.)rlond(i,j)=rlond(i,j)-360.	 
       enddo
      enddo
      write(6,*)"grid x,n=",gridx,gridn
      write(6,*)"rlon x,n=",rlonx,rlonn
      write(6,*)"rlat x,n=",rlatx,rlatn
      !write(6,*)"i    x,n=",imax,imin
      !write(6,*)"j    x,n=",jmax,jmin

      !write(97,'(30f6.1)') grid
      !write(98,'(30f6.2)') rlond
      !write(99,'(30f6.2)') rlatd

      write(6,*)"grid spc(km)"
      do j=1,jl
        do i=1,il
          inumx(i,j)= nint(grid(i,j))
        enddo !i=1,il
      enddo !j=1,jl
      do j=96,49,-1
       write(6,'(48i2)')(inumx(i,j),i=1,48)
      enddo

c     nested model topography (output)
      write(6,*) 'open',luout,fileout
      open(luout,file=fileout,form='formatted',status='unknown')

! netcdf file
      call ncdf_setup("zsm.nc",idnc,3,il,jl,1,20020101,0000)
      call nc2out(grid ,il,jl,1,1.,idnc,"grid","grid","km",0.,5000.)
      call nc2out(rlond,il,jl,1,1.,idnc,"lon","lon","degrees",0.,360.)
      call nc2out(rlatd,il,jl,1,1.,idnc,"lat","lat","degrees",-90.,90.)
      call ncsnc(idnc,ier)

c initialize min,max, and sd arrays
      do j=1,jl
       do i=1,il
        tmin(i,j)=99999.
        tmax(i,j)=-99999.
        tsd(i,j)=0.
        almsk(i,j)=0.
        zss(i,j)=0.
        inum(i,j)=0
       enddo
      enddo
      
      write(6,*) 'inum(1,1)=',inum(0,0)
      
!================================
      if ( do250 ) then
!================================
      write(6,*)"Proccess Australian 250m files"

      do n = 1,37

        open(13,file=ausfile(n),form="unformatted")
        read(13) nx,ny,lons,lats,dl
        lone=lons+(nx-1)*dl
        late=lats+(ny-1)*dl

        okx = .false.
        okx = okx .or. (lons.ge.rlonn .and. lons.le.rlonx)
        okx = okx .or. (lone.ge.rlonn .and. lone.le.rlonx)
        okx = okx .or. (lons.le.rlonn .and. lone.ge.rlonx)
        oky = .false.
        oky = oky .or. (lats.ge.rlatn .and. lats.le.rlatx) ! start lat in area
        oky = oky .or. (late.ge.rlatn .and. late.le.rlatx) !  end  lat in area
        oky = oky .or. (lats.le.rlatn .and. late.ge.rlatx) ! grid encompasses hr reg
        oky = oky .or. (lats.ge.rlatn .and. late.le.rlatx) ! grid within hr reg

        write(6,'("lons,lone=",2f6.1," lats,late=",2f6.1,2l4,2i6,a10)')
     &             lons,lone,lats,late,okx,oky,nx,ny,ausfile(n)

        if ( okx.and.oky ) then
          write(6,*)"call read250(nx,ny,lons,lats,dl,debug,idia,jdia)"
          call read250(nx,ny,lons,lats,dl,debug,idia,jdia)
        endif

      enddo ! n

      write(6,*)"inum"
      do j=96,49,-1
       write(6,'(48i2)')(inum(i,j),i=1,48)
      enddo

      write(6,*)"almsk"
      do j=1,jl
        do i=1,il
          if ( inum(i,j) .gt. 0 ) then
             inumx(i,j)=nint(real(almsk(i,j))/real(inum(i,j)))
          else
             inumx(i,j)=0
          endif
        enddo !i=1,il
      enddo !j=1,jl
      do j=96,49,-1
       write(6,'(48i2)')(inumx(i,j),i=1,48)
      enddo

      write(6,*)"zss/10"
      do j=1,jl
        do i=1,il
          if ( inum(i,j) .gt. 0 ) then
             inumx(i,j)=nint(zss(i,j)/(10.*real(inum(i,j))))
          else
             inumx(i,j)=0
          endif
        enddo !i=1,il
      enddo !j=1,jl
      do j=96,49,-1
       write(6,'(48i2)')(inumx(i,j),i=1,48)
      enddo


      write(6,*)"Done with Australian 250m files"

!================================
      endif ! ( do250 ) then
!================================

!================================
      write(6,*)"################## rlat x,n=",rlatx,rlatn
      do1km = do1km .and. ( rlatx .gt. rlatn )
      write(6,*)"###################### do1km=",do1km
      if ( do1km ) then
!================================
      write(6,*)"Processing DEM files"
      dl = 1.0/120.0
      dlh = 1.0/240.0
      nx=4800
      ny=6000
      do i = 1,9
       lons=-180.+real(i-1)*40.
       lone=lons+real(nx-1)*dl
       ew="W"
       if ( nint(lons).gt.0 ) ew="E"

       do j = 1,3
        lats=90.-real(j-1)*50.
        late=lats-real(ny-1)*dl 
        ns="N"
        if ( nint(lats).lt.0 ) ns="S"

        write(file,'(a1,i3.3,a1,i2.2,".DEM")')
     &               ew,abs(nint(lons)),ns,abs(nint(lats))

        okx = .false.
        okx = okx .or. (lons.gt.rlonn .and. lons.lt.rlonx)
        okx = okx .or. (lone.gt.rlonn .and. lone.lt.rlonx)
        okx = okx .or. (lons.lt.rlonn .and. lone.gt.rlonx)
        okx = okx .or. (lons.lt.-179. .and. rlonx.gt. 180.) ! fix for dateline data
        okx = okx .or. (lons.gt. 179. .and. rlonn.lt.-180.) ! fix for dateline data
        oky = .false.
        oky = oky .or. (lats.gt.rlatn .and. lats.lt.rlatx)
        oky = oky .or. (late.gt.rlatn .and. late.lt.rlatx)
        oky = oky .or. (lats.gt.rlatx .and. late.lt.rlatn)

        write(6,'("lons,lone=",2f6.1," lats,late=",2f6.1,a12,2l4)')
     &             lons,lone,lats,late,file,okx,oky

        if ( okx.and.oky ) then
          write(6,*)"call read1km(file,nx,ny,lons,lats,debug,idia,jdia)"
! 9 sept 2004
          clons=lons+dlh
          clats=lats-dlh
          call read1km(file,nx,ny,clons,clats,debug,idia,jdia)
          !call read1km(file,nx,ny,lons,lats,debug,idia,jdia)
! 9 sept 2004
        endif

       enddo ! j

      enddo ! i

      nx=7200
      ny=3600
      do i = 1,6
        lons=-180.+real(i-1)*60.
        lone=lons+real(nx-1)*dl
        ew="W"
        if ( nint(lons).gt.0 ) ew="E"
        lats=-60.
        late=lats-real(ny-1)*dl 
        ns="S"

        write(file,'(a1,i3.3,a1,i2.2,".DEM")')
     &               ew,abs(nint(lons)),ns,abs(nint(lats))

        okx = .false.
        okx = okx .or. (lons.gt.rlonn .and. lons.lt.rlonx)
        okx = okx .or. (lone.gt.rlonn .and. lone.lt.rlonx)
        okx = okx .or. (lons.lt.rlonn .and. lone.gt.rlonx)
        oky = .false.
        oky = oky .or. (lats.gt.rlatn .and. lats.lt.rlatx)
        oky = oky .or. (late.gt.rlatn .and. late.lt.rlatx)
        oky = oky .or. (lats.gt.rlatx .and. late.lt.rlatn)

        write(6,'("lons,lone=",2f6.1," lats,late=",2f6.1,a12,2l4)')
     &             lons,lone,lats,late,file,okx,oky

        if ( okx.and.oky ) then
          write(6,*)"call read1km(file,nx,ny,lons,lats,debug,idia,jdia)"
! 9 sept 2004
          clons=lons+dlh
          clats=lats-dlh
          call read1km(file,nx,ny,clons,clats,debug,idia,jdia)
          !call read1km(file,nx,ny,lons,lats,debug,idia,jdia)
! 9 sept 2004
        endif

      enddo ! i

      write(6,*)"after do1km"

      write(6,*)"inum"
      do j=96,49,-1
       write(6,'(48i2)')(inum(i,j),i=1,48)
      enddo

!     write(6,*)"almsk"
!     do j=96,49,-1
!      write(6,'(48i2)')(nint(almsk(i,j)/real(inum(i,j))),i=1,48)
!     enddo
!     write(6,*)"zss"
!     do j=96,49,-1
!      write(6,'(48i2)')(nint(zss(i,j)/10/real(inum(i,j))),i=1,48)
!     enddo

      write(6,*)"almsk"
      do j=1,jl
        do i=1,il
          if ( inum(i,j) .gt. 0 ) then
             inumx(i,j)=nint(almsk(i,j)/real(inum(i,j)))
          else
             inumx(i,j)=0
          endif
        enddo !i=1,il
      enddo !j=1,jl
      do j=96,49,-1
       write(6,'(48i2)')(inumx(i,j),i=1,48)
      enddo

      write(6,*)"zss/10"
      do j=1,jl
        do i=1,il
          if ( inum(i,j) .gt. 0 ) then
             inumx(i,j)= nint(zss(i,j)/10/real(inum(i,j)))
          else
             inumx(i,j)=0
          endif
        enddo !i=1,il
      enddo !j=1,jl
      do j=96,49,-1
       write(6,'(48i2)')(inumx(i,j),i=1,48)
      enddo

      write(6,*)"Done with DEM files"

!================================
      endif ! ( do1km ) then
!================================

!***********************************************************************

      write(6,*)"Now calling read10km"
      call read10km(debug,do1km)

!================================

      write(6,*)"after read10km"
      write(6,*)"inum"
      do j=96,49,-1
       write(6,'(48i2)')(inum(i,j),i=1,48)
      enddo
      write(6,*)"almsk"
      do j=96,49,-1
       write(6,'(48i2)')(nint(almsk(i,j)/real(inum(i,j))),i=1,48)
      enddo
      write(6,*)"zss"
      do j=96,49,-1
       write(6,'(48i2)')(nint(zss(i,j)/10/real(inum(i,j))),i=1,48)
      enddo



      ! PATCH !!!!
      Do i=1,il
        Do j=1,jl
        !  If ((zss(i,j)*real(inum(i,j))).LT.0.) almsk(i,j)=0.

          aglon=rlond(i,j)
          aglat=rlatd(i,j)
          If ((aglon.GT.135.).AND.(aglon.LT.140.)) then
            If ((aglat.GT.-32.).AND.(aglat.LT.-25.)) then
              almsk(i,j)=real(inum(i,j))
            End If
          End If
        End Do
      End Do


      write(6,*)"now compute grid average values"
      numzer=0
      imin=9999
      jmin=9999
      imax=0
      jmax=0

      do j=1,jl
       do i=1,il
          !write(6,*)i,j,inum(i,j),zss(i,j)
          if ( inum(i,j).eq.0 ) then
! no input points found, assume sea point
             numzer=numzer+1
             imin=min(i,imin)
             imax=max(i,imax)
             jmin=min(j,jmin)
             jmax=max(j,jmax)
             zss(i,j) = 0.
             rmsk(i,j)= 0.
             tsd(i,j)= 0.
             tmin(i,j)= 0.
             tmax(i,j)= 0.
          else
! set land-sea mask (rmsk)
             almsk(i,j) = almsk(i,j)/real(inum(i,j))
             zss(i,j) = zss(i,j)/real(inum(i,j))
	     
! if less than half of pnts are ocean, assume it is land pnt
             if ( almsk(i,j).le. 0.5 ) then
! ocean point
                rmsk(i,j) = 0.
                zss(i,j)=0. ! jjk added 18-8-2004
	  	    tsd(i,j)=0. !mjt 12-5-05
		        tmax(i,j)=0. !mjt 12-5-05
		        tmin(i,j)=0. !mjt 12-5-05
             else ! almsk
! land point
                rmsk(i,j) = 1.
                zss(i,j)=max(.1,zss(i,j))
             endif ! almsk
! compute sd
             tsd(i,j)=sqrt(abs(tsd(i,j)/real(inum(i,j))-zss(i,j)**2)) ! abs for rounding?
!            tsd(i,j)=sqrt(tsd(i,j)/real(inum(i,j))-zss(i,j)**2)
          end if

       enddo ! i
      enddo ! j

2     if(numzer.gt.0)then
!       this checks for non-data points, and inserts a neighbouring value
        print *,'**** numzer= ',numzer
        numzer=0
        do j=1,jl
         do i=1,il
          iq=i+(j-1)*il
          inumx(i,j)=inum(i,j)
          if(inum(i,j).eq.0)then
            iq2=0
            if(inum(in(iq),1).ne.0)then
              iq2=in(iq)
            endif
            if(inum(ie(iq),1).ne.0)then
              iq2=ie(iq)
            endif
            if(inum(iw(iq),1).ne.0)then
              iq2=iw(iq)
            endif
            if(inum(is(iq),1).ne.0)then
              iq2=is(iq)
            endif
            if(iq2.ne.0)then
              zss(i,j)=zss(iq2,1)
              rmsk(i,j)=rmsk(iq2,1)
              tsd(i,j)=tsd(iq2,1)
              almsk(i,j)=almsk(iq2,1)
              inumx(i,j)=1
            else
              numzer=numzer+1
            endif    !  (iq2.ne.0)
          endif      !  (inum(i,j).eq.0)
         enddo
        enddo
        do j=1,jl
         do i=1,il
          inum(i,j)=inumx(i,j)
         enddo ! i
        enddo ! j
      endif         ! (numzer.gt.0)
      if(numzer.gt.0)go to 2

      Write(6,*) 'final rmsk'
      Do j=96,49,-1
        Write(6,'(48i2)')(nint(rmsk(i,j)),i=1,48)
      End Do
      Write(6,*) 'final zs/10'
      Do j=96,49,-1
        Write(6,'(48i2)')(nint(zss(i,j)/10.),i=1,48)
      End Do

      write(luout,'(i3,i4,2f7.2,f6.3,f8.0,''  orog-mask-var'')')
     &                           il,jl,rlong0,rlat0,schmidt,ds
!      write(luout,*)il,jl,rlong0,rlat0,schmidt,ds,"orog-mask-var"

      call nc2out(zss  ,il,jl,1,1.,idnc,"zs","zs","m",-100.,30000.)
      call nc2out(rmsk ,il,jl,1,1.,idnc,"lsm","lsm","none",-1.,1.)
      call nc2out(tsd  ,il,jl,1,1.,idnc,"tsd","tsd","m",0.,30000.)
      call nc2out(tmax ,il,jl,1,1.,idnc,"zmax","zmax","m",-100.,30000.)
      call nc2out(tmin ,il,jl,1,1.,idnc,"zmin","zmin","m",-100.,30000.)

c     write out g*zs(il,jl) to formatted file
      do j=1,jl
        do i=1,il
          zss(i,j)=9.80616*zss(i,j)
        end do ! i=1,il
      end do ! j=1,jl
      ilout=min(il,30)
      write (formout,'(''('',i3,''f7.0)'')')ilout   !  i.e. (<il>f7.0)
      write(luout,formout) zss

c     write out land/sea mask to formatted file
      write (formout,'(''('',i3,''f4.1)'')')ilout   !  i.e. (<il>f4.1)
      write(luout,formout) rmsk

c     write out std.dev of top. to formatted file
      write (formout,'(''('',i3,''f6.0)'')')ilout   !  i.e. (<il>f6.0)
      write(luout,formout) tsd
      print *,'zss, rmsk and tsd written to unit',luout,' for il=',il

      do j=1,jl
        do i=1,il
          zss(i,j)=real(inum(i,j))
        end do ! i=1,il
      end do ! j=1,jl
      call nc2out(zss,il,jl,1,1.,idnc,"inum","inum","none",0.,2000.)
      call ncsnc(idnc,ier)

      stop
      end
c=======================================================================