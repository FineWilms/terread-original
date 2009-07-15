      program global_terrain

      parameter (rtd=180.0/3.14159265 )

      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'map.h'   ! em
      include 'xyzinfo.h'   ! rlat,rlong
      include 'indices.h'
      include 'parm.h'   ! to pass rlong0,rlat0,schmidt  to setxyz
      logical olam

      include 'rwork.h'

      logical debug, do1km, okx, oky
      character*60 fileout
      character*9 formout

      integer nx,ny
      integer i,j
      character*1 ns,ew
      character*11 file
      real lons,lats
      real lone,late

      data debug/.true./
      data do1km/.true./
      data idia/24/
      data jdia/72/

      namelist / topnml / ds, du, tanl, rnml, stl1, stl2, debug
     &  ,luout, fileout, olam, wbd, sbd, dlon, dlat
     &  ,idia, jdia ,rlong0, rlat0, schmidt
     &  ,do1km

      open ( unit=5,file='top.nml',status='unknown' )
      read ( 5,topnml, end=5 )
 5    write ( 6,topnml )

!================================
      dl = 1.0/120.0
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

        write(6,*)"call read1km(file,nx,ny,lons,lats,debug)"

        call read1km(file,nx,ny,lons,lats,debug)

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

        write(demfile,'(a1,i3.3,a1,i2.2,".DEM")')
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
     &             lons,lone,lats,late,demfile,okx,oky

        if ( okx.and.oky ) then
          write(6,*)"call read1km(demfile,nx,ny,lons,lats,debug)"
          call read1km(demfile,nx,ny,lons,lats,debug)
        endif

      enddo ! i

      stop
      end
c=======================================================================
      subroutine read1km(file,nx,ny,lons,lats,debug)

      character*11 file
      integer*2 zz(7200)
      real zzi(7200)
      real zzo(7200)
      real dl
      integer nrecl
      real lons,lats
      parameter ( zmin=-100 )

      logical debug,ok

      include 'newmpar.h'
      include 'rwork.h'

      nrecl = nx*2
      dl = 1.0/120.0

      write(6,*)"read1km file,nx,ny=",file,nx,ny

      write(6,*)"file=",file," recl=",nrecl
!     open(21,file=file,access='direct',form='unformatted',recl=nrecl)
      open(CONVERT='BIG_ENDIAN',unit=21,file=file,access='direct'
     &            ,form='unformatted',recl=nrecl)

!#######################################################################
      do j=1,ny
!#######################################################################

        aglat=(lats)-real(j-1)*dl
        write(6,*)j,aglat

        write(6,*)"lat,latn,latx,lons=",aglat,rlatn,rlatx,lons

         read(21,rec=j) (zz(i),i=1,nx)

!#######################################################################
         do i=1,nx
!#######################################################################

           n=i+(nj+1-j)*nx
           aglon=lons+real(i-1)*dl

!-----------------------------------------------------------------------
           ok=(aglat.gt.rlatn.and.aglat.lt.rlatx)
           ok=ok .and. (aglon.gt.rlonn.and.aglon.lt.rlonx)
           if ( ok ) then
           !write(6,*)"i,lon,lonn,lonx=",i,aglon,rlonn,rlonx
!-----------------------------------------------------------------------

! compute model grid i,j
           call latltoij(aglon,aglat,alci,alcj,nface)  ! con-cubic
           lci = nint(alci)
           lcj = nint(alcj)
! convert to "double" (i,j) notation
           lcj=lcj+nface*il
           !write(6,*)i,j,aglon,aglat,lci,lcj

! check to make sure within grid dimensions
           if(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then
             !write(6,*)"i,j,zi2(i)=",i,j,zi2(i)

             if( zi2(i).lt.1 ) then
! ocean point
               amask = 0.
               zs=0.
             else  ! zs>=-1000
! land point
               amask = 1.
               zs=real(zi2(i))
             end if  ! zs<-1000

             id1km(lci,lcj)=1
! accumulate number of pnts in grid box
             inum(lci,lcj) = inum(lci,lcj) + 1
! accumulate topog. pnts
             zss(lci,lcj) = zss(lci,lcj) + zs
! accumulate lmask pnts
             almsk(lci,lcj) = almsk(lci,lcj) + amask
! find max/min topog. pnts in grid box
             tmax(lci,lcj)=max(tmax(lci,lcj),zs)
             tmin(lci,lcj)=min(tmin(lci,lcj),zs)
! sum of squares for sd calc.
             tsd(lci,lcj)=tsd(lci,lcj)+zs**2

             if ( debug ) then
                if ( lci.eq.idia .and. lcj.eq.jdia ) then
                  write(6,'("lci,lcj,i,zi2(i),zs,zss(),almsk(), inum()="
     &           ,4i6,3f10.3,i6)') lci,lcj,i,zi2(i),zs,zss(lci,lcj)
     &                       ,almsk(lci,lcj), inum(lci,lcj)
                endif  ! selected points only
             endif  ! debug

           endif ! (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then

!-----------------------------------------------------------------------
           endif ! ( ok ) then
!-----------------------------------------------------------------------

!#######################################################################
         enddo ! i
        endif! (ok)then
!#######################################################################

        if ( mod(j,100).eq.0 .and. ok ) then
          write(6,*)"i,j,lon,lat,lci,lcj=",i,j,aglon,aglat,lci,lcj
        endif

!#######################################################################
      enddo ! j
!#######################################################################

      close(21)

      return
      end ! read1km
c=======================================================================
