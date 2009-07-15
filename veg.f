      program veg

      include 'newmpar.h'

      parameter ( rtd=180.0/3.14159265 )
      parameter ( dl = 1.0/120.0 )
      parameter ( gslon = -180+dl/2. )
      parameter ( gslat = 90-dl/2. )
      parameter ( mx=43200, my=21600 )

      character*13 gfile
      parameter( gfile='gsib2_0ll.img')

      integer nsum(il,jl,0:20) ! 44 with old ones

      integer i,j
      byte zza(mx)
      integer iz(mx)
      character*9 formout
      character*80 vegin,vegout,topo,topout
      character*80 title

      include 'dates.h'   ! to pass ds from setxyz
      include 'indices.h'
      include 'map.h'     ! em
      include 'parm.h'    ! to pass rlong0,rlat0,schmidt  to setxyz
      include 'xyzinfo.h' ! rlat,rlong
      include 'rwork.h'   ! grid,tsd

      integer iarr(il,jl)
      real arr(ifull)
      real zs(ifull),lsm(ifull)
      real waterfrac(ifull),npoints(ifull),nsummax(ifull)

      namelist / vegnml / topo,vegin,vegout,topout,id,jd,ds1km

      data id/37/,jd/81/,ds1km/60./

      open (85,file='veg.nml',status='old')
      read (85,vegnml)
      write (6,vegnml)

      ijd = id+(jd-1)*il

      write(6,*)"read global 1 km data"
      write(*,*) 'gfile=',gfile
      open(20,file=gfile,access='direct',form='unformatted',recl=mx/4)
      jmin=1
      jmax=my
      do jg=jmin,jmax
! read each row of relevant data from 1 km data set (row=jg)
        read(20,rec=jg) (zza(ig),ig=1,mx) ! read in of byte data
        izmax=-9999
        izmin=+9999
        do ig=1,mx
          iz(ig)=zza(ig) ! convert to integer
          if(iz(ig).ge.19)iz(ig)=0
          izmax=max(izmax,iz(ig))
          izmin=min(izmin,iz(ig))
        enddo ! ig=1,mx
!       if(mod(jg,100).eq.0)then
!          write(6,*)"jg=",jg," izmax=",izmax," izmin=",izmin
!       endif! (mod(jg,100).eq.0)then
      enddo ! jg=jmin,jmax

!     if (jg.ne.0)stop

      write(6,*)"usual topo data"
      open (21,file=topo,status='old')
      read (21,'(i3,i4,2f6.2,f6.3,f8.0,a20)')
     &           ili,jli,rlong0,rlat0,schmidt,ds,title
      write(6,*) ili,jli,rlong0,rlat0,schmidt,ds,title
c     read in g*zs(il,jl) to formatted file
      ilout=min(il,30)
      write (formout,'(''(''i3''f7.0)'')')ilout   !  i.e. (<il>f7.0)
      read(21,formout) zs

c     write out land/sea mask to formatted file
      write (formout,'(''(''i3''f4.1)'')')ilout   !  i.e. (<il>f4.1)
      read(21,formout) lsm

c     write out std.dev of top. to formatted file
      write (formout,'(''(''i3''f6.0)'')')ilout   !  i.e. (<il>f6.0)
      read(21,formout) tsd
      print *,'zss, rmsk and tsd written to unit',luout,' for il=',il

! usual veg data and map proj info
      open (22,file=vegin,status='old')
      read (22,'(i3,i4,2f6.2,f6.3,f8.0,a20)')
     &           ili,jli,rlong0,rlat0,schmidt,ds,title
      write(6,*) ili,jli,rlong0,rlat0,schmidt,ds,title

      write (formout,'(1h(,i3,2hi3,1h))')ili   !  i.e. (<il>i3)
      write(6,*)"formout=",formout
      read (22,formout) iarr     ! read old vegetation type

      write(6,*)"old id,jd,zs,lsm,tsd,veg=",id,jd
     &             ,zs(ijd),lsm(ijd),tsd(id,jd),iarr(id,jd)

! open output file (vegout)
      open (23,file=vegout,status='unknown')
      open (24,file=topout,status='unknown')

! now start good stuff
      call setxyz

!     initialize summing array
      do ntype=0,20
       do iq=1,ifull
        nsum(iq,1,ntype)=0
       enddo ! iq
      enddo ! ntype
! preinit. arr array with old data
       do iq=1,ifull
         arr(iq)=iarr(iq,1)
       enddo ! iq

!      do j=1,jl
!        write(6,*)"il/2,j,iarr(il/2)=",il/2,j,iarr(il/2,j)
!      enddo ! iq

!     is = nint((rlons - gslon)/dl) + 1 
!     js = nint((gslat - rlats)/dl) + 1
!     write(6,*)"is=",is," js=",js

      gridx=-1.e35
      rlatx=-1.e35
      rlonx=-1.e35
      gridn=+1.e35
      rlatn=+1.e35
      rlonn=+1.e35

      do j=1,jl
       do i=1,il
         n=i+(j-1)*il
         grid(i,j)=(ds/em(i,j))/1.e3 ! km
         rlond(i,j)=rlong(n)*rtd
         rlatd(i,j)=rlat (n)*rtd
         if(rlond(i,j).gt.180.)rlond(i,j)=rlond(i,j)-360.
         gridx=max(gridx,grid(i,j))
         gridn=min(gridn,grid(i,j))
         !write(6,*)i,j,rlond(i,j),rlatd(i,j)
         if ( grid(i,j).lt.ds1km ) then    ! use 1km when grid spc < 60 km
           rlatx=max(rlatx,rlatd(i,j))
           rlonx=max(rlonx,rlond(i,j))
           rlatn=min(rlatn,rlatd(i,j))
           rlonn=min(rlonn,rlond(i,j))
         endif ! grid < 60 km
       enddo ! il
       !write(6,*)j,rlond(1,j),rlond(il,j)
       !write(6,*)j,rlatd(1,j),rlatd(il,j)
      enddo ! jl

      imin = nint((rlonn - gslon)/dl) + 1
      imax = nint((rlonx - gslon)/dl) + 1
      jmin = nint((gslat - rlatx)/dl) + 1
      jmax = nint((gslat - rlatn)/dl) + 1

      write(6,*)"grid n,x=",gridn,gridx
      write(6,*)"rlon n,x=",rlonn,rlonx
      write(6,*)"rlat n,x=",rlatn,rlatx
      write(6,*)"i    n,x=",imin,imax
      write(6,*)"j    n,x=",jmin,jmax

      write(6,*)"loop over all valid rows in 1 km data set"
      do jg=jmin,jmax

        glat=gslat-(jg-1)*dl

        if (glat.gt.rlatn .and. glat.lt.rlatx) then
          if(mod(jg,10).eq.0)then
            write(6,*)"jg,glat,rlatn,glat,rlatx=",jg,glat,rlatn,glat,rlatx
          endif! (mod(jg,10).eq.0)then

! read each row of relevant data from 1 km data set (row=jg)
          read(20,rec=jg) (zza(ig),ig=1,mx)

          do ig = 1 , mx
!           if(zza(ig).gt.12) write(6,*)"Big input value ",ig,jg,zza(ig)
          enddo
          !write(6,*) jg,glat,rlatn,rlatx
          !write(6,'(100i2)') (zza(ig),ig=1,mx,500)

          do ig=1,mx

            glon = gslon + (ig-1)*dl
            if (glon.gt.rlonn .and. glon.lt.rlonx) then

              iz(ig)=zza(ig)
              if(iz(ig).ge.19)iz(ig)=0
! compute model grid i,j
              call latltoij(glon,glat,alci,alcj,nface)  ! con-cubic/octagon
              lci = nint(alci)
              lcj = nint(alcj)
! convert to "double" (i,j) notation
              lcj=lcj+nface*il

! check to make sure within grid dimensions
              if(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then
                nsum(lci,lcj,iz(ig))=nsum(lci,lcj,iz(ig))+1
              endif!(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then

              if(lci.eq.id.and.lcj.eq.jd)then
                write(6,'("ig=",i6," jg=",i6," glon=",f7.2," glat=",f6.2
     &           ," i=",i4," j=",i4," iz=",i4," nsum=",i6)')
     &            ig,jg,glon,glat,lci,lcj,iz(ig),nsum(lci,lcj,iz(ig))
!                write(6,*)"ig=",ig," jg=",jg," glon=",glon," glat=",glat
!    &           ," i=",lci," j=",lcj," iz=",iz(ig)," nsum="
!    &           ,nsum(lci,lcj,iz(ig))
              endif

            endif ! glon

          enddo ! ig

        endif ! glat

      enddo ! jg
! finished loop over all valid rows in 1 km data set

      write(6,*)"(nsum(id,jd,nn),nn=1,20)=",(nsum(id,jd,nn),nn=0,20)

! close global datafile
      close(20)

      write(6,*)"old veg for high res panel (j=1:96)"
      do j=96,1,-1
       write(6,'(49i2)')j,(iarr(i,j),i=1,48)
      enddo
      write(6,'(49i2)')j,(i,i=1,48)

      write(6,*)"old lsm for high res panel (j=1:96)"
      do j=96,1,-1
       write(6,'(49i2)')j,(nint(lsm(i+(j-1)*il)),i=1,48)
      enddo
      write(6,'(49i2)')j,(i,i=1,48)

      write(6,*)"nsum"
      do j=96,1,-1
       write(6,'(49i2)')j,(nsum(i,j,12),i=1,48)
      enddo
      write(6,'(49i2)')j,(i,i=1,48)

      write(6,*)"now update vegie by most popular one ifull=",ifull

!-----------------------------------------------------------------------
! loop over all model grid points
      do iq=1,ifull
!-----------------------------------------------------------------------

       if ( lsm(iq).gt..5 .and. iarr(iq,1).eq.0 ) then
         write(6,*)"Input problem iq=",iq," lsm=",lsm(iq)
     &                           ," iarr=",iarr(iq,1)
         stop
       endif ! ( lsm(iq).lt..5 .and. iarr(i,1).eq.0 ) then

       ntypmax=0
       nsummax(iq)=0
       npoints(iq)=0
       nwater=0
       waterfrac(iq) = 0.
       if ( lsm(iq) .gt. .1 ) waterfrac(iq) = 1.

       do ntype=0,20
        !if(ntype.eq.12) write(6,*)iq,ntype,nsum(iq,1,ntype),nsummax(iq)
! sum up number of water input data points
        if(nsum(iq,1,0).gt.0)nwater=nwater+nsum(iq,1,0)
! pick veg type with largest number of input values
        if(nsum(iq,1,ntype).gt.nsummax(iq))then
          nsummax(iq)=nsum(iq,1,ntype)
          ntypmax=ntype
        endif ! (nsum(iq,1,ntype).gt.nsummax(iq))then
! sum up all 1 km data inputs for this grid box
        npoints(iq)=npoints(iq)+nsum(iq,1,ntype)
       enddo ! ntype

! only for points which have some input data!!!
       if ( npoints(iq) .gt. 1 ) then

         waterfrac(iq) = real(nwater)/real(npoints(iq))

! fix land-sea mask
         if ( ntypmax .lt. 1 ) lsm(iq)=0

         if ( ntypmax .gt. 0 ) then
!          print *,'iq,iarr(iq,1),ntypmax,nsummax(iq) ',
!    &            iq,iarr(iq,1),ntypmax,nsummax(iq)

! here we specify new veg type using high-res data
           iarr(iq,1) = ntypmax+31

           if(iarr(iq,1).gt.43)then
             write(6,*)"BIG IARR iq=",iq," iarr=",iarr(iq,1)
           endif
           arr(iq) = iarr(iq,1)
           lsm(iq)=1
         else ! ntypmax
! water point
           iarr(iq,1)=0
           arr(iq) = iarr(iq,1)
           lsm(iq)=0
         endif ! ntypmax

       endif ! ( npoints(iq) .gt. 0 ) then

       if ( lsm(iq).gt..5 .and. iarr(iq,1).eq.0 ) then
         write(6,*)"Output problem iq=",iq," lsm=",lsm(iq)
     &                           ," iarr=",iarr(iq,1)
       endif ! ( lsm(iq).lt..5 .and. iarr(i,1).eq.0 ) then

!-----------------------------------------------------------------------
      enddo   ! iq loop
!-----------------------------------------------------------------------

      write(6,*)"new id,jd,zs,lsm,tsd,veg=",id,jd
     &             ,zs(ijd),lsm(ijd),tsd(id,jd),iarr(id,jd)

!     do j=1,jl
!        write(6,*)"il/2,j,iarr(il/2)=",il/2,j,iarr(il/2,j)
!     enddo ! iq

      write(6,*)"new veg for high res panel (j=1:96)"
      do j=96,1,-1
       write(6,'(49i2)')j,(iarr(i,j),i=1,48)
      enddo
      write(6,'(49i2)')j,(i,i=1,48)

      write(6,*)"new lsm for high res panel (j=1:96)"
      do j=96,1,-1
       write(6,'(49i2)')j,(nint(lsm(i+(j-1)*il)),i=1,48)
      enddo
      write(6,'(49i2)')j,(i,i=1,48)

!     write(6,*)"veg for high res panel (j=49:96)"
!     do j=96,1,-1
!      write(6,'(49i2)')j,(iarr(i,j),i=1,48)
!     enddo
!     write(6,'(49i2)')j,(i,i=1,48)

      write(23,'(i3,i4,2f6.2,f6.3,f8.0,''  veg type'')')
     &            il,jl,rlong0,rlat0,schmidt,ds
      write(23,formout) iarr     ! write vegetation type
      close(23)

      write(24,'(i3,i4,2f6.2,f6.3,f8.0,'' orog-mask-var'')')
     &            il,jl,rlong0,rlat0,schmidt,ds
      ilout=min(il,30)
      write (formout,'(''(''i3''f7.0)'')')ilout   !  i.e. (<il>f7.0)
      write(24,formout) zs

c     write out land/sea mask to formatted file
      write (formout,'(''(''i3''f4.1)'')')ilout   !  i.e. (<il>f4.1)
      write(24,formout) lsm

c     write out std.dev of top. to formatted file
      write (formout,'(''(''i3''f6.0)'')')ilout   !  i.e. (<il>f6.0)
      write(24,formout) tsd
      print *,'zss, rmsk and tsd written for il=',il


      call ncdf_setup("veg.nc",idnc,3,il,jl,1,20020101,0000)
      call nc2out(grid ,il,jl,1,1.,idnc,"grid","grid","km",0.,5000.)
      call nc2out(rlond,il,jl,1,1.,idnc,"lon","lon","degrees",0.,360.)
      call nc2out(rlatd,il,jl,1,1.,idnc,"lat","lat","degrees",-90.,90.)
      call nc2out(arr  ,il,jl,1,1.,idnc,"veg","veg","none",-1.,50.)
      call nc2out(zs   ,il,jl,1,1.,idnc,"zs ","zs ","m",-1.,50000.)
      call nc2out(lsm  ,il,jl,1,1.,idnc,"lsm","lsm","none",-1.,1.)
      call nc2out(waterfrac,il,jl,1,1.,idnc,"wfrc","wfrc","none",0.,1.)
      call nc2out(npoints,il,jl,1,1.,idnc,"npts","npts","none",0.,6000.)
      call nc2out(nsummax,il,jl,1,1.,idnc,"nmax","nmax","none",0.,6000.)
      call ncclos(idnc,ier)

      stop
      end
c=======================================================================
