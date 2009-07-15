c=======================================================================
      subroutine read1km(file,nx,ny,lons,lats,debug,idia,jdia)

      character*11 file
      integer*2 izs(7200)
      real dl
      integer nrecl
      real lons,lats
      parameter ( zmin=-100 )

      logical debug,ok

      include 'newmpar.h'
      include 'rwork.h' ! rmsk,inum,inumx,zss,almsk,tmax,tmin,tsd,rlatd,rlond,grid,id1km,rlonx,rlonn,rlatx,rlatn

      write(6,*)"read1km file,nx,ny=",file,nx,ny

      dl = 1.0/120.0

! on sol-as
!      nrecl = nx*2
!      open(21,file=file,access='direct',form='unformatted',recl=nrecl)

! on jjk-as
      nrecl = nx/2
      open(CONVERT='BIG_ENDIAN',unit=21,file=file,access='direct'
     &            ,form='unformatted',recl=nrecl)

!#######################################################################
      do j=1,ny
!#######################################################################

        aglat=(lats)-real(j-1)*dl
        !write(6,*)j,aglat
        aglata=aglat

        ok = (aglat.gt.rlatn.and.aglat.lt.rlatx)

        !if ( ok ) write(6,*)"lat,latn,latx,lons=",aglat,rlatn,rlatx,lons

         read(21,rec=j) (izs(i),i=1,nx)
         !rmax=-1.e35
         !rmin=+1.e35
         !do i=1,nx
         !  rmax=max(rmax,real(izs(i)))
         !  rmin=min(rmin,real(izs(i)))
         !enddo !i=1,nx
         !if(rmax.gt.-999.)write(6,*)"j,rmax,rmin=",j,rmax,rmin
        izsmaxj=-99999
        izsminj=+99999

       if ( ok ) then

!#######################################################################
         do i=1,nx
!#######################################################################

           n=i+(nj+1-j)*nx
           aglon=lons+real(i-1)*dl
!!! fix up for fiji latitude
   !        if (aglon > 170. . and . aglon < 180. .and. 
   !  .         aglata > -20. .and. aglata < -10. ) then
   !           aglat=aglata + .06
   !        else
              aglat=aglata
   !        endif

!-----------------------------------------------------------------------
          izsmaxj=max(izsmaxj,izs(i))
          izsminj=min(izsminj,izs(i))

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

!          if ( inum(lci,lcj) .eq. 0 ) then ! this would only take first input point

! check to make sure within grid dimensions
           if(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then
             !write(6,*)"i,j,izs(i)=",i,j,izs(i)

             if( izs(i).lt.1 ) then
! ocean point zs < 1
               amask = 0.
               zs=0.
             else  ! zs>=-1000
! land point zs >= 1
               amask = 1.
               zs=real(izs(i))
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
                  write(6,'("lci,lcj,i,izs(i),zs,zss(),almsk(), inum()="
     &           ,4i6,3f10.3,i6)') lci,lcj,i,izs(i),zs,zss(lci,lcj)
     &                       ,almsk(lci,lcj), inum(lci,lcj)
                endif  ! selected points only
             endif  ! debug

           endif ! (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then

!          endif ! ( inum(i,j) .eq. 0 ) then

!-----------------------------------------------------------------------
           endif ! ( ok ) then
!-----------------------------------------------------------------------

!#######################################################################

         enddo ! i

!#######################################################################

        if ( mod(j,100).eq.0 ) then
          write(6,'("i,j,aglon,aglat,lci,lcj=",2i8,2f10.2,4i8)')
     &               i,j,aglon,aglat,lci,lcj,izsmaxj,izsminj
        endif

!#######################################################################

        endif! ( ok ) then

!#######################################################################
      enddo ! j
!#######################################################################

      close(21)

      return
      end ! read1km
c=======================================================================
