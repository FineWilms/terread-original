module rwork

private
public rmsk,inum,inumx,zss,almsk,tmax,tmin,tsd,rlatd,rlond,grid,id1km,rlonx,rlonn,rlatx,rlatn, &
       rworkalloc,rworkdealloc

real, save :: rlonx,rlonn,rlatx,rlatn
real, dimension(:,:), allocatable, save :: rmsk,zss,almsk,tmax,tmin,tsd,rlatd,rlond,grid
integer, dimension(:,:), allocatable, save :: inum,inumx,id1km

contains

subroutine rworkalloc(il)

implicit none

integer, intent(in) :: il
integer jl

jl=6*il

allocate(rmsk(il,jl),zss(il,jl),almsk(il,jl),tmax(il,jl),tmin(il,jl),tsd(il,jl))
allocate(rlatd(il,jl),rlond(il,jl),grid(il,jl),inum(il,jl),inumx(il,jl),id1km(il,jl))

return
end subroutine rworkalloc

subroutine rworkdealloc

implicit none

deallocate(rmsk,zss,almsk,tmax,tmin,tsd)
deallocate(rlatd,rlond,grid,inum,inumx,id1km)

return
end subroutine rworkdealloc

end module rwork