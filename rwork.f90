! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
module rwork

private
public rmsk,inum,inumx,zss,almsk,tmax,tmin,tsd,rlatd,rlond,grid,id1km,rlonx,rlonn,rlatx,rlatn, &
       dum,rworkalloc,rworkdealloc

real, save :: rlonx,rlonn,rlatx,rlatn
real, dimension(:,:), allocatable, save :: rmsk,zss,almsk,tmax,tmin,tsd,rlatd,rlond,grid,dum
integer, dimension(:,:), allocatable, save :: inum,inumx,id1km

contains

subroutine rworkalloc(il)

implicit none

integer, intent(in) :: il
integer jl

jl=6*il

allocate(rmsk(il,jl),zss(il,jl),almsk(il,jl),tmax(il,jl),tmin(il,jl),tsd(il,jl))
allocate(rlatd(il,jl),rlond(il,jl),grid(il,jl),inum(il,jl),inumx(il,jl),id1km(il,jl))
allocate(dum(il,jl))

return
end subroutine rworkalloc

subroutine rworkdealloc

implicit none

deallocate(rmsk,zss,almsk,tmax,tmin,tsd)
deallocate(rlatd,rlond,grid,inum,inumx,id1km)
deallocate(dum)

return
end subroutine rworkdealloc

end module rwork