subroutine FillBond(na, iza,conn,bondneed,fakebmx)
implicit none

integer  :: na,i,j,icycle,ineed,bondadd
integer  ::  startp,chainp
integer   ::   iza(na),coord(na),bondneed(na),needlist(na),possible(na),aval(na)
integer   ::   conn(na,na),fakebmx(na,na)
integer   ::   neighbor(4,na)

coord= 0
fakebmx = 0

do i =1 ,na
    if (iza(i) > 18) then
        bondneed(i) =0
        cycle
    end if
    bondneed(i) =min(abs(iza(i)-10), abs(iza(i)-2),abs(iza(i)-18))
end do

!write(*,*) bondneed


do i = 1,na
    do j = i+1, na
        if (conn(i,j) >0 .and. bondneed(j) >0 ) then
            coord(i) =coord(i)+1
            !neighbor(coord(i),i) =j
            coord(j) =coord(j)+1
            !neighbor(coord(j),j) = i
            fakebmx(i,j) = fakebmx(i,j)+1
            fakebmx(j,i) = fakebmx(j,i)+1
        end if
    end do
    bondneed(i) = bondneed(i)-coord(i)
end do


do icycle = 1,10
    !write(*,*) bondneed
    ineed = 0
    possible = 0
    needlist = 0
    aval =1
    !write(*,*) possible
    do i =1,na
        if (bondneed(i) <=0) cycle
        ineed =ineed +1
        needlist(ineed) =i
        do j =i+1,na
            if (conn(i,j) >0 .and. bondneed(j) >0 ) then
                possible(i) =possible(i)+1
                possible(j) = possible(j)+1
            end if
        end do
    end do

    !write(*,*) 'possible'
    !write(*,*) possible
    !write(*,*) 'needlist'
    !write(*,*) needlist


    ! sort by possible neighbor
    do i = 1,ineed-1
        do j =i+1 ,ineed
            if (possible(needlist(i)) > possible(needlist(j))) then
                call change_atom(needlist(i),needlist(j))
            end if
        end do
    end do


    ! sort by chain
    chainp = ineed
    do i = 1,ineed
        if (possible(needlist(i))> 1) then
            chainp =i
            exit
        end if
    end do

    do i = 1,ineed-1
        if (i >= chainp) chainp =chainp +1
        startp = chainp
        do j = startp ,ineed
            if (conn(needlist(i), needlist(j))> 0) then
                call change_atom(needlist(chainp),needlist(j))
                chainp =chainp +1
            end if
        end do
    end do

    !write(*,*) 'needlist'
    !write(*,*) needlist


    bondadd =0
    do i = 1, ineed
        do j = i+1,ineed
            if (conn(needlist(i),needlist(j))>0 .and. aval(needlist(j)) >0 .and. aval(needlist(i)) >0 ) then
                bondneed(needlist(i)) = bondneed(needlist(i)) -1
                bondneed(needlist(j)) = bondneed(needlist(j)) -1
                aval(needlist(i)) = 0
                aval(needlist(j)) = 0
                bondadd = bondadd +1
                fakebmx(needlist(i),needlist(j)) = fakebmx(needlist(i),needlist(j))+1
                fakebmx(needlist(j),needlist(i)) = fakebmx(needlist(j),needlist(i))+1

                exit
            end if
        end do
    end do

    if (ineed <= 1 .or.  bondadd  ==0 )  exit
end do

end subroutine


subroutine FillBondSurface(na, iza,conn,bondneed,surface, fakebmx)
implicit none

integer  :: na,i,j,icycle,ineed,bondadd,itail
integer  ::  startp,chainp
integer   ::  iza(na),coord(na),bondneed(na),needlist(na),possible(na),aval(na)
integer   ::  surface(na)
integer   ::   conn(na,na), fakebmx(na,na)
integer   ::   neighbor(4,na)

coord= 0
do i =1 ,na
    if (iza(i) > 18) then
        bondneed(i) =0
        cycle
    end if
    bondneed(i) =min(abs(iza(i)-10), abs(iza(i)-2),abs(iza(i)-18))
end do

!write(*,*) bondneed


do i = 1,na
    do j = i+1, na
        if (conn(i,j) >0 .and. bondneed(j) >0 ) then
            coord(i) =coord(i)+1
            !neighbor(coord(i),i) =j
            coord(j) =coord(j)+1
            !neighbor(coord(j),j) = i
            fakebmx(i,j) = fakebmx(i,j)+1
            fakebmx(j,i) = fakebmx(j,i)+1

        end if
    end do
    bondneed(i) = bondneed(i)-coord(i)
end do
    

do icycle = 1,9999
    !write(*,*) bondneed
    ineed = 0
    possible = 0
    needlist = 0
    aval =1
    bondadd = 0

    !write(*,*) possible
    do i =1,na 
        if (bondneed(i) <=0) cycle
        ineed =ineed +1
        needlist(ineed) =i
        do j =i+1,na
            if (conn(i,j) >0 .and. bondneed(j) >0 ) then
                possible(i) =possible(i)+1
                possible(j) = possible(j)+1
            end if
        end do
    end do
    
   !write(*,*) 'wwwwwwwwwww'
   ! write(*,*) 'possible'
   ! write(*,*) possible
   ! write(*,*) 'needlist'
   ! write(*,*) needlist

    ! mv surface atom to tail
    itail = ineed
    do i =ineed,1,-1
        if (surface(needlist(i)) == 1) then
            call change_atom(needlist(i),needlist(itail))
            itail = itail -1
            !write(*,*) itail
        end if
    end do


    ! sort by possible neighbor
    do i = 1,ineed-1
        do j =i+1 ,ineed
            if (possible(needlist(i)) > possible(needlist(j)) .and. (surface(needlist(i))+surface(needlist(j)))/=1) then 
                call change_atom(needlist(i),needlist(j))
            end if
        end do
    end do

    
    ! sort by chain
    chainp = ineed
    do i = 1,ineed
        if (possible(needlist(i))> 1) then
            chainp =i
            exit
        end if
    end do

    !write(*,*)  'here?'

    !write (*,*) chainp
    !write (*,*) needlist
    if (chainp /= 1 ) then
        do i = 1, chainp -1
            if (surface(needlist(i)) ==1) cycle
            do j = i+1,ineed
                if (conn(needlist(i),needlist(j))>0 .and. aval(needlist(j)) >0 .and. aval(needlist(i)) >0 ) then
                    bondneed(needlist(i)) = bondneed(needlist(i)) -1
                    bondneed(needlist(j)) = bondneed(needlist(j)) -1
                    aval(needlist(i)) = 0
                    aval(needlist(j)) = 0
                    bondadd = bondadd +1
                    fakebmx(needlist(i),needlist(j)) =fakebmx(needlist(i),needlist(j))+1
                    fakebmx(needlist(j),needlist(i)) =fakebmx(needlist(j),needlist(i))+1
                    exit
                end if
            end do
        end do
        if (bondadd > 0)   cycle
    end if


    do i = 1,ineed-1
        if (i >= chainp) chainp =chainp +1
        startp = chainp
        do j = startp ,ineed
            if (conn(needlist(i), needlist(j))> 0 .and. surface(needlist(j))==0) then
                call change_atom(needlist(chainp),needlist(j))
                chainp =chainp +1
            end if
        end do
    end do

    !write(*,*) 'needlist'
    !write(*,*) needlist

    
    do i = 1, ineed
        if (surface(needlist(i)) ==1) cycle
        do j = i+1,ineed
            if (conn(needlist(i),needlist(j))>0 .and. aval(needlist(j)) >0 .and. aval(needlist(i)) >0 ) then
                bondneed(needlist(i)) = bondneed(needlist(i)) -1
                bondneed(needlist(j)) = bondneed(needlist(j)) -1
                aval(needlist(i)) = 0
                aval(needlist(j)) = 0
                bondadd = bondadd +1
                fakebmx(needlist(i),needlist(j)) = fakebmx(needlist(i),needlist(j))+1
                fakebmx(needlist(j),needlist(i)) = fakebmx(needlist(j),needlist(i))+1
                exit
            end if
        end do
    end do

    !write(*,*)    'test I'
    if (bondadd ==0 ) then

        !write(*,*)    'the part I want'
        if (ineed <= 1 ) then 
            exit
        else
            do i = 1, ineed
                do j = i+1,ineed
                    if (conn(needlist(i),needlist(j))>0 .and. aval(needlist(j)) >0 .and. aval(needlist(i)) >0 ) then
                        bondneed(needlist(i)) = bondneed(needlist(i)) -1
                        bondneed(needlist(j)) = bondneed(needlist(j)) -1
                        aval(needlist(i)) = 0
                        aval(needlist(j)) = 0
                        bondadd = bondadd +1
                        fakebmx(needlist(i),needlist(j)) = fakebmx(needlist(i),needlist(j))+1
                        fakebmx(needlist(j),needlist(i)) = fakebmx(needlist(j),needlist(i))+1
                        exit
                    end if
                end do
            end do
            if ( bondadd  ==0 )  exit
        end if
    end if
end do

!for CO
!do i = 1,na
!    if (bondneed(i) > 0 .and. iza(i)== 6 .and. coord(i) == 1) then
!        do j = 1,na 
!            if (conn(i,j) >0 .and. iza(j)==8 ) then
!                fakebmx(i,j) = fakebmx(i,j)+1
!                fakebmx(j,i) = fakebmx(j,i)+1
!                bondneed(i) =0
!                bondneed(j) =0
!            end if
!        end do
!    end if
!end do
 

end subroutine

subroutine JudgeBond(na,iza,xa,cell,minstr,fakebmx,bondneed)
implicit none

integer ::  na,i
integer ::  iza(na),bondneed(na),surface(na)
integer ::  conn(na,na),fakebmx(na,na)
double precision  :: cell(3,3),xa(3,na)
logical  :: minstr,Lsurface
!write(*,*)  'test start'
minstr = .true.

    call get_conn(na,iza,xa,conn,cell)
!write(*,*)  'connbmx finish'

    call FillBond(na,iza,conn,bondneed,fakebmx)


!write(*,*)  bondneed
do i = 1,na
    if (bondneed(i) /= 0) then
        minstr= .false.
        exit
    end if
end do

end subroutine



subroutine JudgeBondSurface(na,iza,xa,cell,minstr,fakebmx, bondneed,surface)
implicit none

integer ::  na,i
integer ::  iza(na),bondneed(na),surface(na)
integer ::  conn(na,na),fakebmx(na,na)
double precision  :: cell(3,3),xa(3,na)
logical  :: minstr,Lsurface
!write(*,*)  'test start'
minstr = .true.

    call ScreenSurface(na,iza,xa,cell,surface)
    call get_conn(na,iza,xa,conn,cell)
!write(*,*)  'connbmx finish'

    call FillBondSurface(na,iza,conn,bondneed,surface,fakebmx)

!write(*,*) surface
!write(*,*)  bondneed
do i =1,na
    if (bondneed(i)> 0) then
        if(surface(i) == 1)  bondneed(i) =0
!    else if(bondneed(i)== 0) then
!        surface(i)= 0
    end if
end do

!write(*,*)  bondneed
do i = 1,na
    if (bondneed(i) /= 0) then
        minstr= .false.
        exit
    end if
end do

end subroutine

subroutine ScreenSurface(na,iza,xa,cell,surface)
implicit none
double precision :: totbond,dist,co1,co2
integer   ::  na,i,j
integer   ::  iza(na),surface(na)
double precision  :: cell(3,3),xa(3,na)
!logical  :: surface

surface = 0
do i =1 ,na
    if (iza(i)> 18) cycle
    do j =1,na
        if (iza(j) > 18) then
            call get_dist(na,dist,i,j,xa,cell)
            !write(*,*) dist
    
            call species_radius(iza(i),co1)
            call species_radius(iza(j),co2)
            totbond = co1 + co2
    
            !write(*,*) totbond
            if(dist < totbond + 0.4d0) then
               surface(i) =1
               exit
            end if
        end if
    end do
end do
end subroutine


subroutine change_atom(n1,n2)
implicit none
integer :: t, n1, n2

t = n1
n1 = n2
n2 = t

return
end subroutine change_atom

