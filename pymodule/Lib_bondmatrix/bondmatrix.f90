!
!  calculate bond matrix
!
!  17/5/26 created
!
!  Sincerely thanks to Sida Huang for these works.
!  The development is time-consuming and exhausting, 
!      any donation is greatly appreciated (Alipay account: goodspeedsida@gmail.com). 
!
!  Copyright Fudan University 2017
!
!  Sida Huang, Fudan University, March 2017
!
  subroutine hRadicalOrder (na, iza, bmatrix, radical, fa, lat)
  implicit none
  integer                   :: na, maxbond, allbond, countatm
  integer, dimension(na)    :: iza, radical, bondorder, coord
  integer, dimension(na,na) :: bmatrix
  integer                   :: i, iatm, junc(2), nc, jatm
  double precision          :: fa(3,na), lat, theta

  call calcoordnumber (na, bmatrix, coord, bondorder)
  do iatm = 1, na 
    select case(iza(iatm))
      case(1)
        radical(iatm) = abs(coord(iatm) - 1)
      case(6)
        maxbond = maxval(bmatrix(iatm,:))
        if (maxbond == 1) then
          radical(iatm) = abs(coord(iatm) - 4)
        else if (maxbond == 2) then
          if (coord(iatm)==2) then
            nc = 0
            junc =(/0,0/)
            do i = 1, na
              if(bmatrix(iatm,i) == 2) then
                nc = nc + 1
                junc(nc) = i
              end if  
            end do
            if (junc(2) == 0) then
             radical(iatm) = abs(coord(iatm) - 3)
            else
                call hcalbondangle (fa(:,iatm), fa(:,junc(1)), fa(:,junc(2)), lat, theta)
            
                if(theta > 150.0) then
                  radical(iatm) = 0
                else
                  radical(iatm) = abs(coord(iatm) - 3)
                end if
            end if
          else if (coord(iatm)==3) then
            if (bondorder(iatm) == 4) then
              radical(iatm) = 0
            else if (bondorder(iatm) == 5) then
              countatm = 0
              do jatm = 1, na
                if (bmatrix(iatm,jatm)>0 .and. bondorder(jatm)==5) countatm = countatm + 1
              end do
              radical(iatm) = abs(countatm-2)
            else
              radical(iatm) = abs(bondorder(iatm) - 4)
            end if
          else
            radical(iatm) = abs(coord(iatm) - 3)
          end if
        else if (maxbond == 3) then
          radical(iatm) = abs(coord(iatm) - 2)
        else 
            radical(iatm) =4
        end if
      case(7)
        maxbond = maxval(bmatrix(iatm,:))
        if (maxbond == 1) then
          radical(iatm) = abs(coord(iatm) - 3)
        else if (maxbond == 2) then
          radical(iatm) = abs(coord(iatm) - 2)
        else if (maxbond == 3) then
          radical(iatm) = abs(coord(iatm) - 1)
        else
            radical(iatm)= 3
        end if
      case (8)
        maxbond = maxval(bmatrix(iatm,:))
        if (maxbond == 1) then
          radical(iatm) = abs(coord(iatm) - 2)
        else if (maxbond == 2) then
          radical(iatm) = abs(coord(iatm) - 1)
        else
          radical(iatm)= 2
        end if

    end select
  end do

!  write(*,*) coord
!  write(*,*) '-------------------------'
!  write(*,*) radical
!  write(*,*) '-------------------------'
!  write(*,*) iza
!  write(*,*) '-------------------------'
!  write(*,*) bondorder
  

  return 
  end subroutine hRadicalOrder

  subroutine calcoordnumber (na, bmatrix, coord, bondorder)
  implicit none
  integer :: na, i, j
  integer :: coord(na), bmatrix(na,na), bondorder(na)

  do i = 1, na
    coord(i) = 0
    bondorder(i) = sum(bmatrix(i,:))
    do j = 1, na
      if (bmatrix(i,j) > 0) coord(i) = coord(i) + 1
    end do
  end do
 
  return
  end subroutine calcoordnumber

  subroutine hsegmentMol (na, fa, iza, lat, group, bmatrix)
  implicit none
  integer                           :: gid
  integer                           :: na
  integer                           :: iatm
  integer, dimension(na)            :: iza
  integer, dimension(na)            :: group
  integer, dimension(na,na)         :: bmatrix
  double precision, dimension(3,na) :: fa
  double precision, dimension(3,3)  :: lat

  !call hbondmatrix( na, fa, iza, lat, bmatrix)

  gid   = 1
  group = 0

  do iatm = 1, na
    if (sum(bmatrix(:,iatm)) == 0) then
!
!     single atom as a group
!
      group(iatm) = gid
      gid = gid + 1
    else 
!
!     find connecting atom in group
!
      if (group(iatm) == 0) then
        call hconnectSeg (na, bmatrix, iatm, group, gid)
        gid = gid + 1
      end if
    end if

  end do 


  end subroutine hsegmentMol

  subroutine hsegmentMol2 (na, fa, iza, lat, group, bmatrix)
  implicit none
  integer                           :: gid, atmid
  integer                           :: na
  integer                           :: iatm
  integer, dimension(na)            :: iza
  integer, dimension(na)            :: group
  integer, dimension(na,na)         :: bmatrix
  double precision, dimension(3,na) :: fa
  double precision, dimension(3,3)  :: lat

  gid = 1  

  do iatm = 1, na
    call hjudge_atom_inGroup (na, fa, iza, lat, group, bmatrix, gid, iatm)
  end do

  return 
  end subroutine hsegmentMol2

  recursive subroutine hjudge_atom_inGroup (na, fa, iza, lat, group, bmatrix, gid ,atmid)
  implicit none
  integer                           :: gid, atmid
  integer                           :: na
  integer                           :: iatm
  integer, dimension(na)            :: iza
  integer, dimension(na)            :: group
  integer, dimension(na,na)         :: bmatrix
  double precision, dimension(3,na) :: fa
  double precision, dimension(3,3)  :: lat
  
  if ( group(atmid) == 0 ) then
    ! init
    do iatm = 1, na
      if (bmatrix(iatm, atmid) > 0 .and. group(iatm) > 0) then
        group(atmid) = group(iatm); exit 
      end if
    end do

    group(atmid) = gid
    gid = gid + 1
  end if

    do iatm = 1, na
      if (bmatrix(iatm, atmid) > 0 .and. group(iatm)/=group(atmid)) then
        group(iatm) = group(atmid)
        call hjudge_atom_inGroup (na, fa, iza, lat, group, bmatrix, gid, iatm)
      end if
    end do

  end subroutine hjudge_atom_inGroup

  recursive subroutine hconnectSeg (na, bmatrix, startatom, group, gid)
  implicit none
  integer                   :: i
  integer                   :: na
  integer                   :: gid
  integer                   :: startatom
  integer, dimension(na,na) :: bmatrix
  integer, dimension(na)    :: group
  logical                   :: flag
 
  do i = 1, na
    if (bmatrix(startatom,i) > 0 .and. group(i) == 0 ) then
      group(i) = gid
      call hconnectSeg (na, bmatrix, i, group, gid)
    end if
  end do

  return
  end subroutine hconnectSeg


  subroutine hbondmatrix (na, fa, iza, lat, bmatrix)
  implicit none
  integer                            :: i
  integer                            :: j, maxbond, idmax
  integer                            :: na
  integer,          dimension(na)    :: iza, conjugate, ctemp
  double precision                   :: d
  double precision, dimension(3,na)  :: fa
  double precision, dimension(3)     :: delt
  double precision, dimension(3,3)   :: lat
  integer,          dimension(na,na) :: bmatrix

  bmatrix = 0  

  do i = 1, na-1
    do j = i+1, na
      if (iza(i) > 18 .or. iza(j) > 18) then
        bmatrix(i,j) = 0; bmatrix(j,i) = 0
        cycle
      end if

      delt = fa(:,i) - fa(:,j)    
      delt = delt - dnint(delt)
      delt = matmul(lat, delt)
      d = dsqrt(sum(delt**2.0d0))

      if (iza(i) <= iza(j)) then
        call bond_order(iza(i), iza(j), d, bmatrix(i,j))
        bmatrix(j,i) = bmatrix(i,j)
      else
        call bond_order(iza(j), iza(i), d, bmatrix(i,j))
        bmatrix(j,i) = bmatrix(i,j)
      end if

    end do
  end do
!
! judge whether conjugation
!
  !conjugate = 0
  !do i = 1, na-3
  !  maxbond = maxval(bmatrix(i,:))
  !  if (maxbond > 2) then
  !    ctemp = 0
  !    do while (.true.)    
  !      if

  !    end do
  !  end if
  !end do

  return
  end subroutine hbondmatrix

  subroutine bond_order(n1, n2, d, b_o)
!********************
! remark : n1 <= n2 *
!********************
  implicit none
  integer          :: n1, n2, b_o
  double precision :: d
  logical, external :: bsta

  b_o = 0
  select case(n1)
  !
  ! judge hydrogen
  !
    case(1)
      select case(n2)
        case(1)
          if (d < 0.85d0 ) b_o = 1
        case(5)
          if (d < 1.3d0 ) b_o = 1
        case(6)
          if (d < 1.2d0 ) b_o = 1
        case(7)
          if (d < 1.15d0 ) b_o = 1
        case(8)
          if (d < 1.1d0 ) b_o = 1
        case(9)
          if (d < 1.05d0 ) b_o = 1
        case(14)
          if (d < 1.58d0 ) b_o = 1
        case(15)
          if (d < 1.55d0 ) b_o = 1
        case(16)
          if (d < 1.45d0 ) b_o = 1
        case(17)
          if (d < 1.4d0 ) b_o = 1
        case(32)
          if (d < 1.64d0 ) b_o = 1
        case(34)
          if (d < 1.56d0 ) b_o = 1
        case(35)
          if (d < 1.45d0 ) b_o = 1
        case(50)
          if (d < 1.8d0 ) b_o = 1
        case(52)
          if (d < 1.8d0 ) b_o = 1
        case(53)
          if (d < 1.71d0 ) b_o = 1
      end select
  !
  ! judge boron
  !
    case(5)
      select case(n2)
        case(5)
          if (d < 2.0d0 )  b_o = 1
        case(6)
          if (d < 1.66d0 ) b_o = 1
        case(17)
          if (d < 1.85d0 ) b_o = 1
      end select
  !
  ! judge Carbon
  !
    case(6)
      select case(n2)
        case(6)
          if (d < 1.25d0 ) then
            b_o = 3
          else if ( d < 1.43d0) then
            b_o = 2
          else if ( d < 1.64) then
            b_o = 1
          end if
        case(7)
          if (d < 1.26d0 ) then
            b_o = 3
          else if ( d < 1.40d0) then
            b_o = 2
          else if ( d < 1.6) then
            b_o = 1
          end if
        case(8)
          if (d<1.15  ) then
            b_o = 3
          else if ( d < 1.3d0) then
            b_o = 2
          else if ( d < 1.55d0) then
            b_o = 1
          end if
        case(9)
          if (d < 1.48d0 ) b_o = 1
        case(14)
          if (d < 1.96d0 ) b_o = 1
        case(15)
          if (d < 1.97d0 ) b_o = 1
        case(16)
          if (d < 1.7d0 ) then
            b_o = 2
          else if ( d < 1.92d0) then
            b_o = 1
          end if
        case(17)
          if (d < 2.1d0 ) b_o = 1
        case(32)
          if (d < 2.05d0 ) b_o = 1
        case(35)
          if (d < 2.04d0 ) b_o = 1
        case(53)
          if (d < 2.25d0 ) b_o = 1
      end select
  !
  ! judge Nitrogen
  !
    case(7)
      select case(n2)
        case(7)
          if (d < 1.2d0 ) then
            b_o = 3
          else if ( d < 1.35d0) then
            b_o = 2
          else if ( d < 1.55) then
            b_o = 1
          end if
        case(8)
          if (d < 1.15d0 ) then
            b_o = 3
          else if ( d < 1.3d0) then
            b_o = 2
          else if ( d < 1.56d0) then
            b_o = 1
          end if
      end select
  !
  !  judge Oxygen
  !
    case(8)
      select case(n2)
        case(8)
          if (d < 1.3d0 ) then
            b_o = 2
          else if ( d < 1.58d0) then
            b_o = 1
          end if
        case(14)
          if (d < 1.73d0 ) b_o = 1
        case(15)
          if (d < 1.48d0 ) then
            b_o = 2
          else if ( d < 1.73d0) then
            b_o = 1
          end if
        case(16)
          if (d < 1.53d0 ) b_o = 1
      end select

  end select

  return
  end subroutine bond_order

  function bsta(d, sd)
  implicit none
  logical :: bsta
  double precision :: d, sd
  double precision :: tol = 0.1d0

  bsta = .false.
  if(d>(sd-tol) .and. d<(sd+tol)) bsta = .true.

  return
  end function bsta
