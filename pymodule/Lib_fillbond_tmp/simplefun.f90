
      subroutine get_dist(na,dist,atom1,atom2,xa,cell,coord)

      implicit none

      integer :: na,atom1,atom2
      double precision,dimension(3,na)  :: xa
!     double precision, allocatable, save :: xa2(:,:),xacna(:,:)
      double precision  xa2(3,2),xacna(3,2)
      double precision cell(3,3),celli(3,3),dist,r,str(3)
      double precision, optional :: coord(3)
      integer                           :: i,j,i1,i2,i3,k
         call reci_latt(cell, celli)
!        do i = 1,na
!           xacna(:,i) = matmul(transpose(celli),xa(:,i))
!           xacna(:,i)= modulo(xacna(:,i)+1000.0d0,1.0d0)
!           do j = 1,3
!             xa(j,i) = cell(j,1) * (xacna(1,i)) + &
!                       cell(j,2) * (xacna(2,i)) + &
!                       cell(j,3) * (xacna(3,i))                   ! xa upper angle cartessian
!           enddo
!        enddo
         xacna(:,1)= matmul(transpose(celli),xa(:,atom1))
         xacna(:,2)= matmul(transpose(celli),xa(:,atom2))
         xacna(:,1)= modulo(xacna(:,1)+1000.0d0,1.0d0)
         xacna(:,2)= modulo(xacna(:,2)+1000.0d0,1.0d0)
         xa2(:,1)= matmul(cell,xacna(:,1))

         dist=1000000d0
         do i1=-1,1
          do i2=-1,1
           do i3=-1,1

            do j = 1,3
              xa2(j,2) = cell(j,1) * (xacna(1,2)+i1) + &
                         cell(j,2) * (xacna(2,2)+i2) + &
                         cell(j,3) * (xacna(3,2)+i3)                   ! xa upper angle cartessian
            enddo
              r= dsqrt(sum((xa2(:,1)-xa2(:,2))**2))
              if(r<dist) then
                dist=r
                !write(*,*)  'reci done'
                !if(present(coord)) coord=xa2(:,2)
                !write(*,*)  'f done'
              endif
            enddo
           enddo
         enddo

      end subroutine

      subroutine get_conn(na,iza,xa,conn,cell)
      implicit none

      integer :: na
      integer,dimension(na)             :: iza
      double precision,dimension(3,na)  :: xa
      integer,dimension(na,na)          :: conn
      integer                           :: i,j
      double precision  :: dist,co1,co2,totbond,cell(3,3)
      
      double precision                  :: minOHDist 
      integer                           :: minOHBond, counter
     
      conn = 0
      do i = 1,na
          do j = 1,i-1
              ! dist = dsqrt(sum((xa(:,i)-xa(:,j))**2))*0.529d0
              call get_dist(na,dist,i,j,xa,cell) ! periodic system
              call fastbond(iza(i),iza(j),totbond)
              if(totbond < 0.05d0) then
                  call species_radius(iza(i),co1)
                  call species_radius(iza(j),co2)
                  totbond = co1 + co2
              end if

! define rules for counting bond
              if(dist < totbond + 0.2d0) then
                 if(iza(i)<18 .and. iza(j) <18 ) then

!             if( (dist < totbond + 0.25d0 &
!                    .and. iza(i)/=1 .and. iza(j)/=1   &
!                    .and. ((iza(i)<18 .and. iza(j) <18) &
!                    .or. (iza(i)==35 .or. iza(j)==35))) & ! count all
!                    non-H bonding, C-C,C-O
!                .or. (dist < totbond + 0.25d0  &
!                    .and. iza(i)==1 .and. iza(j)==1 ) &! count H-H bond
!                .or. (dist < totbond + 0.25d0   &
!                    .and. (iza(i)==6 .or. iza(j)==6) )  ) then ! count
!                    all C-H bond

                 ! write(para%ioout,*) "connection",i,j,dist
                  conn(i,j) = 1
                  conn(j,i) = 1
              endif
             endif
!              write(para%ioout,*) i,j,dist
          end do
      end do

      ! syf
      ! in the above O-H bond length is as high as 1.35 A is allowed
      ! to detect unnormal HO*~H~OCH3* Hydrogen bond.
      ! now remove unnecessary Hydrogen bond (remove the longer one in O~H~O)
      ! eg. only the shortest O-H bond for H atom with more than 1 bond 
      ! and O-H bond longer than 1.147A(normal cut off) for H atom with only 1
      ! bond is preserved
      do i = 1,na
          if (iza(i)==1) then
              minOHBond = -1
              minOHDist = 999999d0
              counter = 0
              do j = 1, na
                  if ( iza(j)==8 .and. conn(i,j)==1 ) then
                      counter = counter + 1
                      call get_dist(na,dist,i,j,xa,cell)                  
                      if (dist < minOHDist) then
                          minOHDist = dist
                          minOHBond = j
                      endif
                  endif
              end do
              
              if (counter > 1) then
                  do j = 1, na
                      if (iza(j)==8 .and. j /= minOHBond) then
                          conn(i,j)=0
                          conn(j,i)=0
                      endif
                  end do
              endif
              if (counter == 1 .and. minOHDist > 1.147) then
                  conn(i,minOHBond) = 0
                  conn(minOHBond,i) = 0
              endif
          endif
      end do 

      !write(*,*)  conn
      !write(*,*) "==========================================="
      end subroutine

      subroutine fastbond(zi,zj,bondlen)

      double precision,dimension(53,53),save :: bondmap
      integer ::  zi,zj,i,j
      double precision :: bondlen
      logical ,save :: init=.true.
      
      if(init) then
        bondmap=0d0
        bondmap(1,1) = 0.850
        bondmap(1,5) = 1.084
        bondmap(1,6) = 1.084
        bondmap(1,7) = 1.001
        bondmap(1,8) = 1.150
        bondmap(1,9) = 0.920
        bondmap(1,14) = 1.48
        bondmap(1,15) = 1.415
        bondmap(1,16) = 1.326
        bondmap(1,17) = 1.28
        bondmap(1,35) = 1.41
        bondmap(1,53) = 1.60
        bondmap(5,5) = 1.70
        bondmap(5,8) = 1.512
        bondmap(6,6) = 1.512
        bondmap(6,7) = 1.439
        bondmap(6,8) = 1.46
        bondmap(6,9) = 1.353
        bondmap(6,14) = 1.86
        bondmap(6,15) = 1.84
        bondmap(6,16) = 1.812
        bondmap(6,17) = 1.781
        bondmap(6,35) = 1.94
        bondmap(6,53) = 2.16
        bondmap(7,7) = 1.283
        bondmap(7,8) = 1.333
        bondmap(7,9) = 1.36
        bondmap(7,14) = 1.74
        bondmap(7,15) = 1.65
        bondmap(7,16) = 1.674
        bondmap(7,17) = 1.75
        bondmap(7,35) = 1.90
        bondmap(7,53) = 2.10
        bondmap(8,8) = 1.45
        bondmap(8,9) = 1.42
        bondmap(8,14) = 1.63
        bondmap(8,15) = 1.66
        bondmap(8,16) = 1.470
        bondmap(8,17) = 1.70
        bondmap(8,22) = 2.1   ! Ti-O
        bondmap(8,35) = 1.85
        bondmap(8,53) = 2.05
        bondmap(9,14) = 1.57
        bondmap(9,15) = 1.54
        bondmap(9,16) = 1.55
        bondmap(14,8) = 2.0
        bondmap(16,35) = 2.24
        bondmap(16,53) = 2.40
        bondmap(17,17) = 1.99
        bondmap(35,35) = 2.28
        bondmap(53,53) = 2.67
        init = .false.
      end if

      if(zi > 53 .or. zj > 53) then
        bondlen = 0.0d0
      else
        if(zi<zj) then
            bondlen = bondmap(zi,zj)
        else
            bondlen = bondmap(zj,zi)
        end if
      end if
      !write(*,*) zi,zj, bondlen
      end subroutine


      SUBROUTINE reci_latt(A,B)

!  CALCULATES RECIPROCAL LATTICE VECTORS. THEIR PRODUCT WITH DIRECT
!  LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1

      integer ::  I
      DOUBLE PRECISION A(3,3),B(3,3), c , ci
      B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
      B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
      B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
      B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
      B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
      B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      C=1.D0
      DO I=1,3
         CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
         B(1,I)=B(1,I)*CI
         B(2,I)=B(2,I)*CI
         B(3,I)=B(3,I)*CI
      ENDDO
      END SUBROUTINE reci_latt

    subroutine species_radius(atom,atomicRadius)
      integer::atom
      double precision atomicRadius

      !atomicRadius = 0.0d0
      select case(atom)
        case(1)
          atomicRadius=0.310d0
        case(2)
          atomicRadius=0.280d0
        case(3)
          atomicRadius=1.280d0
        case(4)
          atomicRadius=0.960d0
        case(5)
          atomicRadius=0.850d0
        case(6)
          atomicRadius=0.760d0
        case(7)
          atomicRadius=0.710d0
        case(8)
          atomicRadius=0.660d0
        case(9)
          atomicRadius=0.570d0
        case(10)
          atomicRadius=0.580d0
        case(11)
          atomicRadius=1.660d0
        case(12)
          atomicRadius=1.410d0
        case(13)
          atomicRadius=1.210d0
        case(14)
          atomicRadius=1.110d0
        case(15)
          atomicRadius=1.070d0
        case(16)
          atomicRadius=1.050d0
        case(17)
          atomicRadius=1.020d0
        case(18)
          atomicRadius=1.060d0
        case(19)
          atomicRadius=2.030d0
        case(20)
          atomicRadius=1.760d0
        case(21)
          atomicRadius=1.700d0
        case(22)
          atomicRadius=1.600d0
        case(23)
          atomicRadius=1.530d0
        case(24)
          atomicRadius=1.390d0
        case(25)
          atomicRadius=1.390d0
        case(26)
          atomicRadius=1.320d0
        case(27)
          atomicRadius=1.260d0
        case(28)
          atomicRadius=1.240d0
        case(29)
          atomicRadius=1.320d0
        case(30)
          atomicRadius=1.220d0
        case(31)
          atomicRadius=1.220d0
        case(32)
          atomicRadius=1.200d0
        case(33)
          atomicRadius=1.190d0
        case(34)
          atomicRadius=1.200d0
        case(35)
          atomicRadius=1.200d0
        case(36)
          atomicRadius=1.160d0
        case(37)
          atomicRadius=2.200d0
        case(38)
          atomicRadius=1.950d0
        case(39)
          atomicRadius=1.900d0
        case(40)
          atomicRadius=1.750d0
        case(41)
          atomicRadius=1.640d0
        case(42)
          atomicRadius=1.540d0
        case(43)
          atomicRadius=1.470d0
        case(44)
          atomicRadius=1.460d0
        case(45)
          atomicRadius=1.420d0
        case(46)
          atomicRadius=1.390d0
        case(47)
          atomicRadius=1.450d0
        case(48)
          atomicRadius=1.440d0
        case(49)
          atomicRadius=1.420d0
        case(50)
          atomicRadius=1.390d0
        case(51)
          atomicRadius=1.390d0
        case(52)
          atomicRadius=1.380d0
        case(53)
          atomicRadius=1.390d0
        case(54)
          atomicRadius=1.400d0
        case(55)
          atomicRadius=2.440d0
        case(56)
          atomicRadius=2.150d0
        case(57)
          atomicRadius=2.070d0
        case(58)
          atomicRadius=2.040d0
        case(59)
          atomicRadius=2.030d0
        case(60)
          atomicRadius=2.010d0
        case(61)
          atomicRadius=1.990d0
        case(62)
          atomicRadius=1.980d0
        case(63)
          atomicRadius=1.980d0
        case(64)
          atomicRadius=1.960d0
        case(65)
          atomicRadius=1.940d0
        case(66)
          atomicRadius=1.920d0
        case(67)
          atomicRadius=1.920d0
        case(68)
          atomicRadius=1.890d0
        case(69)
          atomicRadius=1.900d0
        case(70)
          atomicRadius=1.870d0
        case(71)
          atomicRadius=1.870d0
        case(72)
          atomicRadius=1.750d0
        case(73)
          atomicRadius=1.700d0
        case(74)
          atomicRadius=1.620d0
        case(75)
          atomicRadius=1.510d0
        case(76)
          atomicRadius=1.440d0
        case(77)
          atomicRadius=1.410d0
        case(78)
          atomicRadius=1.360d0
        case(79)
          atomicRadius=1.360d0
        case(80)
          atomicRadius=1.320d0
        case(81)
          atomicRadius=1.450d0
        case(82)
          atomicRadius=1.460d0
        case(83)
          atomicRadius=1.480d0
        case(84)
          atomicRadius=1.400d0
        case(85)
          atomicRadius=1.500d0
        case(86)
          atomicRadius=1.500d0
        case(87)
          atomicRadius=2.600d0
        case(88)
          atomicRadius=2.210d0
        case(89)
          atomicRadius=2.150d0
        case(90)
          atomicRadius=2.060d0
        case(91)
          atomicRadius=2.000d0
        case(92)
          atomicRadius=1.960d0
        case(93)
          atomicRadius=1.900d0
        case(94)
          atomicRadius=1.870d0
        case(95)
          atomicRadius=1.800d0
        case(96)
          atomicRadius=1.690d0
        case default
          atomicRadius=2.0d0
      end select

    end subroutine species_radius


