!
!  process array 
!
!  17/6/26 created
!
!  Sincerely thanks to Sida Huang for these works.
!  The development is time-consuming and exhausting, 
!      any donation is greatly appreciated (Alipay account: goodspeedsida@gmail.com). 
!
!  Copyright Fudan University 2017
!
!  Sida Huang, Fudan University, March 2017
!
  subroutine hcount_int(n, arr, item, times)
  implicit none
  integer :: n, times, item, i
  integer, dimension(n) :: arr

  times = 0
  do i = 1, n
    if (arr(i) == item) times = times + 1
  end do
 
  return
  end subroutine hcount_int
