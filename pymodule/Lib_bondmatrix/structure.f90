

  subroutine hcalbondangle (fi, fj, fk, lat, theta)
  implicit none
  double precision, dimension(3) :: fi, fj, fk
  double precision               :: lat(3,3)
  double precision               :: theta, costheta 
  double precision, dimension(3) :: rij, rik, rjk
  double precision               :: dij, dik, djk

    !write(*,*) fk
  rij = fj - fi; rik = fk - fi
  !write(*,*)  rij,rik
  rij = rij - dnint(rij); rik = rik - dnint(rik)! modified rjk to rik
  rjk = rik - rij
     !write(*,*)  rij,rik,rjk

  dij = dsqrt(sum(matmul (lat, rij)**2.0))
  dik = dsqrt(sum(matmul (lat, rik)**2.0))
  djk = dsqrt(sum(matmul (lat, rjk)**2.0))

  !write(*,*) dij,dik,djk
  costheta = (dij**2.0 + dik**2.0 - djk**2.0) / (2.0 * dij * dik)
  theta    = dacos(costheta) / 3.1415926 * 180.0  

  return
  end subroutine hcalbondangle 
