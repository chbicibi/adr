program pareto
  implicit none

  integer        :: nobj, nind, i, j, k, l, idum, ndv, np, ngene, ig, igene(500000)
  real(8)        :: obj(3, 500000), obj_arc(3, 500000), a
  real(8)        :: val(0:5, 500000), dum, hvol, max_spread, gdist, spacing
  real(8)        :: obj_max(3), obj_min(3)
  integer        :: np2(1000), pmemc(500000, 5)
  character(128) :: filename

  open(10, file = 'ga-result.txt')
    read(10,*) ndv, nobj, nind

    do i = 1, nind
      read(10,*) igene(i), (obj(k, i), k = 1, nobj), val(0, i), (pmemc(i, k), k = 1, ndv), (val(k, i), k = 1, ndv)
      print *, 'igene', igene(i), i, (pmemc(i, k), k = 1, ndv)
      obj(2, i) = -obj(2, i)
    end do
  close(10)

  ngene = igene(nind)

! ----------------------------------------------------------
! -------calculation of pareto------------------------------

  open(12, file = 'pareto-param.dat')
    if (nobj == 2) then
!................ Final ranking result
!        np=1
!        do l=1,nobj
!         obj_max(l)=-1000000.0
!         obj_min(l)=1000000.0
!        end do
!!
!        do i=1,nind     !shikakuno motos
!        k=0
!          do j=1,nind     !taishou no bangou
!             if (obj(1,j) > obj(1,i)-1d20 .and. obj(1,j) < obj(1,i)
!     &.and.  obj(2,j) > obj(2,i)-1d20 .and. obj(2,j) < obj(2,i)) then
!             k=k+1
!             end if
!          end do
!         if (k = 0) then
!          write(12,*)np,igene(i),obj(1,i),obj(2,i),(pmemc(i,j),j=1,ndv)
!     &         ,(val(i,j),j=0,ndv)
!             np=np+1
!          do l=1,nobj
!           if (obj_max(l) < obj(l,i))obj_max(l)=obj(l,i)
!           if (obj_min(l) > obj(l,i))obj_min(l)=obj(l,i)
!          end do
!         end if
!         end do
!
!        do l=1,nobj
!         print*, "objmax,min",obj_max(l),obj_min(l)
!        end do
!.................. Rank history

      open(111, file = 'metrics.txt')
        write(111, *) 'Generation  Maximum_spread  Spacing  Hypervolume'

        do ig = 1, ngene
          write(filename, '(a,i4.4,a)') 'pareto-hist', ig, '.csv'

          open(11, file = filename)
            np = 1
            do i = 1, nind ! history
              if (igene(i) <= ig) then
                k = 0

                do j = 1, nind ! taishou no bangou
                  if (igene(j) <= ig) then
                    if (obj(1, j) >= obj(1, i) - 1d20 .and. obj(1, j) <= obj(1, i) .and. &
                        obj(2, j) >= obj(2, i) - 1d20 .and. obj(2, j) <= obj(2, i)) then
                      k = k + 1
                    end if
                  end if
                end do

                if (k == 0) then
                ! if (k <= 1) then
                  write(11, *) ig, np, igene(i), obj(1, i), -1. * obj(2, i), (pmemc(i, j), j = 1, ndv), (val(j, i), j = 0, ndv)

                  do l = 1, nobj ! ...For metrics' calculation
                    ! np2(ig) = np
                    ! obj_arc2(l,np,ig)=obj(l,i)
                    ! obj_arc(l,np)=(obj(l,i)-obj_min(l))/(obj_max(l)-obj_min(l))-1.0
                    obj_arc(l, np) = obj(l, i)
                  end do
                  np2(ig) = np

                  np = np + 1
                end if
              end if
            end do
          close(11)
          np = np - 1

          ! call mspread(nobj, np, obj_arc, max_spread)
          call HV(nobj, np, obj_arc, hvol)
          ! call space(nobj, np, obj_arc, spacing)
          write(111, *) ig, max_spread, spacing, hvol
        end do
      close(111)

      ! call generalDist(nobj, ngene, np2, obj_arc2)

    elseif (nobj == 3) then
      np = 1
      do i = 1, nind ! shikakuno moto
        k = 0
        do j = 1, nind ! taishou no bangou
          if (obj(1, j) > obj(1, i) - 1d20 .and. obj(1, j) < obj(1, i) .and. &
              obj(2, j) > obj(2, i) - 1d20 .and. obj(2, j) < obj(2, i) .and. &
              obj(3, j) > obj(3, i) - 1d20 .and. obj(3, j) < obj(3, i)) then
            k = k + 1
          end if
        end do
        if (k == 0) then
          write(11, *) np, igene(i), obj(1, i), obj(2, i), obj(3, i)
          write(12, *) np, igene(i), (pmemc(i, j), j = 1, ndv), obj(1, i), obj(2, i), obj(3, i)
          np = np + 1
        end if
      end do
    end if
  close(12)
  write(*, *)'wrote pareto.dat'
  stop
end program pareto

subroutine HV(nobj, np, obj_arc, hvol)
  implicit none

  integer :: nobj, np
  integer :: i, j, k, l, idum
  real(8) :: obj_arc(3, 100000), vol(100000)
  real(8) :: obj_max(3), hvol
  real(8) :: objlength, objlength_temp
  integer :: i_volmax

  do i = 1, nobj
    obj_max(i) = 0.0
  end do

  do j = 1, np
    vol(j) = 1.0
    do i = 1, nobj
      objlength = obj_max(i) - obj_arc(i,j)
      if (i /= nobj) then
        do k = 1, np
          if (obj_arc(i, j) < obj_arc(i, k)) then
            objlength_temp = obj_arc(i, k) - obj_arc(i, j)
            if (objlength_temp < objlength) objlength = objlength_temp
          end if
        end do
      end if
      vol(j) = vol(j) * objlength
    end do
  end do

  hvol = 0.0
  do j = 1, np
    hvol = hvol + vol(j)
  end do
  return
end subroutine HV
