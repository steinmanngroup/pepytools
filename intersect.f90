subroutine intersect(ids, nids, idsadd, nc1, c1, nc2, c2)
    implicit none
    integer ids
    integer nids
    integer nc1
    integer nc2
    double precision c1
    double precision c2
    dimension ids(nids,2)
    dimension c1(nc1,3)
    dimension c2(nc2,3)
    integer idsadd
!f2py intent(out) :: idsadd
!f2py intent(out) :: ids
!f2py intent(in) :: nids
!f2py intent(in) :: c1
!f2py intent(in) :: nc1
!f2py intent(in) :: c2
!f2py intent(in) :: nc2

    integer n, i, j
    double precision dr, R2, EPS2
    dimension dr(3)

    ! Subroutine declarations
    double precision dot3
    logical inlist

    EPS2 = 1.0D-4
    ids = -1

    ! we take the longer one as the reference
    n = nc1
    if (nc2 .gt. nc1) then
        n = nc2
    endif

    ! linear search first
    idsadd = 1
    do i = 1, n
        dr = c1(i,:) - c2(i,:)
        R2 = dot3(dr, dr)
        if (R2 .lt. EPS2) then
            ids(idsadd,1) = i
            ids(idsadd,2) = i
            idsadd = idsadd +1
        endif
    enddo

!    ! stuff might be moved around a little because
!    ! atoms are deleted
!    do i = 1, n
!        if (inlist(i, nids, ids(:,1))) then
!            goto 90
!        endif
!
!        do j = -1, 1
!            if (j.eq.0) then
!                goto 80
!            endif
!            dr = c1(i,:) - c2(i+j,:)
!            R2 = dot3(dr, dr)
!            if (R2 .lt. EPS2) then
!                ids(idsadd,1) = i
!                ids(idsadd,2) = i+j
!                idsadd = idsadd + 1
!                goto 90
!            endif
!80          continue
!        enddo
!
!90     continue
!    enddo

    ! now we just do a pair-wise comparison for the
    ! rest of the coordinates
    do i = 1, n
        if (inlist(i, nids, ids(:,1))) then
            goto 100
        endif

        do j = 1, n
            if (inlist(j, nids, ids(:,2))) then
                goto 95
            endif

            dr = c1(i,:) - c2(j,:)
            R2 = dot3(dr,dr)
            if (R2 .lt. EPS2) then
                 ids(idsadd,1) = i
                 ids(idsadd,2) = j
                 idsadd = idsadd +1
                 goto 100
            end if
 95         continue
        enddo
 100    continue
    enddo

    ! return python-style indexing
    idsadd = idsadd - 1
    ids = ids - 1

    return
end subroutine

function dot3( U, V )
    implicit none
    double precision U, V
    dimension U(3), V(3)
    double precision dot3
    dot3 = U(1)*V(1) + U(2)*V(2) + U(3)*V(3)
    return
end function

function inlist( ival, n, ilist )
    implicit none
    integer ival, ilist, n, i
    logical inlist
    dimension ilist(n)
    inlist = .false.
    do i=1,n
        if (ilist(i) .eq. ival) then
            inlist = .true.
            exit
        endif
    enddo
    return
end function
