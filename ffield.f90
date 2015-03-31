!> Calculates the static electric field from a set of nsites charges, dipoles
!> and eventually quadrupoles on npols polarizable points.
subroutine fstatic_field(field,nsites,nexclude,npols,coord,hasalpha,exclusion_list,charges,dipoles)

    use omp_lib

    implicit none
    integer nsites
    integer npols
    integer nexclude
    integer hasalpha
    integer exclusion_list
    double precision field
    double precision coord
    double precision charges
    double precision dipoles
    dimension coord(nsites,3)
    dimension hasalpha(nsites)
    dimension field(3*npols)
    dimension exclusion_list(nsites,nexclude)
    dimension charges(nsites)
    dimension dipoles(nsites,3)
!f2py intent(in) :: nsites
!f2py intent(in) :: npols
!f2py intent(in) :: nexclude
!f2py intent(in) :: coord
!f2py intent(in) :: charges
!f2py intent(in) :: dipoles
!f2py intent(in) :: hasalpha
!f2py intent(in) :: exclusion_list
!f2py intent(out) :: field

    ! internal
    integer i, itensor, iexclusion_list
    integer j
    double precision dRij, Rij, uRij
    double precision R1i, R2i, R3i, R5i
    double precision M0, M1
    dimension iexclusion_list(nexclude)
    dimension dRij(3), uRij(3), M0(3), M1(3)
    integer ioffset

    ! Subroutine declarations
    double precision dot3
    logical inlist

    save dRij, Rij, uRij, R1i, R2i, R3i, R5i, M0, M1
    !$OMP THREADPRIVATE(dRij, Rij, uRij, R1i, R2i, R3i, R5i, M0, M1)


    field = 0.0d0
    ! increase by 1 to get fortran style indices
    hasalpha = hasalpha+1

    !$OMP PARALLEL DO PRIVATE(i,j,ioffset,itensor,iexclusion_list) REDUCTION(+:field)
    do i=1,nsites
        itensor = hasalpha(i)
        ioffset = (itensor-1)*3+1
        if (itensor .eq. 0) then
            cycle
        endif
        ! fortran style indices
        iexclusion_list = exclusion_list(i,:) + 1

        do j=1,nsites
            ! never interact with yourself
            if( i .eq. j ) then
                cycle
            endif

            ! are we to skip a point?
            if(charges(j).eq.0) then
                if(abs(sum(dipoles(j,:))).eq.0.0)then
                    cycle
                endif
            endif

            ! never interact with yourself if you were
            ! polarized during generation
            if (inlist(j, nexclude, iexclusion_list)) then
                cycle
            endif
            dRij = coord(j,:) - coord(i,:)
            Rij = sqrt(dot3( dRij, dRij ))
            R1i = 1.0d0 / Rij
            uRij = dRij * R1i
            R2i = R1i * R1i
            R3i = R2i * R1i
            R5i = R3i * R2i

            ! M0: charge, M1: dipole
            M0 = uRij * R2i * charges(j)
            M1 = - dipoles(j,:) * R3i
            M1 = M1 + 3.0d0*dRij*R5i*dot3( dRij, dipoles(j,:) )

            field(ioffset:ioffset+2) = field(ioffset:ioffset+2) + M0
            field(ioffset:ioffset+2) = field(ioffset:ioffset+2) + M1
        enddo
        !ioffset = ioffset + 3
    enddo
    !$OMP END PARALLEL DO
    return
end subroutine fstatic_field

subroutine generate_tt( TT, nsites, nexclude, npols, coord, hasalpha, exclusion_list )
    implicit none
    integer nsites
    integer npols
    integer nexclude
    double precision TT
    double precision coord
    double precision exclusion_list
    integer hasalpha
    dimension TT(3*npols,3*npols)
    dimension coord(nsites,3)
    dimension hasalpha(nsites)
    dimension exclusion_list(nsites,nexclude)
!f2py intent(out) :: TT
!f2py intent(in) :: nsites
!f2py intent(in) :: npols
!f2py intent(in) :: nexclude
!f2py intent(in) :: coordinates
!f2py intent(in) :: hasalpha
!f2py intent(in) :: exclusion_list

    integer i, itensor, iexclusion_list
    integer ii,jj,iii,jjj
    integer j, jtensor
    double precision dRij, Rij, R1i, R2i, R3i, R5i
    dimension dRij(3)
    dimension iexclusion_list(nexclude)

    double precision dot3
    logical inlist

    save dRij, Rij, R1i, R2i, R3i, R5i, ii,jj,iii,jjj,itensor,jtensor
    !$OMP THREADPRIVATE(dRij, Rij, R1i, R2i, R3i, R5i,ii,jj,iii,jjj,itensor,jtensor)

    ! increase by 1 to get fortran style indices
    hasalpha = hasalpha+1

    TT = 0.0d0

    ! so a REDUCTION is apparently not needed here because
    ! the matrix TT is only ever updated in discreet sub blocks.
    ! trying to do REDUCTION crashes the program, perhaps because
    ! of a memory issue
    !$OMP PARALLEL DO PRIVATE (i,j,iexclusion_list)
    do i = 1, nsites
        itensor = hasalpha(i)

        ! fortran style indices
        iexclusion_list = exclusion_list(i,:) + 1
        if (itensor .eq. 0) then
            cycle
        endif

        do j = 1, nsites
            if (j .eq. i) then
                cycle
            endif

            if (inlist(j, nexclude, iexclusion_list)) then
                cycle
            endif

            jtensor = hasalpha(j)
            if (jtensor .eq. 0) then
                cycle
            endif

            dRij = coord(j,:) - coord(i,:)
            Rij = sqrt(dot3( dRij, dRij ))
            R1i = 1.0d0 / Rij
            R2i = R1i * R1i
            R3i = R2i * R1i
            R5i = R3i * R2i

            ! update subblocks of interaction matrix
            do ii=1,3
               do jj=1,3
                   iii = 3*(itensor-1)+ii
                   jjj = 3*(jtensor-1)+jj
                   TT(iii,jjj) = 3.0d0 * dRij(ii) * dRij(jj) * R5i
                   if (ii.eq.jj) then
                       TT(iii,jjj) = TT(iii,jjj) - R3i
                   endif
               enddo
            enddo
        enddo
    enddo
    !$OMP END PARALLEL DO

    return

end subroutine generate_tt

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
