!> Calculates the static electric field from a set of nsites charges, dipoles
!> and eventually quadrupoles on npols polarizable points.
subroutine static_field(field,nsites,nexclude,npols,coord,hasalpha,exclusion_list,charges,dipoles)

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

    !$OMP PARALLEL
    !$OMP DO PRIVATE(i,j,ioffset,itensor,iexclusion_list) & 
    !$OMP REDUCTION(+:field)
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
            M0 = dRij * R3i * charges(j)
            M1 = - dipoles(j,:) * R3i
            M1 = M1 + 3.0d0*dRij*R5i*dot3( dRij, dipoles(j,:) )

            field(ioffset:ioffset+2) = field(ioffset:ioffset+2) + M0
            field(ioffset:ioffset+2) = field(ioffset:ioffset+2) + M1
        enddo
        !ioffset = ioffset + 3
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    return
end subroutine static_field

subroutine interaction_matrix( TT, nsites, nexclude, npols, coord, hasalpha, exclusion_list, alphas, damping, damping_factor )
    implicit none
    integer nsites
    integer npols
    integer nexclude
    double precision TT
    double precision coord
    double precision exclusion_list
    integer hasalpha
    double precision alphas
    logical damping
    double precision damping_factor
    dimension TT(3*npols,3*npols)
    dimension coord(nsites,3)
    dimension hasalpha(nsites)
    dimension exclusion_list(nsites,nexclude)
    dimension alphas(npols)
!f2py intent(out) :: TT
!f2py intent(in) :: nsites
!f2py intent(in) :: npols
!f2py intent(in) :: nexclude
!f2py intent(in) :: coordinates
!f2py intent(in) :: hasalpha
!f2py intent(in) :: exclusion_list
!f2py intent(in) :: alphas

    integer i, itensor, iexclusion_list
    integer ii,jj,iii,jjj
    integer j, jtensor
    double precision dRij, Rij, R1i, R2i, R3i, R5i
    dimension dRij(3)
    dimension iexclusion_list(nexclude)

    double precision dot3
    logical inlist

    ! screening variables
    double precision d6i, FE, FT, factor, temp
    parameter( d6i = 1.0d0 / 6.0d0 )

    save dRij, Rij, R1i, R2i, R3i, R5i, ii,jj,iii,jjj,itensor,jtensor,temp,factor,FT,FE
    !$OMP THREADPRIVATE(dRij, Rij, R1i, R2i, R3i, R5i,ii,jj,iii,jjj,itensor,jtensor,temp,factor,FT,FE)

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

            ! Thole damping
            ! JPC A 102 (1998) 2399 and Mol. Sim. 32 (2006) 471
            ! factor = a * u , where a = 2.1304 (default) and u = R / (alpha_i * alpha_j)**(1/6)
            FE = 1.0d0
            FT = 1.0d0
            factor = 0.0d0

            if (damping) then
                ! we will bail in the event that there are zero polarizabilites
                if (alphas(i) .lt. 0.0001 .or. alphas(j) .lt. 0.0001) then
                    cycle
                endif
                temp = (alphas(i)*alphas(j))**d6i

                factor = damping_factor * Rij / temp
                ! the screening is from appendix A in the mol. sim. paper
                FE = 1.0d0 - (1.0d0 + factor + 0.5d0*factor**2)*exp(-factor)
                FT = FE - (d6i * factor**3)*exp(-factor)
                !print '(A,4I6,3F20.4)', "CSS", i,j, itensor, jtensor, factor, FE, FT
            endif


            ! update subblocks of interaction matrix
            do ii=1,3
               do jj=1,3
                   iii = 3*(itensor-1)+ii
                   jjj = 3*(jtensor-1)+jj
                   TT(iii,jjj) = FT * 3.0d0 * dRij(ii) * dRij(jj) * R5i
                   if (ii.eq.jj) then
                       TT(iii,jjj) = TT(iii,jjj) - FE*R3i
                   endif
               enddo
            enddo
        enddo
    enddo
    !$OMP END PARALLEL DO

    return

end subroutine interaction_matrix

function dot3( U, V )
    implicit none
    double precision U, V
    dimension U(3), V(3)
    double precision dot3
    dot3 = U(1)*V(1) + U(2)*V(2) + U(3)*V(3)
    return
end function dot3

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
end function inlist
