module fea

    !! This module contains procedures that are common to FEM-analysis

    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover

contains

!
!--------------------------------------------------------------------------------------------------
!
    subroutine initial

        !! This subroutine is mainly used to allocate vectors and matrices

        use fedata
        use link1
        use plane42

! Hint for continuum elements:
        integer, dimension(mdim) :: edof
        integer :: e, i, nen , ig


        ! This subroutine computes the number of global equation,
        ! half bandwidth, etc and allocates global arrays.

        ! Calculate number of equations
        neqn = 2*nn
        ! Calculate bw
    bw = 0
    do e = 1, ne

        nen = element(e)%numnode
        do i = 1, nen
            edof(2*i-1) = 2 * element(e)%ix(i) - 1
            edof(2*i)   = 2 * element(e)%ix(i)
        end do
        bw = max(bw, maxval(edof(1:2*nen)) - minval(edof(1:2*nen)) + 1)
    end do
    print *, 'Half-bandwidth = ', bw

        if (.not. banded) then
            allocate (kmat(neqn, neqn))
        else
            allocate(kmat(bw, neqn))
        end if

        allocate (d(neqn),deltaP(neqn))
        allocate (strain(ne, 3), stress(ne, 3))
        allocate (pfinal(neqn), p_ext(neqn), p_int(neqn), res(neqn))
        allocate (p(neqn))


        p = 0.0_wp
        pfinal = 0.0_wp
        p_ext  = 0.0_wp
        p_int  = 0.0_wp
        res    = 0.0_wp
        d = 0.0_wp

        ! Initial stress and strain
        strain = 0
        stress = 0


    ! ==============
    !   PLASTICITY INPUTS
    ! ==============

    mprop%sigy0 = 250e6_wp      
    mprop%Et    = 2e9_wp
    mprop%h     = mprop%young * mprop%Et / (mprop%young - mprop%Et)
    mprop%Ckine = 5e9_wp

    ! Choice of hardening
    ! 1 = isotrope, 2 = kinematic
    mprop%hard_type = 1



        allocate( sigma_gauss(ne, ngauss, nstr) )
        allocate( epsp_gauss(ne, ngauss) )
        allocate( sigy_gauss(ne, ngauss) )
        allocate( dlambda_gauss(ne, ngauss) )
        allocate( eps_gauss(ne, ngauss, nstr) ) 
        allocate(alpha_gauss(ne, ngauss, nstr))
        sigma_gauss   = 0.0_wp
        epsp_gauss    = 0.0_wp
        dlambda_gauss = 0.0_wp
        eps_gauss     = 0.0_wp
        alpha_gauss   = 0.0_wp

        ! Initialiser la contrainte de fluage à sigma_y0
        do e = 1, ne
            do ig = 1, ngauss
                sigy_gauss(e, ig) = mprop( element(e)%mat )%sigy0
            end do
        end do


    end subroutine initial
!
!--------------------------------------------------------------------------------------------------
!
    subroutine displ

        !! This subroutine calculates displacements

        use fedata
        use numeth
        use processor

        integer :: unit_gp
        integer :: e,n, ig
        real(wp), dimension(:), allocatable :: plotval
        real(wp), dimension(ne,3) :: principal_stresses
        real(wp) :: s1, s2, s12, vonm, c, s
        real(wp) :: compliance, stressA, stressB
        real(wp), dimension(:), allocatable :: deltaD
        real(wp), dimension(:), allocatable :: rhs
        integer :: itmax
        real(wp) :: alpha
        character(len=60) :: title_step
        real(wp) :: epsp_global_max


        integer :: it
        real(wp) :: res_norm, res_norm0, tol_rel,tol_abs
        logical :: converged


        ! only add relative tolerance for complex strucutre
        tol_rel = 1.0e-5_wp
        tol_abs = 1.0e-6_wp
        itmax = 1000
        epsp_global_max = 0.0_wp


        ninc=150;
        n1=50
        n2=50
        n3=50

        allocate(deltaD(neqn))
        allocate(rhs(neqn))
        allocate(plotval(ne))
        deltaD = 0.0_wp
        plotval = 0.0_wp
        principal_stresses = 0.0_wp
        rhs    = 0.0_wp

        d = 0.0_wp
        p = 0.0_wp


open(newunit=unit_gp, file='stress_log.dat', status='replace', action='write')
write(unit_gp,'(A)') '# n  e  ig  sigma_x  sigma_y  tau_xy  eps_x  eps_y  gamma_xy  epsp  sigy'



        ! Build load-vector final
        call buildload
        ! Load increment
        deltaP=pfinal/ninc



    ! Incremental loop over load with 3 phases:
    !  1) 0 -> pfinal
    !  2) pfinal -> 0
    !  3) 0 -> 1.5*pfinal
    do n = 1, ninc

        if (n <= n1) then
            ! Phase 1 : 0 -> +0.6  (0 -> +300 MPa)
            alpha = 0.6_wp * real(n, wp) / real(n1, wp)

        else if (n <= n1 + n2) then
            ! Phase 2 : +0.6 -> -1.0  (+300 -> -500 MPa)
            alpha = 0.6_wp + (-1.0_wp - 0.6_wp) * real(n - n1, wp) / real(n2, wp)
            ! = 0.6 - 1.6 * (n-n1)/n2

        else
            ! Phase 3 : -1.0 -> 0  (-500 -> 0)
            alpha = -1.0_wp + (0.0_wp + 1.0_wp) * real(n - n1 - n2, wp) / real(n3, wp)
            ! = -1 + 1 * (n-n1-n2)/n3
        end if

        ! Historique de charge (au cas où tu en aies besoin)
        p_ext = alpha * pfinal
        converged = .false.
        res_norm0 = -1.0_wp

    do it = 1, itmax

        kmat=0.0

        ! Build stiffness matrix
        call buildstiff
        deltaD = 0.0_wp
        call recover(deltaD, .false.)
        res = p_ext - p_int
        p = res
        call enforce                 ! modifie kmat et p pour les appuis
        res = p
        ! 2e) CONVERGENCE ABSOLUE UNIQUEMENT
        res_norm = sqrt( sum( res**2 ) )

        if (it == 1) then
            ! éviter la division par zéro
            res_norm0 = max(res_norm, 1.0e-16_wp)
        end if



        if ( (res_norm / res_norm0) < tol_rel .or. res_norm < tol_abs ) then
            converged = .true.
            exit   ! on sort de la boucle Newton pour ce pas de charge
        end if


! More accurate for single test
!        if ( res_norm < tol_abs ) then
!            converged = .true.
!            exit   ! on sort de la boucle Newton pour ce pas de charge
!        end if

        deltaD = res
        

        if (.not. banded) then
            call factor(kmat)
            call solve(kmat, deltaD)
        else
            call bfactor(kmat)
            call bsolve(kmat, deltaD)
        end if

        ! 2g) Mettre à jour le déplacement total
        d = d + deltaD

        ! 2h) Mettre à jour L'ÉTAT MATÉRIEL (plasticité) avec Δd
        call recover(deltaD, .true.)

        if(it==itmax) then
            print *, 'ERROR in itmax reached without convergence'
        end if
        end do   ! fin boucle Newton

        do e = 1, ne
            if (element(e)%id == 2) then
                plotval(e) = maxval( epsp_gauss(e,1:ngauss) )
            else
                plotval(e) = 0.0_wp
            end if
        end do

        epsp_global_max = max( epsp_global_max, maxval(epsp_gauss) )

            ! Call MATLAB   once every  5 steps
        if (mod(n,5)==0 .or. n==1 .or. n==ninc) then

            write(title_step,'("Plasticity, step ",I0)') n
            call plotmatlabeval(title_step, plotval)
        end if

        do e = 1, ne
            do ig = 1, ngauss
                    write(unit_gp,'(3I5,8ES15.6)') n, e, ig, sigma_gauss(e,ig,1), sigma_gauss(e,ig,2), sigma_gauss(e,ig,3), eps_gauss(e,ig,1),  eps_gauss(e,ig,2),  eps_gauss(e,ig,3), epsp_gauss(e,ig), sigy_gauss(e,ig)
            end do
        end do

        ! file to get the scale of platic strain
        open(33, file="plasticity_scale.m", status="replace")
        write(33,*) "globalmin = 0;"
        write(33,'("globalmax = ",ES15.7,";")') epsp_global_max
        close(33)
        
        ! Output results
        !call output
    end do

        !File matlab
        close(unit_gp)
        call dat_to_matlab('stress_log.dat', 'stress_log.m')

        !call plotmatlabeval('Stresses',plotval)
        !call plotmatlabevec('Principal stresses',principal_stresses(:,1),principal_stresses(:,2),principal_stresses(:,3))

    end subroutine displ
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildload
        !! This subroutine builds the global load vector
        use fedata
        use plane42

        integer :: i, j
! Hint for continuum elements:
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim) :: re
        integer :: nen, e, eface
        real(wp) :: fe, thk

        
        ! Build load vector
        p(1:neqn) = 0


        do i = 1, np
            select case(int(loads(i, 1)))
            case( 1 )
            	! Build nodal load contribution
                pfinal(int((loads(i,2)-1.0_wp)*2 + loads(i,3))) = loads(i,4)

            case( 2 )
            	! Build uniformly distributed surface (pressure) load contribution
            ! Uniform surface pressure load
            e = int(loads(i,2))    ! Element number
            eface = int(loads(i,3))
            fe = loads(i,4)
            thk = mprop(element(e)%mat)%thk

            nen = element(e)%numnode
            do j = 1, nen
                xe(2*j-1) = x(element(e)%ix(j),1)
                xe(2*j  ) = x(element(e)%ix(j),2)
                edof(2*j-1) = 2 * element(e)%ix(j) - 1
                edof(2*j)   = 2 * element(e)%ix(j)
            end do

            call plane42_re(xe, eface, fe, thk, re)

            do j = 1, 2*nen
                pfinal(edof(j)) = pfinal(edof(j)) + re(j)
            end do

            case default
                print *, 'ERROR in fea/buildload'
                print *, 'Load type not known'
                stop
            end select
        end do
    end subroutine buildload
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildstiff

        !! This subroutine builds the global stiffness matrix from
        !! the local element stiffness matrices

        use fedata
        use link1
        use plane42

        integer :: e, i, j, krow, kcol
        integer :: nen
! Hint for system matrix in band form:
!        integer :: irow, icol
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke
! Hint for modal analysis:
!        real(wp), dimension(mdim, mdim) :: me
        real(wp) :: young, area
! Hint for modal analysis and continuum elements:
        real(wp) :: nu, dens, thk, h

        ! Reset stiffness matrix
        if (.not. banded) then
            kmat = 0.0

        else

            kmat = 0.0

        end if

        do e = 1, ne

            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode
            do i = 1, nen
                xe(2*i-1) = x(element(e)%ix(i),1)
                xe(2*i  ) = x(element(e)%ix(i),2)
                 edof(2*i-1) = 2 * element(e)%ix(i) - 1
                 edof(2*i)   = 2 * element(e)%ix(i)
            end do

            ! Gather material properties and find element stiffness matrix
            select case( element(e)%id )
            case( 1 )
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area
                call link1_ke(xe, young, area, ke)

            case( 2 )
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area
                nu = mprop(element(e)%mat)%nu
                thk = mprop(element(e)%mat)%thk
                h = mprop(element(e)%mat)%h

                ! input e for which element plastic - elastic 
                call plane42_ke(e, xe, young, nu, thk, h, ke)
            end select

            ! Assemble into global matrix
            if (.not. banded) then
                do i = 1, 2*nen
                    do j = 1, 2*nen
                        kmat(edof(i), edof(j)) = kmat(edof(i), edof(j)) + ke(i, j)
                    end do
                end do
! Hint: Can you eliminate the loops above by using a different Fortran array syntax?
! Assemble in band form
            else
            do i = 1, 2*nen
                do j = 1, 2*nen
                    ! only fill within the band
                    krow = edof(i)
                    kcol = edof(j)
                    if (kcol <= krow .and. krow - kcol < bw) then
                        kmat(krow - kcol + 1, kcol) = kmat(krow - kcol + 1, kcol) + ke(i,j)
                    end if
                end do
            end do
        end if
    end do
    end subroutine buildstiff
!
!--------------------------------------------------------------------------------------------------
!
    subroutine enforce

        !! This subroutine enforces the support boundary conditions

        use fedata

        integer :: i, idof
        real(wp) :: penal
        integer :: krow, kcol, rowmin

        ! Correct for supports
        if (.not. banded) then
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    p(1:neqn) = p(1:neqn) - kmat(1:neqn, idof) * bound(i, 3)
                    p(idof) = bound(i, 3)
                    kmat(1:neqn, idof) = 0
                    kmat(idof, 1:neqn) = 0
                    kmat(idof, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(idof, idof) = kmat(idof, idof) + penal
                    p(idof) = penal * bound(i, 3)
                end do
            end if
        else
        
        do i = 1, nb

            idof = int(2*(bound(i,1)-1) + bound(i,2))
            ! Zero row (only elements within the band)
            rowmin = max(1, idof-bw+1)
            do kcol = rowmin, idof
                kmat(idof - kcol + 1, kcol) = 0
            end do
        
            ! Zero column (only elements within the band)
            do krow = 1, min(bw, neqn - idof + 1)
                kmat(krow, idof) = 0
            end do
            ! Set diagonal = 1
            kmat(1,idof) = 1
            ! Set RHS
            p(idof) = bound(i,3)
        end do
    end if
    end subroutine enforce
!
!--------------------------------------------------------------------------------------------------
!
subroutine recover(deltaD, update_state)

        !! This subroutine recovers the element stress, element strain,
        !! and nodal reaction forces
            !! This subroutine recovers the element stress, element strain,
        !! and INTERNAL FORCES p_int.
        !! If update_state = .true., it also updates the Gauss-point
        !! plasticity state via gauss_stress_strain.

        use fedata
        use link1
        use plane42

        real(wp), intent(in) :: deltaD(:)
        logical, intent(in)  :: update_state   ! .true. => appel gauss_stress_strain


        integer :: e, i,j, nen
        integer :: ig
        integer :: edof(mdim)
        real(wp), dimension(mdim) :: xe, de, deltaD_e
        real(wp), dimension(mdim, mdim) :: ke
        real(wp) :: young, area
        real(wp):: nu, dens, thk, h
        real(wp), dimension(3) :: estrain, estress
        real(wp) :: sigma_loc(3)
        real(wp) :: bmat(3,8), N(4), detJ
        real(wp) :: xi(2), eta(2), gp_xi, gp_eta, w(2), wt


        ! Reset force vector 
        p_int = 0.0_wp
        xi  = (/ -1.0_wp/sqrt(3.0_wp),  1.0_wp/sqrt(3.0_wp) /)
        eta = xi
        w   = (/ 1.0_wp, 1.0_wp /)              

        do e = 1, ne

            ! Find coordinates etc...
            nen = element(e)%numnode
            do i = 1,nen
                xe(2*i-1) = x(element(e)%ix(i), 1)
                xe(2*i)   = x(element(e)%ix(i), 2)
                edof(2*i-1) = 2 * element(e)%ix(i) - 1
                edof(2*i)   = 2 * element(e)%ix(i)
                de(2*i-1) = d(edof(2*i-1))
                de(2*i)   = d(edof(2*i))
            end do



            ! Find stress and strain
            select case( element(e)%id )
            case( 1 )
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area
                call link1_ke(xe, young, area, ke)
                p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                call link1_ss(xe, de, young, estress, estrain)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain
            case( 2 )

                ! Call subroutine gauss_stress_strain
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area
                nu = mprop(element(e)%mat)%nu
                thk = mprop(element(e)%mat)%thk
                h = mprop(element(e)%mat)%h

                call plane42_ke(e, xe, young, nu, thk, h, ke)

            ! mise à jour de l’état de Gauss SEULEMENT si demandé
            if (update_state) then
                deltaD_e(1:2*nen) = deltaD(edof(1:2*nen))
                call gauss_stress_strain(e, xe, deltaD_e, young, nu, h)
            end if

              ! 2) Forces internes = intégrale B^T * sigma
            do j = 1, 2
                gp_eta = eta(j)
                do i = 1, 2
                    gp_xi = xi(i)
                    wt    = w(i) * w(j)
                    ig    = 2*(j-1) + i

                    call shape(xe, gp_xi, gp_eta, bmat, detJ, N)

                    sigma_loc(1:3) = sigma_gauss(e, ig, 1:3)

                    p_int(edof(1:2*nen)) = p_int(edof(1:2*nen)) + matmul( transpose(bmat), sigma_loc ) * detJ * wt * thk
                end do
            end do



                    ! Element stress = average of Gauss point stresses
                    stress(e,:) = 0.0_wp
                    do ig = 1, 4
                        stress(e,:) = stress(e,:) + sigma_gauss(e,ig,:)
                    end do
                    stress(e,:) = stress(e,:) / 4.0_wp
            end select
        end do
    end subroutine recover

end module fea

