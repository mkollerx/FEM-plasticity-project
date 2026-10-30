module plane42

    !! This module contains subroutines specific to the plane42 element
    !!
    !! The plane42 element has 4 nodes. Each node has 2 degrees-of-freedom,
    !! namely, displacement along the \(x\)- and \(y\)-coordinate directions.
    !!
    !! The nodes are numbered counter-clockwise as indicated below. The
    !! figure also shows how the element edges are labelled. For example,
    !! edge 1 is the edge between element node 1 and 2.
    !!
    !!       N4    E3    N3
    !!          o------o
    !!          |      |
    !!       E4 |      | E2
    !!          |      |
    !!          o------o
    !!       N1    E1    N2
    !!
    !!
    !! `N1` = element node 1, `N2` = element node 2, etc
    !! `E1` = element face 1, `E2` = element face 2, etc

    use types
    implicit none
    save

    private
    public :: plane42_ke, plane42_re, plane42_ss, gauss_stress_strain, shape


contains

subroutine plane42_ke(e, xe, young, nu, thk, h, ke)

        use types
        use fedata
        implicit none

        !! This subroutine constructs the stiffness matrix for
        !! a rectangular 4-noded quad element.
        
        integer,  intent(in)  :: e  
            !! Index of element
        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), intent(in) :: thk
                    !! Poisson's Ratio for this element
        real(wp), intent(in) :: h
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe

            !! Nodal coordinates of this element in undeformed configuration
            !!
            !! * `xe(1:2)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(3:4)` = \((x,y)\)-coordinates of element node 2
            !! * `xe(5:6)` = \((x,y)\)-coordinates of element node 3
            !! * `xe(7:8)` = \((x,y)\)-coordinates of element node 4
            !!
            !! See also [[plane42rect]]
        real(wp), dimension(:,:), intent(out) :: ke
            !! Stiffness matrix

        real(wp) :: cmat(3,3), fact

        real(wp) :: cmat_ep(3,3)
        real(wp) :: sigma_e, F 
        real(wp) :: xi(2), eta(2), w(2)
        real(wp) :: gp_xi, gp_eta, wt
        real(wp) :: bmat(3,8), detJ, N(4)
        real(wp) :: sigma_old(3), sigy_old,vm_arg
        real(wp) :: dFdS(3), CdFdS(3), denom
        integer  :: i, j, ig, a, b
        integer  :: hardtype, mat_id
        real(wp) :: sigy0, Ckine          ! <--- MANQUAIT
        real(wp) :: alpha_old(3)          ! <--- MANQUAIT
        real(wp) :: sigma_rel(3)          ! <--- MANQUAIT

        ! build constitutive matrix (plane stress)
        cmat = 0.0
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2

        ke = 0.0

        ! Gauss points
        xi  = (/ -1.0/sqrt(3.0),  1.0/sqrt(3.0) /)
        eta = xi 
        w   = (/ 1.0, 1.0 /)

        ! Sigma e (von Mises)
        sigma_e=0.0
        ! Define the yield function
        F=0.0

        bmat=0.0
    ! Loop on gauss points

    do j = 1, 2
        gp_eta = eta(j)
        do i = 1, 2
            gp_xi = xi(i)
            wt = w(i)*w(j)
            
            ! global gauss element index
            ig = 2*(j-1) + i
            call shape(xe, gp_xi, gp_eta, bmat, detJ, N)

            ! old stress state
            sigma_old(1:3) = sigma_gauss(e, ig, 1:3)
            sigy_old       = sigy_gauss(e, ig)

            mat_id   = element(e)%mat
            hardtype = mprop(mat_id)%hard_type

            select case (hardtype)

            case(1) ! ISOTROPE
            ! von Mises from old stress state
                sigma_e = sqrt( sigma_old(1)**2 + sigma_old(2)**2 - sigma_old(1)*sigma_old(2) + 3.0_wp*sigma_old(3)**2 )

                F = sigma_e - sigy_old
            if( F < 0) then
                ! Normal element stiffness matrix formulation
                ke = ke + matmul( transpose(bmat), matmul(cmat, bmat) ) * detJ * wt * thk
            else

                dFdS(1) = (2.0_wp*sigma_old(1) - sigma_old(2)) / (2.0_wp*sigma_e)
                dFdS(2) = (2.0_wp*sigma_old(2) - sigma_old(1)) / (2.0_wp*sigma_e)
                dFdS(3) =  6.0_wp*sigma_old(3)                 / (2.0_wp*sigma_e)


                CdFdS = matmul(cmat, dFdS)
                denom = dot_product(dFdS, CdFdS) + h

                do a = 1, 3
                    do b = 1, 3
                        cmat_ep(a,b) = cmat(a,b) - CdFdS(a)*CdFdS(b)/denom
                    end do
                end do

                ! Element stiffness matrix formulation using Cep
                ke = ke + matmul( transpose(bmat), matmul(cmat_ep, bmat) ) * detJ * wt * thk
            end if 

            case (2)
                sigy0 = mprop(mat_id)%sigy0
                Ckine = mprop(mat_id)%Ckine

                alpha_old(1:3) = alpha_gauss(e, ig, 1:3)
                sigma_rel(1:3) = sigma_old(1:3) - alpha_old(1:3)

                vm_arg = sigma_rel(1)**2 + sigma_rel(2)**2 - sigma_rel(1)*sigma_rel(2) + 3.0_wp*sigma_rel(3)**2
                if (vm_arg < 0.0_wp) vm_arg = 0.0_wp
                sigma_e = sqrt(vm_arg)

                F = sigma_e - sigy0

                if (F < 0.0_wp .or. sigma_e < 1.0e-12_wp) then
                    ! élastique
                    ke = ke + matmul( transpose(bmat), matmul(cmat, bmat) ) * detJ * wt * thk
                else
                    dFdS(1) = (2.0_wp*sigma_rel(1) - sigma_rel(2)) / (2.0_wp*sigma_e)
                    dFdS(2) = (2.0_wp*sigma_rel(2) - sigma_rel(1)) / (2.0_wp*sigma_e)
                    dFdS(3) =  6.0_wp*sigma_rel(3)                 / (2.0_wp*sigma_e)

                    CdFdS = matmul(cmat, dFdS)
                    denom = dot_product(dFdS, CdFdS) + (2.0_wp/3.0_wp) * Ckine
                    if (abs(denom) < 1.0e-16_wp) denom = sign(1.0e-16_wp, denom)

                    do a = 1, 3
                        do b = 1, 3
                            cmat_ep(a,b) = cmat(a,b) - CdFdS(a)*CdFdS(b)/denom
                        end do
                    end do

                    ke = ke + matmul( transpose(bmat), matmul(cmat_ep, bmat) ) * detJ * wt * thk
                end if

                case default
                ! fallback élastique
                ke = ke + matmul( transpose(bmat), matmul(cmat, bmat) ) * detJ * wt * thk

            end select

        end do
    end do
    end subroutine plane42_ke
!
!--------------------------------------------------------------------------------------------------
!
subroutine plane42_re(xe, eface, pres, thk, re)
    use types
    implicit none
    real(wp), intent(in) :: xe(:), pres, thk
    integer,  intent(in) :: eface
    real(wp), intent(out) :: re(8)

    integer :: i, a
    real(wp) :: gp(2), w(2), xi, eta, detJ
    real(wp) :: bmat(3,8), N(4)
    real(wp) :: x1, y1, x2, y2
    real(wp) :: tvec(2), nvec(2), Le, tmag
    real(wp) :: Nmat(2,8)

    ! 2-point Gauss sur l’arête
    gp = (/ -1.0_wp/sqrt(3.0_wp),  1.0_wp/sqrt(3.0_wp) /)
    w  = (/ 1.0_wp, 1.0_wp /)
    re = 0.0_wp

    ! --- 1) Géométrie de l’arête (fixe pour tous les points de Gauss) ---
    select case (eface)
    case (1)    ! entre noeud 1 et 2
        x1 = xe(1);  y1 = xe(2)
        x2 = xe(3);  y2 = xe(4)
    case (2)    ! 2–3
        x1 = xe(3);  y1 = xe(4)
        x2 = xe(5);  y2 = xe(6)
    case (3)    ! 3–4
        x1 = xe(5);  y1 = xe(6)
        x2 = xe(7);  y2 = xe(8)
    case (4)    ! 4–1
        x1 = xe(7);  y1 = xe(8)
        x2 = xe(1);  y2 = xe(2)
    end select

    ! vecteur tangent et longueur de l’arête
    tvec(1) = x2 - x1
    tvec(2) = y2 - y1
    Le = sqrt( tvec(1)**2 + tvec(2)**2 )

    ! normale unitaire (en supposant numérotation anti-horaire des noeuds)
    nvec(1) = -tvec(2) / Le
    nvec(2) =  tvec(1) / Le

    ! --- 2) Intégration numérique le long de l’arête ---
    do i = 1, 2

        ! coord. naturelles selon la face
        select case (eface)
        case (1)
            xi  = gp(i)
            eta = -1.0_wp
        case (2)
            xi  =  1.0_wp
            eta = gp(i)
        case (3)
            xi  = gp(i)
            eta =  1.0_wp
        case (4)
            xi  = -1.0_wp
            eta = gp(i)
        end select

        ! récupère N(ξ,η) via ta routine shape
        call shape(xe, xi, eta, bmat, detJ, N)

        ! matrice [N] en forme 2×8
        Nmat = 0.0_wp
        do a = 1, 4
            Nmat(1, 2*a-1) = N(a)   ! composante x
            Nmat(2, 2*a  ) = N(a)   ! composante y
        end do

        ! Jacobien ligne : |∂x/∂s| = Le/2 pour 2-point Gauss sur [-1,1]
        tmag = Le / 2.0_wp

        ! contribution au vecteur de charge : [N]^T * p * n * |J| * w * thk
        re = re + matmul( transpose(Nmat), pres * nvec ) * tmag * w(i) * thk

    end do

end subroutine plane42_re

!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_ss(xe, de, young, nu, estress, estrain)

        !! This subrotuine computes the element stress and strain (The location inside the element
        !! where stress and and strain is evaluated, is defined inside the subroutine).

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), dimension(:), intent(in)  :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), dimension(:), intent(in)  :: de
            !! Displacement field
            !!
            !! * `de(1:2)` = displacement of element node 1 along \(x\)- and \(y\)-axis, respectively
            !! * `de(3:4)` = displacement of element node 2 along \(x\)- and \(y\)-axis, respectively
            !! * etc...
        real(wp), dimension(:), intent(out) :: estress
            !! Stress at a point inside the element
            !!
            !! * `estress(1)` = \(\sigma_{11}\)
            !! * `estress(2)` = \(\sigma_{22}\)
            !! * `estress(3)` = \(\sigma_{12}\)
        real(wp), dimension(:), intent(out) :: estrain
            !! Strain at a point inside the element
            !!
            !! * `estrain(1)` = \(\epsilon_{11}\)
            !! * `estrain(2)` = \(\epsilon_{22}\)
            !! * `estrain(3)` = \(\epsilon_{12}\)
    real(wp) :: cmat(3,3), fact
    real(wp) :: bmat(3,8), N(4), detJ
    real(wp), parameter :: xi0 = 0.0_wp, eta0 = 0.0_wp

        cmat=0.0
        ! build constitutive matrix (plane stress)
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2

        ! Build strain-displacement matrix
        ! Build strain-displacement matrix (B-matrix) for a 4-noded rectangle
        bmat = 0.0

        ! Build B at the center using the same shape() as ke
        call shape(xe, xi0, eta0, bmat, detJ, N)

            ! Strain and stress at center
        estrain = matmul(bmat, de)     
        estress = matmul(cmat, estrain)

        ! Find eps 12 
        estrain(3) = estrain(3) / 2.0 

    end subroutine plane42_ss

subroutine shape(xe, xi, eta, bmat, detJ, N)

    use types
    implicit none

    real(wp), intent(in)  :: xe(:)         ! [x1 y1 x2 y2 x3 y3 x4 y4]
    real(wp), intent(in)  :: xi, eta       ! natural coords (ξ, η)
    real(wp), intent(out) :: bmat(3,8)        ! strain-displacement at (ξ,η)
    real(wp), intent(out) :: detJ          ! |J|
    real(wp), intent(out) :: N(4)          ! shape functions at (ξ,η)

        real(wp) :: dN_dxi(4), dN_deta(4)
    real(wp) :: J(2,2), gamma(2,2),invJ(2,2)
    real(wp) :: dN_dx(4), dN_dy(4)
    real(wp) :: x1,y1,x2,y2,x3,y3,x4,y4
    real(wp) :: invdet

    real(wp) :: Ntil(4,8), gammatil(4,4)

    real(wp) :: L(3,4)                      ! L matrix 

    L = 0.0
    L(3,4) = 0.0
    L(1,1)=1.0
    L(2,4)=1.0
    L(3,2)=1.0
    L(3,3)=1.0

    ! Unpack nodal coord
    x1=xe(1)
    y1=xe(2)
    x2=xe(3)
    y2=xe(4)
    x3=xe(5)
    y3=xe(6)
    x4=xe(7)
    y4=xe(8)



    ! Shape functions
    N(1) = 0.25*(1.0 - xi)*(1.0 - eta)
    N(2) = 0.25*(1.0 + xi)*(1.0 - eta)
    N(3) = 0.25*(1.0 + xi)*(1.0 + eta)
    N(4) = 0.25*(1.0 - xi)*(1.0 + eta)

    ! Derivatives w.r.t. natural coords
    dN_dxi(1)  = -0.25*(1.0 - eta)
    dN_dxi(2)  =  0.25*(1.0 - eta)
    dN_dxi(3)  =  0.25*(1.0 + eta)
    dN_dxi(4)  = -0.25*(1.0 + eta)

    dN_deta(1) = -0.25*(1.0 - xi)
    dN_deta(2) = -0.25*(1.0 + xi)
    dN_deta(3) =  0.25*(1.0 + xi)
    dN_deta(4) =  0.25*(1.0 - xi)

    ! Jacobian J = [dx/dξ  dy/dξ; dx/dη  dy/dη]
    J = 0.0
    J(1,1) = dN_dxi(1)*x1 + dN_dxi(2)*x2 + dN_dxi(3)*x3 + dN_dxi(4)*x4
    J(1,2) = dN_dxi(1)*y1 + dN_dxi(2)*y2 + dN_dxi(3)*y3 + dN_dxi(4)*y4
    J(2,1) = dN_deta(1)*x1 + dN_deta(2)*x2 + dN_deta(3)*x3 + dN_deta(4)*x4
    J(2,2) = dN_deta(1)*y1 + dN_deta(2)*y2 + dN_deta(3)*y3 + dN_deta(4)*y4


    ! Compute determinant of Jacobian matrix
    detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)

    ! Inverse Jacobian
    invJ=0.0
    invdet   = 1.0/detJ
    invJ(1,1)=  J(2,2)
    invJ(1,2)= -J(1,2)
    invJ(2,1)= -J(2,1)
    invJ(2,2)=  J(1,1)

    gamma=0.0
    gamma = invdet * invJ

gammatil=0.0
gammatil(1,1)=gamma(1,1)
gammatil(1,2)=gamma(1,2)
gammatil(2,1)=gamma(2,1)
gammatil(2,2)=gamma(2,2)
gammatil(3,3)=gamma(1,1)
gammatil(3,4)=gamma(1,2)
gammatil(4,3)=gamma(2,1)
gammatil(4,4)=gamma(2,2)

Ntil=0.0
Ntil(1,1)=dN_dxi(1)
Ntil(2,1)=dN_deta(1)
Ntil(3,2)=dN_dxi(1)
Ntil(4,2)=dN_deta(1)
Ntil(1,3)=dN_dxi(2)
Ntil(2,3)=dN_deta(2)
Ntil(3,4)=dN_dxi(2)
Ntil(4,4)=dN_deta(2)
Ntil(1,5)=dN_dxi(3)
Ntil(2,5)=dN_deta(3)
Ntil(3,6)=dN_dxi(3)
Ntil(4,6)=dN_deta(3)
Ntil(1,7)=dN_dxi(4)
Ntil(2,7)=dN_deta(4)
Ntil(3,8)=dN_dxi(4)
Ntil(4,8)=dN_deta(4)

bmat=0.0
bmat = matmul(L, matmul(gammatil, Ntil))

    end subroutine shape



subroutine gauss_stress_strain(e, xe, deltaD_e, young, nu, h)


    use types
    use fedata
    implicit none

    !-------------------------
    ! Arguments
    !-------------------------
    integer, intent(in) :: e
    real(wp), intent(in)    :: xe(:)         ! [x1 y1 x2 y2 x3 y3 x4 y4]
    real(wp), intent(in)    :: deltaD_e(:)     ! element displacement increment (size 8)
    real(wp), intent(in)    :: young         ! Young's modulus
    real(wp), intent(in)    :: nu            ! Poisson's ratio
    real(wp), intent(in)    :: h             ! hardening modulus

    !-------------------------
    ! Locals
    !-------------------------
    real(wp) :: cmat(3,3), fact
    real(wp) :: xi(2), eta(2)
    real(wp) :: gp_xi, gp_eta
    real(wp) :: bmat(3,8), N(4), detJ
    real(wp) :: dStrain(3), dStress(3)
    real(wp) :: sigma_old(3), sigma_trial(3), sigma_new(3)
    real(wp) :: sigma_rel(3)  
    real(wp) :: dFdS(3), CdFdS(3)
    real(wp) :: sigma_e, F, dLambda, denom, vm_arg
    real(wp) :: sigy_old, sigy_new
    real(wp) :: G, sigy0, Ckine
    real(wp) :: alpha_old(3), alpha_new(3)
    integer  :: i, j, ig
    integer  :: mat_id, hardtype

        ! build constitutive matrix (plane stress)
        cmat = 0.0
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2


        xi  = (/ -1.0/sqrt(3.0),  1.0/sqrt(3.0) /)
        eta = xi 




    do j = 1,2
        gp_eta = eta(j)
        do i = 1,2
            gp_xi = xi(i)
            ig    = 2*(j-1) + i

            call shape(xe, gp_xi, gp_eta, bmat, detJ, N)

            ! Strain = B * deltaD
            dStrain = matmul(bmat,deltaD_e)

            ! Cumulative TOTAL strain at Gauss point
            eps_gauss(e, ig, 1:3) = eps_gauss(e, ig, 1:3) + dStrain

                        ! Matériau associé
            mat_id   = element(e)%mat
            hardtype = mprop(mat_id)%hard_type

            ! n-1 state
            sigma_old(1:3) = sigma_gauss(e, ig, 1:3)
            sigy_old       = sigy_gauss(e, ig)

            ! elastic try : σ_trial = σ_old + C Δε
            dStress = matmul(cmat, dStrain)
            sigma_trial = sigma_old + dStress

                        ! Par défaut (au cas où)
            sigma_new  = sigma_trial
            sigy_new   = sigy_old
            dLambda    = 0.0_wp
            alpha_new  = 0.0_wp

            select case (hardtype)

            !====================================================
            !  CASE (1) : ISOTROPIC HARDENING (ton code actuel)
            !====================================================
            case (1)


            ! von Mises
            sigma_e = sqrt( sigma_trial(1)**2 + sigma_trial(2)**2 - sigma_trial(1)*sigma_trial(2)+ 3.0_wp*sigma_trial(3)**2 )

            F = sigma_e - sigy_old

            sigy_new = sigy_old 

            if(F <= 0.0) then ! Elastic 
                sigma_new = sigma_trial
                dLambda   = 0.0_wp
                sigy_new  = sigy_old
            else ! Plastic 

                dFdS(1) = (2.0_wp*sigma_trial(1) - sigma_trial(2)) / (2.0_wp*sigma_e)
                dFdS(2) = (2.0_wp*sigma_trial(2) - sigma_trial(1)) / (2.0_wp*sigma_e)
                dFdS(3) =  6.0_wp*sigma_trial(3)                  / (2.0_wp*sigma_e)
                CdFdS = matmul(cmat, dFdS)
                denom = dot_product(dFdS, CdFdS) + h
                dLambda = dot_product(dFdS, dStress) / denom


                sigma_new = sigma_trial - dLambda * CdFdS
                if (dLambda > 0.0) then
                    sigy_new  = sigy_old  + h*dLambda
                    epsp_gauss(e, ig)   = epsp_gauss(e, ig) + dLambda
                endif

            endif


    !====================================================
            !  CASE (2) : KINEMATIC HARDENING (Prager linéaire)
            !====================================================
                case (2)

                ! Paramètres matériaux
                sigy0 = mprop(mat_id)%sigy0
                Ckine = mprop(mat_id)%Ckine

                ! Backstress précédent
                alpha_old(1:3) = alpha_gauss(e, ig, 1:3)

                ! Contrainte relative
                sigma_rel(1:3) = sigma_trial(1:3) - alpha_old(1:3)

                ! von Mises sur sigma_rel
                vm_arg = sigma_rel(1)**2 + sigma_rel(2)**2 - sigma_rel(1)*sigma_rel(2) + 3.0_wp*sigma_rel(3)**2
                if (vm_arg < 0.0_wp) vm_arg = 0.0_wp
                sigma_e = sqrt(vm_arg)

                F = sigma_e - sigy0   ! rayon constant


                if (F <= sigy0*1e-8) then
                    sigma_new = sigma_trial
                    sigy_new  = sigy0
                    dLambda   = 0.0_wp
                    alpha_new = alpha_old
                else
                    if (sigma_e < 1.0e-12_wp) then
                        sigma_new = sigma_trial
                        sigy_new  = sigy0
                        dLambda   = 0.0_wp
                        alpha_new = alpha_old
                    else
                        dFdS(1) = (2.0_wp*sigma_rel(1) - sigma_rel(2)) / (2.0_wp*sigma_e)
                        dFdS(2) = (2.0_wp*sigma_rel(2) - sigma_rel(1)) / (2.0_wp*sigma_e)
                        dFdS(3) =  6.0_wp*sigma_rel(3)                 / (2.0_wp*sigma_e)

                        CdFdS = matmul(cmat, dFdS)
                        denom = dot_product(dFdS, CdFdS) + (2.0_wp/3.0_wp)*Ckine
                        if (abs(denom) < 1.0e-16_wp) denom = 1.0e-16_wp  ! éviter crash

                        dLambda = dot_product(dFdS, dStress) / denom
                        sigma_new = sigma_trial - dLambda * CdFdS
                        alpha_new = alpha_old + (2.0_wp/3.0_wp)*Ckine*dLambda * dFdS
                        epsp_gauss(e, ig) = epsp_gauss(e, ig) + dLambda
                        sigy_new  = sigy0

                        
                end if
            end if


            !====================================================
            !  AUTRES CAS : on peut traiter comme élastique/iso
            !====================================================
            case default
                ! pour l'instant : élastique
                sigma_new = sigma_trial
                sigy_new  = sigy_old
                dLambda   = 0.0_wp
                alpha_new = 0.0_wp

            end select



            sigma_gauss(e, ig, 1:3) = sigma_new
            sigy_gauss(e, ig)       = sigy_new
            dlambda_gauss(e, ig)    = dLambda
            alpha_gauss(e, ig, 1:3) = alpha_new(1:3)

        end do
    end do

    end subroutine gauss_stress_strain


end module plane42