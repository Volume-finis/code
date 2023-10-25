program Principal_Program
    use mod_precision
    use mod_maillage
    use mod_sortie


    implicit none

    real(PR), dimension(:), allocatable     :: T, Tnp1
    integer, dimension(:), allocatable      :: cl_arete !cl à arrête k si elle est au bord
    real(PR), dimension(:), allocatable     :: aire_maille, d_arete, l_arete
    real(PR), dimension(:,:), allocatable   :: milieu_arete, coord_noeud !tq coord_noeud(i,1) = abscisse noeud i 
    integer, dimension(:,:), allocatable    :: arete_maille, noeud_maille, maille_arete

    real(PR)                                :: tmax, D, cfl, Tinit, Rmax, Rint, Rext, T1, T2, T3
    real(PR)                                :: delta_t, sum_arete, F
    real(PR)                                :: Sol_exact, sum_Sol_exact_carre, X_centre,Y_centre, X_sum, Y_sum, Rayon
    real(PR)                                :: Err_N_Euclidienne, sum
    real(PR)                                :: L_triangle_min

    integer                                 :: nb_aretes, nb_mailles
    integer                                 :: i, j, n, a, M_1, M_2, nmax, m, Cas_Test, Affichage_vtk

    character(len=20)                       :: Ntest

    character(len=30) :: fichier_maillage

    print*, "Afin de résoudre l'approximation en volumes finis de l'équation de la chaleur " 
    print*, "en dimension 2, vous aurez le choix d'effectuer 4 tests sur des structures différentes :" 
    print*, " "
    print*, "- Cas test 1 : Carré", "                  ", "- Cas test 2 : Cercle"
    print*, "- Cas test 3 : Couronne", "               ", "- Cas test 4 : Domaine quelconque"
    print*, "Quels cas choisissez vous d'étudier (entre 1 et 4) : " 
    ! Cas_Test = 4
    read*, Cas_Test

    if (Cas_Test /= 1 .AND. Cas_Test /= 2 .AND. Cas_Test /= 3 .AND. Cas_Test /= 4)then 
        stop "Vous n'avez pas choisi un Cas Test existant"
    end if 

    print*, "Voulez-vous afficher les fichiers .vtk dans votre dossier :"
    print*, "Si oui tappez 1"
    read*, Affichage_vtk        !! permet d'écrire ou non les fichiers .vtk pour les utiliser avec paraview
    ! Affichage_vtk = 0

    write(Ntest,*) Cas_Test

    open(11, file="Donnee_Cas_test"//trim(adjustl(Ntest))//".dat")
        read(11,*) tmax
        read(11,*) D
        read(11,*) cfl
        read(11,*) Tinit
        read(11,*) T1  ! Tg
        read(11,*) T2  ! Td ou Tb c'est la valeur de référence
        if (Cas_Test==1 .OR. Cas_Test ==4) then 
            read(11,*) T3  ! Phi_b
        elseif (Cas_Test==2) then 
            read(11,*) Rmax
        end if
        if (Cas_Test==3) then 
            read(11,*) Rint
            read(11,*) Rext
        end if
        read(11,*) fichier_maillage 
    close(11)

    print*, "--------------------------------"
    print*, "Les paramètres utilisés : "
    print*, "Tmax = ", tmax,", D = ", D, ", cfl = ", cfl, ", Tinit = ", Tinit
    print*,  "Tg = ", T1, ", Tb = ", T2, ", Phi_b = ", T3
    print*, "--------------------------------"

    call maillage(fichier_maillage , nb_mailles, nb_aretes &
        & , coord_noeud, noeud_maille, aire_maille, l_arete, d_arete &
        & , milieu_arete, arete_maille, maille_arete, cl_arete )


    allocate(T(1:nb_mailles), Tnp1(1:nb_mailles))
    
    delta_t = 10._PR

    do i = 1, nb_mailles
        sum_arete = 0._PR 
        do j = 1, 3 
            sum_arete = sum_arete + (l_arete(arete_maille(i,j))*D) / d_arete(arete_maille(i,j))
        end do 
        if (delta_t >= aire_maille(i)/(sum_arete)) then 
            delta_t = aire_maille(i)/(sum_arete)
        end if 

    end do 

    print*, "--------------------------------"
    print*, "delta_t = ", delta_t  !! delta_t =  2.8887438505932329E − 004 pour "cas_1d.mesh"


    T = Tinit
    Tnp1 = Tinit

    nmax = int(tmax / delta_t)  
    call sortie (0, Tnp1, coord_noeud, noeud_maille)

    do n = 1, nmax
        do a = 1, nb_aretes

            M_1 = maille_arete(a,1)
            M_2 = maille_arete(a,2)

            if ( M_2 == 0 ) then !arete est de bord
                if ((cl_arete(a)) == 10) then   !! 7 ou 10
                    F = - D * (T1-T(M_1))/ d_arete(a)
                    Tnp1(M_1) = Tnp1(M_1) - (delta_t/aire_maille(M_1))* F * l_arete(a)
                elseif ((cl_arete(a)) == 11) then !! 6 ou 11
                    F = - D * (T2-T(M_1))/ d_arete(a)
                    Tnp1(M_1) = Tnp1(M_1) - (delta_t/aire_maille(M_1))* F * l_arete(a)
                else ! cl_arete(a)) == 20
                    F = T3 !!Phi_b
                    Tnp1(M_1) = Tnp1(M_1) - (delta_t/aire_maille(M_1))* F * l_arete(a)
                end if 
            
            else !! maille intérieur
                F = -D*(T(M_2)-T(M_1))/ d_arete(a)
                Tnp1(M_1) = Tnp1(M_1) - (delta_t/aire_maille(M_1))* F * l_arete(a)
                Tnp1(M_2) = Tnp1(M_2) + (delta_t/aire_maille(M_2))* F * l_arete(a)
        
            end if 
               
        end do 

        if (Cas_Test == 3) then
            if (n==nmax) then 
                if (Norme_infini(T) / Norme_infini(Tnp1)  < 0.99) then 
                    print*, "Pas en stationnaire, augmentez tmax"  !! permet de calculer l'erreur en 
                else                                               !! stationnaire
                    print*, "Vous êtes en stationnaire, vous pouvez utiliser tmax"
                end if 
            end if
        end if 

        T = Tnp1
    
        if (Affichage_vtk == 1 .AND.  (modulo(n, 20) == 0)) then  !! Augmenter le quotient si vous avez 
            call sortie (n, Tnp1, coord_noeud, noeud_maille)      !! trop de fichiers .vtk et diminuez
        end if                                                    !! le pour un écoulement du temps plus 
                                                                  !! précis. Ex mettre 50 pour Cas = 4 
                                                                  !! et 10 pour les autres cas
    end do 

    Err_N_Euclidienne = 0._PR      !! Cela va nous permettre de calculer la norme 2 relative global
    sum_Sol_exact_carre = 0._PR

    L_triangle_min = 10._PR

    do m = 1, nb_mailles 

        if (L_triangle_min > sqrt(aire_maille(m))) then         !! On prend cette longueur mais ce n'est pas forcément la hauteur min
            L_triangle_min = sqrt(aire_maille(m))               !! il s'agit de la hauteur min pour un triangle équilatéral
        end if

        if (Cas_Test == 1) then
            X_sum = 0._PR
            do j = 1, 3 
                X_sum = X_sum + milieu_arete(arete_maille(m, j),1)  !! j'ai mit 1 pour avoir l'abscisse 
            end do 
            X_centre = X_sum / 3

            sum = 0._PR

            do i = 1,100
                sum = sum + (2/(i*pi))*((-1)**i) * exp(-D*(i*pi)*(i*pi)*nmax * delta_t)*sin(i*pi*X_centre)
            end do 

            Sol_exact = T1 + (T2-T1)*(X_centre + sum)

            
        elseif (Cas_Test == 2) then !le cercle
            X_sum = 0._PR
            Y_sum = 0._PR
            do j = 1, 3 
                X_sum = X_sum + milieu_arete(arete_maille(m, j),1)  !! j'ai mit 1 pour avoir l'abscisse 
                Y_sum = Y_sum + milieu_arete(arete_maille(m, j),2)  !! j'ai mit 2 pour avoir l'ordonnée 
            end do 
            X_centre = X_sum / 3
            Y_centre = Y_sum / 3

            Rayon = sqrt(X_centre*X_centre +Y_centre*Y_centre)
            Sol_exact = T1d_axi(nmax * delta_t,Rayon)

        elseif (Cas_Test == 3) then !la couronne
            X_sum = 0._PR
            Y_sum = 0._PR
            do j = 1, 3 
                X_sum = X_sum + milieu_arete(arete_maille(m, j),1)  !! j'ai mit 1 pour avoir l'abscisse 
                Y_sum = Y_sum + milieu_arete(arete_maille(m, j),2)  !! j'ai mit 2 pour avoir l'ordonnée 
            end do 
            X_centre = X_sum / 3
            Y_centre = Y_sum / 3

            Rayon = sqrt(X_centre*X_centre +Y_centre*Y_centre)

            Sol_exact = T1 + (T2-T1)*log(Rayon/Rint)/log(Rext/Rint)
        
        else !! pour le Cas 4
            stop "Pas de calcul d'erreur pour le Cas test"
        end if

        sum_Sol_exact_carre = sum_Sol_exact_carre + Sol_exact * Sol_exact * aire_maille(m)
        Err_N_Euclidienne = Err_N_Euclidienne + ((Sol_exact - Tnp1(m))**2) * aire_maille(m)


    end do 

    
    Err_N_Euclidienne = Err_N_Euclidienne / (sum_Sol_exact_carre)
    Err_N_Euclidienne = sqrt(Err_N_Euclidienne)

    print*, "Hauteur min des mailles triangulaires = ", L_triangle_min
    print*, "--------------------------------"
    print*, "Erreur Euclidienne / relative est : ", Err_N_Euclidienne !! il s'agit de l'erreur relative global
   
   
   
   
    !! 1er maillage,  
    ! nb_noeuds :         109
    ! nb_mailles         184
    ! delta_t =    2.8887438505932329E-004
    ! Erreur Euclidienne est :    8.6661426550605405E-003


    ! nb_noeuds :         369
    ! nb_mailles         672
    ! delta t = 2.8215729248389001E-005
    ! Erreur Euclidienne est :    3.6121179818330751E-003

    ! Carré fait à la main fin test 1.
    ! nb_noeuds :          98
    ! nb_mailles         162
    ! delta_t =    6.1613231431596456E-004
    ! Erreur Absolue est :   0.89119633683193811     
    ! Erreur Euclidienne est :    5.4937740745659817E-003

    
    deallocate(T, Tnp1)


contains 

function T1d_axi(t,r) result(y)

    real(pr), intent(in) :: t, r
    real(pr) :: y
  
    integer :: n
    real(pr), dimension(1:20) :: z, c
  
    !--- zeros de J0
    z = (/ 2.404825558, 5.520078110, 8.653727913, 11.79153444, 14.93091771 &
         & , 18.07106397, 21.21163663, 24.35247153, 27.49347913, 30.63460647 &
         & , 33.77582021, 36.91709835, 40.05842576, 43.19979171, 46.34118837 &
         & , 49.48260990, 52.62405184, 55.76551076, 58.90698393, 62.04846919/)
  
    !--- coeff cn
    do n=1,20
       c(n) = 2 * (Tinit - T2) / ( z(n) * Bessel_JN(1,z(n)) )
    end do
  
    !--- T exacte
    y  = T2
    do n=1,20
       y = y + c(n) * exp(-z(n)**2 * D/Rmax**2 * t) * Bessel_JN(0,r*z(n)/Rmax)
    end do
  
end function T1d_axi

function Norme_infini(T) result(norme)
    real(PR), dimension(:) ::  T
    real(PR) :: norme
    integer :: i

    norme = 0._PR 
    do i = 1, size(T)
        if (norme < T(i)) then 
            norme = T(i)
        end if 
    end do 

end function Norme_infini






end program Principal_Program