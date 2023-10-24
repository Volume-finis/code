!-------------------------------------------------------------------
! module contenant les subroutines :
!  maillage
!  lecture_maillage_medit
!  connectivite
!  aire
!
!  VERSION
!  	15/10/2022 par Luc Mieussens
!-------------------------------------------------------------------

module mod_maillage

  use mod_precision

  implicit none


contains 

  !-------------------------------------------------------------------
  ! appel a differentes subroutines pour lire le maillage, construire
  ! les informations de connectivite, et calculer des quantites geometriques
  !-------------------------------------------------------------------
  subroutine maillage(fichier_maillage, nb_mailles, nb_aretes                        &
       &        , coord_noeud, noeud_maille, aire_maille, l_arete, d_arete           &
       &        , milieu_arete, arete_maille, maille_arete, cl_arete)

    !--- entree
    character(len=*) :: fichier_maillage

    !--- sorties
    integer, intent(out) :: nb_mailles, nb_aretes
    real(pr), dimension(:), allocatable, intent(out) :: aire_maille, d_arete, l_arete
    real(pr), dimension(:,:), allocatable, intent(out) :: coord_noeud, milieu_arete
    integer, dimension(:,:), allocatable, intent(out) :: arete_maille, noeud_maille, maille_arete
    !-- format msh    character(len=50), dimension(:), allocatable, intent(out) :: cl_arete
    !-- format medit
    integer, dimension(:), allocatable, intent(out) :: cl_arete

    !--- locales
    integer :: nb_noeuds, i, j, nb_aretes_bord
    integer, dimension(:,:), allocatable :: noeud_arete_bord, noeud_arete
    !- pour format msh :     character(len=50), dimension(:), allocatable :: cl_arete_bord
    !- pour format medit :
    integer, dimension(:), allocatable :: cl_arete_bord
    real(pr), dimension(:,:), allocatable :: coord_centre_maille
    real(pr), dimension(2) :: A, B, C, Cg, Cd, M
    real(pr) :: delta, ra2, rb2, rc2, dx, dy
    integer :: tg, td


    !--- lecture du maillage
    call lecture_maillage_medit(fichier_maillage,noeud_maille,coord_noeud,cl_arete_bord,noeud_arete_bord)

    nb_noeuds = size(coord_noeud,1)
    nb_mailles = size(noeud_maille,1)
    nb_aretes_bord = size(cl_arete_bord)

    print*,'nb_noeuds :',nb_noeuds
    print*,'nb_mailles',nb_mailles


    !--- connectivité
    call connectivite(coord_noeud,noeud_maille,noeud_arete,arete_maille,maille_arete)
    nb_aretes = size(noeud_arete,1)


    !--- numeros de cl pour les aretes de bord
    !  il nous faut le code cl d'une arete de bord par son numéro global
    ! boucle aretes i
    !! quand un n° triangle = 0, alors arete de bord
    !! boucle aretes de bord j
!!! quand noeuds = noeuds arete, alors c'est elle
!!! cl_arete(i) = cl_arete_bord(j)
!!!! nb : pas optimal, car cl_arete trop gros tableau
!!! on pourrait donne la cl au triangle de bord

    allocate(cl_arete(nb_aretes))
    do i=1,nb_aretes

       if (maille_arete(i,1)==0.or.maille_arete(i,2)==0) then

          do j=1,nb_aretes_bord

             if ( (noeud_arete_bord(j,1)==noeud_arete(i,1)                            &
                  &              .and.noeud_arete_bord(j,2)==noeud_arete(i,2))        &
                  & .or.(noeud_arete_bord(j,1)==noeud_arete(i,2)                      &
                  &              .and.noeud_arete_bord(j,2)==noeud_arete(i,1))) then

                cl_arete(i) = cl_arete_bord(j)

             end if

          end do

       end if

    end do


    !--- aires des mailles
    allocate(aire_maille(size(noeud_maille,1)))
    call aire(noeud_maille,coord_noeud,aire_maille)


    !--- centres cercles circonscrits
    allocate(coord_centre_maille(nb_mailles,2))
    do i = 1,nb_mailles

       !--- sommets de la maille
       A = coord_noeud(noeud_maille(i,1),:)
       B = coord_noeud(noeud_maille(i,2),:)
       C = coord_noeud(noeud_maille(i,3),:)

       delta = 2 * ( (B(1)-A(1))*(C(2)-A(2)) - (B(2)-A(2))*(C(1)-A(1)) )

       ra2 = sum(A**2) ; rb2 = sum(B**2) ; rc2 = sum(C**2) 
       dx = (rb2-ra2) * (C(2) - A(2)) - (rc2-ra2) * (B(2) - A(2))
       dy = (rb2-ra2) * (C(1) - A(1)) - (rc2-ra2) * (B(1) - A(1)) 

       coord_centre_maille(i,1) = dx / delta
       coord_centre_maille(i,2) = - dy / delta

    end do


    !--- coef aretes : longueur des normales
    allocate(d_arete(nb_aretes))
    do i=1,nb_aretes

       tg = maille_arete(i,1)
       Cg = coord_centre_maille(tg,:)

       td = maille_arete(i,2)

       !--- si pas arete de bord
       if (td/=0) then

          Cd = coord_centre_maille(td,:)

          d_arete(i) = sqrt(sum( (Cd - Cg)**2 ))

       else

          A = coord_noeud(noeud_arete(i,1),:)
          B = coord_noeud(noeud_arete(i,2),:)
          M = (A + B) / 2

          d_arete(i) = sqrt(sum( (M - Cg)**2 ))

       end if

    end do


    !--- coef aretes : coord milieu arete
    allocate(milieu_arete(nb_aretes,2))
    allocate(l_arete(nb_aretes))
    do i=1,nb_aretes

       A = coord_noeud(noeud_arete(i,1),:)
       B = coord_noeud(noeud_arete(i,2),:)
       milieu_arete(i,:) = (A + B) / 2
       l_arete(i) = sqrt(sum((B-A)**2))

    end do


    !--- verification : sortie du maillage au format vtk

    open(unit=1,file='maillage.vtk')

    write(1,'(1A26)') '# vtk DataFile Version 2.0'
    write(1,*) 'test maillage'
    write(1,*) 'ASCII'
    write(1,*) 'DATASET UNSTRUCTURED_GRID'

    write(1,*) 'POINTS',nb_noeuds,' double'
    do i=1,nb_noeuds
       write(1,*) coord_noeud(i,1),coord_noeud(i,2),0.0_pr
    end do

    write(1,*) 'CELLS ',nb_mailles,nb_mailles*4
    do i=1,nb_mailles
       write(1,*) 3,noeud_maille(i,1)-1, noeud_maille(i,2)-1, noeud_maille(i,3)-1
    end do

    write(1,*) 'CELL_TYPES ',nb_mailles
    do i=1,nb_mailles
       write(1,*) 5
    end do

    close(1)

  end subroutine maillage


  !-------------------------------------------------------------------
  ! lecture du maillage sous format .mesh medit MeshVersionFormatted 2
  !-------------------------------------------------------------------
  subroutine lecture_maillage_medit(fich,noeud_maille,coord_noeud,cl_arete_bord,noeud_arete_bord)

    implicit none

    !--- entrees
    character(len=*), intent(in) :: fich

    !--- sorties
    integer, dimension(:,:), allocatable, intent(out) :: noeud_maille
    real(pr), dimension(:,:), allocatable, intent(out) :: coord_noeud
    integer, dimension(:), allocatable, intent(out) :: cl_arete_bord
    integer, dimension(:,:), allocatable, intent(out) :: noeud_arete_bord


    !--- locales    
    integer :: nb_noeuds 
    integer :: i
    integer :: nb_aretes_bord, nb_mailles


    open(unit=10,file=trim(adjustl(fich)))

    !--- 3 premières lignes (format)
    read(10,*)
    read(10,*)
    read(10,*) 

    !--- lecture sommets
    read(10,*)
    read(10,*) nb_noeuds
    allocate(coord_noeud(1:nb_noeuds,1:2))
    do i=1, nb_noeuds

       read(10,*) coord_noeud(i,1), coord_noeud(i,2) ! le 3eme nbre est la coord z, le 4eme est une reference inutilisee ici

    end do

    !--- lecture aretes de bord
    read(10,*)
    read(10,*) nb_aretes_bord
    allocate(noeud_arete_bord(nb_aretes_bord,2))
    allocate(cl_arete_bord(nb_aretes_bord))
    do i=1, nb_aretes_bord

       read(10,*) noeud_arete_bord(i,1), noeud_arete_bord(i,2), cl_arete_bord(i)

    end do

    !--- lecture mailles
    read(10,*)
    read(10,*) nb_mailles
    allocate(noeud_maille(nb_mailles,3))
    do i=1, nb_mailles

       read(10,*) noeud_maille(i,1), noeud_maille(i,2), noeud_maille(i,3) ! le 4eme nbre est le code physique de la maille, utile pour CI, mais pas utilise ici

    end do

  end subroutine lecture_maillage_medit



  !-------------------------------------------------------------------
  ! calcul des informations de connectivite
  !-------------------------------------------------------------------
  subroutine connectivite(p,s,e,ar,trig)

    implicit none

    !--- entrees
    real(pr), dimension(:,:), intent(in) :: p
    integer, dimension(:,:), intent(in) :: s

    !--- sorties
    integer, dimension(:,:), allocatable, intent(out) :: e, ar, trig

    !--- locales
    integer :: i,nb_noeuds,nb_trig,nb_cotes,n1,n2,n3
    integer, dimension(:,:), allocatable :: flag


    !O)- nb noeuds et mailles
    nb_noeuds = size(p,1)
    nb_trig = size(s,1)


    !1)-calcul du nombre d'arêtes nb_cotes
    nb_cotes=0
    allocate(flag(1:nb_noeuds,1:nb_noeuds))      
    ! flag(n1,n2) contient en fin de boucle le numéro de l'arête liant les sommets de numéros n1 et n2
    flag=0
    do i=1,nb_trig
       n1=S(i,1); n2=S(i,2); n3=S(i,3)
       if (flag(n1,n2)==0.or.flag(n2,n1)==0) then ! alors arete non numerotee
          nb_cotes=nb_cotes+1
          flag(n1,n2)=nb_cotes; flag(n2,n1)=nb_cotes
       end if
       if (flag(n2,n3)==0.or.flag(n3,n2)==0) then ! alors arete non numerotee
          nb_cotes=nb_cotes+1
          flag(n2,n3)=nb_cotes; flag(n3,n2)=nb_cotes
       end if
       if (flag(n3,n1)==0.or.flag(n1,n3)==0) then ! alors arete non numerotee
          nb_cotes=nb_cotes+1
          flag(n3,n1)=nb_cotes; flag(n1,n3)=nb_cotes
       end if
    end do

    !2)-remplissage de e
    allocate(e(1:nb_cotes,2))    ! numéros des noeuds extrémités de chaque arête
    do i=1,nb_trig
       n1=S(i,1); n2=S(i,2); n3=S(i,3)
       e(flag(n1,n2),1)=n1; e(flag(n1,n2),2)=n2
       e(flag(n2,n3),1)=n2; e(flag(n2,n3),2)=n3
       e(flag(n3,n1),1)=n3; e(flag(n3,n1),2)=n1
    end do

    !3)-remplissage de ar et trig
    allocate(ar(1:nb_trig,3))    ! numéros des arêtes de chaque triangle
    allocate(trig(1:nb_cotes,2)) ! numéro des triangles ayant l'arête i en commun
    trig=0
    do i=1,nb_trig
       n1=S(i,1); n2=S(i,2); n3=S(i,3)
       !-arete (n1,n2)
       ar(i,1)=flag(n1,n2)
       if (trig(ar(i,1),1)/=0) then
          trig(ar(i,1),2)=i
       else
          trig(ar(i,1),1)=i
       end if
       !-arete (n2,n3)
       ar(i,2)=flag(n2,n3)
       if (trig(ar(i,2),1)/=0) then
          trig(ar(i,2),2)=i
       else
          trig(ar(i,2),1)=i
       end if
       !-arete (n3,n1)
       ar(i,3)=flag(n3,n1)
       if (trig(ar(i,3),1)/=0) then
          trig(ar(i,3),2)=i
       else
          trig(ar(i,3),1)=i
       end if

    end do

  end subroutine connectivite



  !-------------------------------------------------------------------
  ! calcul des aires des mailles
  !-------------------------------------------------------------------
  subroutine aire(noeud_maille,coord_noeud,aire_maille)

    !--- entrees
    integer, dimension(:,:), intent(in) :: noeud_maille
    real(pr), dimension(:,:), intent(in) :: coord_noeud

    !--- sorties
    real(pr), dimension(size(noeud_maille,1)), intent(out) :: aire_maille

    !--- locales
    integer :: i
    real(pr), dimension(2) :: a,b,c


    !--- calcul de l'aire des mailles
    do i = 1,size(noeud_maille,1)

       a = coord_noeud(noeud_maille(i,1),:)
       b = coord_noeud(noeud_maille(i,2),:)
       c = coord_noeud(noeud_maille(i,3),:)

       aire_maille(i) = abs( ((b(1) - a(1))*(c(2) - a(2)))   &
            &               - ((b(2) - a(2))*(c(1) - a(1))) ) / 2
    end do


  end subroutine aire


end module mod_maillage
