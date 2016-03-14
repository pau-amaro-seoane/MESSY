c======================================================================
c     routines d'entree/sortie pour Arbre Binaire d'amas
c======================================================================
c
c----------------------------------------------------------------------
      subroutine EcrireArbreASC(iFich)
c----------------------------------------------------------------------
c     Ecriture ASCII de l'arbre a fin de verifications
c     sous forme de tableau
c      
      include 'EvolAmas_ArbreIO_inc.f'
c
      integer iFich
c
      integer i,iSE
c......................................................................
c
      write(iFich,120)iAdam
 120  format(
     $     '#  TABLEAU DE STRUCTURE D''ARBRE BINAIRE',/,
     $     '# ======================================',/,
     $     '#  iAdam = ',I6 )
      write(iFich,130)
     $     'iNoeud','iSE','iPere','iFilsG','iFilsD','NbDec_G',
     $     'R','M','dA','dB'
 130  format(10('|',A12),/,10('|------------'))
      do i = 1,NbNoeuds
         iSE=iSE_arbre(i)
         write(iFich,140)
     $        i,iSE,iPere_arbre(i),
     $        iFilsG_arbre(i),iFilsD_arbre(i),iNbDecG_arbre(i),
     $        R_SE(iSE),M_SE(iSE),dA_arbre(i),dB_arbre(i)
 140     format(6('|',I12),4('|',E12.4))
      end do
c
      end   
c  
c----------------------------------------------------------------------
      subroutine EcrireArbreASC2(iFich)
c----------------------------------------------------------------------
c     Ecriture ASCII de l'arbre a fin de verifications
c     sous forme de tableau
c      
      include 'EvolAmas_ArbreIO_inc.f'
c
      integer iFich
c
      integer i,iSE
c......................................................................
c
      write(iFich,120)iAdam
 120  format(
     $     '#  TABLEAU DE STRUCTURE D''ARBRE BINAIRE',/,
     $     '# ======================================',/,
     $     '#  iAdam = ',I6 )
      write(iFich,130)
     $     '1:iNoeud','2:iSE','3:iPere','4:iFilsG','5:iFilsD','6:NbDec_G',
     $     '7:R','8:M','9:dA','10:dB'
 130  format(11(' ',A12),/,10('-------------'))
      do i = 1,NbNoeuds
         iSE=iSE_arbre(i)
         write(iFich,140)
     $        i,iSE,iPere_arbre(i),
     $        iFilsG_arbre(i),iFilsD_arbre(i),iNbDecG_arbre(i),
     $        R_SE(iSE),M_SE(iSE),dA_arbre(i),dB_arbre(i)
 140     format(6(' ',I12),5(' ',E12.4))
      end do
c
      end   
c   
#if (_HIER_ARBRE_>0) 
c----------------------------------------------------------------------
      subroutine EcrireArbreASC3(iFich)
c----------------------------------------------------------------------
c     Ecriture ASCII de l'arbre a fin de verifications
c     sous forme de tableau
c     Version indiquant le niveau hierarchique
c      
      include 'EvolAmas_ArbreIO_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c var locales :
c ^^^^^^^^^^^^^
      integer i,Niv_Hier(iDimArbre)
c......................................................................
c
      call Determiner_Hierarchie(Niv_Hier)

      write(iFich,120)NbSE,iAdam
 120  format(
     $     '#  TABLEAU DE STRUCTURE D''ARBRE BINAIRE',/,
     $     '# ======================================',/,
     $     '#  NbSE  = ',I6,/,
     $     '#  iAdam = ',I6,/,
     $     '# 1: iSE 2: iPere 3: iFilsG 4: iFilsD 5: iNbDecG',
     $     ' 6: dA 7: dB 8: R 9: Niv_Hier')
      do i=1,NbNoeuds
         write(iFich,130) iSE_arbre(i),iPere_arbre(i),
     $        iFilsG_arbre(i),iFilsD_arbre(i),iNbDecG_arbre(i),
     $        dA_arbre(i),dB_arbre(i),R_SE(iSE_arbre(i)),Niv_Hier(i)
 130     format(5(I7,' '),3(E12.5,' '),I5)
      end do
c
      end
c
#endif /* if (_HIER_ARBRE_>0) */
c----------------------------------------------------------------------
      subroutine EcrireArbreBIN(iFich)
c----------------------------------------------------------------------
c     Ecriture binaire de l'arbre a fin de verifications
c     sous forme de tableau
c      
      include 'EvolAmas_ArbreIO_inc.f'
c
      integer iFich

      integer i
c
c......................................................................
c
      write(iFich)
     $     NbNoeuds,NbSE,iAdam,B0,
     $     (iNoeud_SE(i)     ,i=1,NbNoeuds),
     $     (iSE_arbre(i)     ,i=1,NbNoeuds),
     $     (dA_arbre(i)      ,i=1,NbNoeuds),
     $     (dB_arbre(i)      ,i=1,NbNoeuds),
     $     (iPere_arbre(i)   ,i=0,NbNoeuds),
     $     (iFilsG_arbre(i)  ,i=0,NbNoeuds),
     $     (iFilsD_arbre(i)  ,i=0,NbNoeuds),
     $     (iNbDecG_arbre(i) ,i=0,NbNoeuds)
c
      end   
c  
c----------------------------------------------------------------------
      subroutine LireArbreASC(iFich)
c----------------------------------------------------------------------
c     Lit une structure d'arbre du type ecrit 
c     par Ecrire_Arbre_Amas_ASC.
c      
      include 'EvolAmas_ArbreIO_inc.f'
c
      integer iFich
c
      integer i
      character*4 car
c......................................................................
c
      NbNoeuds=0
      Nb_Parcours_Arbre = 0
c
c---- d'abord sauter les remarques de debut
c
 1    continue
      read(iFich,'(A)') car
      if (car(1:1).EQ.'#' .OR. car(1:1).EQ.' '
     $     .OR. LEN(car).EQ.0) goto 1

      read(iFich,'(A)') car     ! on saute l en-tete
      
 2    continue                  ! debut de la lecture des tableaux
      read(iFich,140,END=9,ERR=99)
     $        i,iSE_arbre(i),iPere_arbre(i),
     $        iFilsG_arbre(i),iFilsD_arbre(i),iNbDecG_arbre(i),
     $        R_SE(iSE_arbre(i)),M_SE(iSE_arbre(i)),dA_arbre(i),dB_arbre(i)
 140  format(6('|',I12),4('|',E12.4))
      NbNoeuds = NbNoeuds+1
      goto 2

 9    continue                  ! fin de lecture
c
c---- on recherche le noeud-racine
c
      iAdam=0
      do i=1,NbNoeuds
         if (iPere_arbre(i).EQ.0) then
            iAdam = i
            goto 19
         end if
      end do
         
 19   continue
      return
c
c---- En cas d'erreur de lecture
c
 99   write(0,99100)
99100 format(/,'!!! Erreur de lecture de l''arbre binaire !!!')
      stop
c
      end
c 
c----------------------------------------------------------------------
      subroutine LireArbreASC2(iFich)
c----------------------------------------------------------------------
c     Lit une structure d'arbre du type ecrit 
c     par Ecrire_Arbre_Amas_ASC.
c      
      include 'EvolAmas_ArbreIO_inc.f'
c
      integer iFich
c
      integer i
      character*4 car
c......................................................................
c
      NbNoeuds=0
      Nb_Parcours_Arbre = 0
c
c---- d'abord sauter les remarques de debut
c
 1    continue
      read(iFich,'(A)') car
      if (car(1:1).EQ.'#' .OR. car(1:1).EQ.' '
     $     .OR. LEN(car).EQ.0) goto 1

      read(iFich,'(A)') car     ! on saute l en-tete
      
 2    continue                  ! debut de la lecture des tableaux
      read(iFich,140,END=9,ERR=99)
     $     i,iSE_arbre(i),iPere_arbre(i),
     $     iFilsG_arbre(i),iFilsD_arbre(i),iNbDecG_arbre(i),
     $     R_SE(iSE_arbre(i)),M_SE(iSE_arbre(i)),
     $     dA_arbre(i),dB_arbre(i)
 140  format(6(' ',I12),4(' ',E12.4))
      NbNoeuds = NbNoeuds+1
      goto 2

 9    continue                  ! fin de lecture
c
c---- on recherche le noeud-racine
c
      iAdam=0
      do i=1,NbNoeuds
         if (iPere_arbre(i).EQ.0) then
            iAdam = i
            goto 19
         end if
      end do
         
 19   continue
      return
c
c---- En cas d'erreur de lecture
c
 99   write(0,99100)
99100 format(/,'!!! Erreur de lecture de l''arbre binaire !!!')
      stop
c
      end
c
c----------------------------------------------------------------------
      subroutine LireArbreBIN(iFich)
c----------------------------------------------------------------------
c     Ecriture binaire de l'arbre a fin de verifications
c     sous forme de tableau
c      
      include 'EvolAmas_ArbreIO_inc.f'
c
      integer iFich
c
      integer i
c......................................................................
c
      Nb_Parcours_Arbre = 0
      read(iFich,ERR=99,END=99)
     $     NbNoeuds,NbSE,iAdam,B0,
     $     (iNoeud_SE(i)     ,i=1,NbNoeuds),
     $     (iSE_arbre(i)     ,i=1,NbNoeuds),
     $     (dA_arbre(i)      ,i=1,NbNoeuds),
     $     (dB_arbre(i)      ,i=1,NbNoeuds),
     $     (iPere_arbre(i)   ,i=0,NbNoeuds),
     $     (iFilsG_arbre(i)  ,i=0,NbNoeuds),
     $     (iFilsD_arbre(i)  ,i=0,NbNoeuds)
      return
c
c---- En cas d'erreur de lecture
c
 99   write(0,99100)
99100 format(/,'!!! Erreur de lecture de l''arbre binaire !!!')
      stop
c
c
      end   
c  
c----------------------------------------------------------------------
      subroutine EcritArbre_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_ArbreIO_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich_XDR         ! Numero XDR du fichier (attribue par initxdr)
c
c var locales :
c ^^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichArbre_XDR)
c
      call XDR_Ecrit_i (iFich_XDR,'NbSE',NbSE)
      call XDR_Ecrit_i (iFich_XDR,'iAdam',iAdam)
      call XDR_Ecrit_ti
     $     (iFich_XDR,'iSE_arbre',    NbNoeuds,  iSE_arbre)
      call XDR_Ecrit_ti
     $     (iFich_XDR,'iPere_arbre',  NbNoeuds+1,iPere_arbre)
      call XDR_Ecrit_ti
     $     (iFich_XDR,'iFilsG_arbre', NbNoeuds+1,iFilsG_arbre)
      call XDR_Ecrit_ti
     $     (iFich_XDR,'iFilsD_arbre', NbNoeuds+1,iFilsD_arbre)
      call XDR_Ecrit_ti
     $     (iFich_XDR,'iNbDecG_arbre',NbNoeuds+1,iNbDecG_arbre)
      call XDR_Ecrit_td
     $     (iFich_XDR,'dA_arbre',     NbNoeuds,  dA_arbre)
      call XDR_Ecrit_td
     $     (iFich_XDR,'dB_arbre',     NbNoeuds,  dB_arbre)
cc
      end
c
c----------------------------------------------------------------------
      subroutine LitArbre_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_ArbreIO_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich_XDR         ! Numero XDR du fichier (attribue par initxdr)
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      character*20 NomType,NomVar
      character*80 EnTete
      integer idum, ixdr_err
c
c fonctions :
c ^^^^^^^^^^^
      integer ixdr_read_string
c......................................................................
c
      Nb_Parcours_Arbre = 0
      ixdr_err = ixdr_read_string(iFich_XDR,EnTete)
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichArbre_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' d''arbre XDR valable !!!'
         call Meurt('LitArbre_XDR')
      end if
    
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.LT.0) goto 99

         if (NomVar.EQ.'NbSE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,NbSE)
         else if (NomVar.EQ.'iAdam') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,iAdam)

         else if (NomVar.EQ.'iSE_arbre') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbNoeuds)
            ixdr_err = ixdrimat(iFich_XDR,NbNoeuds,iSE_arbre)
         else if (NomVar.EQ.'iPere_arbre') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbNoeuds)
            NbNoeuds = NbNoeuds-1
            ixdr_err = ixdrimat(iFich_XDR,NbNoeuds+1,iPere_arbre)
         else if (NomVar.EQ.'iFilsG_arbre') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbNoeuds)
            NbNoeuds = NbNoeuds-1
            ixdr_err = ixdrimat(iFich_XDR,NbNoeuds+1,iFilsG_arbre)
         else if (NomVar.EQ.'iFilsD_arbre') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbNoeuds)
            NbNoeuds = NbNoeuds-1
            ixdr_err = ixdrimat(iFich_XDR,NbNoeuds+1,iFilsD_arbre)
         else if (NomVar.EQ.'iNbDecG_arbre') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbNoeuds)
            NbNoeuds = NbNoeuds-1
            ixdr_err = ixdrimat(iFich_XDR,NbNoeuds+1,iNbDecG_arbre)
         else if (NomVar.EQ.'dA_arbre') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbNoeuds)
            ixdr_err = ixdrdmat(iFich_XDR,NbNoeuds,dA_arbre)
         else if (NomVar.EQ.'dB_arbre') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbNoeuds)
            ixdr_err = ixdrdmat(iFich_XDR,NbNoeuds,dB_arbre)

         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitArbre_XDR')
         end if

      end do

 99   continue
      if (NbSE.GT.iDimSE) then
         write(0,*) '!!! NbSE(',NbSE,') > iDimSE(',iDimSE,') !!!'
         call Meurt('LitArbre_XDR')
      else if (NbNoeuds.GT.iDimArbre) then
         write(0,*) '!!! NbNoeuds(',NbNoeuds,') > iDimArbre(',
     $        iDimArbre,') !!!'
         call Meurt('LitArbre_XDR')
      end if
c
      end
c
#if (_ECRIRE_POT_>1)
c
c----------------------------------------------------------------------
      subroutine EcrirePotArbre(iFich)
c----------------------------------------------------------------------
c     
      include 'EvolAmas_ArbreIO_inc.f'
c
      integer iFich
c
      integer iNoeud
      double precision Pot_Amas_Arbre_i
c......................................................................
c
      write(iFich,100)
 100  format('# 1:iSE 2:R 3:Pot')
      do iNoeud=1,NbSE
         write(iFich,110)
     $        iSE_arbre(iNoeud),R_SE(iSE_arbre(iNoeud)),
     $        Pot_Amas_Arbre_i(iNoeud)
 110     format(I5,' ',2(E16.6,' '))
      end do
c
      end
c
#endif
c
#if (_HIER_ARBRE_>0)
c----------------------------------------------------------------------
      subroutine Determiner_Hierarchie(Tab_Hier)
c----------------------------------------------------------------------
c     place dans le tableau Tab_Hier, les niveaux hierarchique
c     de chaque noeud (avec un niveau de 0 pour "Adam" )
c
      include 'EvolAmas_ArbreIO_inc.f'
c
      integer Tab_Hier(iDimArbre)
c
      integer iNoeud, iNiv
c
      integer Pop_Pile
c......................................................................
c
c---- Le point de depart du parcours est la racine...
c
      call Init_Pile
      iNiv   = 0
      iNoeud = iAdam

 1    continue                  ! Debut de la boucle de parcourt
                                ! de l arbre
c
c---- on attribue le niveau hierarchique atteint au noeud courant
c
      Tab_Hier(iNoeud) = iNiv
c
c---- si il y a un fils a droite on le memorise sur la pile (ainsi
c     que son niveau hierarchique)
c
      if (iFilsD_arbre(iNoeud).GT.0) then
         call Push_Pile(iNiv+1)
         call Push_Pile(iFilsD_arbre(iNoeud))
      end if
c
c---- si il y a un fils a gauche, on descend la ligne de ce cote ...
c
      if (iFilsG_arbre(iNoeud).GT.0) then
         iNiv=iNiv+1
         iNoeud=iFilsG_arbre(iNoeud)
      else
c
c---- sinon, on prend le prochain element de la pile et on recommence
c     a moins qu'il n y ait plus rien sur la pile, auquel cas on
c     a termine
c
         iNoeud = Pop_Pile()
         if (iNoeud.LE.0) goto 9
         iNiv   = Pop_Pile()
      end if
      goto 1

 9    continue                  ! Termine !
c
      end 
c
c----------------------------------------------------------------------
      subroutine Ecrire_Hier_Asc(iFich,Tab_Hier)
c----------------------------------------------------------------------
c     ecrit les niveaux de hierarchie (tableau Tab_Hier) dans un
c     fichier ASCII
c
      include 'EvolAmas_ArbreIO_inc.f'
c
      integer iFich,Tab_Hier(iDimArbre)
c
      integer i
c......................................................................
c
      write(iFich,120)
 120  format(
     $     '#  TABLEAU DE HIERARCHIE D''ARBRE BINAIRE',/,
     $     '# ======================================='
     $     )
      write(iFich,130) 'iNoeud','Niveau'
 130  format(2(A8,'|'),/,2('--------|'))
      do i = 1,NbNoeuds
         write(iFich,140) i, Tab_Hier(i)
 140     format(2(I8,'|'))
      end do
c
      end   
c  
c----------------------------------------------------------------------
      subroutine EcrireArbreHier(iFich)
c----------------------------------------------------------------------
c     Ecriture ASCII de l'arbre 
c      
      include 'EvolAmas_ArbreIO_inc.f'
c
      integer Tab_Hier(iDimArbre)
      common /Hier/ Tab_Hier

      integer iFich
c
      integer i,iSE
c......................................................................
c
      write(iFich,120)NbSE,iAdam
 120  format(
     $     '#  TABLEAU DE STRUCTURE D''ARBRE BINAIRE',/,
     $     '# ======================================',/,
     $     '#  NbSE = ',I6,/,
     $     '#  iAdam = ',I6,/
     $     )
      write(iFich,130)
     $     'iNoeud','iSE','iPere','iFilsG','iFilsD','NbDec_G',
     $     'R','M','dA','dB','Niveau'
 130  format(11('|',A12),/,11('|------------'))
      do i = 1,NbNoeuds
         iSE=iSE_arbre(i)
         write(iFich,140)
     $        i,iSE,iPere_arbre(i),
     $        iFilsG_arbre(i),iFilsD_arbre(i),iNbDecG_arbre(i),
     $        R_SE(iSE),M_SE(iSE),dA_arbre(i),dB_arbre(i),
     $        Tab_Hier(i)
 140     format(6('|',I12),4('|',E12.4),('|',I12))
      end do
c
      end   
c  
#endif
c
