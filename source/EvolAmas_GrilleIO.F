c======================================================================
c     routines d'entree-sorties pour grilles d'amas
c======================================================================

c----------------------------------------------------------------------
      subroutine EcrireGrilleASC (iOut)
c----------------------------------------------------------------------
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iOut              ! num d`unite
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer i
c......................................................................
c
      write(iOut,100)N_PG,rNSEc_PG,Nmin_PG,Nmax_PG
 100  format(
     $     '# Profil de densitee grossier pour "EvolAmasXXX"',/,
     $     '# Nb de cellules : ',I6,/,
     $     '# Valeur moy,min,max de SE par cellule : ',F8.2,2(' ',I5),/,
     $     '# 1:Rsup_couche 2:M_couche 3:V_couche 4:T_couche ',
     $     '5:Ttg_couche 6:MR_couche 7:NSE_couche')
      write(iOut,110) Rsup_PG(0),0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0
      do i=1,N_PG
         write(iOut,110) Rsup_PG(i),M_PG(i),
     $        pi4_3*(Rsup_PG(i)**3-Rsup_PG(i-1)**3),
     $        T_PG(i),Tt_PG(i),MR_PG(i),NSE_PG(i)
      end do
 110  format(6(E12.5,' '),I5)
c
      end
c
c----------------------------------------------------------------------
      subroutine EcrireGrilleBIN (iOut)
c----------------------------------------------------------------------
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iOut              ! num d`unite
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer i
c......................................................................
c
      write(iOut) N_PG,
     $     (Rsup_PG(i),i=0,N_PG),(M_PG(i)  ,i=1,N_PG),
     $     (V_PG(i)   ,i=1,N_PG),(T_PG(i)  ,i=1,N_PG),
     $     (MR_PG(i)  ,i=1,N_PG),(NSE_PG(i),i=1,N_PG),
c$$$#if (_CHAUFF_BIN_3C_>0)
c$$$     $     (Y3c_PG(i),i=1,N_PG),
c$$$#endif
     $     Nmin_PG,Nmax_PG
c
      end
c
c----------------------------------------------------------------------
      subroutine EcrireGrilleBIN2 (iOut)
c----------------------------------------------------------------------
c     version "econome" n'ecrivant que les limites
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iOut              ! num d`unite
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer i
c......................................................................
c
      write(iOut) N_PG,(Rsup_PG(i),i=0,N_PG)
c
      end
c
c----------------------------------------------------------------------
      subroutine LireGrilleBIN (iIn)
c----------------------------------------------------------------------
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iIN              ! num d`unite
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer i
c......................................................................
c
      read(iIn) N_PG,
     $     (Rsup_PG(i),i=0,N_PG),(M_PG(i)  ,i=1,N_PG),
     $     (V_PG(i)   ,i=1,N_PG),(T_PG(i)  ,i=1,N_PG),
     $     (MR_PG(i)  ,i=1,N_PG),(NSE_PG(i),i=1,N_PG),
c$$$#if (_CHAUFF_BIN_3C_>0)
c$$$     $     (Y3c_PG(i),i=1,N_PG),
c$$$#endif
     $     Nmin_PG,Nmax_PG
c
      end
c
c----------------------------------------------------------------------
      subroutine LireGrilleBIN2(iIn)
c----------------------------------------------------------------------
c     version "econome" se contentant des limites radiales
c     il faut que l'amas soit deja charge et trie !
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iIN              ! num d`unite
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer i,iSE,iPG,NSEpg
      double precision Rsup,Mpg,Tpg,MRpg
c......................................................................
c
      read(iIn) N_PG,(Rsup_PG(i),i=0,N_PG)
c
c---- il faut a present faire les sommes des diverses quantites
c     dans les cellules
c
      i=1
      do iPG=1,N_PG
         iSE = iListeRay(i)
         Rsup = Rsup_PG(iPG)
         Mpg   = 0.0d0
         Tpg   = 0.0d0
         MRpg  = 0.0d0
c$$$         Ypg   = 0.0d0
         NSEpg = 0
         do while ( R_SE(iSE).LE.Rsup .AND. i.LE.(NbSE-NbSE_perte) )
            Mpg   = Mpg   + M_SE(iSE)
            Tpg   = Tpg   + T_SE(iSE)
            MRpg  = MRpg  + M_SE(iSE)*R_SE(iSE)
            NSEpg = NSEpg + 1
c$$$#if (_CHAUFF_BIN_3C_>0)
c$$$            Ypg   = Ypg  + stellar_mass(iSE)**2/sqrt(2.0d0*T_SE(iSE))**3
c$$$#endif
         end do
         M_PG(iPG)   = Mpg
         T_PG(iPG)   = Tpg
         MR_PG(iPG)  = MRpg
         NSE_PG(iPG) = NSEpg
c$$$         Y3c_PG(iPG) = Ypg
         V_PG(iPG)   = Pi4_3*( Rsup**3 - Rsup_PG(iPG-1)**3 )
         i=i+1
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine EcritGrille_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Grille_inc.f'
      include 'fxdr.inc'
c
c argument :
c ^^^^^^^^^^
      integer iFich_XDR
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichGrille_XDR)
c
      call XDR_Ecrit_i(iFich_XDR,'Nmin_PG',Nmin_PG)
      call XDR_Ecrit_i(iFich_XDR,'Nmax_PG',Nmax_PG)
      call XDR_Ecrit_td(iFich_XDR,'M_PG',N_PG,M_PG)
      call XDR_Ecrit_td(iFich_XDR,'MR_PG',N_PG,MR_PG)
      call XDR_Ecrit_td(iFich_XDR,'T_PG',N_PG,T_PG)
      call XDR_Ecrit_td(iFich_XDR,'Rsup_PG',N_PG+1,Rsup_PG)
      call XDR_Ecrit_ti(iFich_XDR,'NSE_PG',N_PG,NSE_PG)
c
      end
c
c----------------------------------------------------------------------
      subroutine EcritGrille_XDR2(iFich_XDR,
     $     N_Grille,
     $     Rsup_Grille,M_Grille,MR_Grille,T_Grille,Tt_Grille,NSE_Grille)
c----------------------------------------------------------------------
      include 'EvolAmas_Grille_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich_XDR
      integer N_Grille,
     $     NSE_Grille(N_Grille)
      double precision
     $     Rsup_Grille(0:N_Grille),M_Grille(N_Grille),
     $     MR_Grille(N_Grille),T_Grille(N_Grille),Tt_Grille(N_Grille)
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichGrille_XDR)
c
      call XDR_Ecrit_td(iFich_XDR,'M_PG',N_Grille,M_Grille)
      call XDR_Ecrit_td(iFich_XDR,'MR_PG',N_Grille,MR_Grille)
      call XDR_Ecrit_td(iFich_XDR,'T_PG',N_Grille,T_Grille)
      call XDR_Ecrit_td(iFich_XDR,'Tt_PG',N_Grille,Tt_Grille)
      call XDR_Ecrit_td(iFich_XDR,'Rsup_PG',N_Grille+1,Rsup_Grille)
      call XDR_Ecrit_ti(iFich_XDR,'NSE_PG',N_Grille,NSE_Grille)
c
      end
c
c----------------------------------------------------------------------
      subroutine LitGrille_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Grille_inc.f'
      include 'fxdr.inc'
c
c argument :
c ^^^^^^^^^^
      integer iFich_XDR
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      character*20 NomType,NomVar
      character*80 EnTete
      integer idum, ixdr_err
c
c fonction :
c ^^^^^^^^^^
      integer ixdr_read_string
c......................................................................
c
      ixdr_err = ixdr_read_string(iFich_XDR,EnTete)
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichGrille_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' de grille XDR valable !!!'
         call Meurt('LitGrille_XDR')
      end if
    
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.LT.0) goto 99

         if (NomVar.EQ.'M_PG') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,N_PG)
            ixdr_err = ixdrdmat(iFich_XDR,N_PG,M_PG)
         else if (NomVar.EQ.'MR_PG') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,N_PG)
            ixdr_err = ixdrdmat(iFich_XDR,N_PG,MR_PG)
         else if (NomVar.EQ.'T_PG') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,N_PG)
            ixdr_err = ixdrdmat(iFich_XDR,N_PG,T_PG)
         else if (NomVar.EQ.'Tt_PG') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,N_PG)
            ixdr_err = ixdrdmat(iFich_XDR,N_PG,Tt_PG)
         else if (NomVar.EQ.'Rsup_PG') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,N_PG)
            N_PG=N_PG-1
            ixdr_err = ixdrdmat(iFich_XDR,N_PG+1,Rsup_PG)
         else if (NomVar.EQ.'NSE_PG') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,N_PG)
            ixdr_err = ixdrimat(iFich_XDR,N_PG,NSE_PG)
         else if (NomVar.EQ.'Nmin_PG') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,Nmin_PG)
         else if (NomVar.EQ.'Nmax_PG') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,Nmax_PG)

         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitGrille_XDR')
         end if

      end do

 99   continue
      if (N_PG.GT.iDimPG) then
         write(0,*) '!!! N_PG(',N_PG,') > iDimPG(',iDimPG,') !!!'
         call Meurt('LitGrille_XDR')
      end if
c
      end
c
