c======================================================================
c     routines pour lire et ecrire des fichiers d'amas en employant
c     le "format binaire universel" XDR
c======================================================================

c----------------------------------------------------------------------
      subroutine EcritAmas_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'AmasIO_XDR_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich_XDR         ! Numero XDR du fichier (attribue par initxdr)
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichAmas_XDR)
c
      call XDR_Ecrit_td(iFich_XDR,'T_SE',NbSE,T_SE)
      call XDR_Ecrit_td(iFich_XDR,'J_SE',NbSE,J_SE)
      call XDR_Ecrit_td(iFich_XDR,'M_SE',NbSE+1,    M_SE) ! pour des raisons historiques, ce tableau a un element 0
      call XDR_Ecrit_td(iFich_XDR,'R_SE',max(NbNoeuds,NbSE)+1,R_SE) ! pour des raisons historiques, ce tableau a un element 0
c
      end
c
c----------------------------------------------------------------------
      integer function NbSE_FichAmasXDR(NomFich)
c----------------------------------------------------------------------
      implicit none
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      character*(*) NomFich
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      character*80 Chaine
      integer iXDR,ixdr_err
      integer N
      logical exist
c
c fonctions :
c ^^^^^^^^^^^
      integer ixdr_read_string
c......................................................................
c
      inquire(file=NomFich,exist=exist)
      if (.NOT.exist) then
         NbSE_FichAmasXDR = -1
         return
      end if

      iXDR = initxdr(NomFich,'r',.TRUE.)
      ixdr_err = ixdr_read_string(iXDR,Chaine) ! lecture de l'en-tete (on s'en fout)
      ixdr_err = ixdr_read_string(iXDR,Chaine)
      if (Chaine.EQ.'T_SE' .OR. Chaine.EQ.'J_SE') then
         ixdr_err = ixdr_read_string(iXDR,Chaine) ! lecture du nom du type de variable (on s'en fout)
         ixdr_err = ixdrint(iXDR,N)
         NbSE_FichAmasXDR = N 
      else
         NbSE_FichAmasXDR = -2
      end if
      ixdr_err = ixdrclose(iXDR)

      end
c
c----------------------------------------------------------------------
      subroutine LitAmas_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'AmasIO_XDR_inc.f'
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
      integer N,iNoeud,ixdr_err
c
c fonctions :
c ^^^^^^^^^^^
      integer ixdr_read_string
c......................................................................
c
      ixdr_err = ixdr_read_string(iFich_XDR,EnTete)

      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichAmas_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     +        ' d''amas XDR valable !!!'
         call Meurt('LitAmas_XDR')
      end if
    
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         !call EnleverBlancs(NomVar)
         if (ixdr_err.LT.0 .OR. NomVar.EQ.'') goto 99
         !if (ixdr_err.LT.0) goto 99
         ixdr_err = ixdr_read_string(iFich_XDR,NomType)

         if (NomType.NE.NomType_td) then
            write(0,*) '!!! Variable du mauvais type dans',
     $           ' le fichier d''amas XDR (',NomVar,':',NomType,') !!!'
            call Meurt('LitAmas_XDR')
         end if

         if (NomVar.EQ.'T_SE') then
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            if (NbSE.GT.iDimSE) then
               write(0,*) '!!! NbSE(',NbSE,') > iDimSE(',iDimSE,') !!!'
               call Meurt('LitAmas_XDR')
            end if
            ixdr_err = ixdrdmat(iFich_XDR,NbSE,T_SE)
         else if (NomVar.EQ.'J_SE') then
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            if (NbSE.GT.iDimSE) then
               write(0,*) '!!! NbSE(',NbSE,') > iDimSE(',iDimSE,') !!!'
               call Meurt('LitAmas_XDR')
            end if
            ixdr_err = ixdrdmat(iFich_XDR,NbSE,J_SE)
         else if (NomVar.EQ.'M_SE') then
            ixdr_err = ixdrint(iFich_XDR,N)
            if (N.GT.iDimArbre+1) then
               write(0,*) '!!! NbNoeuds(',N-1,') > iDimArbre(',
     $              iDimArbre,') !!!'
               call Meurt('LitAmas_XDR')
            end if
            ixdr_err = ixdrdmat(iFich_XDR,N,M_SE)
         else if (NomVar.EQ.'R_SE') then
            ixdr_err = ixdrint(iFich_XDR,NbNoeuds)
            if (NbNoeuds.GT.iDimArbre+1) then
               write(0,*) '!!! NbNoeuds(',N-1,') > iDimArbre(',
     $              iDimArbre,') !!!'
               call Meurt('LitAmas_XDR')
            end if
            ixdr_err = ixdrdmat(iFich_XDR,NbNoeuds,R_SE)
            NbNoeuds = NbNoeuds-1 ! pour des raisons historiques, ce tableau a un element 0
         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitAmas_XDR')
         end if

      end do

 99   continue
                                ! on met a 0 la masse des noeuds vides
      do iNoeud=NbSE+1,NbNoeuds
         M_SE(iNoeud) = 0.0d0
      end do
                                ! on determine certaines grandeurs globales 
      call DetGrandGlobAmas     ! indispensables
c     
      end
c
c----------------------------------------------------------------------
      subroutine LitAmas_UNIV(NomFich) ! Lit un fichier d'amas au format sun ou xdr
c----------------------------------------------------------------------
      include 'AmasIO_XDR_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      character*(*) NomFich
c
c constante :
c ^^^^^^^^^^^
      integer iFich
      parameter (iFich=55)
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer iXDR,l,i, ixdr_err
      double precision dum
c
c fonctions :
c ^^^^^^^^^^^
      integer LongChaine
c......................................................................
c
      l=LongChaine(NomFich)
      if (NomFich(l-3:l).EQ.'.xdr') then
         iXDR = initxdr(NomFich,'r',.TRUE.)
         call LitAmas_XDR(iXDR)
         ixdr_err = ixdrclose(iXDR)
      else
         open(iFich,file=NomFich,status='old',form='unformatted')
         read(iFich)
     $        NbSE,NbNoeuds,M_amas,dum,
     $        (T_SE(i),i=1,NbSE),(J_SE(i),i=1,NbSE),
     $        (M_SE(i),i=1,NbNoeuds),(R_SE(i),i=1,NbNoeuds)
         R_SE(0) = 0.0
         close(iFich)
         call DetGrandGlobAmas
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine DetGrandGlobAmas
c----------------------------------------------------------------------
      include 'AmasIO_XDR_inc.f'
      include 'EvolAmas_Pertes_inc.f' ! comptabilisation des SE perdues
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i
c......................................................................
c
      M_amas         = 0.0d0
      B0             = 0.0d0
      NbSE_perte     = 0
      NbSE_evap      = 0
      NbSE_dest_coll = 0
      NbSE_dest_dechir = 0
      NbSE_dest_disp = 0
      NbSE_dest_captGW = 0
      do i=1,NbSE
         if (R_SE(i).LT.R_Perte) then ! SE "normale"
            M_amas = M_amas + M_SE(i)
            B0 = B0 + M_SE(i)/R_SE(i) 
         else                         ! SE "perdue"
            NbSE_perte = NbSE_perte + 1
            if (R_SE(i).LT.R_dest_coll) then ! Evaporee
               NbSE_evap = NbSE_evap + 1
            else if  (R_SE(i).LT.R_dechir) then ! detruite par collision
               NbSE_dest_coll = NbSE_dest_coll + 1
            else if  (R_SE(i).LT.R_disp) then ! dechiree par le TN
               NbSE_dest_dechir = NbSE_dest_dechir + 1
            else if  (R_SE(i).LT.R_captGW) then ! disparue a travers l'horizon
               NbSE_dest_disp = NbSE_dest_disp + 1
            else                                ! captee sur orbite a forte emissions d'ondes grav.
               NbSE_dest_captGW = NbSE_dest_captGW + 1
            end if
         end if
      end do
c
      end
c
