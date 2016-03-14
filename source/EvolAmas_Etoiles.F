c======================================================================
c     routines liees a la "gestion" des proprietes stellaires dans
c     amas a la Henon
c======================================================================

#if (_ETOILES_IO_OLD_>0)
c----------------------------------------------------------------------
      subroutine LireEtoilesBIN_old(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
c
c var locales : 
c ^^^^^^^^^^^^^
      integer i
      double precision dummy_NbEt
c......................................................................
c
c---- lecture de l'enregistrement
c
      read(iFich,end=999,err=999) NbSE,dummy_NbEt,
     $     (Met_SE(i),i=1,NbSE),Coef_ConvR,(Ret_SE(i),i=1,NbSE)
c
c     Coef_ConvR permet de passer des unites de distance en rayon solaire au unites
c     de distance de l'amas. Ce coefficient, qui peut etre choisi librement, fixe
c     donc les dimensions physiques de l'amas ce qui n'a d'importance que lorsque
c     la taille finie des etoiles entre en jeu
c
c     dummy_NbEt est inutile mais est place la pour garder la compatibilite avec le format,
c     initialement defini avec le nb d'etoiles
c        
      call DetGrandStellGlob

      return
c
 999  continue
      write(0,*) '!!! Erreur lors de la lecture du fichier',
     $     ' de donnees stellaires !!!'
      call Meurt('LireEtoilesBIN')
c
      end
c
#endif /* (_ETOILES_IO_OLD_>0) */
c----------------------------------------------------------------------
      subroutine LireEtoilesBIN(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
c
c var locales : 
c ^^^^^^^^^^^^^
      integer i
c......................................................................
c
c---- lecture de l'enregistrement
c
      read(iFich,end=999,err=999) NbSE,
     $     (Met_SE(i),i=1,NbSE),(DNet_SE(i),i=1,NbSE),
     $     (iTet_SE(i),i=1,NbSE)
c        
      return
c
 999  continue
      write(0,*) '!!! Erreur lors de la lecture du fichier',
     $     ' de donnees stellaires !!!'
      call Meurt('LireEtoilesBIN')
c
      end
c
#if (_ETOILES_IO_OLD_>0)
c----------------------------------------------------------------------
      subroutine EcrireEtoilesBIN_old(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
c
c var locales : 
c ^^^^^^^^^^^^^
      integer i
      double precision dummy_NbEt
      data dummy_NbEt /0.0d0/
c......................................................................
c
      call DetGrandStellGlob

      write(iFich) NbSE,dummy_NbEt,
     $     (Met_SE(i),i=1,NbSE),Coef_ConvR,(Ret_SE(i),i=1,NbSE)
c
      end
c
#endif /* (_ETOILES_IO_OLD_>0) */
c----------------------------------------------------------------------
      subroutine EcrireEtoilesBIN(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
c
c var locales : 
c ^^^^^^^^^^^^^
      integer i
c......................................................................
c
      write(iFich) NbSE,
     $     (Met_SE(i),i=1,NbSE),(DNet_SE(i),i=1,NbSE),
     $     (iTet_SE(i),i=1,NbSE)
c
      end
c
#if (_ETOILES_IO_OLD_>0)
c----------------------------------------------------------------------
      subroutine EcrireEtoilesASC_old(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
c
c var locales : 
c ^^^^^^^^^^^^^
      integer i
c......................................................................
c
      write(iFich,100) NbSE,NEt_par_SE
 100  format(
     $     '# Donnes Stellaires d''amas',/,
     $     '# NbSE = ',I6,/,
     $     '# NbEt_par_SE = ',E12.4,/,
     $     '# 1: i 2: Met_SE 3: Ret_SE')
      do i=1,NbSE
         write(iFich,110) i,Met_SE(i),Ret_SE(i)
 110     format(I6,2(' ',E12.4))
      end do
c
      end
c
#endif /* (_ETOILES_IO_OLD_>0) */
c----------------------------------------------------------------------
      subroutine EcrireEtoilesASC(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
c
c var locales : 
c ^^^^^^^^^^^^^
      integer i
c......................................................................
c
      write(iFich,100) NbSE,NEt_par_SE
 100  format(
     $     '# Donnes Stellaires d''amas',/,
     $     '# NbSE = ',I6,/,
     $     '# NbEt_par_SE = ',E12.4,/,
     $     '# 1: i 2: Met_SE 3: DNet_SE 4: iTet_SE')
      do i=1,NbSE
         write(iFich,110) i,Met_SE(i),DNet_SE(i),iTet_SE(i)
 110     format(I7,2(' ',E12.4),' ',I3)
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine EcrireEtoilesRDB(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
c
c var locales : 
c ^^^^^^^^^^^^^
      integer i
      character*1 t_            ! tabulation
c......................................................................
c
      t_ = CHAR(9)
c
c--- ecriture en-tete
c
      write(iFich,100) 'Met_SE',t_,'DNet_SE',t_,'iTet_SE'
      write(iFich,100) '------',t_,'-------',t_,'-------'
 100  format(2(A,A1),A)
      do i=1,NbSE
         write(iFich,110) Met_SE(i),t_,DNet_SE(i),t_,iTet_SE(i)
 110     format(2(E15.8,A1),I3)
      end do
c
      end
c
#if (_ETOILES_IO_OLD_>0)
c----------------------------------------------------------------------
      subroutine EcritEtoiles_XDR_old(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich_XDR         ! Numero XDR du fichier (attribue par initxdr)
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      call DetGrandStellGlob

      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichEtoiles_XDR)

      call XDR_Ecrit_td(iFich_XDR,'Met_SE',NbSE,Met_SE)
      call XDR_Ecrit_td(iFich_XDR,'Ret_SE',NbSE,Ret_SE)
c
      end
c
#endif /* (_ETOILES_IO_OLD_>0) */
#if (_ETOILES_IO_OLD_>0)
c----------------------------------------------------------------------
      subroutine LitEtoiles_XDR_old(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
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
      double precision fdum
c
c fonction :
c ^^^^^^^^^^
      integer ixdr_read_string
c......................................................................
c
      ixdr_err = ixdr_read_string(iFich_XDR,EnTete)
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichEtoiles_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' d''etoiles XDR valable !!!'
         call Meurt('LitTemps_XDR')
      end if
    
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.NE.0) goto 99

         if (NomVar.EQ.'Net_par_SE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Net_par_SE)
         else if (NomVar.EQ.'Coef_ConvR') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Coef_ConvR)
         else if (NomVar.EQ.'Coef_ConvM') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Coef_ConvM)
         else if (NomVar.EQ.'Met_SE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            ixdr_err = ixdrdmat(iFich_XDR,NbSE,Met_SE)
         else if (NomVar.EQ.'Ret_SE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            ixdr_err = ixdrdmat(iFich_XDR,NbSE,Ret_SE)
         else if (NomVar.EQ.'rNbEtoiles') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,fdum)
         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitEtoiles_XDR')
         end if

      end do

 99   continue
      if (NbSE.GT.iDimSE) then
         write(0,*) '!!! NbSE(',NbSE,') > iDimSE(',iDimSE,') !!!'
         call Meurt('LitEtoilesE_XDR')
      end if
      call DetGrandStellGlob
c
      end
c
#endif /* (_ETOILES_IO_OLD_>0) */
c----------------------------------------------------------------------
      subroutine EcritEtoiles_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich_XDR         ! Numero XDR du fichier (attribue par initxdr)
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err, iType(iDimSE), i
c......................................................................
c
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichEtoiles_XDR)

      call XDR_Ecrit_td(iFich_XDR,'Met_SE',NbSE,Met_SE)
      call XDR_Ecrit_td(iFich_XDR,'DNet_SE',NbSE,DNet_SE)
                                ! passage integer*1 -> integer (fxdr ne connait pas le type integer*1)
      do i=1,NbSE
         iType(i) = iTet_SE(i)
      end do
      call XDR_Ecrit_ti(iFich_XDR,'iTet_SE',NbSE,iType)
c
      end
c
c----------------------------------------------------------------------
      subroutine LitEtoiles_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
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
      integer ixdr_err, iType(iDimSE), i
c
c fonction :
c ^^^^^^^^^^
      integer ixdr_read_string
c......................................................................
c
      ixdr_err = ixdr_read_string(iFich_XDR,EnTete)
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichEtoiles_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' d''etoiles XDR valable !!!'
         call Meurt('LitTemps_XDR')
      end if
    
      do while (.TRUE.)
         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.NE.0) goto 99
         if (NomVar.EQ.'Met_SE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            ixdr_err = ixdrdmat(iFich_XDR,NbSE,Met_SE)
         else if (NomVar.EQ.'DNet_SE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            ixdr_err = ixdrdmat(iFich_XDR,NbSE,DNet_SE)
         else if (NomVar.EQ.'iTet_SE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            ixdr_err = ixdrimat(iFich_XDR,NbSE,iType)
                                ! passage integer -> integer*1 (fxdr ne connait pas le type integer*1)
            do i=1,NbSE
               iTet_SE(i) = iType(i)
            end do
         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitEtoiles_XDR')
         end if

      end do

 99   continue
      if (NbSE.GT.iDimSE) then
         write(0,*) '!!! NbSE(',NbSE,') > iDimSE(',iDimSE,') !!!'
         call Meurt('LitEtoilesE_XDR')
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine LitEtoiles_XDR_tol(iFich_XDR)
c----------------------------------------------------------------------
c     version tolerante admettant les variables inconnues
c
      include 'EvolAmas_Etoiles_inc.f'
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
      integer ixdr_err, iType(iDimSE), i
c
c fonction :
c ^^^^^^^^^^
      integer ixdr_read_string
c......................................................................
c
      ixdr_err = ixdr_read_string(iFich_XDR,EnTete)
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichEtoiles_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' d''etoiles XDR valable !!!'
         call Meurt('LitTemps_XDR')
      end if
    
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.NE.0) goto 99
         if (NomVar.EQ.'Met_SE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            ixdr_err = ixdrdmat(iFich_XDR,NbSE,Met_SE)
         else if (NomVar.EQ.'DNet_SE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            ixdr_err = ixdrdmat(iFich_XDR,NbSE,DNet_SE)
         else if (NomVar.EQ.'iTet_SE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            ixdr_err = ixdrimat(iFich_XDR,NbSE,iType)
                                ! passage integer -> integer*1 (fxdr ne connait pas le type integer*1)
            do i=1,NbSE
               iTet_SE(i) = iType(i)
            end do
         end if

      end do

 99   continue
      if (NbSE.GT.iDimSE) then
         write(0,*) '!!! NbSE(',NbSE,') > iDimSE(',iDimSE,') !!!'
         call Meurt('LitEtoilesE_XDR')
      end if
c
      end
c

#if (_ETOILES_IO_OLD_>0)
c----------------------------------------------------------------------
      subroutine DetGrandStellGlob_old
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_Param_common.f'
c
c var locales : 
c ^^^^^^^^^^^^^
      integer i
      double precision Mtot_en_Msol,Mtot
c......................................................................
                                ! determination du coefficient de conversion
                                ! unites de masse stellaire -> unites de masse de l'amas
      Mtot = 0.0d0
      Mtot_en_Msol = 0.0d0
      if (rNbEtoiles.GT.0.0d0) then
         Net_par_SE = rNbEtoiles/NbSE
      else
         rNbEtoiles = NbSE*Net_par_SE
      end if
      do i=1,NbSE
         Mtot = Mtot + M_SE(i)
         Mtot_en_Msol =  Mtot_en_Msol + Net_par_SE*Met_SE(i)
      end do
      if (Mtot.EQ.0.0d0) Mtot=1.0d0
      Coef_ConvM = Mtot/Mtot_en_Msol 
c
      end
c
#endif /* (_ETOILES_IO_OLD_>0) */
