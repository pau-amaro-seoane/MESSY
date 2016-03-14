c======================================================================
c     routines d'entree/sortie pour les grandeurs de suivi de la
c     conservation de la masse et de l'energie de l'amas
c======================================================================
c
c----------------------------------------------------------------------
      subroutine EcrireConservAsc(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Conserv_IO_inc.f'
      include 'EvolAmas_Param_common.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
c     code genere parle script :
c
c     #! /usr/bin/ksh
c
c     cat EvolAmas_Conserv_inc.f | grep '\$' |\
c     sed 's/,/ /g; s/\$//;s/  */ /g' | tr ' ' '\012' | sort -u |\
c     grep -v '^ *$' |\
c     gawk '{
c         print "       write(iFich,~(A20,\" = \",E14.6)~) ~"$1"~,"$1
c     }' | tr '~' "'"
c
      write(iFich,100)
 100  format('# Grandeurs de suivi de la conservation de la masse',
     $     ' et de l''energie de l''amas')
      write(iFich,110) 'Etot_ini',Etot_ini
      write(iFich,110) 'Mamas_ini',Mamas_ini
      write(iFich,110) 'MSEtot_ini',MSEtot_ini
      write(iFich,110) 'Mtn_ini',Mtn_ini
      write(iFich,110) 'Ttot_ini',Ttot_ini
      write(iFich,110) 'Estell_grav_ini',Estell_grav_ini
      write(iFich,110) 'dEtot_coll',dEtot_coll
      write(iFich,110) 'dEtot_evst',dEtot_evst
      write(iFich,110) 'dEtot_dechir',dEtot_dechir
      write(iFich,110) 'dEtot_disp',dEtot_disp
      write(iFich,110) 'dEtot_captGW',dEtot_disp
      write(iFich,110) 'dEtot_evap',dEtot_evap
      write(iFich,110) 'dEtot_bin3c',dEtot_bin3c
      write(iFich,110) 'dEtot_nk',dEtot_nk
      write(iFich,110) 'dMaccr_coll',dMaccr_coll
      write(iFich,110) 'dMaccr_evst',dMaccr_evst
      write(iFich,110) 'dMaccr_dechir',dMaccr_dechir
      write(iFich,110) 'dMaccr_disp',dMaccr_disp
      write(iFich,110) 'dMaccr_captGW',dMaccr_captGW
      write(iFich,110) 'dMejec_coll',dMejec_coll
      write(iFich,110) 'dMejec_evst',dMejec_evst
      write(iFich,110) 'dMejec_dechir',dMejec_dechir
      write(iFich,110) 'dMejec_evap',dMejec_evap
      write(iFich,110) 'dTtot_coll',dTtot_coll
      write(iFich,110) 'dTtot_evst',dTtot_evst
      write(iFich,110) 'dTtot_dechir',dTtot_dechir
      write(iFich,110) 'dTtot_disp',dTtot_disp
      write(iFich,110) 'dTtot_captGW',dTtot_captGW
      write(iFich,110) 'dTtot_evap',dTtot_evap
      write(iFich,110) 'dTtot_nk',dTtot_nk
      write(iFich,110) 'rNbEtoiles_ini',rNbEtoiles_ini
      write(iFich,110) 'rNbEtoiles',rNbEtoiles
 110  format(A20,' = ',E14.6)
c
      end
c
c----------------------------------------------------------------------
      subroutine EcritCons_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Conserv_IO_inc.f'
      include 'EvolAmas_Param_common.f'
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
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichCons_XDR)
c
c---- code genere par le script :
c
c     #! /usr/bin/ksh
c
c     cat EvolAmas_Conserv_inc.f | grep '\$' |\
c     sed 's/,/ /g; s/\$//;s/  */ /g' | tr ' ' '\012' | sort -u |\
c     grep -v '^ *$' | gawk '{
c         print "      call XDR_Ecrit_d(iFich_XDR,~"$1"~,"$1")"
c     }' | tr '~' "'"
c
      call XDR_Ecrit_d(iFich_XDR,'Ecin_amas',Ecin_amas)
      call XDR_Ecrit_d(iFich_XDR,'Egrav_amas',Egrav_amas)
      call XDR_Ecrit_d(iFich_XDR,'Etot_amas',Etot_amas)
      call XDR_Ecrit_d(iFich_XDR,'Estell_grav_amas',Estell_grav_amas)
      call XDR_Ecrit_d(iFich_XDR,'Etot_ini',Etot_ini)
      call XDR_Ecrit_i(iFich_XDR,'NbSE_ini',NbSE_ini)
      call XDR_Ecrit_d(iFich_XDR,'Mamas_ini',Mamas_ini)
      call XDR_Ecrit_d(iFich_XDR,'MSEtot_ini',MSEtot_ini)
      call XDR_Ecrit_d(iFich_XDR,'Mtn_ini',Mtn_ini)
      call XDR_Ecrit_d(iFich_XDR,'Ttot_ini',Ttot_ini)
      call XDR_Ecrit_d(iFich_XDR,'Estell_grav_ini',Estell_grav_ini)
      call XDR_Ecrit_d(iFich_XDR,'dEtot_bin3c',dEtot_bin3c)
      call XDR_Ecrit_d(iFich_XDR,'dEtot_coll',dEtot_coll)
      call XDR_Ecrit_d(iFich_XDR,'dEtot_evst',dEtot_evst)
      call XDR_Ecrit_d(iFich_XDR,'dEtot_dechir',dEtot_dechir)
      call XDR_Ecrit_d(iFich_XDR,'dEtot_disp',dEtot_disp)
      call XDR_Ecrit_d(iFich_XDR,'dEtot_captGW',dEtot_captGW)
      call XDR_Ecrit_d(iFich_XDR,'dEtot_evap',dEtot_evap)
      call XDR_Ecrit_d(iFich_XDR,'dEtot_nk',dEtot_nk)
      call XDR_Ecrit_d(iFich_XDR,'dMaccr_coll',dMaccr_coll)
      call XDR_Ecrit_d(iFich_XDR,'dMaccr_evst',dMaccr_evst)
      call XDR_Ecrit_d(iFich_XDR,'dMaccr_dechir',dMaccr_dechir)
      call XDR_Ecrit_d(iFich_XDR,'dMaccr_disp',dMaccr_disp)
      call XDR_Ecrit_d(iFich_XDR,'dMaccr_captGW',dMaccr_captGW)
      call XDR_Ecrit_d(iFich_XDR,'dMejec_coll',dMejec_coll)
      call XDR_Ecrit_d(iFich_XDR,'dMejec_evst',dMejec_evst)
      call XDR_Ecrit_d(iFich_XDR,'dMejec_dechir',dMejec_dechir)
      call XDR_Ecrit_d(iFich_XDR,'dMejec_evap',dMejec_evap)
      call XDR_Ecrit_d(iFich_XDR,'dTtot_coll',dTtot_coll)
      call XDR_Ecrit_d(iFich_XDR,'dTtot_evst',dTtot_evst)
      call XDR_Ecrit_d(iFich_XDR,'dTtot_dechir',dTtot_dechir)
      call XDR_Ecrit_d(iFich_XDR,'dTtot_disp',dTtot_disp)
      call XDR_Ecrit_d(iFich_XDR,'dTtot_captGW',dTtot_captGW)
      call XDR_Ecrit_d(iFich_XDR,'dTtot_evap',dTtot_evap)
      call XDR_Ecrit_d(iFich_XDR,'dTtot_nk',dTtot_nk)
      call XDR_Ecrit_d(iFich_XDR,'rNbEtoiles_ini',rNbEtoiles_ini)
      call XDR_Ecrit_d(iFich_XDR,'rNbEtoiles',rNbEtoiles)
c
      end
c
c----------------------------------------------------------------------
      subroutine LitCons_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Conserv_IO_inc.f'
      include 'EvolAmas_Param_common.f'
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
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichCons_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' de conservation XDR valable !!!'
         call Meurt('LitCons_XDR')
      end if
    
      Estell_grav_ini = 1.0d30   ! Valeur-drapeau pour pouvoir verifier 
                                ! si elle a ete lue ou non

      MSEtot_ini = -1.0d0
      rNbEtoiles = -1.0d0

      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.LT.0) goto 99
c
c---- code genere par le script :
c
c     #! /usr/bin/ksh
c
c     cat EvolAmas_Conserv_inc.f | grep '\$' |\
c     sed 's/,/ /g; s/\$//;s/  */ /g' | tr ' ' '\012' | sort -u |\
c     grep -v '^ *$' | gawk '{
c         print "        else if (NomVar.EQ.~"$1"~) then"
c         print "           ixdr_err = ixdr_read_string(iFich_XDR,NomType)"
c         print "           ixdr_err = ixdrint(iFich_XDR,idum)"
c         print "           ixdr_err = ixdrdouble(iFich_XDR,"$1")"
c     }' | tr '~' "'"
c
         if (NomVar.EQ.'Ecin_amas') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Ecin_amas)
         else if (NomVar.EQ.'Egrav_amas') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Egrav_amas)
         else if (NomVar.EQ.'Etot_amas') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Etot_amas)
         else if (NomVar.EQ.'Estell_grav_amas') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Estell_grav_amas)
         else if (NomVar.EQ.'Etot_ini') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Etot_ini)
         else if (NomVar.EQ.'Mamas_ini') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Mamas_ini)
         else if (NomVar.EQ.'MSEtot_ini') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,MSEtot_ini)
         else if (NomVar.EQ.'NbSE_ini') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,NbSE_ini)
         else if (NomVar.EQ.'Mtn_ini') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Mtn_ini)
         else if (NomVar.EQ.'Ttot_ini') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Ttot_ini)
         else if (NomVar.EQ.'Estell_grav_ini') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Estell_grav_ini)
         else if (NomVar.EQ.'dEtot_bin3c') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dEtot_bin3c)
         else if (NomVar.EQ.'dEtot_coll') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dEtot_coll)
         else if (NomVar.EQ.'dEtot_evst') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dEtot_evst)
         else if (NomVar.EQ.'dEtot_dechir') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dEtot_dechir)
         else if (NomVar.EQ.'dEtot_disp') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dEtot_disp)
         else if (NomVar.EQ.'dEtot_captGW') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dEtot_captGW)
         else if (NomVar.EQ.'dEtot_evap') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dEtot_evap)
         else if (NomVar.EQ.'dEtot_nk') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dEtot_nk)
         else if (NomVar.EQ.'dMaccr_coll') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dMaccr_coll)
         else if (NomVar.EQ.'dMaccr_evst') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dMaccr_evst)
         else if (NomVar.EQ.'dMaccr_dechir') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dMaccr_dechir)
         else if (NomVar.EQ.'dMaccr_disp') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dMaccr_disp)
         else if (NomVar.EQ.'dMaccr_captGW') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dMaccr_captGW)
         else if (NomVar.EQ.'dMejec_coll') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dMejec_coll)
         else if (NomVar.EQ.'dMejec_evst') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dMejec_evst)
         else if (NomVar.EQ.'dMejec_dechir') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dMejec_dechir)
         else if (NomVar.EQ.'dMejec_evap') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dMejec_evap)
         else if (NomVar.EQ.'dTtot_coll') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dTtot_coll)
         else if (NomVar.EQ.'dTtot_evst') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dTtot_evst)
         else if (NomVar.EQ.'dTtot_dechir') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dTtot_dechir)
         else if (NomVar.EQ.'dTtot_disp') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dTtot_disp)
         else if (NomVar.EQ.'dTtot_captGW') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dTtot_captGW)
         else if (NomVar.EQ.'dTtot_evap') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dTtot_evap)
         else if (NomVar.EQ.'dTtot_nk') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,dTtot_nk)
         else if (NomVar.EQ.'rNbEtoiles_ini') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,rNbEtoiles_ini)
         else if (NomVar.EQ.'rNbEtoiles') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,rNbEtoiles)

         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitCons_XDR')
         end if

      end do

 99   continue
                                ! si la valeur de l'energie stellaire n'etait pas sauvegardee
                                ! on ne peut faire mieux que d'assumer que c'est la meme que 
                                ! l'energie totale (ceci est exact si il n'y a pas de TN).
      if (Estell_grav_ini.EQ.1.0d30) Estell_grav_ini = Etot_ini-Ttot_ini
      if (MSEtot_ini.LT.0.0d0) MSEtot_ini = 1.0d0
c
      end
c
c----------------------------------------------------------------------
      block data BD_Init_Conserv
c----------------------------------------------------------------------
      include 'EvolAmas_Conserv_IO_inc.f'
c
      data Mamas_ini,MSEtot_ini,Mtn_ini /1.0d0,1.0d0,0.0d0/
      data Etot_ini,Ttot_ini,Estell_grav_ini /-0.25d0,0.25d0,-0.5d0/
c
      end
