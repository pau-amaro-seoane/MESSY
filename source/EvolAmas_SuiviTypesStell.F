c======================================================================
c     routines pour Le suivi des differents types stellaires
c======================================================================
c----------------------------------------------------------------------
      subroutine Lire_FracM_SuiviTypesStell(iIn)
c----------------------------------------------------------------------
c     read the fractional masses for the Lagragian radii, in file iIn
c
      include 'EvolAmas_SuiviTypesStell_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iIn
c
c var locales :
c ^^^^^^^^^^^^^
      integer i
      double precision frac, frac_prec
c......................................................................
c
      N_RayLag_STS = 0
      frac_prec = 0.0d0
      do while(.TRUE.)
         read(iIn,*,ERR=99,END=99) frac
         if (frac.GT.1.0d0 .OR. frac.LE.0.0d0) then
            write(0,*) '!!! Erreur dans le fichier de fractions ',
     $           'lagrangiennes : les fractions doivent etre ',
     $           'dans ]0;1] !!!',frac
            call Meurt('Lire_FracM_RayLag')
         end if
         if (frac.LE.frac_prec) then
            write(0,*) '!!! Erreur dans le fichier de fractions ',
     $           'lagrangiennes : les fractions doivent etre ',
     $           'en ordre croissant !!!'
            call Meurt('Lire_FracM_SuiviTypesStell')
         end if
         if (N_RayLag_STS.EQ.N_RayLag_STS_Max) then
            write(0,*) '!!! Erreur dans le fichier de fractions ',
     $           'lagrangiennes : ',N_RayLag_STS_Max,
     $           ' fractions au maximum !!!'
            call Meurt('Lire_FracM_SuiviTypesStell')
         end if
         N_RayLag_STS=N_RayLag_STS+1
         f_RayLag_STS(N_RayLag_STS)=frac
         frac_prec = frac
      end do
         
 99   continue
      if (N_RayLag_STS.LT.1) then
         write(0,*) '!!! Erreur dans le fichier de fractions ',
     $        'lagrangiennes : aucune fraction specifiee !!!'
         call Meurt('Lire_FracM_SuiviTypesStell')
      end if
         
      do i=N_RayLag_STS+1,N_RayLag_STS_Max
         f_RayLag_STS(i) = -1.0d0
      end do
c
      end
c----------------------------------------------------------------------
      subroutine Init_SuiviTypesStell
c----------------------------------------------------------------------
      include 'EvolAmas_SuiviTypesStell_inc.f'

c
c var locales :
c ^^^^^^^^^^^^^
      integer i,iOut,iType, ifrac, ibase
      logical lExiste
      character*64 TypeFich
      integer Nformat
      character*256 format_str
      double precision form_multiplier
c......................................................................
c
      do iType=1,Nb_TypesEtoile

         dMaccr_evst_STS(iType)   = 0.0d0
         dMaccr_dechir_STS(iType) = 0.0d0
         dMaccr_disp_STS(iType)   = 0.0d0
         dMaccr_captGW_STS(iType) = 0.0d0
         dMaccr_coll_STS(iType)   = 0.0d0
         dMejec_evst_STS(iType)   = 0.0d0
         dMejec_dechir_STS(iType) = 0.0d0
         dMejec_coll_STS(iType)   = 0.0d0
         dMejec_evap_STS(iType)   = 0.0d0

         iOut=iFich_STS(iType) 
         TypeFich=TypeFich_STS(iType)
         call OuvrirFichASC(iOut,TypeFich(1:index(TypeFich,' ')-1),5,
     $        lExiste)

         if (lExiste) then
            write(0,*) ' > Fichier de suivi du Type stellaire ',iType,
     $           ' deja present'
         else
         
                                ! One forms the format (!) for the header 
                                ! so that there is enough information in the name of variables
                                ! to reconstruct the mass fractions
            Nformat = 3
            do i=1,N_RayLag_STS
               Nformat = max( Nformat,
     $             int(abs(log10(f_RayLag_STS(i))))+3 )
            end do
            form_multiplier = 10.0d0**(Nformat-1)
            write(format_str,200) (Nformat,i=1,12)
 200        format(
     $           '(" ",I3,": N",I',I1,'.',I1,
     $           '," ",I3,": M",I',I1,'.',I1,
     $           '," ",I3,": R",I',I1,'.',I1,
     $           '," ",I3,": Trad",I',I1,'.',I1,
     $           '," ",I3,": Ttrs",I',I1,'.',I1,
     $           '," ",I3,": Mavrg",I',I1,'.',I1,',$)')

            write(iOut,100) iType
 100        format(
     $           '# Suivi des SE de type stellaire ',I2,/,
     $           '# Les grandeures "dM" doivent etre accumulees',/,
     $           '# 1: iPas_Evol 2: Nombre 3: Masse',
     $           ' 4: dMaccr_evst 5: dMaccr_dechir 6: dMaccr_disp',
     $           ' 7: dMaccr_captGW 8: dMaccr_coll 9: dMejec_evst',
     $           ' 10: dMejec_dechir 11: dMejec_coll 12: dMejec_evap',
     $           $)
            do i=1,N_RayLag_STS
               ibase= 12 +(i-1)*6
               ifrac= NINT(form_multiplier*f_RayLag_STS(i))
               write(iOut,format_str)ibase+1,ifrac, ibase+2,ifrac,
     $              ibase+3,ifrac, ibase+4,ifrac, ibase+5,ifrac,
     $              ibase+6,ifrac
            end do
            write(iOut,*)
            call flush_MD(iOut)
         end if

      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine SuiviTypesStell
c----------------------------------------------------------------------
c     L'amas doit etre trie dans iListeRay.
c
      include 'EvolAmas_SuiviTypesStell_inc.f'
      include 'EvolAmas_Compteurs_common.f'
c
c var locales :
c ^^^^^^^^^^^^^
      integer lSE, iSE, icouche, iType, iOut, Npart, Ntot
      double precision Msomme, Mlim, Mtot, Ttrs, Ttrs_i, Trad,
     $     Mstell_avg
c
c functions :
c ^^^^^^^^^^^
      integer*1 stellar_generaltype
      double precision stellar_mass
c......................................................................
c
      do iType=1,Nb_TypesEtoile
         M_STS(iType)=0.0d0
         N_STS(iType)=0
      end do
      do lSE=1,NbSE-NbSE_perte
         iSE=iListeRay(lSE)
         iType=stellar_generaltype(iSE)
         M_STS(iType)=M_STS(iType)+M_SE(iSE)
         N_STS(iType)=N_STS(iType)+1
      end do
      
      do iType=1,Nb_TypesEtoile

         iOut =iFich_STS(iType)  
         write(iOut,100) iPas_evol, N_STS(iType), M_STS(iType),
     $        dMaccr_evst_STS(iType),dMaccr_dechir_STS(iType),
     $        dMaccr_disp_STS(iType),dMaccr_captGW_STS(iType),
     $        dMaccr_coll_STS(iType),
     $        dMejec_evst_STS(iType),dMejec_dechir_STS(iType),
     $        dMejec_coll_STS(iType),dMejec_evap_STS(iType)
 100     format(I10,' ',I7,' ',10(E14.5,' '),$)

         Mtot =M_STS(iType)
         Ntot =N_STS(iType)
         Mlim = 0.0d0
         Msomme = 0.0d0
         Trad = 0.0d0
         Ttrs = 0.0d0
         Mstell_avg = 0.0d0
         Npart = 0
         lSE = 0

         if (Mtot.GT.0.0d0) then
            do icouche=1,N_RayLag_STS
               Mlim = min(f_RayLag_STS(icouche),0.9999d0)*Mtot
               do while (Msomme.LT.Mlim .AND. Npart.LT.Ntot)
                  lSE = lSE + 1
                  iSE = iListeRay(lSE)
                  if (stellar_generaltype(iSE).EQ.iType) then
                     Npart = Npart+1
                     Msomme = Msomme + M_SE(iSE)
                     Ttrs_i = 0.5d0*(J_SE(iSE)/R_SE(iSE))**2
                     Trad = Trad +  M_SE(iSE)*(T_SE(iSE)-Ttrs_i)
                     Ttrs = Ttrs + M_SE(iSE)*Ttrs_i
                     Mstell_avg = Mstell_avg + stellar_mass(iSE)
                  end if
               end do
               if (Npart.GT.0) then
                  write(iOut,110) Npart, Msomme, R_SE(iSE), 
     $                 Trad/Msomme, Ttrs/Msomme, Mstell_avg/dble(Npart)
               else
                  write(iOut,110) 0, 0.0d0, -1.0d0, 
     $                 -1.0d0, -1.0d0, -1.0d0
               end if
 110           format(I8,' ',5(E14.5,' '),$)
            end do
         else                   ! Il n'y a plus d'etoiles du type voulu
            do icouche=1,N_RayLag_STS
               write(iOut,110) 0, 0.0d0, -1.0d0, 
     $              -1.0d0, -1.0d0, -1.0d0
            end do
         end if
         write(iOut,*)
         call flush_MD(iOut)

         dMaccr_evst_STS(iType)   = 0.0d0
         dMaccr_dechir_STS(iType) = 0.0d0
         dMaccr_disp_STS(iType)   = 0.0d0
         dMaccr_captGW_STS(iType) = 0.0d0
         dMaccr_coll_STS(iType)   = 0.0d0
         dMejec_evst_STS(iType)   = 0.0d0
         dMejec_dechir_STS(iType) = 0.0d0
         dMejec_coll_STS(iType)   = 0.0d0
         dMejec_evap_STS(iType)   = 0.0d0

      end do                  
c
      end
c
c----------------------------------------------------------------------
      block data BD_Init_SuiviTypesStell
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_SuiviTypesStell_common.f'
      include 'EvolAmas_param.f'
c
c---- Fichiers de suivi
c
      data iFich_STS
     $     /iFich_STS1,iFich_STS2,iFich_STS3,iFich_STS4,iFich_STS5/
      data TypeFich_STS
     $     /'SuiviMS','SuiviRG','SuiviWD','SuiviNS','SuiviBH'/
c
c---- valeurs standard des limites en masses pour les rayons lagrangiens
c
      data N_RayLag_STS /7/
      data f_RayLag_STS /0.02,0.05,0.1,0.25,0.5,0.75,0.9,93*-1.0d0/
c
      end
c----------------------------------------------------------------------
