#define _ECRIRE_TEST_TO_ 0
#define _VERB_TO_ 1
#define _TEST_RTSAFE_ 0
#define _NO_REFINEMENT_ 0
#if (_SANS_TIR_ORB_CORRIGE_==0)
c----------------------------------------------------------------------
      subroutine TirerRorb_SE ( Masse,Ener,Mom_cin,iRang,
     $     Peri,Apo,A_peri,A_apo,iRang_peri,
     $     Rnouveau,Vrad2nouveau )
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      double precision
                                ! EN ENTREE :
     $     Masse,Ener,Mom_cin,Peri,A_peri,Apo,A_apo
      integer
     $     iRang, iRang_peri
      double precision
                                ! EN SORTIE :
     $     Rnouveau,Vrad2nouveau
c
c var locales :
c ^^^^^^^^^^^^^
      logical lErr
c......................................................................
c
#if (_TEST_MODIF_PROP_SE_ACT_>0) /* test 26.4.2001 */
      call test_modif_prop_orb_SE_act (.FALSE.,
     $     Masse,Ener,Mom_cin,1.0d30,1.0d30,Peri,Apo, lErr)
      if (lErr) call Meurt('TirerRorb_SE')
#endif
      call  Realiser_TirOrb(
     $     Masse,Ener,Mom_cin,Peri,Apo,A_peri,A_apo,iRang_peri,
     $     Rnouveau,Vrad2nouveau )
      return
c
      end
c
c----------------------------------------------------------------------
      subroutine Realiser_TirOrb(
     $     Masse,Ener,Mom_cin,Peri,Apo,A_peri,A_apo,iRang_peri,
     $     Rnouveau,Vrad2nouveau )
c----------------------------------------------------------------------
      include 'EvolAmas_TirOrb2_inc.f'
      include 'EvolAmas_Compteurs_common.f'
      include 'EvolAmas_LimExt_common.f'
      include 'Flag_common.f'
      include 'test_modif_prop_orb_SE_act_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision          ! EN ENTREE :
     $     Masse,Ener,Mom_cin,Peri,A_peri,Apo,A_apo
      integer
     $     iRang_peri
      double precision          ! EN SORTIE :
     $     Rnouveau,Vrad2nouveau
c
c constantes et apparentees :
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer N_Raffinements
#if (_NO_REFINEMENT_==0)
      parameter (N_Raffinements=5)
      integer N_lim_VN(N_Raffinements),N_interv(N_Raffinements)
      data N_lim_VN             ! Valeurs limites du nb d'essais apres lesquelles
     $     /10,100,1000,10000,1000000/  ! on raffine l'evaluation de la borne sur P_TH
                                ! en considerant plus de sous-intervalles sur [Peri,Apo]
      data N_interv /1,5,20,50,100/ ! Nb successifs de sous-intervalles de "raffinement"
#else
      parameter (N_Raffinements=1)
      integer N_lim_VN(N_Raffinements),N_interv(N_Raffinements)
      data N_lim_VN             ! Valeurs limites du nb d'essais apres lesquelles
     $     /10000000/  ! on raffine l'evaluation de la borne sur P_TH
                                ! en considerant plus de sous-intervalles sur [Peri,Apo]
      data N_interv /1/ ! Nb successifs de sous-intervalles de "raffinement"
#endif
      integer NVN_verb
      parameter (NVN_verb=50)

      double precision ecc_min,Fact_prec_R,Fact_DistPeriApo_min,
     $     Fact_Vr2_min,coef_rapp_proba
      parameter (
     $     ecc_min = 5.0d-4,    ! eccentricite en dessous de-laquelle on considere l'orbite comme circulaire
     $     Fact_prec_R = 1.0d-5, ! Precision relative de la determination de Rnouveau par Newton_Raphson
     $     Fact_DistPeriApo_min = 1.0d-3, ! delimite les zones "proches" du peri- ou de l'apo-centre
     $     Fact_Vr2_min = 1.0d-5, ! Determine la limites des valeurs de Vr2 "dangereusement petites"
     $     coef_rapp_proba = 0.9999d0 ! coef de "securite" applique au rapport dPdr_reel/borne_dPdR
     $     )
      double precision DeuxSurPi
      parameter (DeuxSurPi=0.63661977236758138243d0)
c
c var locales :
c ^^^^^^^^^^^^^
      integer i_VN,iRang,i_Raff,N_lim,i
      double precision prec_R,DistPeriApo_min,Vr2_min,DeuxEB,
     $     rapport_P_TH,rapport_Vrad,Feff, la,lp, A_,B_, Ikep,dIkep,
     $     arg_sqrt, B_dummy, dB_dummy, B_peri, B_apo
      logical l_interv_borne_multiples
c
c fonctions :
c ^^^^^^^^^^^
      integer iRang_TO
      double precision random, rtsafe, P_TH
      external BdB_1interv_NR, BdB_Ninterv_NR
c
c common particulier :
c ^^^^^^^^^^^^^^^^^^^^
      double precision B_rand
      common /BdB_NR_common/ B_rand
c......................................................................
c
c---- Initialisations
c
      call Init_TirOrb(Masse,Ener,Mom_cin,Peri,Apo)
      prec_R = Fact_prec_R*(Ra_TO-Rp_TO)
      DistPeriApo_min = Fact_DistPeriApo_min*(Ra_TO-Rp_TO)
      Vr2_min = Fact_Vr2_min*A_apo/Apo
      DeuxEB = -J_TO**2/(Ra_TO*Rp_TO)
c
c---- traitons a part le cas d'une orbite quasi-circulaire
c
      if (ecc_TO.LT.ecc_min) then
         Rnouveau = 0.5d0*(Ra_TO+Rp_TO)
#if (_LIMITE_EXT_==1)
         Rnouveau = min(Rnouveau,R_LimExt)
#endif
         call Determiner_Vr2_Rang(Rnouveau,Vrad2nouveau,iRang)
         goto 9
      end if

      i_Raff= 1
      N_lim = N_lim_VN(i_Raff)
      call Init_borne_TirOrb(N_interv(i_Raff))

      i_VN=0

      do while(.TRUE.)

         i_VN = i_VN+1
         if (i_VN.GT.N_lim) then
            if (i_Raff.LT.N_Raffinements) then
               i_Raff = i_Raff+1
               N_lim = N_lim_VN(i_Raff)
#if (_VERB_TO_>=1)
               if (i_Raff.GT.2) then
                  write(0,100) i_VN,i_Raff,N_lim
 100              format( '@@@ Realiser_TirOrb : tirage = ',I4,
     $                 ' Raffinement de la grille ',I2,
     $                 ' Nb de tirages limite : ',I7)
               end if
#endif
               l_interv_borne_multiples = (N_interv(i_Raff).GT.1)
               if (l_interv_borne_multiples) then      
                  call Init_borne_TirOrb(N_interv(i_Raff))   
               end if
#if (_ECRIRE_TEST_TO_>=1)
               if (lFlag) then
                  write(0,*) '###### _ECRIRE_TEST_TO_ ######'
                  B_rand=0.0d0
                  open(11,file='_out_test_TO.asc')
                  write(11,900) iPas_Evol
 900              format('# iPas_Evol =',I10) 
                  call Ecrire_Test_TO(11,5000)
                  close(11)
                  open(11,file='_out_TO.asc')
                  call Ecrire_TO_asc(11)
                  close(11)
                  write(0,*) '###### TERMINE #######'
                  stop
               end if
#endif
            else
                                ! le nb max de rejections a ete atteint
               write(0,*) ('!!! Nb Max de tirage de VN realises !!!')
               call Meurt('Realiser_TirOrb')
            end if
         end if
               
#if (_VERB_TO_>=2)
         if (i_VN.GE.NVN_verb .AND. MOD(i_VN,NVN_verb).EQ.0) then
            write(0,110) i_VN
 110        format( '@@ Realiser_TirOrb : tirage = ',I4)
         end if
#endif
         iCompt_essais_TO = iCompt_essais_TO+1
c
c---- I) on tire un R en fonction d'une fonction-borne superieure
c        a la densite de proba que l'on cherche a appliquer
c     
c     dist. prob ~ P_TH/Vrad
c     borne      ~ borne_P_TH_TO/Vrad_kep
c
         if (l_interv_borne_multiples) then
            B_rand=random()*Bint_TO(N_borne_TO)
#if (_LIMITE_EXT_==1)
            if (Ra_TO.GT.R_LimExt) then
               call Calc_IdI_kepler(ecc_TO,UnSurDR_TO*(R_LimExt-Rp_TO),
     $              Ikep,dIkep)

               B_rand=DeuxSurPi*Ikep*B_rand
               Rnouveau = rtsafe(BdB_Ninterv_NR,Rp_TO,R_LimExt,prec_R)
            else
               Rnouveau = rtsafe(BdB_Ninterv_NR,Rp_TO,Ra_TO,prec_R)
            end if
#else
#   if (_TEST_RTSAFE_>0)
c$$$            call BdB_Ninterv_NR_(Ra_TO,B_apo,dB_dummy)
c$$$            if (Bint_TO(N_borne_TO).GT.B_apo) then
c$$$               write(0,'(A,2(E25.20," "),A)')
c$$$     $              '!!! Bint_TO(N_borne_TO) > B_apo : ',
c$$$     $              Bint_TO(N_borne_TO),B_apo,' !!!'
c$$$               call Meurt('Realiser_TirOrb')
c$$$            end if
            call BdB_Ninterv_NR(Rp_TO,B_peri,dB_dummy)
            call BdB_Ninterv_NR(Ra_TO,B_apo, dB_dummy)
            if ( (B_peri.GT.0.0d0 .AND. B_apo.GT.0.0d0) .OR.
     $           (B_peri.LT.0.0d0 .AND. B_apo.LT.0.0d0) ) then
               write(0,*) '!!! B_peri and B_apo from BdB_Ninterv_NR',
     $              ' have same sign !!!'
               write(0,*) '    Rp_TO,B_peri = ',Rp_TO,B_peri
               write(0,*) '    Ra_TO,B_apo  = ',Ra_TO,B_apo
               open(11,file='_RIP_BdB_Ninterv_NR.asc')
               do i=0,100000
                  Rnouveau = Rp_TO + i*(Ra_TO-Rp_TO)/100000.0d0
                  call BdB_Ninterv_NR(Rnouveau,B_dummy,dB_dummy)
                  write(11,*)Rnouveau,B_dummy,dB_dummy
               end do
               close(11)
               call Meurt('Realiser_TirOrb')
            end if
#   endif
            Rnouveau = rtsafe(BdB_Ninterv_NR,Rp_TO,Ra_TO,prec_R)
#endif
         else
            B_rand=random()
#if (_LIMITE_EXT_==1)
            if (Ra_TO.GT.R_LimExt) then
               call Calc_IdI_kepler(ecc_TO,UnSurDR_TO*(R_LimExt-Rp_TO),
     $              Ikep,dIkep)
               B_rand=DeuxSurPi*Ikep*B_rand
               Rnouveau = rtsafe(BdB_1interv_NR,Rp_TO,R_LimExt,prec_R)
            else
               Rnouveau = rtsafe(BdB_1interv_NR,Rp_TO,Ra_TO,prec_R)
            end if
#else
#   if (_TEST_RTSAFE_>0)
            call BdB_1interv_NR(Rp_TO,B_peri,dB_dummy)
            call BdB_1interv_NR(Ra_TO,B_apo, dB_dummy)
            if ( (B_peri.GT.0.0d0 .AND. B_apo.GT.0.0d0) .OR.
     $           (B_peri.LT.0.0d0 .AND. B_apo.LT.0.0d0) ) then
               write(0,*) '!!! B_peri and B_apo from BdB_1interv_NR',
     $              ' have same sign !!!'
               write(0,*) '    Rp_TO,B_peri = ',Rp_TO,B_peri
               write(0,*) '    Ra_TO,B_apo  = ',Ra_TO,B_apo
               open(11,file='_RIP_BdB_1interv_NR.asc')
               do i=0,100000
                  Rnouveau = Rp_TO + i*(Ra_TO-Rp_TO)/100000.0d0
                  call BdB_1interv_NR(Rnouveau,B_dummy,dB_dummy)
                  write(11,*)Rnouveau,B_dummy,dB_dummy
               end do
               close(11)
               call Meurt('Realiser_TirOrb')
            end if
#   endif

            Rnouveau = rtsafe(BdB_1interv_NR,Rp_TO,Ra_TO,prec_R)
#endif
         end if

         call Determiner_Vr2_Rang(Rnouveau,Vrad2nouveau,iRang)
         if (iRang.LT.iRang_peri+1) iRang=iRang+1
                                ! calcul du rapport des vitesses radiales
                                ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (l_interv_borne_multiples) then
            rapport_P_TH = P_TH(iRang)/borne_P_TH_TO(iRang_TO(Rnouveau))
         else
            rapport_P_TH = P_TH(iRang)/P_TH(iRang_peri)
         end if
                                ! Traitement special a proximite du pericentre
                                ! ou de l'apocentre car Vrad y est tres petit
c$$$XXX         if (Vrad2nouveau.LT.Vr2_min) then

            if ((Rnouveau-Rp_TO).LT.DistPeriApo_min)
     $           then           ! proximite du pericentre
#if (_SURVEILLE_SQRT_ > 0)
               arg_sqrt = (-J_TO**2*(Rp_TO-Ra_TO))/
     $              (2.0d0*Ra_TO*(Rp_TO*(A_peri-FracMasseInt*M_TO)+
     $              J_TO**2))
               if (arg_sqrt.LT.0.0d0) then
                  write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
                  call Meurt('Realiser_TirOrb (a)')
               end if
               rapport_Vrad = sqrt( arg_sqrt )
#else
               rapport_Vrad = sqrt( (-J_TO**2*(Rp_TO-Ra_TO))/
     $              (2.0d0*Ra_TO*(Rp_TO*(A_peri-FracMasseInt*M_TO)+
     $              J_TO**2)) )
#endif
            else if ((Ra_TO-Rnouveau).LT.DistPeriApo_min)
     $              then        ! proximite de l'apocentre
#if (_SURVEILLE_SQRT_ > 0)
               arg_sqrt = (-J_TO**2*(Ra_TO-Rp_TO))/
     $              (2.0d0*Rp_TO*(Ra_TO*(A_apo-FracMasseInt*M_TO)+
     $              J_TO**2))
               if (arg_sqrt.LT.0.0d0) then
                  write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
                  call Meurt('Realiser_TirOrb (b)')
               end if
               rapport_Vrad = sqrt( arg_sqrt )
#else
               rapport_Vrad = sqrt( (-J_TO**2*(Ra_TO-Rp_TO))/
     $              (2.0d0*Rp_TO*(Ra_TO*(A_apo-FracMasseInt*M_TO)+
     $              J_TO**2)) )
#endif
c$$$            else                ! probleme !!!
c$$$               write(0,*) '!!! Vrad2nouveau trop petit : ',
c$$$     $              Vrad2nouveau,' !!!'
c$$$               call Meurt('Realiser_TirOrb')
c$$$            end if
                                ! Loin du pericentre et de l'apocentre
         else
#if (_SURVEILLE_SQRT_ > 0)
            arg_sqrt = DeuxEB*(Rnouveau-Rp_TO)*
     $           (Rnouveau-Ra_TO)/Vrad2nouveau
            if (arg_sqrt.LT.0.0d0) then
               write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
               call Meurt('Realiser_TirOrb (c)')
            end if
            rapport_Vrad = sqrt( arg_sqrt ) / Rnouveau
#else
            rapport_Vrad = sqrt(
     $           DeuxEB*(Rnouveau-Rp_TO)*(Rnouveau-Ra_TO)/Vrad2nouveau
     $           ) / Rnouveau
#endif
         end if

         rapport_P_TH = coef_rapp_proba*rapport_P_TH
         rapport_Vrad = coef_rapp_proba*rapport_Vrad
#if (_SURV_TIRORB_>=1)
         if (rapport_P_TH.GT.1.0d0 .OR. rapport_Vrad.GT.1.0d0) then
            write(0,120) rapport_P_TH, rapport_Vrad, Rnouveau,
     $           Vrad2nouveau, (Apo-Peri)/(Apo+Peri)
 120        format('!!! rapport_P_TH, rapport_Vrad, Rnouveau, ',
     $           'Vrad2nouveau, eccentricite =',5(' ',E12.5),' !!!')
            write(0,*) '!!! ',
     $           (Rnouveau-Rp_TO)/(Ra_TO-Rp_TO),
     $           (Rnouveau-Ra_TO)/(Rp_TO-Ra_TO)
c##############DEBUT
            lp=log10(Peri)
            la=log10(Apo)
            do i=1,10000
               Rnouveau = Peri + 0.000001d0*(Apo-Peri) +
     $              i*0.99999d0*(Apo-Peri)/dble(10001)
               iRang = -1
               call Determiner_Vr2_Rang(Rnouveau,Vrad2nouveau,iRang)
               call Determiner_PotRang_Arbre_r(Rnouveau,A_,B_,iRang)
                                ! calcul du rapport des vitesses radiales

                                ! Traitement special a proximite du pericentre
                                ! ou de l'apocentre car Vrad y est tres petit      
c$$$               if (Vrad2nouveau.LT.Vr2_min) then

                  if (abs(Rnouveau-Rp_TO).LT.DistPeriApo_min)
     $                 then     ! proximite du pericentre
                     rapport_Vrad = sqrt( (-J_TO**2*(Rp_TO-Ra_TO))/
     $                    (2.0d0*Ra_TO*(Rp_TO*(A_peri-FracMasseInt*M_TO)+
     $                    J_TO**2)) )
                  else if (abs(Ra_TO-Rnouveau).LT.DistPeriApo_min)
     $                    then  ! proximite de l'apocentre
#if (_SURVEILLE_SQRT_ > 0)
                     arg_sqrt = (-J_TO**2*(Ra_TO-Rp_TO))/
     $                    (2.0d0*Rp_TO*(Ra_TO*(A_apo-FracMasseInt*M_TO)+
     $                    J_TO**2))
                     if (arg_sqrt.LT.0.0d0) then
                        write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
                        call Meurt('Realiser_TirOrb (d)')
                     end if
                     rapport_Vrad = sqrt(arg_sqrt)
#else
                     rapport_Vrad = sqrt( (-J_TO**2*(Ra_TO-Rp_TO))/
     $                    (2.0d0*Rp_TO*(Ra_TO*(A_apo-FracMasseInt*M_TO)+
     $                    J_TO**2)) )
#endif
c$$$                  else          ! probleme !!!
c$$$                     rapport_Vrad = -1.0d0
c$$$                  end if
                                ! Loin du pericentre et de l'apocentre
               else
#if (_SURVEILLE_SQRT_ > 0)
                  arg_sqrt = DeuxEB*(Rnouveau-Rp_TO)*(Rnouveau-Ra_TO)/
     $                 Vrad2nouveau
                  if (arg_sqrt.LT.0.0d0) then
                     write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
                     call Meurt('Realiser_TirOrb (e)')
                  end if
                  rapport_Vrad = sqrt(arg_sqrt)/ Rnouveau
#else
                  rapport_Vrad = sqrt(
     $                 DeuxEB*(Rnouveau-Rp_TO)*(Rnouveau-Ra_TO)/
     $                    Vrad2nouveau
     $                 ) / Rnouveau
#endif
               end if
               rapport_Vrad = coef_rapp_proba*rapport_Vrad
               write(98,'(5(E20.13," "),I7)') Rnouveau,Vrad2nouveau,
     $              rapport_Vrad,A_,B_,iRang
            end do
            write(0,*)'*********** ',-J_TO**2*(Ra_TO+Rp_TO)/(2.0d0*Ra_TO*Rp_TO)
c##############FIN
               
            call Meurt('Realiser_TirOrb')
         end if
#endif /* if (_SURV_TIRORB_>=1) */
c
c---- II) on conserve ce point avec une probabilite egale
c         au rapport entre la vraie densite de proba et sa borne.
c
         if (random().LT.rapport_P_TH*rapport_Vrad) goto 9

      end do

 9    continue      
c
      end
c
#endif /* (_SANS_TIR_ORB_CORRIGE_==0) */
#if (_SANS_TIR_ORB_CORRIGE_==1)
c----------------------------------------------------------------------
      subroutine TirerRorb_SE ( Masse,Ener,Mom_cin,iRang,
     $     Peri,Apo,A_peri,A_apo,iRang_peri,
     $     Rnouveau,Vrad2nouveau )
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      double precision
                                ! EN ENTREE :
     $     Masse,Ener,Mom_cin,Peri,A_peri,Apo,A_apo
      integer
     $     iRang, iRang_peri
      double precision
                                ! EN SORTIE :
     $     Rnouveau,Vrad2nouveau
c
c var locales :
c ^^^^^^^^^^^^^
      logical lErr
c......................................................................
c
      call  Realiser_TirOrb_direct(
     $     Masse,Ener,Mom_cin,Peri,Apo,A_peri,A_apo,
     $     Rnouveau,Vrad2nouveau )
      return
c
      end
c
#endif /* (_SANS_TIR_ORB_CORRIGE_==1) */
c----------------------------------------------------------------------
      subroutine Realiser_TirOrb_direct(
     $     Masse,Ener,Mom_cin,Peri,Apo,A_peri,A_apo,
     $     Rnouveau,Vrad2nouveau )
c----------------------------------------------------------------------
c     Tirage sans correction pour la probabilite de replacement
c     Version utile pour la duplication des SE
c
      include 'EvolAmas_TirOrb2_inc.f'
      include 'EvolAmas_Compteurs_common.f'
      include 'EvolAmas_LimExt_common.f'
      include 'Flag_common.f'
      include 'test_modif_prop_orb_SE_act_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision          ! EN ENTREE :
     $     Masse,Ener,Mom_cin,Peri,A_peri,Apo,A_apo
      double precision          ! EN SORTIE :
     $     Rnouveau,Vrad2nouveau
c
c constantes et apparentees :
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      double precision ecc_min,Fact_prec_R,Fact_DistPeriApo_min,
     $     Fact_Vr2_min,coef_rapp_proba
      parameter (
     $     ecc_min = 5.0d-4,    ! eccentricite en dessous de-laquelle on considere l'orbite comme circulaire
     $     Fact_prec_R = 1.0d-4, ! Precision relative de la determination de Rnouveau par Newton_Raphson
     $     Fact_DistPeriApo_min = 1.0d-2, ! delimite les zones "proches" du peri- ou de l'apo-centre
     $     Fact_Vr2_min = 1.0d-12, ! Determine la limites des valeurs de Vr2 "dangereusement petites"
     $     coef_rapp_proba = 0.9999d0 ! coef de "securite" applique au rapport dPdr_reel/borne_dPdR
     $     )
      integer i_VNmax
      parameter (i_VNmax=1000)
      double precision DeuxSurPi
      parameter (DeuxSurPi=0.6366197724d0)
c
c var locales :
c ^^^^^^^^^^^^^
      integer i_VN,iRang
      double precision prec_R,DistPeriApo_min,Vr2_min,DeuxEB,
     $     rapport_Vrad,Feff,Ikep,dIkep,
     $     arg_sqrt
c
c fonctions :
c ^^^^^^^^^^^
      double precision random, rtsafe, P_TH
      external BdB_1interv_NR
c
c common particulier :
c ^^^^^^^^^^^^^^^^^^^^
      double precision B_rand
      common /BdB_NR_common/ B_rand
c......................................................................
c
c---- Initialisations
c
      call Init_TirOrb(Masse,Ener,Mom_cin,Peri,Apo)
      prec_R = Fact_prec_R*(Ra_TO-Rp_TO)
      DistPeriApo_min = Fact_DistPeriApo_min*(Ra_TO-Rp_TO)
      Vr2_min = Fact_Vr2_min*A_apo/Apo
      DeuxEB = -J_TO**2/(Ra_TO*Rp_TO)
c
c---- traitons a part le cas d'une orbite quasi-circulaire
c
      if (ecc_TO.LT.ecc_min) then
         Rnouveau = 0.5d0*(Ra_TO+Rp_TO)
#if (_LIMITE_EXT_==1)
         Rnouveau = min(Rnouveau,R_LimExt)
#endif
         call Determiner_Vr2_Rang(Rnouveau,Vrad2nouveau,iRang)
         goto 9
      end if

      do i_VN=0,i_VNmax
c
c---- I) on tire un R en fonction d'une fonction-borne superieure
c        a la densite de proba que l'on cherche a appliquer
c     
c     dist. prob ~ 1/Vrad
c     borne      ~ 1/Vrad_kep
c
         B_rand=random()
#if (_LIMITE_EXT_==1)
         if (Ra_TO.GT.R_LimExt) then
            call Calc_IdI_kepler(ecc_TO,UnSurDR_TO*(R_LimExt-Rp_TO),
     $           Ikep,dIkep)
            B_rand=DeuxSurPi*Ikep*B_rand
            Rnouveau = rtsafe(BdB_1interv_NR,Rp_TO,R_LimExt,prec_R)
         else
            Rnouveau = rtsafe(BdB_1interv_NR,Rp_TO,Ra_TO,prec_R)
         end if
#else 
         Rnouveau = rtsafe(BdB_1interv_NR,Rp_TO,Ra_TO,prec_R)

#endif
         call Determiner_Vr2_Rang(Rnouveau,Vrad2nouveau,iRang)
                                ! calcul du rapport des vitesses radiales
                                ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                ! Traitement special a proximite du pericentre
                                ! ou de l'apocentre car Vrad y est tres petit
         if (Vrad2nouveau.LT.Vr2_min) then

            if (abs(Rnouveau-Rp_TO).LT.DistPeriApo_min)
     $           then           ! proximite du pericentre
               Feff = (2.0d0*J_TO**2/Rp_TO+A_peri-FracMasseInt*M_TO)/
     $              Rp_TO**2    ! Force effective au pericentre
#if (_SURVEILLE_SQRT_ > 0)
               arg_sqrt = (2.0d0*J_TO**2*(Ra_TO-Rnouveau))/
     $              (Ra_TO*Rp_TO*Feff)
               if (arg_sqrt.LT.0.0d0) then
                  write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
                  call Meurt('Realiser_TirOrb_direct (a)')
               end if
               rapport_Vrad = sqrt(arg_sqrt)/Rnouveau
#else
               rapport_Vrad = sqrt( (2.0d0*J_TO**2*(Ra_TO-Rnouveau))/
     $              (Ra_TO*Rp_TO*Feff) ) /Rnouveau
#endif
            else if (abs(Ra_TO-Rnouveau).LT.DistPeriApo_min)
     $              then        ! proximite de l'apocentre
               Feff = (2.0d0*J_TO**2/Ra_TO+A_apo-FracMasseInt*M_TO)/
     $              Ra_TO**2    ! Force effective a l'apocentre
#if (_SURVEILLE_SQRT_ > 0)
               arg_sqrt = (2.0d0*J_TO**2*(Rp_TO-Rnouveau))/
     $              (-Ra_TO*Rp_TO*Feff)
               if (arg_sqrt.LT.0.0d0) then
                  write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
                  call Meurt('Realiser_TirOrb_direct (b)')
               end if
               rapport_Vrad = sqrt(arg_sqrt)/Rnouveau
#else
               rapport_Vrad = sqrt( (2.0d0*J_TO**2*(Rp_TO-Rnouveau))/
     $              (-Ra_TO*Rp_TO*Feff) ) /Rnouveau
#endif
            else                ! probleme !!!
               write(0,*) '!!! Vrad2nouveau trop petit : ',
     $              Vrad2nouveau,' !!!'
               call Meurt('Realiser_TirOrb_direct')
            end if
                                ! Loin du pericentre et de l'apocentre
         else
#if (_SURVEILLE_SQRT_ > 0)
            arg_sqrt = DeuxEB*(Rnouveau-Rp_TO)*(Rnouveau-Ra_TO)/
     $           Vrad2nouveau
            if (arg_sqrt.LT.0.0d0) then
               write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
               call Meurt('Realiser_TirOrb_direct (c)')
            end if
            rapport_Vrad = sqrt(arg_sqrt)/Rnouveau
#else
            rapport_Vrad = sqrt(
     $           DeuxEB*(Rnouveau-Rp_TO)*(Rnouveau-Ra_TO)/Vrad2nouveau
     $           ) / Rnouveau
#endif
         end if
         rapport_Vrad = coef_rapp_proba*rapport_Vrad
#if (_SURV_TIRORB_>=1)
         if (rapport_Vrad.GT.1.0d0) then
            write(0,120) rapport_Vrad, Rnouveau,
     $           Vrad2nouveau, (Apo-Peri)/(Apo+Peri)
 120        format('!!! rapport_Vrad, Rnouveau, ',
     $           'Vrad2nouveau, eccentricite = ',5(E12.5,' '),' !!!')
               
            call Meurt('Realiser_TirOrb_direct')
         end if
#endif
c
c---- II) on conserve ce point avec une probabilite egale
c         au rapport entre la vraie densite de proba et sa borne.
c
         if (random().LT.rapport_Vrad) goto 9
      end do                    ! le nb max de rejections a ete atteint

      write(0,*) ('!!! Nb Max de tirage de VN realises !!!')
      call Meurt('Realiser_TirOrb_direct') 
c
 9    continue

#if (_SURV_TIRORB_>0)
      if (Rnouveau.LT.Peri) then
         write(0,*) '!!! Rnouveau<Peri : ',Rnouveau,Peri,' !!!'
         call Meurt('Realiser_TirOrb_direct')
      end if
      if (Rnouveau.GT.Apo) then
         write(0,*) '!!! Rnouveau>Apo : ',Rnouveau,Apo,' !!!'
         call Meurt('Realiser_TirOrb_direct')
      end if
#endif            
c
      
      end
c
c----------------------------------------------------------------------
      subroutine Init_TirOrb(Masse,Ener,Mom_cin,Peri,Apo)
c----------------------------------------------------------------------
      include 'EvolAmas_TirOrb2_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision Masse,Ener,Mom_cin,Peri,Apo
                                ! Ener,Mom_cin sont les valeurs specifiques !
c......................................................................
c
      M_TO = Masse
      E_TO = Ener
      J_TO = Mom_cin
      Ra_TO = Apo
      Rp_TO = Peri
      ecc_TO = (Apo-Peri)/(Apo+Peri)
      UnSurDR_TO = 1.0d0/(Ra_TO-Rp_TO)
      iCompt_essais_TO = 0
c
      end
c
#if (_SANS_TIR_ORB_CORRIGE_==0)
c----------------------------------------------------------------------
      subroutine Init_borne_TirOrb(N_borne)
c----------------------------------------------------------------------
      include 'EvolAmas_TirOrb2_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer N_borne           ! le nombre de sous-intervalles ou P_TH est evaluee
                                ! afin d'avoir une "bonne" borne pour dP/dR
c
c var locales :
c ^^^^^^^^^^^^^
      double precision A,B,Ri,Ri1,Ikep,dIkep,Ikep_prec,int_borne
      integer iRang,i
c
c fonctions :
c ^^^^^^^^^^^
      double precision P_TH
c......................................................................
c
      if (N_borne.GT.Nmax_borne_TO) then
         write(0,*) '>> N_borne_TO ramene a ',Nmax_borne_TO,
     $        ' dans Init_TirOrb'
      end if
      N_borne_TO = min(Nmax_borne_TO,N_borne)

      dlR_TO = LOG(Ra_TO/Rp_TO)/N_borne_TO

      ! modif 9.6.99 !! call Determiner_PotRang_Arbre_r(Rp_TO,A,B,iRang) ! < inutile ?
      Ikep_prec = 0.0d0
      int_borne = 0.0d0
      Ri = Rp_TO

      do i=1,N_borne_TO
         Ri1 = Ri
         Ri = Rp_TO*exp(i*dlR_TO) ! on opte pour des intervalles equidistants en Log(R)
         call Determiner_PotRang_Arbre_r(Ri1,A,B,iRang)
         borne_P_TH_TO(i) = P_TH(iRang)
         call Calc_IdI_kepler(ecc_TO,UnSurDR_TO*(Ri-Rp_TO),Ikep,dIkep)
         int_borne = int_borne + borne_P_TH_TO(i)*(Ikep-Ikep_prec)
         Bint_TO(i) = int_borne
         Ikep_prec = Ikep
      end do
c
      end
c
#endif /* (_SANS_TIR_ORB_CORRIGE_==0) */
c----------------------------------------------------------------------
      integer function iRang_TO(R)
c----------------------------------------------------------------------
      include 'EvolAmas_TirOrb2_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      !!iRang_TO = min(LOG(R/Rp_TO)/dlR_TO + 1,N_borne_TO)
      iRang_TO = int(min(LOG(R/Rp_TO)/dlR_TO + 1.0d0,dble(N_borne_TO))) !! modif. 2.3.2000
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_IdI_kepler(e,x,I,dI)
c----------------------------------------------------------------------
      implicit none
      include 'CstesMath.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision e,x,I,dI
c
c constantes :
c ^^^^^^^^^^^^
      double precision PiSur2,epseps,eps,Enorme
      parameter (
     $     PiSur2 = rPi/2.0d0,
     $     epseps = 1.0d-13,
     $     eps = 1.0d-9,
     $     Enorme = 1.0d30 )
c
c var locales :
c ^^^^^^^^^^^^^
      double precision UnMoinsX,racX1X,
     $     arg_sqrt
c......................................................................
c     
      UnMoinsX = 1.0d0-x
      if (x.LT.epseps) then
         if (x.LT.-10.0d0*epseps) then
            write(0,*) '!!! Calc_IdI_kepler called with negative x :',
     $           x,' !!!'
            call Meurt('Calc_IdI_kepler')
         end if
         I = 0.0d0
         dI = Enorme
      else if (UnMoinsX.LT.epseps) then
         if (UnMoinsX.LT.-10.0d0*epseps) then
            write(0,*) '!!! Calc_IdI_kepler called with negative 1-x :',
     $           UnMoinsX,' !!!'
            call Meurt('Calc_IdI_kepler')
         end if
         I = PiSur2
         dI = Enorme

      else
      
         racX1X = x*UnMoinsX
         if (racX1X.LT.0.0d0) then
            if (racX1X.GT.-eps) then
               racX1X=0.0d0
            else
               write(0,*) '!!! X*(1-X) = ',racX1X,' !!!'
               call Meurt('Calc_IdI_kepler')
            end if
         else
#if (_SURVEILLE_SQRT_ > 0)
            arg_sqrt = x*UnMoinsX
            if (arg_sqrt.LT.0.0d0) then
               write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
               call Meurt('Calc_IdI_kepler (a)')
            end if
            racX1X = sqrt(arg_sqrt)
#else         
            racX1X = sqrt(x*UnMoinsX)
#endif
         end if
         if (UnMoinsX.LT.eps) then
            I = PiSur2-UnMoinsX -e*racX1X
         else
#if (_SURVEILLE_SQRT_ > 0)
            arg_sqrt = x/UnMoinsX
            if (arg_sqrt.LT.0.0d0) then
               write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
               call Meurt('Calc_IdI_kepler (b)')
            end if
            I = atan(sqrt(arg_sqrt)) - e*racX1X
#else 
            I = atan(sqrt(x/UnMoinsX)) - e*racX1X
#endif
         end if
         if (racX1X.LT.eps) then
            dI = Enorme
         else
            dI = (0.5d0+e*(x-0.5d0))/racX1X
         end if

      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine BdB_1interv_NR(R,B,dB)
c----------------------------------------------------------------------
      include 'EvolAmas_TirOrb2_inc.f' 
      include 'CstesMath.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision R,B,dB
c
c constantes :
c ^^^^^^^^^^^^
      double precision DeuxSurPi
      parameter (DeuxSurPi=2.0d0/rPi)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Ikep,dIkep
c
c common particulier :
c ^^^^^^^^^^^^^^^^^^^^
      double precision B_rand
      common /BdB_NR_common/ B_rand
c......................................................................
c
      call Calc_IdI_kepler(ecc_TO,UnSurDR_TO*(R-Rp_TO),Ikep,dIkep)
      B = DeuxSurPi*Ikep  - B_rand
      dB = DeuxSurPi*dIkep
c
      end
c
c----------------------------------------------------------------------
      subroutine BdB_Ninterv_NR(R,B,dB)
c----------------------------------------------------------------------
      include 'EvolAmas_TirOrb2_inc.f' 
c
c arguments :
c ^^^^^^^^^^^
      double precision R,B,dB
c
c var locales :
c ^^^^^^^^^^^^^
      integer i_TO
      double precision Ikep,dIkep,Ikep_i,dIkep_i,Bi1,Ri
c
c fonction :
c ^^^^^^^^^^
      integer iRang_TO
c
c common particulier :
c ^^^^^^^^^^^^^^^^^^^^
      double precision B_rand
      common /BdB_NR_common/ B_rand
c......................................................................
c
      i_TO = iRang_TO(R)
      if (i_TO.GT.1) then
         Ri = Rp_TO*exp((i_TO-1)*dlR_TO)
         call Calc_IdI_kepler
     $        (ecc_TO,UnSurDR_TO*(Ri-Rp_TO),Ikep_i,dIkep_i)
         Bi1 = Bint_TO(i_TO-1)
      else
         Ikep_i = 0.0d0
         Bi1 = 0.0d0
      end if
      call Calc_IdI_kepler(ecc_TO,UnSurDR_TO*(R-Rp_TO),Ikep,dIkep)
      B = Bi1 + borne_P_TH_TO(i_TO)*(Ikep-Ikep_i) - B_rand
      dB = borne_P_TH_TO(i_TO)*dIkep
c
      end
c
c----------------------------------------------------------------------
      subroutine BdB_Ninterv_NR_(R,B,dB)
c----------------------------------------------------------------------
      include 'EvolAmas_TirOrb2_inc.f' 
c
c arguments :
c ^^^^^^^^^^^
      double precision R,B,dB
c
c var locales :
c ^^^^^^^^^^^^^
      integer i_TO
      double precision Ikep,dIkep,Ikep_i,dIkep_i,Bi1,Ri
c
c fonction :
c ^^^^^^^^^^
      integer iRang_TO
c
c common particulier :
c ^^^^^^^^^^^^^^^^^^^^
      double precision B_rand
      common /BdB_NR_common/ B_rand
c......................................................................
c
      i_TO = iRang_TO(R)
      if (i_TO.GT.1) then
         Ri = Rp_TO*exp((i_TO-1)*dlR_TO)
         call Calc_IdI_kepler
     $        (ecc_TO,UnSurDR_TO*(Ri-Rp_TO),Ikep_i,dIkep_i)
         Bi1 = Bint_TO(i_TO-1)
      else
         Ikep_i = 0.0d0
         Bi1 = 0.0d0
      end if
      call Calc_IdI_kepler(ecc_TO,UnSurDR_TO*(R-Rp_TO),Ikep,dIkep)
      B = Bi1 + borne_P_TH_TO(i_TO)*(Ikep-Ikep_i)
      dB = borne_P_TH_TO(i_TO)*dIkep
      write(89,'(2(I3," "),5(E25.20," "))')N_borne_TO, i_TO, Bi1, Ikep,Ikep_i, borne_P_TH_TO(i_TO),UnSurDR_TO*(R-Rp_TO)
c
      end
c
c----------------------------------------------------------------------
      subroutine Determiner_Vr2_Rang(R,Vr2,iRang)
c----------------------------------------------------------------------
      include 'EvolAmas_TirOrb2_inc.f' 
c
c arguments :
c ^^^^^^^^^^^
      double precision R,Vr2
      integer iRang
c
c var locales :
c ^^^^^^^^^^^^^
      double precision A,B
c......................................................................
c
      call Determiner_PotRang_Arbre_r(R,A,B,iRang)
      Vr2 = 2.0d0*(E_TO-B)-(2.0d0*(A-FracMasseInt*M_TO)+J_TO**2/R)/R
c
      end
c
#if (_SANS_TIR_ORB_CORRIGE_==0)
c----------------------------------------------------------------------
      subroutine Ecrire_Test_TO(iFich,Npoints)
c----------------------------------------------------------------------
      include 'EvolAmas_TirOrb2_inc.f' 
c
c arguments :
c ^^^^^^^^^^^
      integer iFich,Npoints
c
c var locales :
c ^^^^^^^^^^^^^
      double precision R,s,g,DeuxEB,Vr2,Vr2_kep,B,dB
      integer i,iRang,i_TO
c
c fonctions :
c ^^^^^^^^^^^
      double precision P_TH
      integer iRang_TO
c......................................................................
c
      DeuxEB = -J_TO**2/(Ra_TO*Rp_TO)
      write(iFich,100)
 100  format(
     $     '# Fonctions utilisees pour tirage position orbitale',/,
     $     '# 1:s 2:R 3:Vr2 4:Vr2_kep 5:P_TH 6:borne_P_TH 7:g',
     $     ' 8:B 9:dB',
     $     ' 10:iRang 11:i_TO' )
      do i=1,Npoints
         s=dble(i)/dble(Npoints+1)
         R=Rp_TO+(Ra_TO-Rp_TO)*s*s*(3.0d0-2.0d0*s)
         call Determiner_Vr2_Rang(R,Vr2,iRang)
         Vr2_kep = DeuxEB*(R-Rp_TO)*(R-Ra_TO)/R**2
         g = 6.0d0*(Ra_TO-Rp_TO)*s*(1.0d0-s)/sqrt(Vr2)
         i_TO = iRang_TO(R)
         call BdB_Ninterv_NR(R,B,dB)
         write(iFich,110)s,R,Vr2,Vr2_kep,
     $        P_TH(iRang),borne_P_TH_TO(i_TO),g,B,dB,iRang,i_TO
 110     format(9(E14.6,' '),I7,' ',I5)
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_TO_asc(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_TirOrb2_inc.f' 
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
c
c var locales :
c ^^^^^^^^^^^^^
      integer i
      double precision Ri
c......................................................................
c
      write(iFich,100) Rp_TO,Ra_TO
 100  format(
     $     '# Borne sup pour tirage position orbitale',/,
     $     '# Peri,Apo = ',2(E14.6,' '),/,
     $     '# 1:R_inf 2:borne_P_TH_TO 3:Bint_TO' )
      do i=1,N_borne_TO
         Ri = Rp_TO*exp((i-1)*dlR_TO)
         write(iFich,110) Ri,borne_P_TH_TO(i),Bint_TO(i)
 110     format(3(E14.6,' '))
      end do
c
      end
c
#endif /* (_SANS_TIR_ORB_CORRIGE_==0) */
