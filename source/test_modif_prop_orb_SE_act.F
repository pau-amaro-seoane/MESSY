c----------------------------------------------------------------------
      subroutine test_modif_prop_orb_SE_act (lInit,
     $     Masse,Ener,Mom_cin,Ecin,Pot,Peri,Apo, lChgt)
c----------------------------------------------------------------------
      implicit none
      include 'test_modif_prop_orb_SE_act_inc.f' 
c
c arguments :
c ^^^^^^^^^^^
      double precision Masse,Ener,Mom_cin,Ecin,Pot,Peri,Apo
      logical lInit, lChgt
c......................................................................

      lChgt=.FALSE.

      if (lInit) then
         M_ref=Masse
         E_ref=Ener
         J_ref=Mom_cin
         T_ref=Ecin
         Pot_ref=Pot
         Peri_ref=Peri
         Apo_ref=Apo
         return
      end if

      if (Masse.LT.1.0d30.AND.Masse.NE.M_ref) then
         write(0,100) '!!! la masse a change : ',M_ref,Masse
         lChgt=.TRUE.
      end if
      if (Ener.LT.1.0d30.AND.abs(Ener-E_ref).GT.-1.0d-9*E_ref) then
         write(0,100) '!!! l''energie a change : ',E_ref,Ener
         lChgt=.TRUE.
      else
         E_ref = Ener
      end if
      if (Mom_cin.LT.1.0d30.AND.Mom_cin.NE.J_ref) then
         write(0,100) '!!! le moment cin. a change : ',J_ref,Mom_cin
         lChgt=.TRUE.
      end if
      if (Ecin.LT.1.0d30.AND.Ecin.NE.T_ref) then
         write(0,100) '!!! l''energie cin. a changee : ',T_ref,Ecin
         lChgt=.TRUE.
      end if
      if (Pot.LT.1.0d30.AND.Pot.NE.Pot_ref) then
         write(0,100) '!!! le potentiel a change : ',Pot_ref,Pot
         lChgt=.TRUE.
      end if
      if (Peri.LT.1.0d30.AND.Peri.NE.Peri_ref) then
         write(0,100) '!!! le pericentre a change : ',Peri_ref,Peri
         lChgt=.TRUE.
      end if
      if (Apo.LT.1.0d30.AND.Apo.NE.Apo_ref) then
         write(0,100) '!!! l''apocentre a change : ',Apo_ref,Apo
         lChgt=.TRUE.
      end if
 100  format(A,E20.11,' -> ',E20.11)
      
c
      end
