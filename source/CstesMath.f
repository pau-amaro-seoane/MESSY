c----------------------------------------------------------------------
c     constantes mathematiques pour codes divers
c----------------------------------------------------------------------
c
      double precision
     $     rPi, r2Pi, r4Pi, r4PiSur3,
     $     rRacDeux, rRacTrois,
     $     rRacDeux1, rRacTrois1,
     $     rEnorme, Minuscule
      parameter (
     $     rPi       = 3.141592653589793116d0,
     $     r2Pi      = 2.0d0*rPi, r4Pi = 4.0d0*rPi,
     $     r4PiSur3  = 4.0d0/3.0d0*rPi,
     $     rRacDeux  = .14142135623730951455d1,
     $     rRacTrois = .17320508075688771932d1,
     $     rRacDeux1 = 1.0d0/rRacDeux, rRacTrois1 = 1.0d0/rRacTrois,
     $     rEnorme  = 1.0d30, Minuscule = 1.0d-30
     $     )
