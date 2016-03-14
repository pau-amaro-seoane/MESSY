
      double precision Fact_Vrel_Lim
      parameter (Fact_Vrel_Lim=1.0d30) 

      double precision Nu_Max_Merger ! Maximum Vrel_infty/Vesc for mergers
                                     ! Used to limit domain of validity of lambda_crit_merger
      parameter (Nu_Max_Merger=2.0d0)

      double precision Mdebris_coll_max ! Collision products less massive than this are not considered star anymore
                                        ! Their mass is added to the escaping gas (or to the merged object)
      parameter (Mdebris_coll_max = 0.01d0)

      double precision Mmerg_coll_max ! Collision products more massive than this are not allowed
      parameter (Mmerg_coll_max = 1.0d6)
