#!/usr/bin/ksh

gawk 'BEGIN{OFS="\t"} ++i==1 {for(j=1;j<=NF;j++){$j=j ": " $j}} {print}' | tr '\011' ' ' | gawk '++i==1 {$0="# " $0} i==2 {getline} i>1 {$0="  " $0}  {print}'
