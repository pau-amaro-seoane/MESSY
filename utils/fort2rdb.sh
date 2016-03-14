#!/usr/bin/ksh
# transforme un fichier fortran avec en-tete style "# 1: X 2: Y 3: Z"
# en fichier /rdb

# version se passant de fichier provisoire

gawk '
    BEGIN{
	ligne="#";
	while (match(ligne,/^#/)) {
	    en_tete=ligne;
	    getline;
	    ligne=$0
	}
	sub(/^# */,"",en_tete);
	gsub(/ *[0-9][0-9]* *: */," ",en_tete);
	sub(/  *$/,"",en_tete);
	sub(/^  */,"",en_tete);
	gsub(/  */," ",en_tete);
	gsub(/ /,"\t",en_tete);
	print en_tete;
	gsub(/[a-zA-Z_0-9]/,"-",en_tete);
	print en_tete;
	Nb=NF
	gsub(/  */,"#",$0); sub(/^#/,"",$0); sub(/#$/,"",$0); gsub(/#/,"\t",$0); print
    }
    NF==Nb { gsub(/  */,"#",$0); sub(/^#/,"",$0); sub(/#$/,"",$0); gsub(/#/,"\t",$0); print }
'
