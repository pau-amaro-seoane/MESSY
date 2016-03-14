#!/usr/bin/env zsh

# produit la liste des fichiers fortran correspondant a un fichier "makefile"
# presente en entree standard

tmp1=tmp1_$$

if [[ -n $1 ]]; then
    Makefile=$1
    ln -s $Makefile $tmp1
else
    Makefile="-"
    cat > $tmp1
fi

List=( $(for name in $(
    {
	cat $tmp1 | grep 'NOM *=' | cut -d"=" -f2
	cat $tmp1 | tr '\011' ' ' | sed 's/  */ /g' |\
	tr ' ' '\012' | sort -u | grep '^[a-zA-Z][a-zA-Z0-9_]*\.[oF]$' | sed 's/\..$//' 
    } | sort -u
); do
  ls -1 $name.F 2> /dev/null
done) )

print ">> Makefile : $Makefile <<"
grep '^\-D' < $tmp1 | sed 's/^\-D//; s/=/ /; s/ *\\$//'

grep -l '^#define ' $List |\
while read file; do
    print ">> Sourcefile : $file <<"
    grep '^#define ' $file | sed 's/ *\/\*//' | gawk '{print $2,$3}'
done



rm $tmp1

