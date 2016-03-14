#!/usr/bin/perl -w

#use strict;

use Getopt::Long;

my %h;
my $pas;
my %h_pas;
my $var;
my @elements;
my $val;
my $i;
my $non_def=1e30;

# options sur la ligne de commande
$Getopt::Long::autoabbrev = 0;
&GetOptions("non_def|nd|n=s") || exit 1;

$non_def=$opt_non_def if ($opt_non_def);

# lecture de toutes les donnees dans un hash of hashes
# ====================================================
while (<>) {
  if ($_=~/^\| .+:.*/) {
    chop;
    s/\| +//;
    s/ *: */:/g;
    s/  */ /g;
    @elements=split /[: ]/;
    for($i=0;$i<=$#elements;$i+=2) {
      if ($elements[$i] =~ /iPas_Evol/) {
	$pas=$elements[$i+1];
	$h_pas{$pas}=1;
      }
      $h{$elements[$i]}{$pas} = $elements[$i+1];
    }
  }
}

# ecriture au format /rdb
# =======================

my @variables = keys %h;
my @pas = sort {$a<=>$b} keys %h_pas;

# en-tete
my $ligne="";
foreach $var (@variables) {
  $ligne.="$var\t";
}
chop $ligne;
print "$ligne\n";
$ligne =~ s/[a-zA-Z0-9_]/-/g;
print "$ligne\n";

# table
foreach $pas (@pas) {
  $ligne="";
  foreach $var (@variables) {
    $h{$var}{$pas}=$non_def unless defined($h{$var}{$pas}); 
    $ligne.="$h{$var}{$pas}\t";
  }
  chop $ligne;
  print "$ligne\n";
}
