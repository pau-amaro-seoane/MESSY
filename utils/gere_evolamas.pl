#! /usr/bin/perl

# script pour lancer des simulations Monte Carlo d'amas sur "isis", Berne, 02.02.2001
#
# emploi type : lance_evolamas.pl --continue --repertoire=PC_839 --option i_fin=500000000

use Getopt::Long;
Getopt::Long::Configure("default");

use sigtrap qw(die normal-signals);

my $executable="./evolamas";
my $input_file="input_EvolAmas";
my $autre_repertoire="";
my %defines = ();
my $nice=19;
my $attente=10;

# Options sur la ligne de commande

GetOptions(
	   "define=s" => \%defines,
	   "option=s" => \%defines,
	   "O=s"      => \%defines,
	   "continuer|c|C"   => \$continuer,
	   "blanc|a_blanc|B" => \$a_blanc,
	   "verbeux"         => \$verbeux,
	   "forcer|F"        => \$forcer,
	   "executable|X=s"  => \$executable,
	   "input|i=s"       => \$input_file,
	   "repertoire|D=s"  => \$autre_repertoire,
	   "ecrit_input|P"   => \$filtre_input,
	   "random|R"        => \$random,
	   "silencieux|T"    => \$se_taire,
	   "output|F=s"      => \$output_file,
	   "attente=i"       => \$attente,
	   "nice|renice=i"   => \$nice
	  ) || die;

# on se place dans le bon repertoire

if ($autre_repertoire) { 
  chdir($autre_repertoire) or 
    die("impossible de se deplacer dans $autre_repertoire") 
  };

# Lecture des options du fichier "input_EvolAmas"

my %h_options;
open INPUT, $input_file;
while (<INPUT>) {
  if (/^(\w+)/) { $NomOpt=$1 };
  if (/^\s+([\w\.+-]+)/) { $ValOpt=$1; $h_options{$NomOpt}=$ValOpt };
}

# options (re)definies sur la ligne de commande

foreach $NomOpt (keys %defines) {
  $h_options{$NomOpt}=$defines{$NomOpt};
}

# Faut il creer une graine pour le generateur de nb aleatoire?

if ($random) {$h_options{'iRand_Seed'}=-1};

# Determination des valeurs des parametres pour continuer une simulation commencee

if ($continuer) {
  open LISTE, "ls -1 *%AMAS.* | grep -v RIP |";
  $NumMax=-1;
  while (<LISTE>) {
    if (/^(\w+)%([0-9]+)%AMAS.(xdr|bin)$/ && ($2+0 > $NumMax)) {
      $NumMax=$2+0; $NomSimul=$1;
    }
  }
  $h_options{"i_ini"}=$NumMax;
  $h_options{"NomSimul"}=$NomSimul;
}

# Base name of (most) files
my $NomSimul= '_EvolAmas' ;
$NomSimul = $h_options{'NomSimul'} if ($h_options{'NomSimul'});

# Read data from log file
my $FichLog = $NomSimul.'%%Log.asc';

if ($continuer) {
  my %quant_val_hash;
  if (open LOG, $FichLog) {
    while (<LOG>) {
      if ($_=~/^\| .+:.*/) {
	chop;
	s/\| +//;
	s/ *: */:/g;
	s/  */ /g;
	@elements=split /[: ]/;
	for($i=0;$i<=$#elements;$i+=2) {
	  $quant_val_hash{$elements[$i]} = $elements[$i+1];
	}
      }
    }
    close LOG;
    #foreach $quantity (sort keys %quant_val_hash) {
    #  print $quantity,'=',$quant_val_hash{$quantity},"\n";
    #}

    # Update average stellar mass
    $h_options{MasseEtoileDef} = $quant_val_hash{M_amas_en_Msol}/$quant_val_hash{Net_subsist};

    # Update mass of central object
    $h_options{M_TN_ini_def} = $quant_val_hash{M_TN};
  }
}

# Ecriture des options en sortie standard si demande

if ($filtre_input) {
  foreach $NomOpt (sort keys %h_options) {
    print $NomOpt,"\n    ",$h_options{$NomOpt},"\n";
  }
  exit;
}

# Lancement de la simulation
#============================

# Tester que la simulation n'a pas deja ete demarree

if (`ls *%Log.asc 2>/dev/null` and (not ($forcer or $continuer))) {
  print STDERR "!!! Il y a deja une simulation entamee dans ce repertoire !!!\n";
  print STDERR "!!! employer option -F ou -C                              !!!\n";
  die;
}
    
open INPUT_FILE, "> _input_" || die("!!! Impossible d'ouvrir fichier '_input_' !!!");
foreach $NomOpt (sort keys %h_options) {
  print INPUT_FILE $NomOpt,"\n    ",$h_options{$NomOpt},"\n";
}

$Ligne_commande = "$executable -f _input_";
$Ligne_commande = $Ligne_commande . ' > /dev/null 2>&1 ' if ($se_taire);
$Ligne_commande = $Ligne_commande . ' > ' . $output_file . ' 2>&1 ' if ($output_file);
$Ligne_commande = "nice -$nice " . $Ligne_commande if ($nice>0);

# on lance la commande en tache de fond...
if (not $a_blanc) {
  if ( -f '_PID_' ) {
    unlink '_PID_' or 
      die "!!! Impossible de supprimer le fichier '_PID_' !!!";
  }
  system($Ligne_commande." &")==0 and print "> Commande ".$Ligne_commande." lancee \n" or
    die ("!!! Erreur a l'execution de la commande '$Ligne_commande' !!!");

  # on recupere le pid de EvolAmas
  my $max=100;
  my $i=0;
  until ( defined  $pid_EvolAmas ) {
    die '!!! Impossible de recuperer le pid de "EvolAmas" !!!' if (++$i>$max);
    sleep 1;
    if (open INPUT, '< _PID_') {
      $pid_EvolAmas=<INPUT>;
      chop $pid_EvolAmas;
      close INPUT;
    }
  }
} else {
  print "Commande a executer : '$Ligne_commande'\n";
}

# Gestion des sauvegardes
#-------------------------

my $Fich_DemSauv='DemSauv.asc';
my $Fich_LogSauv=$NomSimul.'%%LogSauv.asc';

# Liste des sauvegardes
my %h_i_sauv = ();    # numero
my %h_cond_sauv = (); # condition a realiser
my %h_com_sauv = ();  # commentaire
my %h_real_sauv = (); # flag indiquant les sauvegardes realisees;

# on lit le fichier de journal des sauvegardes realisees 

if (not $forcer and (open INPUT, $Fich_LogSauv)) {
  while (<INPUT>) {
    if (/^id_sauvegarde/) {
      chop;
      ($champ, $valeur) = split /\t/, $_;
      $h_real_sauv{$valeur}=1;
      ($verbeux) and print "**> Sauvegarde ",$valeur," deja realisee\n";
     }
  }
  close INPUT;
}
 
#----------------------------------------------------------------------
# sous-routine : on lit le fichier de demandes de sauvegardes (DemSauv.asc)
#----------------------------------------------------------------------

# Example of the content of such file :

# id_sauvegarde  one_snapshot
# condition   M_amas < 0.9
#
# id_sauvegarde  another_snapshot
# condition   M_amas < 0.5
#
# id_sauvegarde  still_another
# condition   Tps_en_yr > 1e6
#
# id_sauvegarde  STOP
# condition   Tps_en_yr > 1e7


sub lire_conditions {
  ($verbeux) and print "**> Lecture des conditions de sauvegardes...\n";
  if (open INPUT, "< ".$Fich_DemSauv) {
    my $num=0;
    while (<INPUT>) {
      chop;
      next if (/^ *$/ or /^ *\#/);
      ($champ, $valeur) = split /\t/, $_, 2;
    SWITCH: {
	if ($champ eq 'id_sauvegarde') {
	  $id_sauv=$valeur; $h_i_sauv{$id_sauv}=++$num;
	  last SWITCH;
	};
	if ($champ eq 'condition')     {
	  #$valeur =~ s/([a-zA-Z][a-zA-Z0-9_]+)/\$\1/g;
	  $h_cond_sauv{$id_sauv}=$valeur;
	  last SWITCH;
	};
	if ($champ eq 'commentaire')   {
	  $h_com_sauv{$id_sauv}=$valeur;
	  last SWITCH;
	};
      }
    }
    close INPUT;
  }
}

# premiere lecture des conditions
&lire_conditions;

#======================================================================
# boucle principale : test periodique des conditions de sauvegarde
#======================================================================

%h_valeur=();

#----------------------------------------------------------------------
# sous-routine : on lit l'etat de la simulation dans le fichier de Log
#----------------------------------------------------------------------
sub lire_etat {
  open INPUT, 'tail -n 100 '.$FichLog.' |' or
    die "!!! impossible de lire dans ".$FichLog." !!!";
  while (<INPUT>) {
    if ($_=~/^\| .+:.*/) {
      chop;
      s/\| +//;
      s/ *: */:/g;
      s/  */ /g;
      @elements=split /[: ]/;
      for($i=0;$i<=$#elements;$i+=2) {
	$val=$elements[$i+1];
	if ( $val =~ /^ *[+-]?\d+\.?\d*([eEdD][+-]?\d+)? *$/ ) { # valeur numerique ou non ?
	  $val =~ s/[dD]/e/;
	} else {
	  $val='"'.$val.'"';
	}
	$h_valeur{$elements[$i]} = $val;
      }
    }
  }
  return # la suite est inutile...
  @variables = keys %h_valeur;
  $eval_str='';
  foreach $variable (@variables) {
    ($verbeux) and print "**> valeur de ".$variable." : ".$h_valeur{$variable}."\n";
    $eval_str .= '$'.$variable.'='.$h_valeur{$variable}.'; '; #'
  }
  eval $eval_str;
}

#----------------------------------------------------------------------
# sous-routine: determine si il faut demander une sauvegarde, en
# fonction des conditions de sauvegarde, evaluees pour les valeurs
# actuelles des grandeurs specifiees dans le fichier "Log"
#----------------------------------------------------------------------

sub test_dem_sauvegarde {
  foreach $id ( keys %h_cond_sauv ) {
    unless ($h_real_sauv{$id}) {
      $condition=$h_cond_sauv{$id};
      foreach $variable ( keys %h_valeur ) {
	$condition =~ s/$variable/$h_valeur{$variable}/g;
      }
      ($verbeux) and print "**> condition evaluee (",$id,"): ",$condition,"\n";
      if (eval $condition) {
	($verbeux) and print "**> condition de sauvegarde ".$id." : ".$h_cond_sauv{$id}." remplie\n";
	$h_real_sauv{$id}=1;
	return $id;
      }
    }
  }
  return 0;
}

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------

#my $attente=60; # delai d'attente en secondes
while ( kill 0 => $pid_EvolAmas ) { # on verifie que le processus est tjs la
  sleep $attente;
  if (-f '_LIRE_COND_SAUV_') {
    &lire_conditions;
    unlink '_LIRE_COND_SAUV_' or 
      die "!!! Impossible de supprimer le fichier '_LIRE_COND_SAUV_' !!!";
  }
  # Lecture de l'etat de la simulation
  &lire_etat;
  # Faut-il faire une sauvegarde
  if (not (-f '_SAUV_DEMANDEE_') and ($id_sauv = &test_dem_sauvegarde())) {
    open OUTPUT, "> _SAUV_DEMANDEE_" || die("!!! Impossible de creer fichier '_SAUV_DEMANDEE_' !!!");
    print OUTPUT $id_sauv,"\n",$h_cond_sauv{$id_sauv},"\n";
    print OUTPUT $h_com_sauv{$id_sauv},"\n" if ($h_com_sauv{$id_sauv});
    close OUTPUT;
  } 
}
print ">> 'EvolAmas' s'est apparemment termine.\n";


END { # Partie executee en quittant le script
  if ( $pid_EvolAmas ) {
    if ( kill 0 => $pid_EvolAmas ) {
      ($verbeux) and print "**>  On tue EvolAmas\n";
      kill KILL => $pid_EvolAmas;
    }
  }
}
__END__


# Follows an exemple of _EvolAmas%%DemSauv.asc
# apparently tabs must be present after "id_sauvegarde" and "id_sauvegarde"

id_sauvegarde   1
condition       Tps_en_yr > 1.5e10

id_sauvegarde   2
condition       Tps_en_yr > 2e10

id_sauvegarde   STOP
condition       Tps_en_yr > 1e11


