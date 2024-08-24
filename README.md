Introduction:
=============

This is a Hénon-style Monte-Carlo code to simulate spherical stellar clusters and galactic nuclei; i.e. Monte Carlo simulations for stellar dynamics.

It is described and used in the following papers (amongst others): Freitag & Benz 01; Freitag & Benz 02; Freitag, Rasio & Baumgardt 06; Freitag, Gürkan & Rasio 06; Freitag, Amaro-Seoane & Kalogera 06).

Beware that the code is not exceptionally user-friendly. Check by yourself in the picture on the left the mental status of Marc when he was debugging the code.

In particular comments and variable names harmoniously combine (Marc’s kind of) French and (some sort of) English. Sorry. But some people managed to run it.

I have updated and added information to the original documentation files. I am also (trying to) keeping updated the code and fixing some things/expanding others.
How to get the code and license.

Get the code:
=============

The code can be cloned from github like this:

```
git clone git@github.com:pau-amaro-seoane/MESSY.git
```

When you expand it, you should get the following directories:

    -codes : the beast’s den (the beast might be curled up into a tar.gz file)
    -initial_conditions : some simple initial condition files (cluster structures in .xdr format)
    -parameter_files : some parameter files (a run requires initial condition and parameter files, see doc)
    -some_results : some results (duh!) so you get an idea of what you should get
    -utils : (very) few scripts to manipulate output and a script to drive the evolution (gere_evolamas.pl)

Furthermore, if you want to switch on stellar collisions and wish to use Marc’s tabulated data from his SPH simulations, you’ll need the big file in my SPH section of the methods page (11M; it has to be gunzipped for use by `ME(SSY)**2`). If you use this data, please be sure to read the relevant parts in Freitag & Benz 02; Freitag & Benz 05 and Freitag, Rasio & Baumgardt 06 (and citing these sources will alow you to put the blame on the right people!).

Finally, as Marc says, and I can actually corroborate his statement, “for the adventurous or the downright foolish, you can try to hurt yourself” with a big spiky tar ball containing all his SM routines to deal with cluster simulations or “crash your computer using some more of his gothic shell/perl/python scripts (196K) to do random things” with `ME(SSY)**2` output (and possibly input).

License:
========

The license is:

```
    Copyright (c) 2008, Marc Dewi Freitag & Pau Amaro Seoane

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the “Software”), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    Note that, as stated in the documentation, the current version of the
    software requires the use of routines from the “Numerical Recipes in
    Fortran 77” (http://www.nrbook.com/a/bookfpdf.php) but is not
    distributed with them as they are not covered by the present
    license. Users of the software should obtain their own copy of these
    routines or replace them by other routines realising the same tasks.

    THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
```

Instructions:
==============

In the next paragraphs related to the code you will find an exhaustive documentation on how to install and run it, and how to read the output. The initial documentation was based on a pdf document that Marc wrote some years ago, but I am expanding it and updating it regularly. Hence, if you have the old pdf and code, I recommend you to check the potential differences with this site if you plan on using the code for something serious.
References

`ME(SSY)**2` stands for “Monte-carlo Experiments with Spherically SYmmetric Stellar SYstems”… and this name actually my suggestion. The name also conveys something of the programming style Marc used, as he himself says. Note that `ME(SSY)**2` is just the “commercial name” but, traditionally, since Marc started developing this code in 1995 or 1996, the executable is called evolamas (something like “evolve cluster” in French) and the main source file is EvolAmas.F.

The physical and algorithmic principles underlying `ME(SSY)**2` are described in:

    -Freitag, M. & Benz, W. 2001, A New Monte Carlo Code for Star Cluster Simulations: I. Relaxation, A&A, 375, 711 (http://uk.arxiv.org/abs/astro-ph/0102139, http://dx.doi.org/10.1051/0004-6361:20010706)
    -Freitag, M. & Benz, W. 2002, A New Monte Carlo Code for Star Cluster Simulations: II. Central Black Hole and Stellar Collisions, A&A, 394, 345 (http://uk.arxiv.org/abs/astro-ph/0204292, http://dx.doi.org/10.1051/0004-6361:20021142)

The code is based on the scheme pioneered by Michel Hénon in the 70s to follow the evolution of globular clusters

    Hénon 1971
    Hénon 1971b
    Hénon 1973
    Hénon 1975

*Numerical Recipes routines*

`ME(SSY)**2` should be free and open. It is distributed under an MIT license. Unfortunately, it relies on the use of some routines taken from ”The numerical recipes in Fortran” which you can only use if you have bought them (http://www.nr.com/, http://www.numerical-recipes.com/com/storefront.html). As of writing these lines, I am trying to find open and free alternatives to these routines. For the time being, they are in files whose names and with NR.f or NR.F, which are not distributed with `ME(SSY)**2`. The routines used are:

```
covsrt dawson dfridr erfcc gammln gaussj

hunt indexx lfit locate midinf midpnt

midsql odeint polint qromb qromo qsimp

rkck rkqs rtsafe selip shell spline

trapzd zriddr
```

The have all been adapted to the use of double precision floating-point number (instead of real). Furthermore, `ME(SSY)**2` also uses a modified verion of `rtsafe`, `rtsafe_rel` which has the following interface:

```
double precision function rtsafe_rel(funcd,x1,x2,xacc_rel)
```

`xacc_rel` is the required relative accuracy of the root determination. I.e. you need to use a test of the following sort:

```
if (abs(dx).lt.(xacc_rel*abs(rtsafe_rel)+1.0d-20)) then
```

at the appropriate place.

Very quick overview
====================

`ME(SSY)**2` evolves star cluster models. The code is based on the following set of core assumptions:

    -Spherical symmetry (hence no rotation, flattening, triaxiality, etc).
    -Dynamical equilibrium (no violent relaxation, etc).
    -Diffusive relaxation. The 2-body relaxation is treated in the Chandrasekhar approximation by assuming that it amounts to the effects of a large number of uncorrelated, small-angle 2-body encounters, each of which is a Keplerian hyperbolic deflection.

The physics included is:

    -Self-gravity of the cluster.
    -2-body diffusive relaxation.
    -Simple stellar evolution. The stars do not evolve on the MS. At the end of their MS lifetime, they turn into compact remnants.
    -Stellar collisions. One can inter- or extrapolate from a large database of SPH simulations.
    -Large-angle 2-body encounters (“kicks”). Those are supposed to be negligible in most circumstances but it was easy to include into the code.
    -Tidal truncation of the cluster (condition on apocentre distance).
    -A central massive object (generally thought to be a (I)MBH). It is considered fixed at the centre. 
    
    It interacts with the stellar system through:

    -Contribution to potential.
    -Tidal disruptions.
    -Direct plunges through horizon.
    -“Extreme mass-ratio inspiral” due to emission of gravitational waves. This process might not be treated accurately enough in the present implementation.

`ME(SSY)**2` does not include:

    -Giant stars (or any detailed stellar evolution).
    -Binaries.
    -Detailed treatment of tidal effects such as variable tidal field, delayed evaporation (see the paper of Fukushige and Heggie 2000, and also Baumgardt 2001), disk and bulge shocking.
    -Effects of gas (contribution to potential, accretion onto stars or the MBH, effects of gas removal from young clusters).

Units
=====

In most places (in particular for input and output), the code uses so called N-body units (see this paper of Hénon, 1971, but also the paper of Heggie and Mathieu of 1986). In the code we define the unit system such that the constant of gravity is G=1, the total stellar mass is initially M_{\rm cl}(0)=1, and the total initial stellar gravitational energy (not accounting for the contribution of the MBH to the potential) is -1/2 (Freitag and Benz 2001, Freitag, Amaro-Seoane and Kalogera 2006). As a time unit, we use the Fokker-Planck time T_{\rm FP} which is connected to the N-body time unit T_{\rm NB} through T_{\rm FP} = (N_{\star}(0)/\ln\Lambda) T_{\rm NB} were N_{\star}(0) is the initial number of stars and \ln\Lambda=\ln(\gamma N_{\star}(0)) the Coulomb logarithm (see e.g. Binney and Tremaine 2008 or Spitzer 1987). In the code we prefer to use T_{\rm FP} rather than T_{\rm NB} because the former is a relaxation time while the latter is a dynamical time. Individual stellar masses and radii (important for stellar evolution and collisions) are in M_{\odot} and R_{\odot}, respectively.

Compiling the code
====================

The code is written in more or less standard fortran-77. It is spread over a large number of files (about 153). The main file is EvolAmas.F and most other files have a name starting with EvolAmas_. The extension .f is for plain fortran files while .F files contain pre-processor directive allowing to customise the code at compile time. Most pre-processor variables are set in the make file EvolAmas.mke. The file Machine_Dependent.F contains all subroutines that might cause a problem when porting from one architecture to another or from a compiler to another. Their name end with _MD.

In general, files whose names end with _inc.f are include files, i.e., they are “incorporated” into other files through the “include” command. Names ending with _param.f and _common.f contain definitions of parameters (aka constants) and common blocks that are included in several places.

You will need the GNU tools. In particular the GNU make and the the GNU fortran compiler g77 or the more recent gfortran. I was able to compile the code on Solaris, various GNU/Linuxes PCs, Mac OSX and even on the Raspberry pi. Long time ago, I was successful using the Solaris fortran compiler, the Portland group compiler (on PCs), and the Intel compiler (on PCs). But I now stick to gnu compiler and EvolAmas.mke probably only works for g77 (if at all!).

Note that the make file might not work with non-GNU versions of make. I had to install gmake to “make” it work on probably the best OS around, OpenBSD ;). and that it needs the Z shell zsh to be installed (its path is set via the variable SHELL at the beginning). Another reason why you need zsh is that it is used by the scripts mke2list_preproc_defines.sh and CreateFortranPrintCodeForListPreprocDefines.sh called by the makefile. These scripts are called to write some fortran routine which will print out the values of most (all?) precompiler variables in CaracCode.asc when `ME(SSY)**2` is run.

If you are lucky, compiling `ME(SSY)**2` is as easy as

```
%make -f EvolAmas.mke
```

where `%` indicates the prompt. This should create the executable `__evolamas__` which you can happily rename to, say, `evolamas`.

However, many aspects of what the code will do, in particular the included physics, have to be decided before compilation as they are controlled by pre-processor variables. You can download as pdf here where I list the pre-processor variables that are set in EvolAmas.mke through the -D switch of the compiler. For instance, the variable _ECRIT_ETOILES_ is set to 1 in through -D_ECRIT_ETOILES_=1. Not all pre-processor variables are set in EvolAmas.mke. Some are (re)set in various .F files using the pre-processor “#define” directive.

As a general rule, only the pre-processor variables with a leading “*” in the table of the link should be tinkered with. Note that many of them determine what kind of information is output.

See the next items of the menu bar at the top of my web site for more information about `ME(SSY)**2`’s rich output.

Before compiling, you should also set the parameter `iDimSE` in `Param_Amas.f` to a value large enough to accommodate your number of particles. `iDimSE` is the size of arrays containing the particles’ data; it can be larger than the number of particles you will use (this might or might not result in more memory usage than required).

Input files and parameters
===========================

Over the years, `ME(SSY)**2` has grown quite complex but there has always been an effort to try to preserve compatibility with the data and parameter files of the previous versions. As a results, the input and output data are scattered in more files (and using more formats) that is really comfortable.

In this text we have a list of input and control files and explanations. These are files containing initial conditions (the structure of the cluster at t = 0) or parameters to set the physics or affecting the numerics. Note that only a very few are compulsory for `ME(SSY)**2` to run. The others can be used to override the default physics and numerics.

In the same list you will find listed files that can be used for some (limited) real-time control of the run, such as requesting a “snapshot” (saving all particle data to disk) or an early termination. Most importantly, there is a Perl script, gere_evolamas.pl which is driver for the fortran executable evolamas.

*The main parameter file: `input_EvolAmas`*

This is the most important file to set nearly all physical and numerical parameters. If one uses gere_evolamas.pl to drive the run, the script will look for a file named input_EvolAmas containing parameters and comments, and format it for evolamas. In this case one can include blank lines and anything starting with # is considered a comment. If gere_evolamas.pl is not used, the format is more strict (see below) and the preferred way to pass the parameters is by using the command:

```
%./evolamas -f my_parameter_file
```

where my_parameter_file is your parameter file. For simplicity, we will assume it is called input_EvolAmas. The stricter format for `input_EvolAmas` accepted by evolamas is something like that:

```
Param_name_1 # Some comments if you like
Param_value_1
Param_name_2
Param_value_2 # More comments if you please
Param_name_3
Param_value_3
```

This will give the value `Param_value_1` to the variable `Param_name_1` and so on. Note that one or more white space must appear before the value of a parameter and none before its name. In principle # is not required before a comment as only the first “word” of each line will be read anyway! evolamas will complain and stop if it doesn’t recognise the name of a parameter. If you are using gere_evolamas.pl, you can change the value of parameters on the command line:

```
%./gere_evolamas.pl --option Param_name_1=Param_value_1 --option Param_name_2=Param_value_2
```

If you use the switch `--Print` with `gere_evolamas.pl`, it will just output on stdout a cleaned-up version of the parameter file but not run evolamas. This is useful for systems where you cannot drive the code with such a script (such as a PC cluster using some queue submission software) but you still want gere_evolamas.pl to prepare a nice parameter file for you.

It should be noted that the code will try to start with the snapshot files (`MyRun%xxxxxxxxxx%AMAS.xdr`, etc) corresponding to the step number given by the parameter i ini (value of “xxxxxxxxxx”). Since one step correspond to the modification of just one pair of particles. Therefore, the step number can reach billions .

Hence, if you want to (re)start from t = 0 set i_ini to 0. If you want to continue a simulation interrupted at some step n (with existing corresponding snapshot files), set it to n. gere_evolamas.pl will not allow you to run a simulation in a directory where there are already snapshot files for a step number larger than 0 unless you use the switch `--Force` to force restart from scratch or `--Continue` to continue from the last snapshot in which case you don’t have to set i_ini manually.

Note that it is important that the number of digits in “xxxxxxxxxx” be exactly 10 (ten, dix, sepuluh). The code is too dumb to find the input file(s) if they do not follow this convention.

Initial condition files
========================

The file `MyRun%0000000000%AMAS.xdr` is one of the very few which is compulsory it contains the initial cluster structure, particle by particle. This data is the mass M, radius R, (specific) kinetic energy E and modulus of (specific) angular momentum J for each particle. These data are stored in the arrays M_SE, R_SE,T_SE and J_SE in a common bloc (defined in VarSE.f). As a note, SE comes from “Super-Etoile” the French for “super star” as it is how Hénon used to call particles to stress the fact that they may represent more than one star each.

In most cases, you will also need a file `MyRun%0000000000%ETOILES.xdr` containing the basic stellar properties of the stars of each particle, its mass M_{*} (in M_{\odot}), type and “date of birth” (in yr, since t = 0). The mass information is redundant with that in `MyRun%0000000000%AMAS.xdr` since M_{*} \propto M but when more detailled stellar evolution will be included, M_* will probably refer for the “zero-age” mass while M will track mass loss. The date of birth is generally when the star acquirred its present type but may be modified in collisions to account for “rejuvanation” (mixing in MS stars). It is used to determine the current (effective) evolutionary age. There are currently for stellar types: 1, 3, 4, 5 for MS, WD, NS and BH; no giant phase!.

These data are stored in the arrays `Met_SE`, (stellar mass) `DNet_SE`, (birth date), `iTet_SE` (stellar type) in a common block (defined in `VarSE.f`). The numerical type of `iTet_SE` is byte. All other arrays are `real*8` (aka double precision).

Typical usable number of particles range from 10^5 to 10^7 . 10^3 - 10^5 can be used for test runs. It is important to insist on the fact that the number of stars is completely independent of the number of particles (unlike in other current Hénon-type MC codes). The number of stars is given by the parameter `rNbEtoiles` in `input_EvolAmas`.

The easiest way to get initial condition files is to ask me or Marc (although Marc will probably not have time these days). The next easiest way is to use some of these programs to generate initial conditions:

`Creer_Amas_FdeE.F` can create `AMAS.xdr` according to a variety of theoeritcal models: Plummer, King, Wooley & Dickens (not tested), Isochrone, Dehnen with gamma = 0, 1, 2, and stellar polytropes (see Dehnen 1993). No central object can be included. Which model will be generated (as well as what type of parameters have to be entered) depends on the value of the pre-processor variable `_TYPE_` (see source file). When compiled for Plummer models, and called “create_plummer”, the command-line syntax is

```
%create_plummer -N 1000000 -Name MyPlummer [-iseed 763367]
```

to create a file `MyPlummer.xdr` containing one million particles. The option `-iseed` can be used to impose a different seed for the random number generator. When compiled for King models, and called “create_king”, the command-line syntax is

```
%create_king -N 1000000 -Name MyKing -W 8 [-iseed 763367]
```

were `-W` is used to give the value of the dimensionless central potential W_0 . In this case a tiny file `MyPlummer_MAREE.xdr` is also created to indicate the value of the tidal radius.

`Creer_Amas_EtaTN.F` can create “eta-models” (see Dehnen 1993, and also Tremaine et al 1994) with or without a central mass. This includes Hernquist and Jaffe models. The syntax is

```
%create_eta -N 1000000 -Name MyEta -eta 1.5 -mu 0.03 [-iseed 763367]
```

where `-eta` if for the parameter \eta (the density at small radii is n \propto R^{\eta -3}) and -mu is for the parameter \mu, the mass of the central object (in units of the total stellar mass, so the total mass of the system is 1+ \mu). If \mu > 0, a file `MyPlummer_TN.xdr` is created to contain the value of the mass of the central object.

These models will have no mass spectrum (and no ETOILES.xdr file). To introduce a mass spectrum and create the corresponding `ETOILES.xdr` file, one can use AppliquerFM3.F which allows the user to apply any initial mass function that can be represented as a piecewise power-law (see Appendix A2 of Freitag and Benz 2002).

Explanations about usage are given in the source code. For a Kroupa IMF, you would do

```
%appliquerfm3 -m "0.01,0.08,0.5,150" -a "0.3,1.3,2.3" MyModel MyModel_WithMF
```

so `-m` gives the list of limit-masses, -a the list of exponents (2.35 is Salpeter), `MyModel` the base name of the input model without mass function and `MyModel_WithMF` the base name of the output model. At the minimum files `MyModel_WithMF.xdr` and `MyModel_WithMF_ETOILES.xdr` will be created.

The radii might be rescaled a bit to enforce strict virialisation and keep with the definition of N−body units, therefore a file `MyModel_WithMF_MAREE.xdr` can be created to contain the value of the rescaled tidal truncation radius.

To start a run, create a new directory. Inside copy your initial condition files with the following change of names:

```
MyModel.xdr → MyRun%0000000000%AMAS.xdr
```

and

```
MyModel_ETOILES.xdr → MyRun%0000000000%ETOILES.xdr,
```

```
MyModel_MAREE.xdr → MyRun%0000000000%MAREE.xdr,
```

etc if it applies. Note that you need exactely 11 (eleven) “0” in the filenames. If you use gere_evolamas.pl, you also need to put the executable `evolamas` and `input_EvolAmas` into this directory.

Output files
=============

`ME(SSY)**2` can produce a large variety of output files. 

Many of them are for specialised application or testing purpose and only a few are of general use. FHR, the formats are xdr, to preserve the full numerical precision of some data which might be crucial for restart. By the way, saving 8-byte real data in ASCII with 15 decimal positions also preserves the full precision. Once compressed with, e.g., gzip, such files are nearly as compact as their xdr counterpart. I am planning on getting rid completely of the xdr dependence of `ME(SSY)**2`, but it might take a while.

There are also ASCII files for readability for data which are most useful for monitoring the run and plotting results. Most files are incremental, with information appended to the file as the simulation proceeds. The important exception are “snapshot files” that contain (in principle) the whole information about all particles and can be used for restart or detailed analysis. In this case a different set of files is written each time a snapshot is dumped, with file name indicating the step number, for instance `MyRun%00001200000%AMAS.xdr` for the mass and orbit information of all particles at step 1200000.

The most important output file is the “master log file” MyRun%%Log.asc. It contains a wealth of information about the evolution of the simulation, from a numerical and physical point of view. The same information is also sent to stderr to entertain the user. The format is supposed more human- than computer-friendly and a script, `Log2rdb.pl`, can be used to turn this file into a flat table format. The output format is “rdb” with a two line header and tabs used to separate columns. It can be turned into the usual ASCII format we use for output from fortran code with another script: rdb2fort.sh.

The information contained in `MyRun%%Log.asc` is presented now.

For the time, I recommend to use the variable `Tps_Ufokpl` which gives it in Fokker-Planck units. In most other output file, no time information is given, only the step number (`iPas_Evol` or similar name) or, if a time is given, it might not be in the proper units! Therefore it is safer to use `MyRun%%Log.asc` to derive the `iPas_Evol → Tps_Ufokpl` relation and used it determine time from step number.

Most other ASCII files have a rather straightforward flat format. The contain data in columns and a header to explain what the data is. The header are the first lines, starting with #. The last of these lines generally looks like

```
# 1: iPas_Evol 2: this_var 3: that_var 4: this_other_quantity
```

which means that the first column contains (the successive values of) iPas_Evol, the second contains this_var, etc. 

`MyRun%%RayLag.asc` is an important output file whose header has a slightly different structure. The last header line looks like

```
# 1: iPas_Evol 2: NbSE_liees 3: M_liee 4: Tps_amas || FracMasse : 5: R0001 (0.001000) \
6: R0010 (0.010000) 7: R0020 (0.020000)
```

which indicates that the 5th column contains the radius of the sphere containing 0.1% of the total mass, the 6th column the radius for 1% and so on. These quantities are known as “Lagrange radii” and a sphere containing a given fraction of the (remaining) cluster mass is a “Lagrange sphere”. Similarly, the file MyRun%%Segr.asc has the following header

```
# 1: iPas_Evol 2: Tps_amas || Masse stellaire moyenne : 3 : m00050 (M<=0.0050) \
4 : m00100 (M<=0.0100) 5 : m00500 (M<=0.0500)
```

indicating that the third column contain the average stellar mass (in M_{\odot}) for all particles inside the 0.5% Lagrange sphere, etc.

If the code detects a error situation it will try and write a lot of information (such as a last snapshot) in files with name starting with `_RIP_` (e.g., `_RIP_%00051278532%AMAS.xdr`) before it exits to help post-mortem investigations into the origin of the problem.

For more information
====================

Please visit

<a href="http://astro-gr.org/monte-carlo-simulations-for-stellar-dynamics">http://astro-gr.org/monte-carlo-simulations-for-stellar-dynamics</a>

Pau Amaro Seoane, Berlin 12/03/2016
