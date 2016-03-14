#!/bin/zsh

# Remove the snapshots corresponding to a time larger than some value (2e20 yrs by default)
# from all the MESSY simulations whose directories are given as parameters on the command line.
# The value of the time (in years) is given via the "-t" option (before the names of directories).

# The snapshots are not actually deleted but put in a subdirectory _TRASH_  in the
# directory of the simulation

LC_ALL=C
export LC_ALL
setopt EXTENDED_GLOB

time_end_yr=2e10

while [[ $1 = -* ]]; do
    case $1 in
        -t ) time_end_yr=$2; shift 2;;
	-(d|dry-run) ) dry_run=1; shift 1;;
        *       )
              print -u2 "!!! What exactly do you mean with $1, you magnificient weirdo? !!!"
              exit 1 ;;
    esac
done

dir_list=( $@ )

base_dir=$(pwd)
for dir in $dir_list; do
    print -u2 "> dealing with $dir"
    cd $dir

    # Using log file, find first step number corresponding to time larger than 
    # requested end time

    LogFile=$(echo *%%Log.asc(|.gz)(.N)) 
       # the (.N) means 'set option NULL_GLOB' for this pattern 
       # so that the shell does not produce an error if file is not found. 
       # We'll just have an empty variable.

    if [[ -z $LogFile ]]; then
	print -u2 "!!! Can not find log file in $dir !!!"
    else
	BaseName=${LogFile%%"%%"*}
	#print -u2 " > BaseName : $BaseName"
	step_end=$(
	    gunzip --force --to-stdout $LogFile |\
	    gawk -v time_end=$time_end_yr '
               # Build arrays of step numbers and times
               /^\| *iPas_Evol / {iPas_Evol=$NF} 
               / Tps_en_yr / {Tps_en_yr=$NF; arr_i[i]=iPas_Evol; arr_t[i]=Tps_en_yr; i++}
              END{ 
                 # Now analyse the arrays to find the step number 
                 # (while allowing for restarts)
                 N=i
                 j=1e30
                 step_found=-1
                 time_found=-1.0
                 for (i=N-1; i>=0; i--) {
                    if (arr_i[i]<j) {
                       if (arr_t[i]>time_end) {
                          step_found=arr_i[i]
                          time_found=arr_t[i]
                       } else {
                          break
                       }
                       j=arr_i[i]
                    }
                 }
                 print step_found
             }'
	)

	# Now move snapshots corresponding to times later than time_end_yr (except the 
        # first one with t > time_end_yr) to the _TRASH_ subdirectory.

        if (( step_end < 1 )); then
	    print -u2 "!!! The simulation in $dir doesn't seem to have even reached $time_end_yr yrs !!!"
	else
	    unset delete
	    if [[ -z $dry_run ]]; then
		[[ -d _TRASH_ ]] || mkdir _TRASH_ || {
		    print -u2 "!!! Could not creat directory _TRASH_ in $dir !!!"
		    exit 1
		}
		print -u2 " Moving files to delete to _TRASH_ subdirectory of $dir"
	    fi
	    for file in $(ls $BaseName%([0-9])##%AMAS.xdr | sort); do
	        step=$(( ${${file%"%"*}#*"%"} ))
	        [[ -n $delete ]] && {
		    [[ -n $dry_run ]] && {
			print -u2 " > Should delete files " ${file%AMAS.xdr}*
		    } || {
			mv ${file%AMAS.xdr}* _TRASH_
                    }
                }
                (( step > step_end )) && delete=1
            done
	fi
    fi
    cd $base_dir
done

