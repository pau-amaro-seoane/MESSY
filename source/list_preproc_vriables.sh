#!/usr/bin/env zsh
for file in $(grep -l '^\#define ' *.F); do; echo; echo "=== $file ==="; cat $file | gawk '/^\#define / {print $2,$3}'; done
