#!/usr/bin/env bash
#
# This shell script takes a list of files or directories and reports whether
# they are on disk or have been truncated.
#
# !! It is useful only for HSM managed filesystems !!
#
#
# Usage: onDisk.sh [-q|-n] <list of files | directories>
#
# Description: 
# The decision is based on the difference reported between the size in bytes
# and the size in kilobytes as reported by 'find'. Sizes in kilobytes are based
# on the number of blocks effectively used by the file on the filesystem, and 
# therefore reflects truncation state of the file.
#
# Author:      Volker Flegel
#
# Date:  19.04.2011
#
# History:
#   - 19.04.2011: First production release.

# Get program name
prgName=`basename "$0"`

################################################################################
# Helper functions
Usage () {
    echo "$prgName - Determine if HSM files are on disk or truncated.

Usage: $prgName [-t] [-d] [-q] [-n] <list of files | directories>
  <files>       List of filenames to be checked.
  <directory>   List of directories to be checked, all files in all 
                subdirectories are also checked.
                You can mix filenames and directories on the command line.

 Options:
   -q          Quiet. No output is made on stdout.
               The program's exit code is set according to:
               $prgName returns 0 (shell 'true') if all files are on disk.
               $prgName returns 1 (shell 'false') if at least one file is truncated.
   -n          Numeric output.
               $prgName prints the filename followed by 0 if it is on disk,
               or followed by 1 if the file is truncated.
               Normal output prints 'on_disk' or 'truncated'.
   -t          Output only files that are truncated.
               Files on disk are not displayed.
   -d          Output only files that are on disk.
               Truncated files are not displayed.

 Note:
   By using both options '-n' and '-q' together $prgName will display only the
   number of truncated files found.

 Exit codes:
   $prgName always returns 0 (shell 'true') if all files are on disk.
   $prgName always returns 1 (shell 'false') if at least one file is truncated.
  
Warning: This script is useful only if used on an HSM managed filesystem." >&2
    exit 1
}

# Msg(): Write messages to stderr
Msg () {
  for MsgLine
  do
    echo -e "$prgName: $MsgLine" >&2
  done
}

# Fatal(): Write messages to stderr and exit
Fatal () { Msg "$@" "Exiting..."; exit 1; } 

################################################################################

################################################################################
# main() 
################################################################################

# Set default values
FINDCMD=`which find`
# Find options (strange writing due to $IFS modification)
FINDOPT[0]="-xdev"
FINDOPT[1]="-type"
FINDOPT[2]="f"
FINDOPT[3]="-printf"
FINDOPT[4]="%p\\n\$(( %s / 1024 - %k ))\\n"
QUIET=0
NUMERIC=0
TRUNC=1
DISK=1

# Set variables
p_ondisk="on_disk"
p_trunc="truncated"
nb_trunc=0
fname=0

################################################################################
# Parse command line arguments
while [ $# -gt 0 ]
do
  case "$1" in
  -q)     quiet=1;;                 # Quiet, no output
  -n)     numeric=1;;               # Numeric output
  -t)     disk=0;;                  # Display only truncated files
  -d)     trunc=0;;                 # Display only files on disk
  --)     shift; break;;
  -h)     Usage;;
  -*)     Usage;;
  *)      break;;	                  # Expect: Filenames / Dirnames
  esac
  shift
done

################################################################################
# Check command line arguments
# We still need one command line option: file or directory names
[[ $# -lt 1 ]] && Usage

# Set behaviour based on options
: ${quiet:=$QUIET}          # Set default or CLI verbosity
: ${numeric:=$NUMERIC}      # Set default or CLI output format
: ${trunc:=$TRUNC}          # Set default or CLI truncated file display
: ${disk:=$DISK}            # Set default or CLI file on disk display

# Test for conflicting options
if [[ $trunc -eq 0 && $disk -eq 0 ]] 
then
  Fatal "Please give only one of options '-t' or '-d'."
fi
if [[ $quiet -ne 0 && $trunc -eq 0 ]] || [[ $quiet -ne 0 && $disk -eq 0 ]]
then
  Fatal "Option '-q' conflicts with option(s) '-t' or  '-d'."
fi
#if [[ $quiet -ne 0 && $disk -eq 0 ]]
#then
#  Fatal "Option '-q' conflicts with option '-d'."
#fi
# Set output text
if [[ $numeric -ne 0 ]]
then
  p_ondisk="0"
  p_trunc="1"
fi

################################################################################
# Loop over output from the 'find' command.
#  The 1st line of output from 'find' should be the filename,
#  the 2nd line should be '$(( <size_b> / 1024 - <size_kb> ))'
OLDIFS=$IFS
IFS=$'\n'

for line in `$FINDCMD "$@" ${FINDOPT[@]}`
do
  if [[ $fname -eq 0 ]]
  then
    # Output line is the filename
    filename="$line"
    fname=1
  else
    # Output line is a bash Arithmetic expression
    eval cpt_size=$line
    fname=0
    # If the size is greater than 0 it should be considered truncated
    # If option 'quiet' just exit as soon as a truncated file is found
    if [[ $quiet -ne 0 && $numeric -ne 0 && $cpt_size -gt 0 ]]
    then
      let nb_trunc++
      continue
    fi
    [[ $quiet -ne 0 && $cpt_size -gt 0 ]] && exit 1
    [[ $quiet -ne 0 ]] && continue

    # Not quiet: print output
    #echo -e -n "$filename:\t"
    [[ $cpt_size -gt 0 && $trunc -ne 0 ]] && echo -e "$filename:\t$p_trunc" && let nb_trunc++
    [[ $cpt_size -le 0 && $disk -ne 0 ]] && echo -e "$filename:\t$p_ondisk"
  fi
done

# If option '-q' _and_ '-n' are specified, print the number of truncated files
if [[ $quiet -ne 0 && $numeric -ne 0 ]] 
then
  echo $nb_trunc
fi

# Set correct exit code (1 if files were truncated, 0 else)
[[ $nb_trunc -gt 0 ]] && exit 1
exit 0


