#!/bin/bash

# this script merge event.lhe files produced in parallel runs by 
# removing the first 30 lines and the lasts lines and then 
# concatenate them.

usage() 
{
    echo " 
Usage:
    ${0} output.file input.dirs ...
    Merge event.lhe file in input.dirs into one output lhe file."
}

if [ $# -lt 2 ]; then
    echo "Short of optioins"
    usage
    exit 1
fi

OUTPUT="$1"
if [ -a $1 ]; then
    echo "The spedified output file already exist, do you want to override it? [y/n]"
    read yesno
    case yesno in
	[yY]*)
	    echo "override eixsting file"
	    ;;
	[nN]*)
	    echo "output file exist, please specify another one"
	    exit 1
	    ;;
	*)
	    echo " Error input."
	    exit 1
	    ;;
    esac
fi

shift
while [ $# -gt 0 ]; do
    if ! [  -d $1 ] || ! [ -f "$1/event.lhe" ]; then
        echo "Input option $1 is not a dir. or file event.lhe doesn't exist within this dir.  ignore it"
	shift
	continue
    fi
    
    cat $1/event.lhe | sed '1,30 d' | sed '$ d' >> ${OUTPUT}
    shift
done

