#/bin/bash
# stript_LHE.sh

# stript fake events ( events not containing 13 lines per event )

ARGS_1=1	# one arguments, only input file
ARGS_2=2	# two arguments, input and output file
E_BADARGS=34	# wrong argument number
E_BADINPUTFILE=35   # input file not LHE file ( suffix .lhe)
E_OUTPUTEXIST=36
E_BADCHOICE=37

usage()
{
    echo    "Usage: stript_LHE.sh input.lhe [output.lhe]"
    echo -e "\tif output file not set, default output is input_strip.lhe"
}

# check parameters
if [ $# -eq "$ARGS_1" ]
then
    INPUT="$1"
    OUTPUT=${INPUT%.lhe}
    OUTPUT="$OUTPUT"_strip.lhe
    echo "Only one argument, read as input lhe file, set output file name to $OUTPUT"
elif [ $# -eq "$ARGS_2" ]
then
    INPUT="$1"
    OUTPUT="$2"
    echo "input LHE file is $INPUT, output file is $OUTPUT"
else
    usage 
    exit $E_BADARGS
fi

# check input LHE file ( require suffix with .lhe )
if [ "${#INPUT}" -lt 4 ] || [ "${INPUT: -4}" != ".lhe" ]
then
    echo "Irregular input LHE file"
    exit $E_BADINPUTFILE
fi

# check output file
if [ -f "${OUTPUT}" ]
then
    echo "$OUTPUT exists, would you like to cover it: [y/n]"
    read -t 60 yesno  # timeout 60s
    case $yesno in
	y)
	    echo "Erase existing $OUTPUT, and create new one."
	    rm $OUTPUT
	    ;;
	n)
	    echo "You don't want to cover old $OUTPUT, exit"
	    exit $E_OUTPUTEXIST
	    ;;
	*)
	    echo "Undefined choice, exit"
	    exit $E_BADCHOICE
	    ;;
    esac
fi


awk -v RS='</event>\n' '{ if(NF == 137) print $0"</event>"}' $INPUT > $OUTPUT 
