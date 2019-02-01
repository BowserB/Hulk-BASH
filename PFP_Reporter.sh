#!/bin/bash

#This is to take the .ALLRESULTS.txt files from the ProteinFinderPipeline and convert them into binary reports

#set default variables
OUTPUT_DEFAULT='NOTSET'
LIST_DEFAULT='~/brad/U12proteinSearch/HMMdonePipeList.txt'
PSIDQSCORE_DEFAULT='P'
POSITIVE_DEFAULT='+'
NEGATIVE_DEFAULT='-'

#set flags
while getopts :w:o:l:hsbrp:y:n: option
do
case "${option}"
in
h)
echo "This script takes the scoring output from ProteinFinderPipe (.ALLRESULTS.txt) and makes it a binary report"
echo "	Usage:"
echo "	./PFP_Reporter.sh -w <WORKING DIRECTORY>"
echo "	where the WORKING DIRECTORY is the location of the ALLRESULTS.txt files"
echo " 	if a list of files is used, specify with the -d flag"
echo "	Available flags:"
echo "	-o 	Output directory"
echo "	-l	List of Organisms/ Directories in txt file"
echo -e "	-b	HMM \033[4mB\033[0mypass	Candidates disqualified due to failing HMM criteria are accepted as positive"
echo -e "	-s	\033[4mS\033[0mize bypass	Candidates disqualified due to failing Size criteria are accepted as positive"
echo -e "	-r	\033[4mR\033[0mecip bypass	Candidates disqualified due to failing Recip criteria are accepted as positive"
echo "	-y	set symbol used to indicate positive in report (default = "+")"
echo "	-n	set symbol used to indicate negative in report (default = "-")"
exit 0;;
w) WORKING_DIR=${OPTARG};;
o) OUTPUT=${OPTARG};;
l) LIST=${OPTARG};;
b) HMMBYPASS=HB;;
s) SIZEBYPASS=SB;;
r) RECIPBYPASS=RB;;
p) PSIDQSCORE=${OPTARG};;
y) POSITIVE=${OPTARG};;
n) NEGATIVE=${OPTARG};;
esac
done

#setting default values if variable is undefined

: ${OUTPUT=$OUTPUT_DEFAULT}
: ${LIST=$LIST_DEFAULT}
: ${PSIDQSCORE=$PSIDQSCORE_DEFAULT}
: ${POSITIVE=$POSITIVE_DEFAULT}
: ${NEGATIVE=$NEGATIVE_DEFAULT}

#This is the list of organisms that reports have been generated
while read ORGANISM
do
#set up folders
if [ "$OUTPUT" = "NOTSET" ]
then
	OUTPUT=$WORKING_DIR
fi

mkdir $OUTPUT/$ORGANISM/REPORTGEN
cd $WORKING_DIR/$ORGANISM

	while read ALLRESULTS
	do
	if [ "$ALLRESULTS" = "++" ]
	then
		echo "$POSITIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
	elif [ "$ALLRESULTS" = "+" ]
	then
		echo "$POSITIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
	elif [ "$ALLRESULTS" = "H" ]
	then
		if [ "$HMMBYPASS" = "HB" ]
		then
			echo "$POSITIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
		else
			echo "$NEGATIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
		fi
	elif [ "$ALLRESULTS" = "mH" ]
    then
        if [ "$HMMBYPASS" = "HB" ]
        then
            echo "$POSITIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
        else
            echo "$NEGATIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
        fi
	elif [ "$ALLRESULTS" = "S" ]
	then
		if [ "$SIZEBYPASS" = "SB" ]
		then
			echo "$POSITIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
        else
			echo "$NEGATIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
		fi
	elif [ "$ALLRESULTS" = "mS" ]
    then
        if [ "$SIZEBYPASS" = "SB" ]
        then
      		echo "POSITIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
        else
            echo "$NEGATIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
        fi
	elif [ "$ALLRESULTS" = "R" ]
	then
		if [ "$RECIPBYPASS" = "RB" ]
		then
			echo "$POSITIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
        else
			echo "$NEGATIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
		fi
 	elif [ "$ALLRESULTS" = "mR" ]
    then
        if [ "$RECIPBYPASS" = "RB" ]
        then
            echo "$POSITIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
        else
            echo "$NEGATIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
        fi
	elif [ "$ALLRESULTS" = "$PSIDQSCORE" ]
	then
		echo "$NEGATIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
	else
		echo "$NEGATIVE" >> $OUTPUT/$ORGANISM/REPORTGEN/$ORGANISM.BINARYREPORT.txt
	fi
	done < $ORGANISM.ALLRESULTS.txt

done < $LIST 
