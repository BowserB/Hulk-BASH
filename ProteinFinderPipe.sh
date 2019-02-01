#!/bin/bash

#This is setting up the script using user defined flags
#THINGS TO DO: figure out directory assignment/ defaults.  Have a flag to save sequences as fasta or not.  Have a flag for results report?
#i=psIresults l=proteinList

#Set default variables
OUTPUT_DEFAULT="NOTSET"
ORGANISM_LIST_DEFAULT=~/brad/U12proteinSearch/HMMdonePipeList.txt
PROTEIN_LIST_DEFAULT=~/brad/SCRATCH/psiOf6List.txt
PBRESULTS_DEFAULT=~/brad/U12proteinSearch/SCOTTpsiBLAST/COPY
PROTEOME_DIR_DEFAULT=~/brad/U12proteinSearch/SCOTTpsiBLAST
REFERENCE_ORG_DEFAULT=GCF_000001405.38_GRCh38.p12_protein.faa
REFERENCE_DEFAULT=~/brad/U12proteinSearch/SCOTTpsiBLAST/COPY/GCF_000001405.38_GRCh38.p12_protein.faa
LENGTHofREF_DIR_DEFAULT=~/brad/U12proteinSearch/SCOTTpsiBLAST/COPY/GCF_000001405.38_GRCh38.p12_protein.faa/domtblout
TOTHITS_DEFAULT="20"
REDUCEDHITS_DEFAULT="10"
ALLOWED_SIZE_RANGE_DEFAULT=".25"
HMM_RESULTS_DIR_DEFAULT=~/brad/U12proteinSearch/SCOTTHMM/RESULTS
DOMAIN_LISTS_DIR_DEFAULT=~/brad/SPLICEOSOMEproteins/DomainSets
SAVEFASTA_DEFAULT="NO"
RECIPLIMIT_DEFAULT="2"

#Set flags
while getopts :hw:o:l:i:p:c:t:bsrfn:u:a:m:g:d:e:j: option
do
case "${option}"
in
h)
echo "This script analyzes data from psiBLASTs, HMMsearches, reciprocal BLASTp's, and Size matching to identify proteins"
echo "	Usage:"
echo "	./ProteinFinderPipe.sh -w <WORKING DIRECTORY>"
echo "	where the WORKING DIRECTORY is the location of the psiBLAST results files (output format 6 is required)"
echo " 	For HMMsearch results to be parsed, HMMsearch output must be in domtblout format and in the directory set with the -m flag."
echo "	Available flags:"
echo "	-o 	Output directory (default is WORKING_DIRECTORY"
echo "	-l	List of Organisms/ Directories in txt file"
echo -e "	-i	Location of ps\033[4mi\033[0mBLAST results (in -outputfmt 6)"
echo -g "	-g	Organism list (text file)."
echo -e "	-m	H\033[4mM\033[0mM Results directory (in domtblout format)."
echo "	-p	Directory containing the Proteomes of the query organisms."
echo -e	"	-c	Location of Referen\033[4mc\033[0me organism proteome (for Recip pBLASTing)."
echo -e "	-e	R\033[4me\033[0mference organism proteome filename"
echo -e "	-t	Leng\033[4mt\033[0mh of Reference proteins directory (for size comparison)."
echo "	-n	Number of total hits pulled from psiBLAST results (default is 20)"
echo "	-u	Number of red\033[4mu\033[0mced hits after de-duplication (default is 10)."
echo "	-a	Allowed size range for comparison (percentage in decimal format - default is .25)."
echo -e "	-b	HMM \033[4mB\033[0mypass	Candidates disqualified due to failing HMM criteria are accepted as positive"
echo -e "	-s	\033[4mS\033[0mize bypass	Candidates disqualified due to failing Size criteria are accepted as positive"
echo -e "	-r	\033[4mR\033[0mecip bypass	Candidates disqualified due to failing Recip criteria are accepted as positive"
echo -e "	-j	Limit for accepting a reciprocal hit (ranking) against Reference proteome."
exit 0;;

\?) echo "Invalid flag, use -h for help"
exit 1;;

w) WORKING_DIR=${OPTARG};;
o) OUTPUT=${OPTARG};;
l) PROTEIN_LIST=${OPTARG};;
i) PBRESULTS=${OPTARG};;
p) PROTEOME_DIR=${OPTARG};;
c) REFERENCE=${OPTARG};;
t) LENGTHofREF_DIR=${OPTARG};;
b) HMMBYPASS=HB;;
s) SIZEBYPASS=SB;;
r) RECIPBYPASS=RB;;
f) SAVEFASTA=SAVE;;
n) TOTHITS=${OPTARG};;
u) REDUCEDHITS=${OPTARG};;
a) ALLOWED_SIZE_RANGE=${OPTARG};;
m) HMM_RESULTS_DIR=${OPTARG};;
g) ORGANISM_LIST=${OPTARG};;
d) DOMAIN_LISTS_DIR=${OPTARG};;
e) REFERENCE_ORG=${OPTARG};;
j) RECIPLIMIT=${OPTARG};;

esac
done

#Setting default values if variable is undefined

: ${OUTPUT=$OUTPUT_DEFAULT}
: ${PROTEIN_LIST=$PROTEIN_LIST_DEFAULT}
: ${PBRESULTS=$PBRESULTS_DEFAULT}
: ${PROTEOME_DIR=$PROTEOME_DIR_DEFAULT}
: ${REFERENCE=$REFERENCE_DEFAULT}
: ${LENGTHofREF_DIR=$LENGTHofREF_DIR_DEFAULT}
: ${TOTHITS=$TOTHITS_DEFAULT}
: ${REDUCEDHITS=$REDUCEDHITS_DEFAULT}
: ${ALLOWED_SIZE_RANGE=$ALLOWED_SIZE_RANGE_DEFAULT}
: ${HMM_RESULTS_DIR=$HMM_RESULTS_DIR_DEFAULT}
: ${ORGANISM_LIST=$ORGANISM_LIST_DEFAULT}
: ${DOMAIN_LISTS_DIR=$DOMAIN_LISTS_DIR_DEFAULT}
: ${REFERENCE_ORG=$REFERENCE_ORG_DEFAULT}
: ${SAVEFASTA=$SAVEFASTA_DEFAULT}
: ${RECIPLIMIT=$RECIPLIMIT_DEFAULT}

#this is the list of organisms
while read ORGLIST; do
#setting up the directories to store the results, temp files, and sequence files.

if [ $OUTPUT = "NOTSET" ]
then
	OUTPUT=$WORKING_DIR
fi

mkdir $OUTPUT/"$ORGLIST"
cd $OUTPUT/"$ORGLIST"
mkdir $OUTPUT/"$ORGLIST"/RECIP
if [ "$SAVEFASTA" = "SAVE" ]
then
mkdir $OUTPUT/"$ORGLIST"/PROTEINS
fi
mkdir $OUTPUT/"$ORGLIST"/CATEGORIES
paste "$PROTEIN_LIST" > FINALOUTTEMP.txt

#this is performed on psiBLAST outputfmt 6 results. List file is in ~/brad/SCRATCH/psiOf6List.txt
	while read psiRESULTS; do

#this evaluates if there is a psiBLAST hit or no and records a value on the running report. If there are no psiBLAST hits, the script exits. 

		if [ $(ls -l $PBRESULTS/$ORGLIST/$psiRESULTS | tr -s ' ' | cut -f5 -d ' ') = '0' ];
		then
			echo 0 > $psiRESULTS.categories.temp
			paste -s <(echo "P") >> $ORGLIST.ALLRESULTS.txt
			paste -s <(echo "-") >> $ORGLIST.REPORT.txt
			paste -s <(echo "0") >> $ORGLIST.HMMRESULTS.txt
continue
		else echo 1 > $psiRESULTS.categories.temp
		fi
#this parses the psiBLASToutput and pulls out the seqID
		cut -f2 $PBRESULTS/$ORGLIST/$psiRESULTS > zzSeqID2.temp;

#This only takes the top 20 (now VAR)  hits from the psiBLAST results
		head -n$TOTHITS zzSeqID2.temp > zzSeqID.temp

#this prevents the 'Search has converged' message from being in results and also sort -u removes duplicates
			if grep -q 'Search' zzSeqID.temp ; then head -n-2 zzSeqID.temp | sort -u > zzSeqID2.temp; else sort -u zzSeqID.temp > zzSeqID2.temp; fi

#NEW this cuts the hits down to 10
		head -n$REDUCEDHITS zzSeqID2.temp > zzSeqID.temp

#this pulls the list id from the title of the query
			HMMID=$(echo $psiRESULTS | cut -f3 -d. );

#this pulls the protein sequence from the proteome. HITS is from zzSeqID list.
			while read HITS; do
				sed -n '/'$HITS'/,/>/p' $PROTEOME_DIR/$ORGLIST | sed \$d > $OUTPUT/$ORGLIST/$HITS.temp;

#NEW this compares the size of the CANDIDIATE from the REF (human). This also stores the query for the recipBLAST
				if grep -q 'NP' <(head -n1 $REFERENCE/$psiRESULTS);
					then REFID=$(head -n1 $REFERENCE/$psiRESULTS | cut -f2);
					else REFID=$(head -n2 $REFERENCE/$psiRESULTS | grep 'NP' | cut -f2);
				fi
# NEW this is continuation of above, defining the size of the REFLENGTH and the CANDIDATE from the HMMresults
#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.1
				 sed -n '/'$REFID'/,/>/p' $PROTEOME_DIR/$REFERENCE_ORG | sed \$d >$WORKING_DIR/$ORGLIST/REFSEQ.temp

			REFLENGTH=$(tail -n+2 $WORKING_DIR/$ORGLIST/REFSEQ.temp | tr -d '\n' | wc -c )
#	REFLENGTH=$(grep $REFID $LENGTHofREF_DIR/*.txt | head -n1 | tr -s ' ' | cut -f3 -d ' ')

# NEW this is using the actual protein sequence to get the length of the candidate.  This eliminates problem above if there are no found domians
			CANDIDATE=$(tail -n+2 $WORKING_DIR/$ORGLIST/$HITS.temp | tr -d '\n' | wc -c )
#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.1.5
#Comparing the sizes
COMPARE=$(echo "(sqrt (( $REFLENGTH - $CANDIDATE )^2)/$REFLENGTH)"| bc -l)
    if (( $(echo "$COMPARE <= $ALLOWED_SIZE_RANGE" |bc -l) ));
        then
            echo $HITS >> $psiRESULTS.RecipReady
            cp $WORKING_DIR/$ORGLIST/$HITS.temp  $WORKING_DIR/$ORGLIST/$HITS.candidate.temp
			cp $psiRESULTS.categories.temp $psiRESULTS.$HITS.categories.txt
			echo 2 >> $psiRESULTS.$HITS.categories.txt
# ADDING size report
			paste -s <(echo $psiRESULTS.$HITS) <(echo $REFLENGTH) <(echo $CANDIDATE) <(echo $COMPARE) >> $ORGLIST.SizeReport.txt
#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.2


#RecipBlast with Binary outs. HITID created from RecipReady list above and is a subset of the HITS
blastp -query $WORKING_DIR/$ORGLIST/$HITS.candidate.temp -db $PROTEOME_DIR/$REFERENCE_ORG -evalue=0.00001 -outfmt 6 -out $OUTPUT/$ORGLIST/RECIP/$HITS.BLAST.txt

	if grep $REFID <(head -n$RECIPLIMIT RECIP/$HITS.BLAST.txt);
		then
			if [ "$SAVEFASTA" = "SAVE" ]
				then
					mv $WORKING_DIR/$ORGLIST/$HITS.candidate.temp  $OUTPUT/$ORGLIST/PROTEINS/$psiRESULTS.faa
				fi
			echo $HMMID.$HITS >> RECIPEDLIST.txt
			echo 4 >> $psiRESULTS.$HITS.categories.txt
		else
			echo 0 >> $psiRESULTS.$HITS.categories.txt
fi
#This else is reporting the size not being in the range tested above
	else
		cp $psiRESULTS.categories.temp $psiRESULTS.$HITS.categories.txt
		echo 0 >> $psiRESULTS.$HITS.categories.txt
#For those that are not the right size, then open it up to reciping ALL.
blastp -query $WORKING_DIR/$ORGLIST/$HITS.temp -db $PROTEOME_DIR/$REFERENCE_ORG -evalue=0.00001 -outfmt 6 -out $OUTPUT/$ORGLIST/RECIP/$HITS.BLAST.txt

	if grep $REFID <(head -n$RECIPLIMIT RECIP/$HITS.BLAST.txt);
		then
			if [ "$SAVEFASTA" = "SAVE" ]
				then
					mv $WORKING_DIR/$ORGLIST/$HITS.candidate.faa  $OUTPUT/$ORGLIST/PROTEINS/$psiRESULTS.faa
			fi
			echo $HMMID.$HITS >> RECIPEDLIST.txt
			echo 4 >> $psiRESULTS.$HITS.categories.txt
#ADDING size report
            paste -s <(echo $psiRESULTS.$HITS) <(echo $REFLENGTH) <(echo $CANDIDATE) <(echo $COMPARE) >> $ORGLIST.SizeReport.txt
		else
			echo 0 >> $psiRESULTS.$HITS.categories.txt
fi
fi

#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.3

# this searches the HMM outs with HMM from HMMID list, temp file appended to be pasted from later.
#Added desigination of HMM variable for later summary

					while read HMMLIST; do
				if grep $HITS $HMM_RESULTS_DIR/$ORGLIST/domtblout/$HMMLIST.txt
					then
#						echo 8 >> $psiRESULTS.$HITS.categories.txt
						echo "1"
					else
						echo "0"
#						echo 0 >> $psiRESULTS.$HITS.categories.txt
				fi | tail -n1 >> $psiRESULTS.$HITS.wc.temp


#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.4

# this is testing if there are hits for each of the expected HMMs on the list

if [ $(wc -l <$DOMAIN_LISTS_DIR/$HMMID.txt) = $(paste -s -d+ $psiRESULTS.$HITS.wc.temp | bc) ];
then
		if [ "$SAVEFASTA" = "SAVE" ]
			then
				mv $WORKING_DIR/$ORGLIST/$HITS.temp  $OUTPUT/$ORGLIST/$HMMID.$HITS.faa
			fi
		echo $HITS >> $OUTPUT/$ORGLIST/$HMMID.txt;
		echo $HMMID.$HITS >> $OUTPUT/$ORGLIST/$psiRESULTS.txt;
		echo 8 >> $psiRESULTS.$HITS.categories.txt
else
		echo 0 >> $psiRESULTS.$HITS.categories.txt
#		rm $HITS.temp
fi

#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.5

				done < $DOMAIN_LISTS_DIR/$HMMID.txt

#This is where we can generate binary figures - add a header, etc UNDER CONSTRUCTION
paste -s <(echo $HMMLIST) $psiRESULTS.$HITS.wc.temp  >> $psiRESULTS.listhits

rm $psiRESULTS.$HITS.wc.temp

#NEW this is checking each criteria and making a final report ADD OTHER VALUES (1,2,4,8,16) for specific criteria
### MODIFIED FOR NEW SCORING/REPORTING BELOW
touch $psiRESULTS.$HITS.categories.txt
if [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -ge '15' ]
then
	paste -s <(echo "8") >> $psiRESULTS.ALLRESULTS.temp
	echo "$HITS		$psiRESULTS" >> $ORGLIST.PROTEIN_LIST.txt
	echo $HITS >> $psiRESULTS.LIST.txt
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '7' ]
then
	paste -s <(echo "7") >> $psiRESULTS.ALLRESULTS.temp
	echo "$HITS		$psiRESULTS" >> HMMmissing.txt
	if [ "$HMMBYPASS" = 'HB' ]
	then
		echo "$HITS     $psiRESULTS		HMM" >> $ORGLIST.PROTEIN_LIST.txt
		echo "$HITS		HMM" >> $psiRESULTS.LIST.txt
		echo "$psiRESULTS" >> $ORGLIST.HMM_INDIVIDUAL_CHECK.txt
	fi
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '13' ]
then
	paste -s <(echo "6") >> $psiRESULTS.ALLRESULTS.temp
	echo "$HITS		$psiRESULTS" >> WrongSize.txt
	if [ "$SIZEBYPASS" = 'SB' ]
	then
		echo "$HITS		$psiRESULTS		SIZE" >> $ORGLIST.PROTEIN_LIST.txt
		echo "$HITS		SIZE" >> $psiRESULTS.LIST.txt
		paste -s <(echo $psiRESULTS.$HITS) <(echo $REFLENGTH) <(echo $CANDIDATE) <(echo $COMPARE) >> BYPASSED_SIZES.txt
	fi
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '9' ]
then
	paste -s <(echo "5") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '5' ]
then
	paste -s <(echo "4") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '3' ]
then
	paste -s <(echo "3") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '1' ]
then
	paste -s <(echo "1") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '11' ]
then
	paste -s <(echo "2") >> $psiRESULTS.ALLRESULTS.temp
	echo "$HITS		$psiRESULTS" >> NotRecip.txt
	if [ "$RECIPBYPASS" = "RB" ]
	then
		echo "$HITS		$psiRESULTS		RECIP" >> $ORGLIST.PROTEIN_LIST.txt
		echo "$HITS		RECIP" >> $psiRESULTS.LIST.txt
	fi
else
paste -s <(echo "0") >> $psiRESULTS.ALLRESULTS.temp
fi

		done < $WORKING_DIR/$ORGLIST/zzSeqID.temp

#Testing ALLRESULTS.temp
#cat $psiRESULTS.ALLRESULTS.temp

#this is final parsing and making result lits per organism
		touch $OUTPUT/$ORGLIST/$psiRESULTS.txt
				if [ $(wc -l <$OUTPUT/$ORGLIST/$psiRESULTS.txt) -ge '1' ];
					then paste -s <(echo "1") >> $ORGLIST.HMMRESULTS.txt
					else paste -s <(echo "0") >> $ORGLIST.HMMRESULTS.txt
				fi
# NEW this is checking each criteria and making a final report

if [ $(grep "8" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then
	paste -s <(echo "++") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "+") >> $ORGLIST.REPORT.txt
elif [ $(grep "8" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then
	paste -s <(echo "+") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "+") >> $ORGLIST.REPORT.txt
elif [ $(grep "6" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then
	paste -s <(echo "mS") >> $ORGLIST.ALLRESULTS.txt
	if [ "$SIZEBYPASS" = 'SB' ]
    then
        paste -s <(echo "+") >> $ORGLIST.REPORT.txt
     else
        paste -s <(echo "-") >> $ORGLIST.REPORT.txt
	fi
elif [ $(grep "6" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then
	paste -s <(echo "S") >> $ORGLIST.ALLRESULTS.txt
	if [ "$SIZEBYPASS" = 'SB' ]
    then
        paste -s <(echo "+") >> $ORGLIST.REPORT.txt
     else
        paste -s <(echo "-") >> $ORGLIST.REPORT.txt
	fi
elif [ $(grep "2" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then
    paste -s <(echo "mR") >> $ORGLIST.ALLRESULTS.txt
    if [ "$RECIPBYPASS" = 'RB' ]
    then
        paste -s <(echo "+") >> $ORGLIST.REPORT.txt
    else
        paste -s <(echo "-") >> $ORGLIST.REPORT.txt
    fi
elif [ $(grep "2" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then
    paste -s <(echo "R") >> $ORGLIST.ALLRESULTS.txt
    if [ "$RECIPBYPASS" = 'RB' ]
    then
        paste -s <(echo "+") >> $ORGLIST.REPORT.txt
    else
        paste -s <(echo "-") >> $ORGLIST.REPORT.txt
    fi
elif [ $(grep "7" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then
    paste -s <(echo "mH") >> $ORGLIST.ALLRESULTS.txt
    echo "$psiRESULTS" >> HMMConservedCheck.txt
    if [ "$HMMBYPASS" = 'HB' ]
    then
        paste -s <(echo "+") >> $ORGLIST.REPORT.txt
    else
        paste -s <(echo "-") >> $ORGLIST.REPORT.txt
    fi
elif [ $(grep "7" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then
    paste -s <(echo "H") >> $ORGLIST.ALLRESULTS.txt
    echo "$psiRESULTS" >> HMMConservedCheck.txt
    if [ "$HMMBYPASS" = 'HB' ]
    then
        paste -s <(echo "+") >> $ORGLIST.REPORT.txt
     else
        paste -s <(echo "-") >> $ORGLIST.REPORT.txt
    fi
elif [ $(grep "5" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then
	paste -s <(echo "mSR") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "-") >> $ORGLIST.REPORT.txt
elif [ $(grep "5" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then
	paste -s <(echo "SR") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "-") >> $ORGLIST.REPORT.txt
elif [ $(grep "4" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then
	paste -s <(echo "mSH") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "-") >> $ORGLIST.REPORT.txt
elif [ $(grep "4" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then
	paste -s <(echo "SH") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "-") >> $ORGLIST.REPORT.txt
elif [ $(grep "3" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then
	paste -s <(echo "mRH") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "-") >> $ORGLIST.REPORT.txt
elif [ $(grep "3" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then
	paste -s <(echo "RH") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "-") >> $ORGLIST.REPORT.txt
elif [ $(grep "1" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then
	paste -s <(echo "mSRH") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "-") >> $ORGLIST.REPORT.txt
elif [ $(grep "1" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then
	paste -s <(echo "SRH") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "-") >> $ORGLIST.REPORT.txt
else
	paste -s <(echo "x") >> $ORGLIST.ALLRESULTS.txt
	paste -s <(echo "x") >> $ORGLIST.REPORT.txt
fi
		done < $PROTEIN_LIST

#this pulls the results out and pastes it to a final output of all organisms in the PipeList.txt
					paste  $OUTPUT/$ORGLIST/FINALOUTTEMP.txt $OUTPUT/$ORGLIST/$ORGLIST.REPORT.txt >> FINALOUTYOUT.woo
#clean up files
rm *.temp
mv *.categories.txt CATEGORIES/
done < $ORGANISM_LIST 

#if cat $HMMID.$HITS.faa; then paste -s <(echo $psiRESULTS) <(echo "1") >
#clean up files

#fix output report (.listhits)ls
