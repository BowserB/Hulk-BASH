#!/bin/bash

#this is the list of organisms
while read ORGLIST; do

mkdir ~/brad/U12proteinSearch/RESULTS/$ORGLIST
cd ~/brad/U12proteinSearch/RESULTS/$ORGLIST
mkdir ~/brad/U12proteinSearch/RESULTS/$ORGLIST/RECIP
mkdir ~/brad/U12proteinSearch/RESULTS/$ORGLIST/PROTEINS
paste ~/brad/SCRATCH/psiOf6List.txt > FINALOUTTEMP.txt

#this is performed on psiBLAST outputfmt 6 results. List file is in ~/brad/SCRATCH/psiOf6List.txt
	while read psiRESULTS; do

#NEW this evaluates if there is a psiBLAST hit or not and stores this information as variable PSI

		if [ $(ls -l ~/brad/U12proteinSearch/SCOTTpsiBLAST/COPY/$ORGLIST/$psiRESULTS | tr -s ' ' | cut -f5 -d ' ') = '0' ];
		then
			echo 0 > $psiRESULTS.categories.temp
			paste -s <(echo "0") >> $ORGLIST.ALLRESULTS.txt
			paste -s <(echo "0") >> $ORGLIST.RESULTS.txt
continue
		else echo 1 > $psiRESULTS.categories.temp
		fi
#this parses the psiBLASToutput and pulls out the seqID
		cut -f2 ~/brad/U12proteinSearch/SCOTTpsiBLAST/COPY/$ORGLIST/$psiRESULTS > zzSeqID2.temp;

#NEW this only takes the top 20 hits
		head -n20 zzSeqID2.temp > zzSeqID.temp

#this prevents the 'Search has converged' message from being in results and also sort -u removes duplicates
			if grep -q 'Search' zzSeqID.temp ; then head -n-2 zzSeqID.temp | sort -u > zzSeqID2.temp; else sort -u zzSeqID.temp > zzSeqID2.temp; fi

#NEW this cuts the hits down to 10
		head -n10 zzSeqID2.temp > zzSeqID.temp

#this pulls the list id from the title of the query
			HMMID=$(echo $psiRESULTS | cut -f3 -d. );

#this pulls the protein sequence from the proteome. HITS is from zzSeqID list.
			while read HITS; do
				sed -n '/'$HITS'/,/>/p' ~/brad/U12proteinSearch/SCOTTpsiBLAST/$ORGLIST | sed \$d > ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.temp;

#NEW this compares the size of the CANDIDIATE from the REF (human). This also stores the query for the recipBLAST
				if grep -q 'NP' <(head -n1 ~/brad/U12proteinSearch/SCOTTpsiBLAST/COPY/GCF_000001405.38_GRCh38.p12_protein.faa/$psiRESULTS);
					then REFID=$(head -n1 ~/brad/U12proteinSearch/SCOTTpsiBLAST/COPY/GCF_000001405.38_GRCh38.p12_protein.faa/$psiRESULTS | cut -f2);
					else REFID=$(head -n2 ~/brad/U12proteinSearch/SCOTTpsiBLAST/COPY/GCF_000001405.38_GRCh38.p12_protein.faa/$psiRESULTS | grep 'NP' | cut -f2);
				fi
# NEW this is continuation of above, defining the size of the REFLENGTH and the CANDIDATE from the HMMresults
#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.1

	REFLENGTH=$(grep $REFID ~/brad/U12proteinSearch/SCOTTHMM/RESULTS/GCF_000001405.38_GRCh38.p12_protein.faa/domtblout/*.txt | head -n1 | tr -s ' ' | cut -f3 -d ' ')
#NEW, removing syntax error for calculating size when there are no candidate hits in HMM output
#	if grep -q $HITS ~/brad/U12proteinSearch/SCOTTHMM/RESULTS/$ORGLIST/domtblout/*.hmm.txt;
#		then
#			CANDIDATE=$(grep $HITS ~/brad/U12proteinSearch/SCOTTHMM/RESULTS/$ORGLIST/domtblout/*.hmm.txt | head -n1 | tr -s ' ' | cut -f3 -d' ')
#		else
#			CANDIDATE=0
#		fi
# NEW this is using the actual protein sequence to get the length of the candidate.  This eliminates problem above if there are no found domians
			CANDIDATE=$(tail -n+2 ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.temp | tr -d '\n' | wc -c )
#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.1.5
#NEW Comparing the sizes
COMPARE=$(echo "(sqrt (( $REFLENGTH - $CANDIDATE )^2)/$REFLENGTH)"| bc -l)
    if (( $(echo "$COMPARE <= .25" |bc -l) ));
        then
            echo $HITS >> $psiRESULTS.RecipReady
            cp ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.temp  ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.candidate.faa
			cp $psiRESULTS.categories.temp $psiRESULTS.$HITS.categories.txt
			echo 2 >> $psiRESULTS.$HITS.categories.txt
##NEW 1.17 ADDING size report
			paste -s <(echo $psiRESULTS.$HITS) <(echo $REFLENGTH) <(echo $CANDIDATE) <(echo $COMPARE) >> $ORGLIST.SizeReport.txt
#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.2

#testing putting the recipblast for loop in this for loop
# fi
#NEW RecipBlast with Binary outs. HITID created from RecipReady list above and is a subset of the HITS
blastp -query ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.candidate.faa -db ~/brad/U12proteinSearch/SCOTTpsiBLAST/GCF_000001405.38_GRCh38.p12_protein.faa -evalue=0.00001 -outfmt 6 -out ~/brad/U12proteinSearch/RESULTS/$ORGLIST/RECIP/$HITS.BLAST.txt

	if grep $REFID <(head -n2 RECIP/$HITS.BLAST.txt);
		then
			mv ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.candidate.faa  ~/brad/U12proteinSearch/RESULTS/$ORGLIST/PROTEINS/$psiRESULTS.faa
			echo $HMMID.$HITS >> RECIPEDLIST.txt
			echo 4 >> $psiRESULTS.$HITS.categories.txt
		else
			echo 0 >> $psiRESULTS.$HITS.categories.txt
fi
#done < ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$psiRESULTS.RecipReady
##This else is reporting the size not being in the range tested above
	else
		cp $psiRESULTS.categories.temp $psiRESULTS.$HITS.categories.txt
		echo 0 >> $psiRESULTS.$HITS.categories.txt
##NEW 1.17 - for those that are not the right size, then open it up to reciping ALL.
blastp -query ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.temp -db ~/brad/U12proteinSearch/SCOTTpsiBLAST/GCF_000001405.38_GRCh38.p12_protein.faa -evalue=0.00001 -outfmt 6 -out ~/brad/U12proteinSearch/RESULTS/$ORGLIST/RECIP/$HITS.BLAST.txt

	if grep $REFID <(head -n2 RECIP/$HITS.BLAST.txt);
		then
			mv ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.candidate.faa  ~/brad/U12proteinSearch/RESULTS/$ORGLIST/PROTEINS/$psiRESULTS.faa
			echo $HMMID.$HITS >> RECIPEDLIST.txt
			echo 4 >> $psiRESULTS.$HITS.categories.txt
##NEW 1.17 ADDING size report
            paste -s <(echo $psiRESULTS.$HITS) <(echo $REFLENGTH) <(echo $CANDIDATE) <(echo $COMPARE) >> $ORGLIST.SizeReport.txt
		else
			echo 0 >> $psiRESULTS.$HITS.categories.txt
fi
fi

#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.3

# this searches the HMM outs with HMM from HMMID list, temp file appended to be pasted from later.
#NEW added desigination of HMM variable for later summary

					while read HMMLIST; do
				if grep $HITS ~/brad/U12proteinSearch/SCOTTHMM/RESULTS/$ORGLIST/domtblout/$HMMLIST.txt
					then
#						echo 8 >> $psiRESULTS.$HITS.categories.txt
						echo "1"
					else
						echo "0"
#						echo 0 >> $psiRESULTS.$HITS.categories.txt
				fi | tail -n1 >> $psiRESULTS.$HITS.wc.temp

# NEW this is continuation of above, defining the size of the REFLENGTH and the CANDIDATE from the HMMresults MOVED ABOVE BUT SLIGHTLY MODIFIED
# REFLENGTH=$(grep $REFID ~/brad/U12proteinSearch/SCOTTHMM/RESULTS/GCF_000001405.38_GRCh38.p12_protein.faa/domtblout/*.txt | head -n1 | tr -s ' ' | cut -f3 -d' ')
              #  CANDIDATE=$(grep $HITS ~/brad/U12proteinSearch/SCOTTHMM/RESULTS/$ORGLIST/domtblout/$HMMLIST.hmm.txt | tr -s ' ' | cut -f3 -d' ')

#NEW Comparing the sizes MOVED ABOVE
#COMPARE=$(echo "sqrt (( $REFLENGTH - $CANDIDATE )^2)/$REFLENGTH")| bc -l)
#	if (( $(echo "COMPARE <= .25" |bc -l) )); 
#		then
#			echo $HITS >> $psiRESULTS.RecipReady
#			cp ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.temp ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.faa
#	fi


#			hmmsearch -o ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HMMID.$HITS.$HMMLIST.txt --domtblout ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HMMID.$HITS.$HMMLIST.dmtbl.txt  ~/brad/HMM/HMMs/AllHMMs/$HMMLIST ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.temp ;
 #	          	if [ $(wc -l <~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HMMID.$HITS.$HMMLIST.dmtbl.txt) -le '13' ];then echo "0" >> $psiRESULTS.$HITS.wc.temp
#						else echo "1" >> $psiRESULTS.$HITS.wc.temp ; fi

#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.4

# this is testing if there are hits for each of the expected HMMs on the list

if [ $(wc -l <~/brad/SPLICEOSOMEproteins/DomainSets/$HMMID.txt) = $(paste -s -d+ $psiRESULTS.$HITS.wc.temp | bc) ];
then
		mv ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HITS.temp  ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HMMID.$HITS.faa;
		echo $HITS >> ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$HMMID.txt;
		echo $HMMID.$HITS >> ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$psiRESULTS.txt;
		echo 8 >> $psiRESULTS.$HITS.categories.txt
else
		echo 0 >> $psiRESULTS.$HITS.categories.txt
#		rm $HITS.temp
fi

#testing for errors in output
#echo $ORGLIST.$psiRESULTS.$HITS.5

# if cat $HMMID.$HITS.faa; then paste -s <(echo $psiRESULTS) <(echo "1") >> $psiRESULTS.ALLHITS.txt; fi
#	sort -u $psiRESULTS.ALLHITS.txt >> $psiRESULTS.ALLdeDUP.txt


	#	if grep '$HMMID.$HITS.faa' $psiRESULTS.ALL; then paste -s <(echo $psiRESULTS) <(echo "1") >> $psiRESULTS.RESULTS.txt
	#		else paste -s <(echo $psiRESULTS) <(echo "0") >> $psiRESULTS.RESULTS.txt; fi

				done < ~/brad/SPLICEOSOMEproteins/DomainSets/$HMMID.txt

#This is where we can generate binary figures - add a header, etc UNDER CONSTRUCTION
paste -s <(echo $HMMLIST) $psiRESULTS.$HITS.wc.temp  >> $psiRESULTS.listhits

rm $psiRESULTS.$HITS.wc.temp

#NEW this is checking each criteria and making a final report ADD OTHER VALUES (1,2,4,8,16) for specific criteria
### MODIFIED FOR NEW SCORING/REPORTING BELOW
touch $psiRESULTS.$HITS.categories.txt
if [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -ge '15' ]
then paste -s <(echo "8") >> $psiRESULTS.ALLRESULTS.temp
#else paste -s <(echo "0") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '7' ]
then paste -s <(echo "7") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '13' ]
then paste -s <(echo "6") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '9' ]
then paste -s <(echo "5") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '5' ]
then paste -s <(echo "4") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '3' ]
then paste -s <(echo "3") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '1' ]
then paste -s <(echo "1") >> $psiRESULTS.ALLRESULTS.temp
elif [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -eq '11' ]
then paste -s <(echo "2") >> $psiRESULTS.ALLRESULTS.temp
else
paste -s <(echo "0") >> $psiRESULTS.ALLRESULTS.temp
fi
#if [ $(paste -s -d+ $psiRESULTS.$HITS.categories.txt| bc) -lt '7' ]
#then paste -s <(echo "0") >> $psiRESULTS.ALLRESULTS.temp
#this is making sure no other combinations exist, but will add zeros
#else paste -s <(echo "0") >> $psiRESULTS.ALLRESULTS.temp
#fi
		done < ~/brad/U12proteinSearch/RESULTS/$ORGLIST/zzSeqID.temp

#Testing ALLRESULTS.temp
#cat $psiRESULTS.ALLRESULTS.temp

#this is final parsing and making result lits per organism
		touch ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$psiRESULTS.txt
				if [ $(wc -l <~/brad/U12proteinSearch/RESULTS/$ORGLIST/$psiRESULTS.txt) -ge '1' ];
					then paste -s <(echo "1") >> $ORGLIST.RESULTS.txt
					else paste -s <(echo "0") >> $ORGLIST.RESULTS.txt
				fi
# NEW this is checking each criteria and making a final report ADD OTHER VALUES (1,2,4,8,16) for specific criteria
#touch $psiRESULTS.ALLRESULTS.temp
###THIS WAS DITCHED FOR MORE SPECIFIC SCORING/REPORTING BELOW
#if [ $(paste -s -d+ $psiRESULTS.ALLRESULTS.temp| bc) -gt '1' ]
#then paste -s <(echo "2") >> $ORGLIST.ALLRESULTS.txt
#fi
#if [ $(paste -s -d+ $psiRESULTS.ALLRESULTS.temp| bc) -eq '1' ]
#then paste -s <(echo "1") >> $ORGLIST.ALLRESULTS.txt
#fi
#if [ $(paste -s -d+ $psiRESULTS.ALLRESULTS.temp| bc) -eq '0' ]
#then paste -s <(echo "0") >> $ORGLIST.ALLRESULTS.txt
#fi

if [ $(grep "8" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then paste -s <(echo "++") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "8" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then paste -s <(echo "+") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "7" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then paste -s <(echo "mH") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "7" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then paste -s <(echo "H") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "6" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then paste -s <(echo "mS") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "6" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then paste -s <(echo "S") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "5" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then paste -s <(echo "mSR") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "5" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then paste -s <(echo "SR") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "4" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then paste -s <(echo "mSH") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "4" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then paste -s <(echo "SH") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "3" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then paste -s <(echo "mRH") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "3" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then paste -s <(echo "RH") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "2" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then paste -s <(echo "mR") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "2" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then paste -s <(echo "R") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "1" $psiRESULTS.ALLRESULTS.temp | wc -l) -gt '1' ]
then paste -s <(echo "mSRH") >> $ORGLIST.ALLRESULTS.txt
elif [ $(grep "1" $psiRESULTS.ALLRESULTS.temp | wc -l) -eq '1' ]
then paste -s <(echo "SRH") >> $ORGLIST.ALLRESULTS.txt
else
	paste -s <(echo "x") >> $ORGLIST.ALLRESULTS.txt
fi
		done < ~/brad/SCRATCH/psiOf6List.txt

#this pulls the results out and pastes it to a final output of all organisms in the PipeList.txt
					paste  ~/brad/U12proteinSearch/RESULTS/$ORGLIST/FINALOUTTEMP.txt ~/brad/U12proteinSearch/RESULTS/$ORGLIST/$ORGLIST.RESULTS.txt >> FINALOUTYOUT.woo
#clean up files
rm *.temp
done < ~/brad/U12proteinSearch/HMMdonePipeList.txt

#if cat $HMMID.$HITS.faa; then paste -s <(echo $psiRESULTS) <(echo "1") >
#clean up files
#rm *.temp
#rm *.txt
 # add RECIP BLAST later
#fix output report (.listhits)ls
# format
#fix clean up
