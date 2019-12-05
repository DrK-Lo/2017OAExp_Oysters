DIR="/shared_lab/20180226_RNAseq_2017OAExp/DNAm" # set your path
cd $DIR

INPUT="processed_samples/05_countSummary"
RSCRIPT="scripts/R/"
OUTPUT="processed_samples/07_brmsSummary"

# Code not error proof - the difference between START and END should be greater than then number of cores
START=1
END=100
NCORES=10
range=$(($END-$START+1)) # The range is inclusive (START=1 and END=10 includes both loci 1 and 10)
echo $range
lociPerCore=$(($range/$NCORES))
echo $lociPerCore
start=$(seq $START $lociPerCore $END)
echo $start

for i in $start
do
  	echo -e "\n\n"
        end=$((i+$lociPerCore-1))

        if [ $end -gt $END ]
        then
            	end=$END
        fi

    ### Call R script
    echo "Running R scripts for lines" echo 'start row' $i 'end row' $end
    echo ${RSCRIPT}07_diffMethylation_brms.R
    
    #Rscript --vanilla ${RSCRIPT}07_diffMethylation_brms.R ${i} $end > ${i}"_R.out" 2> ${i}"_R.error" & echo $!

    # the above line, using the $! at the end, will run the script in the background on a single core
    # in the next iteration, it will run the script in the background on another core for a different set of rows
    # do not wait until the last background process is finished

done


