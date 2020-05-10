# splits a multi-fasta file into many constituent files

while read line

do

    if [[ ${line:0:1} == '>' ]]

    then

        outfile=${line#>}.fa

        echo $line > $outfile

    else

        echo $line >> $outfile

    fi

done < 15664.fasta
