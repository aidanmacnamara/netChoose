while read -r gwas
do 
	GS_file=${genesetdir}/${gset}.gmt
	GWAS_file=${gwasdir}/gwas_results_${gwas}.txt
	O=$POUT/${gwas}/${gset}

	if [ -e ${O}/*sum.txt ]; then
	  	 echo "Folder $O exists, moving on...**************************************************************"
	else
		echo "Processing $O ..."
		cd $PASCAL

		mkdir -p $O
		sbatch --export=ALL ${WD}/slurm_pascal.sh $GWAS_file ${GS_file} ${O} 50000  50000  10000 sum
	fi 

        if [ -e ${O}_100up/*sum.txt ]; then
                 echo "Folder ${O}_100up exists, moving on...**************************************************************"
        else
                echo "Processing ${O}_100up ..."
                cd $PASCAL

                mkdir -p ${O}_100up
                sbatch --export=ALL ${WD}/slurm_pascal.sh $GWAS_file ${GS_file} ${O}_100up 100000  50000  10000 sum
        fi

        if [ -e ${O}_max/*sum.txt ]; then
                 echo "Folder ${O}_max exists, moving on...**************************************************************"
        else
                echo "Processing ${O}_max ..."
                cd $PASCAL

		mkdir -p ${O}_max
		sbatch --export=ALL ${WD}/slurm_pascal.sh $GWAS_file ${GS_file} ${O}_max 50000  50000  10000 max
	fi

	cd $POUT

done < gwas_65traits.txt
