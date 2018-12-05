###Single cell analysis - Nagano data


```
#Cliq size of >= 5
#Converting normal bed and sorting to get uniq tads
#Calculating pixel size for each tads
#Reformatiing back so that it can be used as a dictionary

awk '$3>=5' ES_cliquestat.txt | cut -f3- -d" " | tr "," "\n" | sed 's/ //g' | tr '[:|-]' '\t' | \
sort -k1,1 -k2,2n -k3,3n -u | tr "\r" " " | \
awk '{printf "%s:%s-%s\t%.0f\n", $1, $2, $3, ($3-$2)/50000}' | \
sed 's/ //g' > tads_in_cliques_atleast5_pixel_size.txt
```

```
#Cliq size of >= 5
#Read the tads_in_cliques_atleast5_pixel_size.txt as dictionary
#Create combination of interactions

awk '$3>=5' ES_cliquestat.txt | python make_clique_bed_files.py tads_in_cliques_atleast5_pixel_size.txt
```

```
#To shuffle the cliqur bedpe files which were created above.

python shuffle_bedpe.py chr8_12_5.bed ~/utils/mm10/mm10.chrom.sizes.sorted
```

```
#Create a dict to liftover single cell data from mm9 to mm10
bedtools makewindows -g mm9.chrom.sizes -w 50000 | awk '{printf "%s\t%s\t%s\t%s|%s|%s\n", $1,$2,$3,$1,$2,$3}'> mm9_50k_windows_to_liftover.bed


```

```
#in server
#Remove header from the bedpe files
#In both haploid and diploid

sed -i '1d' *.csv

#Launch converting from mm9 to mm10
#Used liftingOver.sh script
#Content below

k=100
(
for i in `ls *.csv`
do
  ((m=m%k)); ((m++==0)) && wait
  echo $i;
  j=${i%%.csv};
  echo $j;
  python convert_mm9_to_mm10_bedpe_single_cell.py $i mm10_50k_liftedOver_from_mm9.bed > ${j}.bedpe &
done
)


```

```
#Decidied to run the pipeline for the top 10 scell files.

wc -l *.bedpe | sort -k1,1nr | head -11 | grep scell

     #2104399 scell.nextera2.NXT_2047.bedpe
     #492572 scell.nextera2.1CDES_p10_H9.bedpe
     #386853 scell.nextera2.1CDS2_287.bedpe
     #367763 scell.nextera2.1CDS1_726.bedpe
     #347441 scell.nextera2.1CDS1_673.bedpe
     #334783 scell.nextera2.1CDS2_411.bedpe
     #303754 scell.nextera2.1CDS2_675.bedpe
     #297246 scell.nextera2.1CDX4_242.bedpe
     #292991 scell.nextera2.1CDS2_461.bedpe
     #286604 scell.nextera2.1CDX4_347.bedpe

#Modify the script scell_pipeline

#In server
nohup bash scell_clique_pipeline_mod.sh scell.nextera2.NXT_2047.bedpe > scell.nextera2.NXT_2047.stat &
nohup bash scell_clique_pipeline_mod.sh scell.nextera2.1CDES_p10_H9.bedpe > scell.nextera2.1CDES_p10_H9.stat  &
nohup bash scell_clique_pipeline_mod.sh scell.nextera2.1CDS2_287.bedpe > scell.nextera2.1CDS2_287.stat  &
nohup bash scell_clique_pipeline_mod.sh scell.nextera2.1CDS1_726.bedpe > scell.nextera2.1CDS1_726.stat &
nohup bash scell_clique_pipeline_mod.sh scell.nextera2.1CDS1_673.bedpe > scell.nextera2.1CDS1_673.stat  &
nohup bash scell_clique_pipeline_mod.sh scell.nextera2.1CDS2_411.bedpe > scell.nextera2.1CDS2_411.stat  &
nohup bash scell_clique_pipeline_mod.sh scell.nextera2.1CDS2_675.bedpe > scell.nextera2.1CDS2_675.stat  &
nohup bash scell_clique_pipeline_mod.sh scell.nextera2.1CDX4_242.bedpe > scell.nextera2.1CDX4_242.stat  &
nohup bash scell_clique_pipeline_mod.sh scell.nextera2.1CDS2_461.bedpe > scell.nextera2.1CDS2_461.stat  &
nohup bash scell_clique_pipeline_mod.sh scell.nextera2.1CDX4_347.bedpe > scell.nextera2.1CDX4_347.stat  &


#Content of scell_clique_pipeline_mod.sh

iterator=`seq 1 500`
bedpe=$1
bedpe_nm=${bedpe%%.bedpe}

for j in `ls chr*.bed`
do
        cliq_nm=${j%%.bed}
        binary_count_intersect=`pairToPair -a $bedpe -b $j | wc -l`
        perm_fivehun=`for n in ${iterator[@]}; do python shuffle_bedpe.py $j mm10.chrom.sizes.sorted | sort -k1,1 -k2,2n | pairToPair -a $bedpe -b - | wc -l; done | awk '{s+=$1}END{print s/500}'`
        echo -e "$bedpe_nm \t $cliq_nm \t $binary_count_intersect \t $perm_fivehun"
done

```

