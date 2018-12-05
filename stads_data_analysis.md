####Working on reprogramming data

```
l=(`seq 1 19` X)

for i in ${l[@]}; do chrom=chr${i}; awk -v c=$chrom '$1==c' GSM2535454_D2_rep1_raw_50kb.tsv | cut -f2,3,4 | awk '$1!=$2' | sort -k1,1n  > GSM2535454_D2/${chrom}_50kb.RAWobserved & done

```

```
#In server lamina
#In each folders

l=(`seq 1 19` X)

for i in ${l[@]};
do
arma_auto.sh chr${i}
done

```

```
#In mac

bedops --merge GSM2535450_B_tads.bed GSM2535452_Ba_tads.bed GSM2535454_D2_tads.bed GSM2535456_D4_tads.bed GSM2535458_D6_tads.bed GSM2535460_D8_tads.bed GSM2535462_ES_tads.bed > all_days_merged_tads.bed

sort GSM25354* | uniq -c | sort -n | awk '$1>=6' | awk '{print $2 "\t" $3 "\t" $4 }' | sort -k1,1 -k2,2n | awk '{if($3-$2>150000) print $0}' > tads_robust_armatus.bed

sort GSM25354* | uniq -c | sort -n | awk '$1>=3 && $1<6' | awk '{print $2 "\t" $3 "\t" $4 }' | sort -k1,1 -k2,2n | awk '{if($3-$2>150000) print $0}' | mergeBed -i - -d -1 > tads_alright_armatus.bed

bedops --difference all_days_merged_tads.bed \
tads_robust_armatus.bed | bedops --difference - \
tads_alright_armatus.bed | awk '{if($3-$2>150000) print $0}' > tads_complex_1.bed

awk '{print $0 "\trobust" }' tads_robust_armatus.bed > tads_robust_armatus.mod.bed

awk '{print $0 "\talright" }' tads_alright_armatus.bed > tads_alright_armatus.mod.bed

awk '{print $0 "\tcomplex" }' tads_complex_1.bed > tads_complex_1.mod.bed

cat tads_robust_armatus.mod.bed tads_alright_armatus.mod.bed tads_complex_1.mod.bed | sort -k1,1 -k2,2n > tads_labelled.bed

complementBed -i tads_labelled.bed -g ~/utils/mm10/mm10.chrom.sizes.sorted | awk '{print $0 "\tnot_tad"}' > mm10_fills_bw_tads_armatus.bed

cat tads_labelled.bed mm10_fills_bw_tads_armatus.bed | sort -k1,1 -k2,2n | awk '{print $0 "\t1000\t.\t" $2 "\t" $3}' | awk '{if($4=="not_tad") print $0"\t0,0,0"; else if($4=="robust") print $0"\t0,255,0"; else if($4=="alright") print $0"\t255,255,0"; else print $0"\t255,0,0"}' > genome_seg_tads_armatus.bed

cut -f1,2,3,4 genome_seg_tads_armatus.bed > genome_seg_tads_armatus.no_col.bed

awk '$4=="robust" || $4=="alright" || $4=="not_tad"' genome_seg_tads_armatus.no_col.bed > tads_robust_plus_alright_notad.bed

bedops --merge tads_alright_armatus.bed tads_robust_armatus.bed > tads_robust_plus_alright.bed

bedtools complement -i tads_robust_plus_alright.bed -g ~/utils/mm10/mm10.chrom.sizes.sorted | awk '$3-$2>150000' > complement_of_robust_alright.bed

sort -k1,1 -k2,2n GSM25354* | bedtools intersect -a - -b complement_of_robust_alright.bed | cut -f1,2 | uniq -c | awk '{if($1>2) print $2"\t"$3"\t"$3+1000"\t"$1"\tleft"}' > tads_complex_left_25.bed

sort -k1,1 -k2,2n GSM25354* | bedtools intersect -a - -b complement_of_robust_alright.bed | cut -f1,3 | uniq -c | awk '{if($1>2) print $2"\t"$3-1000"\t"$3"\t"$1"\tright"}' > tads_complex_right_25.bed

bedtools merge -i tads_complex_left_25.bed -d 150000 | awk '{print $1"\t"$2"\t"$2+1000"\tleft"}' > tads_complex_left_25_reduced.bed

sort -k1,1 -k2,2n tads_complex_right_25.bed | bedtools merge -i - -d 150000 | awk '{print $1"\t"$3-1000"\t"$3"\tright"}' > tads_complex_right_25_reduced.bed

cat tads_complex_left_25_reduced.bed tads_complex_right_25_reduced.bed tads_robust_plus_alright_notad.bed | sort -k1,1 -k2,2n > tads_to_process.bed

python ~/XX_ASC/reduce_further.py tads_to_process.bed | cut -f1,2,3,4 | sort -k1,1 -k2,2n | python ~/XX_ASC/reduce_further_further.py > genome_seg_tads_armatus.intermediate.bed

awk '$3-$2>150000' genome_seg_tads_armatus.intermediate.bed > genome_seg_tads_armatus.intermediate.no_small.bed

bedtools merge -i genome_seg_tads_armatus.intermediate.no_small.bed -c 4 -o collapse -d -1 | cut -d',' -f1 > genome_seg_tads_armatus.intermediate.no_small.no_overlap.bed

complementBed -i genome_seg_tads_armatus.intermediate.no_small.no_overlap.bed -g ~/utils/mm10/mm10.chrom.sizes.sorted | awk '{print $0 "\tnot_tad"}' > mm10_fills_bw_tads_armatus.bed

cat genome_seg_tads_armatus.intermediate.no_small.no_overlap.bed mm10_fills_bw_tads_armatus.bed | sort -k1,1 -k2,2n | awk '{print $0 "\t1000\t.\t" $2 "\t" $3}' | awk '{if($4=="not_tad") print $0"\t0,0,0"; else if($4=="robust") print $0"\t0,255,0"; else if($4=="alright") print $0"\t255,255,0"; else print $0"\t255,0,0"}' > genome_seg_tads_armatus_reprog.perf.bed

cut -f1,2,3,4 genome_seg_tads_armatus_reprog.perf.bed > reprogramming_25_5.domains

```

```
for i in ${l[@]}
do
out=chr${i}_reprogramming_25_5.domains
q=chr${i}
echo $out
awk -v chr=$q '$1==chr' reprogramming_25_5.domains > ../$out
done

k=("GSM2535450_B" "GSM2535452_Ba" "GSM2535454_D2" "GSM2535456_D4" "GSM2535458_D6" "GSM2535460_D8" "GSM2535462_ES")

l=(`seq 1 19` X)

for i in ${k[@]}
do
for j in ${l[@]}
do
out=${i}/chr${j}_${i}.50k.bedpe 
echo $out
in=${i}/chr${j}_50k.RAWobserved 
bash make_NCHG_input.sh ../chr${j}_reprogramming_25_5.domains ${in} chr${j} > ${out}
done
done

for i in ${k[@]}
do
cat ${i}/chr*.bedpe > ${i}/${i}.50kb.bedpe
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen | pairToBed -a ${i}/${i}.50kb.bedpe -b stdin -type neither > ${i}/${i}.50kb.no_cen.bedpe
~/Downloads/NCHG_hic/NCHG -m 40000 -p ${i}/${i}.50kb.no_cen.bedpe > ${i}/${i}.50kb.no_cen.NCHG.out
done

for i in ${k[@]}
do
echo $i
python ~/XX_ASC/NCHG_fdr_obs_by_exp_calc.py ${i}/${i}.50kb.no_cen.NCHG.out fdr_bh 10 0.01 > ${i}/${i}.50kb.no_cen.NCHG.sig
done

for i in ${k[@]}
do
bash make_gtrack.sh ${i}/${i}.50kb.no_cen.NCHG.sig ../tads/reprogramming_25_5.domains ${i}/${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01.gtrack
done

for i in ${k[@]}
do
bash make_gtrack_incl_lad.sh ${i}/${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01.gtrack GSE17051_cLAD_regions_mm10.bed ${i}/${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack
done

k=("GSM2535450_B" "GSM2535452_Ba" "GSM2535454_D2" "GSM2535456_D4" "GSM2535458_D6" "GSM2535460_D8" "GSM2535462_ES")

l=(`seq 1 19` X)

for i in ${k[@]}
do
for j in ${l[@]}
do
python ~/XX_ASC/getCliques.py ${i}/${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack chr${j} 3 2> ${i}/chr${j}.${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.clicstat > ${i}/chr${j}.${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.clicnum
done
done

for i in ${k[@]}
do
cat ${i}/*.${i}.*.clicnum | sort -k1,1 -k2,2n > ${i}/${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.clicnum
done

for i in ${k[@]}
do
awk '$6==1' ${i}/${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack | cut -f3 > ${i}/${i}.lads.ids
done

#Clicnums ids - for *.clicnum8 -> clic greater and equal to 8

clic_cut=(1 2 3 4 5 6 7)

for i in ${k[@]} 
do 
echo $i;
for c in ${clic_cut[@]}
do
awk -v cut=$c '{if($4==cut) print $1":"$2"-"$3}' ${i}/${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.clicnum > ${i}/${i}.clicnum${c}.ids
done
awk '{if($4>=8) print $1":"$2"-"$3}' ${i}/${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.clicnum > ${i}/${i}.clicnum8.ids
done

clic_cut=(1 2 3 4 5 6 7 8)
for i in ${k[@]}
do
echo $i
for c in ${clic_cut[@]}
do
python setop_lad.py ${i}/${i}.clicnum${c}.ids ${i}/${i}.lads.ids
done > ${i}/${i}.lads.po
done


for i in ${k[@]}
do
cat ${i}/*.${i}.*.clicnum | sort -k1,1 -k2,2n > ${i}/${i}.50kb.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.clicnum
done


for i in ${k[@]}
do
echo -n $i " "
cat ${i}/*.${i}.*.clicstat  | grep CLICSTAT | sort -k1,1 -k2,2n | awk '$3>=3' | wc -l
done
``` 

```
k=("GSM2535450_B" "GSM2535452_Ba" "GSM2535454_D2" "GSM2535456_D4" "GSM2535458_D6" "GSM2535460_D8" "GSM2535462_ES")

for i in ${k[@]}
do
 awk '$2!=$3' ${i}_rep1_raw_50kb.tsv | awk '{print $1 "\t" $2 "\t" $2+50000 "\t" $1 "\t" $3 "\t" $3+50000 "\t" $4}' > ${i}_rep1_raw_50kb.bedpe &
done

for i in ${k[@]}
do
 python bedpe_to_ijv.py ${i}_rep1_raw_50kb.bedpe mm10_50kb_ord.bed > ${i}_ijv.bed
done


#FOR GC content
~/ngs_tools/bigWigToBedGraph gc5Base.bw gc5Base.bg
awk '{ if ($1 ~ /^chr/) { print $1"\t"$2"\t"$3"\tid-"NR"\t"$4; } }' gc5Base.bg > gc5Base.bed
bedtools map -a mm10_50kb_windows.bed -b /Users/tmali/utils/mm10/gc5Base.bed -c 5 -o mean > /Users/tmali/utils/mm10/gc50kbBase.bed

#For gc content
samples=("B" "Ba" "D2" "D4" "D6" "D8" "ES")
l=(`seq 1 19` X)

for i in ${samples[@]}
do
  for j in ${l[@]}
  do
    awk -v c=$j '$1=="chr"c' /Users/tmali/utils/mm10/gc50kbBase.bed | cut -f5 | paste ${i}_chr${j}_compartments.txt - > tmp
    awk '$4<0' tmp > neg_tmp
    awk '$4>0 && $4!="NA"' tmp > pos_tmp
    neg_val=`awk '{s+=$5}END{print s/NR}' neg_tmp`
    pos_val=`awk '{s+=$5}END{print s/NR}' pos_tmp`
    
    echo ${i}_chr${j}_compartments Neg $neg_val Pos $pos_val
    if (( $(echo "$neg_val > $pos_val" | bc -l) ))
    then
      awk '{if($4<0) print $1 "\t" $2 "\t" $3 "\t" $4 "\tA"; if($4>0 && $4!="NA") print $1 "\t" $2 "\t" $3 "\t" $4 "\tB"; if($4=="NA") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $4 }' tmp > ${i}_${j}_comped.txt
    else
      awk '{if($4>0 && $4!="NA") print $1 "\t" $2 "\t" $3 "\t" $4 "\tA"; if($4<0) print $1 "\t" $2 "\t" $3 "\t" $4 "\tB"; if($4=="NA" || $4==0) print $1 "\t" $2 "\t" $3 "\t" $4 "\tNA"}' tmp > ${i}_${j}_comped.txt    
    fi
    
    rm tmp neg_tmp pos_tmp
  done 
done

samples=("B" "Ba" "D2" "D4" "D6" "D8" "ES")

for i in ${samples[@]}
do
  cat ${i}_*_comped.txt | sort -k1,1 -k2,2n > ${i}_comped.txt
done

rm *_*_comped.txt

for i in ${k[@]}
do
  s=${i#*_}
  B=`cat ${i}/*.${i}.*.clicstat  | grep CLICSTAT | sort -k1,1 -k2,2n | awk '$3>=3' | cut -f2 -d '#' | tr ',' '\n' | sort -u | sed 's/[:-]/ /g' | sort -k1,1 -k2,2n | awk '{print $1 "\t" ($2+$3)/2 "\t" (($2+$3)/2+1)}' | bedtools intersect -a - -b ${s}_comped.txt -wao | awk '$8=="B"' | wc -l`
  A=`cat ${i}/*.${i}.*.clicstat  | grep CLICSTAT | sort -k1,1 -k2,2n | awk '$3>=3' | cut -f2 -d '#' | tr ',' '\n' | sort -u | sed 's/[:-]/ /g' | sort -k1,1 -k2,2n | awk '{print $1 "\t" ($2+$3)/2 "\t" (($2+$3)/2+1)}' | bedtools intersect -a - -b ${s}_comped.txt -wao | awk '$8=="A"' | wc -l`
  total=`cat ${i}/*.${i}.*.clicstat  | grep CLICSTAT | sort -k1,1 -k2,2n | awk '$3>=3' | cut -f2 -d '#' | tr ',' '\n' | sort -u | sed 's/[:-]/ /g' | sort -k1,1 -k2,2n | awk '{print $1 "\t" ($2+$3)/2 "\t" (($2+$3)/2+1)}' | bedtools intersect -a - -b ${s}_comped.txt -wao | wc -l`
  
  B_prop=$(bc -l <<< "scale=2;$B/$total")
  A_prop=$(bc -l <<< "scale=2;$A/$total")
  
  echo $s $A_prop $B_prop  
done > compartment_status_of_tad_proportion_in_clique.txt



```

