Processing tads - All old tads + day3 new rep1 tads called using armatus 2.1 with hiclib.

```
bedops --merge day0_rep1_tads.bed day0_rep2_tads.bed \
day1_rep1_tads.bed day1_rep2_tads.bed day2_rep1_tads.bed \
day2_rep2_tads.bed day3_newrep1_tads.bed day3_rep2_tads.bed > \
all_days_merged_tads.bed 
```

Robust TADs defined by Armatus

These TADs should be present in at least 75% across all the days

```
sort *day*_*rep*_tads.bed | uniq -c | sort -n | \
awk '$1>=9' | awk '{print $2 "\t" $3 "\t" $4 }' | \
sort -k1,1 -k2,2n | awk '{if($3-$2>100000) print $0}' > tads_robust_armatus_an.bed

```

Alright TADs defined by Armatus

These TADs should be present between 50% and 75% across all the days

```
sort *day*_*rep*_tads.bed | uniq -c | sort -n | awk '$1>=6 && $1<9' | awk '{print $2 "\t" $3 "\t" $4 }' | sort -k1,1 -k2,2n | awk '{if($3-$2>100000) print $0}' | mergeBed -i - -d -1 > tads_alright_armatus_an.bed
```

Complex tads - Tentative

```
bedops --difference all_days_merged_tads.bed \
tads_robust_armatus_an.bed | bedops --difference - \
tads_alright_armatus_an.bed  | \
awk '{if($3-$2>100000) print $0}' > tads_complex_1_an.bed
```

Labelling

```
awk '{print $0 "\trobust" }' tads_robust_armatus_an.bed > tads_robust_armatus_an.mod.bed
awk '{print $0 "\talright" }' tads_alright_armatus_an.bed > tads_alright_armatus_an.mod.bed
awk '{print $0 "\tcomplex" }' tads_complex_1_an.bed > tads_complex_1_an.mod.bed
```

Merging and colouring - TENTATIVE

```
cat tads_robust_armatus_an.mod.bed tads_alright_armatus_an.mod.bed tads_complex_1_an.mod.bed | sort -k1,1 -k2,2n > tads_labelled_an.bed

complementBed -i tads_labelled_an.bed -g ~/utils/windowed_hg19/hg19.chrom.sizes | awk '{print $0 "\tnot_tad"}' > hg19_fills_bw_tads_armatus_an.bed

cat tads_labelled_an.bed hg19_fills_bw_tads_armatus_an.bed | sort -k1,1 -k2,2n | awk '{print $0 "\t1000\t.\t" $2 "\t" $3}' | awk '{if($4=="not_tad") print $0"\t0,0,0"; else if($4=="robust") print $0"\t0,255,0"; else if($4=="alright") print $0"\t255,255,0"; else print $0"\t255,0,0"}' > genome_seg_tads_armatus_asc_neuro.bed

cut -f1,2,3,4 genome_seg_tads_armatus_asc_neuro.bed > genome_seg_tads_armatus_asc_neuro.no_col.bed

awk '$4=="robust" || $4=="alright" || $4=="not_tad"' genome_seg_tads_armatus_asc_neuro.no_col.bed > tads_robust_plus_alright_notad.bed
```

Final TAD def based on Armatus

```
bedops --merge tads_alright_armatus_an.bed tads_robust_armatus_an.bed > tads_robust_plus_alright.bed

bedtools complement -i tads_robust_plus_alright.bed -g ~/utils/windowed_hg19/hg19.chrom.sizes | awk '$3-$2>100000' > complement_of_robust_alright.bed

sort -k1,1 -k2,2n *day*_*rep*_tads.bed | bedtools intersect -a - -b complement_of_robust_alright.bed | cut -f1,2 | uniq -c | awk '{if($1>3) print $2"\t"$3"\t"$3+1000"\t"$1"\tleft"}' > tads_complex_left_an_25.bed

sort -k1,1 -k2,2n *day*_*rep*_tads.bed | bedtools intersect -a - -b complement_of_robust_alright.bed | cut -f1,3 | uniq -c | awk '{if($1>3) print $2"\t"$3-1000"\t"$3"\t"$1"\tright"}' > tads_complex_right_an_25.bed

bedtools merge -i tads_complex_left_an_25.bed -d 120000 | awk '{print $1"\t"$2"\t"$2+1000"\tleft"}' > tads_complex_left_an_25_reduced.bed

sort -k1,1 -k2,2n tads_complex_right_an_25.bed | bedtools merge -i - -d 120000 | awk '{print $1"\t"$3-1000"\t"$3"\tright"}' > tads_complex_right_an_25_reduced.bed

cat tads_complex_left_an_25_reduced.bed tads_complex_right_an_25_reduced.bed tads_robust_plus_alright_notad.bed | sort -k1,1 -k2,2n > tads_to_process.bed

python ~/XX_ASC/reduce_further.py tads_to_process.bed | cut -f1,2,3,4 | sort -k1,1 -k2,2n | python ~/XX_ASC/reduce_further_further.py > genome_seg_tads_armatus_asc_neuro.intermediate.bed

awk '$3-$2>120000' genome_seg_tads_armatus_asc_neuro.intermediate.bed > genome_seg_tads_armatus_asc_neuro.intermediate.no_small.bed

#The below line deals with the overlap between not_tad and complex sometimes arises near the unmappable region in very few chromosome.. Do not panic, this shit changed nothing in the statistics
bedtools merge -i genome_seg_tads_armatus_asc_neuro.intermediate.no_small.bed -c 4 -o collapse -d -1 | cut -d',' -f1 > genome_seg_tads_armatus_asc_neuro.intermediate.no_small.no_overlap.bed

complementBed -i genome_seg_tads_armatus_asc_neuro.intermediate.no_small.no_overlap.bed -g ~/utils/windowed_hg19/hg19.chrom.sizes | awk '{print $0 "\tnot_tad"}' > hg19_fills_bw_tads_armatus_an.bed

cat genome_seg_tads_armatus_asc_neuro.intermediate.no_small.no_overlap.bed hg19_fills_bw_tads_armatus_an.bed | sort -k1,1 -k2,2n | awk '{print $0 "\t1000\t.\t" $2 "\t" $3}' | awk '{if($4=="not_tad") print $0"\t0,0,0"; else if($4=="robust") print $0"\t0,255,0"; else if($4=="alright") print $0"\t255,255,0"; else print $0"\t255,0,0"}' > genome_seg_tads_armatus_asc_neuro.perf.bed

```

```
cut -f1,2,3 genome_seg_tads_armatus_asc_neuro.perf.bed > ../ASC_neuro_final_3_5_18_hiclib.domains

l=(`seq 1 22` X)

out=chr${i}_ASC_neuro_final_3_5_18_hiclib.domains
q=chr${i} 
echo $out
awk -v chr=$q '$1==chr' ../ASC_neuro_final_3_5_18_hiclib.domains > $out
done
```

day3-newrep1

```
l=(`seq 1 22` X)

for i in ${l[@]}
do
  out=chr${i}_day3_newrep1.40k.bedpe
  echo $out
  in=day3_newrep1-40k.${i}.${i}.RAWobserved
  bash make_NCHG_input.sh chr${i}_ASC_neuro_final_3_5_18_hiclib.domains ${in} chr${i} > ${out}
done

cat *.bedpe > day3_newrep1.40k.bedpe

curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen | pairToBed -a day3_newrep1.40k.bedpe -b stdin -type neither > day3_newrep1.40k.no_cen.bedpe

~/Downloads/NCHG_hic/NCHG -m 40000 -p day3_newrep1.40k.no_cen.bedpe > day3_newrep1.40k.no_cen.NCHG.out

python ~/XX_ASC/NCHG_fdr_obs_by_exp_calc.py day3_newrep1.40k.no_cen.NCHG.out fdr_bh 10 0.01 > day3_newrep1.40k.no_cen.NCHG.sig

cut -f1,2,3,4 hiclib_analysis_with_newd3/arma_tads_g_0.3/genome_seg_tads_armatus_asc_neuro.perf.bed > hiclib_analysis_with_newd3/ASC_neuro_final_3_5_18_hiclib_w_cat.domains

bash make_gtrack.sh hiclib_analysis_with_newd3/RAWobserved/day3_newrep1.40k.no_cen.NCHG.sig hiclib_analysis_with_newd3/ASC_neuro_final_3_5_18_hiclib_w_cat.domains hiclib_analysis_with_newd3/RAWobserved/day3_newrep1.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/RAWobserved/day3_newrep1.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_d3.bed hiclib_analysis_with_newd3/RAWobserved/day3_newrep1.40k.no_cen.NCHG.obs_by_exp10.0.01_w_LADs.gtrack

l=(`seq 1 22` X)

for i in ${l[@]}
do
  python getCliques.py hiclib_analysis_with_newd3/RAWobserved/day3_newrep1.40k.no_cen.NCHG.obs_by_exp10.0.01_w_LADs.gtrack chr${i} 3 2> hiclib_analysis_with_newd3/RAWobserved/day3_newrep1.40k.no_cen.NCHG.obs_by_exp10.0.01_w_LADs.clicstat | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1000\t.\t"$2"\t"$3}' | awk '{if($4==1) print $0"\t255,255,255"; else if($4==2) print $0"\t190,190,190"; else if($4==3) print $0"\t255,198,198"; else if($4==4) print $0"\t255,169,169"; else if($4==5) print $0"\t255,141,141"; else if($4==6) print $0"\t255,113,113"; else if($4==7) print $0"\t255,84,84"; else if($4==8) print $0"\t255,56,56"; else if($4==9) print $0"\t255,28,28"; else if($4>=10) print $0"\t255,0,0"}' > hiclib_analysis_with_newd3/RAWobserved/chr${i}_day3_newrep1.40k.no_cen.NCHG.obs_by_exp10.0.01_w_LADs.clicnum.bed
done
```

```
k=("day0_rep1" "day0_rep2" "day1_rep1" "day1_rep2" "day2_rep1" "day2_rep2" "day3_rep2" "neuro_day1_rep1" "neuro_day1_rep2" "neuro_day3_rep1" "neuro_day3_rep2")
l=(`seq 1 22` X)

for i in ${k[@]}
do
for j in ${l[@]}
do
out=${i}/chr${j}_${i}.40k.bedpe 
echo $out
in=${i}/chr${j}_${i}.all.${j}.${j}.40k.RAWobserved 
bash make_NCHG_input.sh chr${j}_ASC_neuro_final_3_5_18_hiclib.domains ${in} chr${j} > ${out}
done
done 

for i in ${k[@]}
do
cat ${i}/chr*.bedpe > ${i}/${i}.40k.bedpe
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen | pairToBed -a ${i}/${i}.40k.bedpe -b stdin -type neither > ${i}/${i}.40k.no_cen.bedpe
~/Downloads/NCHG_hic/NCHG -m 40000 -p ${i}/${i}.40k.no_cen.bedpe > ${i}/${i}.40k.no_cen.NCHG.out
done

for i in ${k[@]}
do
echo $i
python ~/XX_ASC/NCHG_fdr_obs_by_exp_calc.py ${i}/${i}.40k.no_cen.NCHG.out fdr_bh 10 0.01 > ${i}/${i}.40k.no_cen.NCHG.sig
done

cd ..

for i in ${k[@]}
do
bash make_gtrack.sh hiclib_analysis_with_newd3/${i}/${i}.40k.no_cen.NCHG.sig hiclib_analysis_with_newd3/ASC_neuro_final_3_5_18_hiclib_w_cat.domains hiclib_analysis_with_newd3/${i}/${i}.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack
done

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/day0_rep1/day0_rep1.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_d0.bed hiclib_analysis_with_newd3/day0_rep1/day0_rep1.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/day0_rep2/day0_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_d0.bed hiclib_analysis_with_newd3/day0_rep2/day0_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/day1_rep1/day1_rep1.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_d1.bed hiclib_analysis_with_newd3/day1_rep1/day1_rep1.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/day1_rep2/day1_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_d1.bed hiclib_analysis_with_newd3/day1_rep2/day1_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/day2_rep1/day2_rep1.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_dm2.bed hiclib_analysis_with_newd3/day2_rep1/day2_rep1.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/day2_rep2/day2_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_dm2.bed hiclib_analysis_with_newd3/day2_rep2/day2_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/day3_rep2/day3_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_d3.bed hiclib_analysis_with_newd3/day3_rep2/day3_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/neuro_day1_rep1/neuro_day1_rep1.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_neuro_d1.bed hiclib_analysis_with_newd3/neuro_day1_rep1/neuro_day1_rep1.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/neuro_day1_rep2/neuro_day1_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/neuro_lad_d1.bed hiclib_analysis_with_newd3/neuro_day1_rep2/neuro_day1_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/neuro_day3_rep1/neuro_day3_rep1.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_neuro_d3.bed hiclib_analysis_with_newd3/neuro_day3_rep1/neuro_day3_rep1.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

bash make_gtrack_incl_lad.sh hiclib_analysis_with_newd3/neuro_day3_rep2/neuro_day3_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01.gtrack lads/lad_neuro_d3.bed hiclib_analysis_with_newd3/neuro_day3_rep2/neuro_day3_rep2.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack

for j in ${k[@]}
do
echo ${j}
for i in ${l[@]}
do
  python getCliques.py hiclib_analysis_with_newd3/${j}/${j}.40k.no_cen.NCHG.obs_by_exp10.0.01._w_LADs.gtrack chr${i} 3 2> hiclib_analysis_with_newd3/${j}/${j}.40k.no_cen.NCHG.obs_by_exp10.0.01_w_LADs.clicstat | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1000\t.\t"$2"\t"$3}' | awk '{if($4==1) print $0"\t255,255,255"; else if($4==2) print $0"\t190,190,190"; else if($4==3) print $0"\t255,198,198"; else if($4==4) print $0"\t255,169,169"; else if($4==5) print $0"\t255,141,141"; else if($4==6) print $0"\t255,113,113"; else if($4==7) print $0"\t255,84,84"; else if($4==8) print $0"\t255,56,56"; else if($4==9) print $0"\t255,28,28"; else if($4>=10) print $0"\t255,0,0"}' > hiclib_analysis_with_newd3/${j}/${j}_chr${i}.40k.no_cen.NCHG.obs_by_exp10.0.01_w_LADs.clicnum.bed
done
done


for j in ${k[@]}
do
cat hiclib_analysis_with_newd3/${j}/${j}_*.40k.no_cen.NCHG.obs_by_exp10.0.01_w_LADs.clicnum.bed > hiclib_analysis_with_newd3/${j}/${j}.40k.no_cen.NCHG.obs_by_exp10.0.01_w_LADs.clicnum.bed
done
```


