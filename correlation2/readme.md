# Data prepare
```
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/gff3/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.58.gff3.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/gff3/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.58.gff3.gz
gunzip *gz
```

# AnchorWave WGDI MCScanX using the same blast file(MCScanX remove gene_name prefix "gene:" due to MCScanX character) and jcvi using default lastz pairwise alignment  
# AnchorWave pro alignment_length=0 bitscore_minimum=100
```
gffread -g Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -y sb.p.fa Sorghum_bicolor.Sorghum_bicolor_NCBIv3.58.gff3 -S
gffread -g Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa -y zea.p.fa Zea_mays.Zm-B73-REFERENCE-NAM-5.0.58.gff3 -S
python3 longestPeps.py -g Zea_mays.Zm-B73-REFERENCE-NAM-5.0.58.gff3 -f Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa -p zea.p.fa -o maize.protein.fa
python3 longestPeps.py -g Sorghum_bicolor.Sorghum_bicolor_NCBIv3.58.gff3 -f Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -p sb.p.fa -o sorghum.protein.fa
diamond makedb --in maize.protein.fa --db maize
diamond blastp --db maize -q sorghum.protein.fa -k 5 -e 1e-10 -o sorghum.maize.blastp
python3 combineBlastAndStrandInformation.py -r Zea_mays.Zm-B73-REFERENCE-NAM-5.0.58.gff3 -q Sorghum_bicolor.Sorghum_bicolor_NCBIv3.58.gff3 -b sorghum.maize.blastp -o sorghum.maize.table
anchorwave pro -i sorghum.maize.table -o sorghum.maize.collinearity -R 1 -Q 2
```
# Get WGDI collinearity file
```
python 01.py Zea_mays.Zm-B73-REFERENCE-NAM-5.0.58.gff3 zm.old.gff
python 01.py Sorghum_bicolor.Sorghum_bicolor_NCBIv3.58.gff3 sb.old.gff
python 02.py sb.old.gff sb.gff sb.lens
python 02.py zm.old.gff zm.gff zm.lens
cp sorghum.maize.blastp WGDI/
mv zm.old.gff sb.old.gff sb.gff zm.gff sb.lens zm.lens WGDI/
conda activate wgdi
cd WGDI
wgdi -icl icl.conf
```

# get MCScanX input file, removing blast file gene_name gene: prefix
```
cd ../
python MCScanX_input.py Sorghum_bicolor.Sorghum_bicolor_NCBIv3.58.gff3 Zea_mays.Zm-B73-REFERENCE-NAM-5.0.58.gff3 zm_sb.gff
python remove_genecolon.py sorghum.maize.blastp zm_sb.blast
mv zm_sb.gff zm_sb.blast MCScanX/zm_sb
cd MCScanX/
MCScanX zm_sb/zm_sb
```
# get longest cds, you should modify cds file name 
```
conda deactivate
cd ../jcvi
python main.py pre_col
python main.py pre_ks
```
# get feature:gene bed format file by jcvi
# jcvi maybe by mcscan(python version) get anchors file
# quota 1:2 this is differnet from anchorwave, but essentially the same and get quota 1:2 anchors file
# get pdf figure
```
python -m jcvi.formats.gff bed Zea_mays.Zm-B73-REFERENCE-NAM-5.0.58.gff3 -o zm.bed
python -m jcvi.formats.gff bed Sorghum_bicolor.Sorghum_bicolor_NCBIv3.58.gff3 -o sb.bed
python -m jcvi.compara.catalog ortholog sb zm --notex
python -m jcvi.compara.quota sb.zm.anchors --qbed sb.bed --sbed zm.bed --quota 1:2 --screen
python -m jcvi.compara.quota sb.zm.lifted.anchors --qbed sb.bed --sbed zm.bed --quota 1:2 --screen
python -m jcvi.graphics.dotplot --notex sb.zm.lifted.1x2.anchors
python -m jcvi.graphics.dotplot --notex sb.zm.1x2.anchors
```
# you should make some folder
# result/four_software/all
# Please move four collinearity file to correlation2 folder then doing correlation analysis in order to compute correlation in different condition
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################            collinearity gene pair correlation   ###############################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
```
time python main.py row_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -e all_special_tpm.xlsx -M zm_sb.collinearity -J sb.zm.1x2.anchors -co result/four_software/gene_pair/raw/correlation_pearson_tpm.csv -m pearson -b F
time python main.py process -i result/four_software/gene_pair/raw/correlation_pearson_tpm.csv -o result/four_software/gene_pair/correlation_pearson_tpm.R.csv
```

```
time python main.py row_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -e all_special_tpm.xlsx -M zm_sb.collinearity -J sb.zm.1x2.anchors -co result/four_software/gene_pair/raw/correlation_spearman_tpm.csv -m spearman -b F
time python main.py process -i result/four_software/gene_pair/raw/correlation_spearman_tpm.csv -o result/four_software/gene_pair/correlation_spearman_tpm_R.csv
```

```
time python main.py row_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -e all_special_tpm.xlsx -M zm_sb.collinearity -J sb.zm.1x2.anchors -co result/four_software/gene_pair/raw/correlation_pearson_log.csv -m pearson -b T
time python main.py process -i result/four_software/gene_pair/raw/correlation_pearson_log.csv -o result/four_software/gene_pair/correlation_pearson_log_R.csv
```

```
time python main.py row_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -e all_special_tpm.xlsx -M zm_sb.collinearity -J sb.zm.1x2.anchors -co result/four_software/gene_pair/raw/correlation_spearman_log.csv -m spearman -b T
time python main.py process -i result/four_software/gene_pair/raw/correlation_spearman_log.csv -o result/four_software/gene_pair/correlation_spearman_log_R.csv
```
# block don't drop_puplicates()  I think block_30 maybe useful
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################          all block compute  correlation   #####################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/all_block/raw/block_R1-d1_pearson_correlation_log.csv -m pearson -c R1-d1 -b T -B 2
time python main.py process -i result/four_software/all_block/raw/block_R1-d1_pearson_correlation_log.csv -o result/four_software/all_block/block_R1-d1_pearson_correlation_log.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/all_block/raw/block_R1-d1_spearman_correlation_log.csv -m spearman -c R1-d1 -b T -B 2
time python main.py process -i result/four_software/all_block/raw/block_R1-d1_spearman_correlation_log.csv -o result/four_software/all_block/block_R1-d1_spearman_correlation_log.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/all_block/raw/block_R1-d1_pearson_correlation_tpm.csv -m pearson -c R1-d1 -b F -B 2
time python main.py process -i result/four_software/all_block/raw/block_R1-d1_pearson_correlation_tpm.csv -o result/four_software/all_block/block_R1-d1_pearson_correlation_tpm.R.csv
```

```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/all_block/raw/block_R1-d1_spearman_correlation_tpm.csv -m spearman -c R1-d1 -b F -B 2
time python main.py process -i result/four_software/all_block/raw/block_R1-d1_spearman_correlation_tpm.csv -o result/four_software/all_block/block_R1-d1_spearman_correlation_tpm.R.csv
```
########################################################################################################################################################
########################################################################################################################################################
########################################################          block length minumim 10              #################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_10/raw/block_R1-d1_pearson_correlation_log.csv -m pearson -c R1-d1 -B 10
time python main.py process -i result/four_software/block_10/raw/block_R1-d1_pearson_correlation_log.csv -o result/four_software/block_10/block_R1-d1_pearson_correlation_log.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_10/raw/block_R1-d1_spearman_correlation_log.csv -m spearman -c R1-d1 -B 10
time python main.py process -i result/four_software/block_10/raw/block_R1-d1_spearman_correlation_log.csv -o result/four_software/block_10/block_R1-d1_spearman_correlation_log.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_10/raw/block_R1-d1_pearson_correlation_tpm.csv -m pearson -c R1-d1 -B 10 -b F
time python main.py process -i result/four_software/block_10/raw/block_R1-d1_pearson_correlation_tpm.csv -o result/four_software/block_10/block_R1-d1_pearson_correlation_tpm.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_10/raw/block_R1-d1_spearman_correlation_tpm.csv -m spearman -c R1-d1 -B 10 -b F
time python main.py process -i result/four_software/block_10/raw/block_R1-d1_spearman_correlation_tpm.csv -o result/four_software/block_10/block_R1-d1_spearman_correlation_tpm.R.csv
```
########################################################################################################################################################
########################################################################################################################################################
########################################################          block minimum length 30       ########################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30/raw/block_R1-d1_pearson_correlation_log.csv -m pearson -c R1-d1 -B 30
time python main.py process -i result/four_software/block_30/raw/block_R1-d1_pearson_correlation_log.csv -o result/four_software/block_30/block_R1-d1_pearson_correlation_log.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30/raw/block_R1-d1_spearman_correlation_log.csv -m spearman -c R1-d1 -B 30
time python main.py process -i result/four_software/block_30/raw/block_R1-d1_spearman_correlation_log.csv -o result/four_software/block_30/block_R1-d1_spearman_correlation_log.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30/raw/block_R1-d1_pearson_correlation_tpm.csv -m pearson -c R1-d1 -B 30 -b F
time python main.py process -i result/four_software/block_30/raw/block_R1-d1_pearson_correlation_tpm.csv -o result/four_software/block_30/block_R1-d1_pearson_correlation_tpm.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30/raw/block_R1-d1_spearman_correlation_tpm.csv -m spearman -c R1-d1 -B 30 -b F
time python main.py process -i result/four_software/block_30/raw/block_R1-d1_spearman_correlation_tpm.csv -o result/four_software/block_30/block_R1-d1_spearman_correlation_tpm.R.csv
```
########################################################################################################################################################
########################################################################################################################################################
########################################################         random query_to_ref block minimum length 30       #####################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30_random/raw/block_R1-d1_pearson_correlation_log.csv -m pearson -c R1-d1 -B 30 -r T
time python main.py process -i result/four_software/block_30_random/raw/block_R1-d1_pearson_correlation_log.csv -o result/four_software/block_30_random/block_R1-d1_pearson_correlation_log.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30_random/raw/block_R1-d1_spearman_correlation_log.csv -m spearman -c R1-d1 -B 30 -r T
time python main.py process -i result/four_software/block_30_random/raw/block_R1-d1_spearman_correlation_log.csv -o result/four_software/block_30_random/block_R1-d1_spearman_correlation_log.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30_random/raw/block_R1-d1_pearson_correlation_tpm.csv -m pearson -c R1-d1 -B 30 -b F -r T
time python main.py process -i result/four_software/block_30_random/raw/block_R1-d1_pearson_correlation_tpm.csv -o result/four_software/block_30_random/block_R1-d1_pearson_correlation_tpm.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30_random/raw/block_R1-d1_spearman_correlation_tpm.csv -m spearman -c R1-d1 -B 30 -b F -r T
time python main.py process -i result/four_software/block_30_random/raw/block_R1-d1_spearman_correlation_tpm.csv -o result/four_software/block_30_random/block_R1-d1_spearman_correlation_tpm.R.csv
```

########################################################################################################################################################
########################################################################################################################################################
########################################################         two block block minimum length 30       #####################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30_two/raw/block_R1-d1_pearson_correlation_log.csv -m pearson -c R1-d1 -B 30 -t T
time python main.py process -i result/four_software/block_30_two/raw/block_R1-d1_pearson_correlation_log.csv -o result/four_software/block_30_two/block_R1-d1_pearson_correlation_log.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30_two/raw/block_R1-d1_spearman_correlation_log.csv -m spearman -c R1-d1 -B 30 -t T
time python main.py process -i result/four_software/block_30_two/raw/block_R1-d1_spearman_correlation_log.csv -o result/four_software/block_30_two/block_R1-d1_spearman_correlation_log.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30_two/raw/block_R1-d1_pearson_correlation_tpm.csv -m pearson -c R1-d1 -B 30 -b F -t T
time python main.py process -i result/four_software/block_30_two/raw/block_R1-d1_pearson_correlation_tpm.csv -o result/four_software/block_30_two/block_R1-d1_pearson_correlation_tpm.R.csv
```
```
python main.py block_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -e all_special_tpm.xlsx -co result/four_software/block_30_two/raw/block_R1-d1_spearman_correlation_tpm.csv -m spearman -c R1-d1 -B 30 -b F -t T
time python main.py process -i result/four_software/block_30_two/raw/block_R1-d1_spearman_correlation_tpm.csv -o result/four_software/block_30_two/block_R1-d1_spearman_correlation_tpm.R.csv
```
###
```
bash merge.sh
```
###
```
bash merge1.sh
```





-
