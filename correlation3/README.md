# Data Prepare(Sorghum bicolor gff is different from corralation directory/sorghum.gff )
```
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Sorghum_bicolor/latest_assembly_versions/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Sorghum_bicolor/latest_assembly_versions/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.gz
gunzip *gz
```

# Delete error annotation(strand position is null)
```
python delete_error_annoation.py -g GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff -o GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.new.gff

```
# Get AnchorWave collinearity file
# AnchorWave WGDI MCScanX using the same blast file(MCScanX remove gene_name prefix "gene-" due to MCScanX character) and jcvi using default lastz pairwise alignment  
# AnchorWave pro alignment_length=0 bitscore_minimum=100

```
gffread -g GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna -y sb.p.fa GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.new.gff -S
gffread -g Zm-B73-REFERENCE-NAM-5.0.fa -y zea.p.fa Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -S
python3 longestPeps.py -g Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -f Zm-B73-REFERENCE-NAM-5.0.fa -p zea.p.fa -o maize.protein.fa
python3 longestPeps.py -g GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.new.gff -f GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna -p sb.p.fa -o sorghum.protein.fa
diamond makedb --in maize.protein.fa --db maize
diamond blastp --db maize -q sorghum.protein.fa -k 5 -e 1e-10 -o sorghum.maize.blastp
python3 combineBlastAndStrandInformation.py -r Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -q GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.new.gff -b sorghum.maize.blastp -o sorghum.maize.table
anchorwave pro -i sorghum.maize.table -o sorghum.maize.collinearity -R 1 -Q 2
```

# Get WGDI collinearity 
```
python 01.py Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 zm.old.gff
python 01.py GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.new.gff sb.old.gff
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
python MCScanX_input.py GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.new.gff Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 zm_sb.gff
mv sorghum.maize.blastp zm_sb.blast
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
python -m jcvi.formats.gff bed Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o zm.bed
python -m jcvi.formats.gff bed GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.new.gff -o sb.bed
python -m jcvi.compara.catalog ortholog sb zm --notex
python -m jcvi.compara.quota sb.zm.anchors --qbed sb.bed --sbed zm.bed --quota 1:2 --screen
python -m jcvi.compara.quota sb.zm.lifted.anchors --qbed sb.bed --sbed zm.bed --quota 1:2 --screen
python -m jcvi.graphics.dotplot --notex sb.zm.lifted.1x2.anchors
python -m jcvi.graphics.dotplot --notex sb.zm.1x2.anchors
```
# collinearity gene pair correlation
```
python main.py row_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -t tassel_cluster_ave.csv -e ear_cluster_ave.csv -s sorghum_cluster_ave.csv -m spearman -b T -co result/gene_pair/raw/spearman_log.csv
python main.py row_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -t tassel_cluster_ave.csv -e ear_cluster_ave.csv -s sorghum_cluster_ave.csv -m spearman -b F -co result/gene_pair/raw/spearman_tpm.csv
```
```
python main.py row_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -t tassel_cluster_ave.csv -e  ear_cluster_ave.csv -s sorghum_cluster_ave.csv -m pearson -b T -co result/gene_pair/raw/pearson_log.csv
python main.py row_corr -A sorghum.maize.collinearity -W zm.sb.collinearity -M zm_sb.collinearity -J sb.zm.1x2.anchors -t tassel_cluster_ave.csv -e ear_cluster_ave.csv -s sorghum_cluster_ave.csv -m pearson -b F -co result/gene_pair/raw/pearson_tpm.csv
```
# transform csv format file to R Input boxplot_violin
```
python main.py process_boxplot_violin -i result/gene_pair/raw/spearman_log.csv -o1 result/gene_pair/log_spearman/o1.csv -o2 result/gene_pair/log_spearman/o2.csv -o3 result/gene_pair/log_spearman/o3.csv
python main.py process_boxplot_violin -i result/gene_pair/raw/spearman_tpm.csv -o1 result/gene_pair/tpm_spearman/o1.csv -o2 result/gene_pair/tpm_spearman/o2.csv -o3 result/gene_pair/tpm_spearman/o3.csv
python main.py process_boxplot_violin -i result/gene_pair/raw/pearson_log.csv -o1 result/gene_pair/log_pearson/o1.csv -o2 result/gene_pair/log_pearson/o2.csv -o3 result/gene_pair/log_pearson/o3.csv
python main.py process_boxplot_violin -i result/gene_pair/raw/pearson_tpm.csv -o1 result/gene_pair/tpm_pearson/o1.csv -o2 result/gene_pair/tpm_pearson/o2.csv -o3 result/gene_pair/tpm_pearson/o3.csv
```
# transform csv format file to R Input venn
```
bash four_multiple_six.sh
```
# transform csv format file to R Input five column
```
bash three_multiple_four.sh
```
