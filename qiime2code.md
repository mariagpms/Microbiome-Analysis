# Qiime2 commands
Our aim is to provide a step by step guide on the qiime2 commands that need to be executed to obtain the necessary files for the R analysis. There are several tutorials on [qiime2 docs](https://docs.qiime2.org/2024.5/) that have been of great use.<br/>
If you are not familiar with linux commands be aware that when giving the file paths '/' must be used, which is different to '\\'. <br/>
## Manifest
The .txt manifest file must be created, assigning to each sample the barcodes and the file path of the forward and reverse sequence.
## Conda
Conda needs to be activated.<br/>
```console
conda activate qiime2-amplicon-2024.2
```
## Import files
The files will be imported into qiime2.<br/>
```console
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /your_path/manifest.txt \
--output-path /your_path/demux-paired-end.qza \ 
--input-format PairedEndFastqManifestPhred33V2
```
## Demultiplexation
The sequences are demultiplexed to check the quality of the reads.<br/>
```console
qiime demux summarize \
--i-data /your_path/demux-paired-end.qza \
--o-visualization /your_path/demux-summary.qzv 
```
## Denoising
The denoising process using the [DADA2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/) algorithm is done. This means that the barcodes, primers and low quality sequences are deleted. <br/>
```console
qiime dada2 denoise-paired \
--i-demultiplexed-seqs /your_path/demux-paired-end.qza \
--p-trim-left-f 6 \
--p-trim-left-r 7 \
--p-trunc-len-f 240 \
--p-trunc-len-r 240 \
--o-representative-sequences /your_path/rep_seqs.qza \
--o-table /your_path/pet-table.qza \
--o-denoising-stats /your_path/denoising_stats.qza
```
The numeric values that are shown above indicate the positions of trimming that depend on the length of the primers and barcodes as well as the quality of the sequences.<br/>
A table of statistics can be generated from the denoising process.<br/>
```console
qiime metadata tabulate \
--m-input-file /your_path/denoising_stats.qza \
--o-visualization /your_path/denoising_stats.qzv
```
## Metadata
A .tsv file in which the metadata is associated to each sample needs to be created.<br/>
## Feature table
A table with the absolute frequences of the representative sequences will be generated.<br/>
```console
qiime feature-table summarize \
--i-table pet-table.qza \
--o-visualization pet-table.qzv \
--m-sample-metadata-file /your_path/metadata.tsv
```
The pet-table.qza will be imported so that a .tsv file can be exported containing the preprocessed otu_mat that will then be used in the R code.<br/>
```console
qiime tools export \
--input-path /your_path/pet-table.qza \
--output-path /your_path/export

biom convert \
--input-fp /your_path/feature-table.biom \
--output-fp /your_path/feature-table.csv --to-tsv
```
A table containing the relative frequences can also be obtained.<br/>
```console
qiime feature-table relative-frequency \
--i-table pet-table.qza \
--o-relative-frequency-table relative-frequency-table.qza

qiime metadata tabulate \
--m-input-file relative-frequency-table.qza \
--o-visualization relative-frequency-summary.qzv
```
## Alignment
The process of aligning the sequences is done using the [MAFFT](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC135756/) algorithm.<br/>
```console
qiime alignment mafft \
--i-sequences rep_seqs.qza \
--o-alignment aligned_rep_seqs.qza

qiime alignment mask \
--i-alignment aligned_rep_seqs.qza \
--o-masked-alignment masked_aligned_rep_seqs.qza

qiime tools export \
--input-path masked_aligned_rep_seqs.qza \
--output-path masked_aligned_rep_seqs
```
## Phylogenetic analysis
### Alpha diversity
It can be done using qiime2 but we have done it in R, for doing it in qiime2 just look for the tutorials mentioned before.
### Beta diversity
First, a tree must be generated. Then, different diversity meassures can be calculated.<br/>
```console
qiime phylogeny fasttree \
--i-alignment masked_aligned_rep_seqs.qza \
--o-tree unrooted_tree

qiime tools export \
--input-path unrooted_tree.qza \
--output-path exported_tree

qiime phylogeny midpoint-root \
--i-tree unrooted_tree.qza \
--o-rooted-tree rooted_tree

qiime tools export \
--input-path rooted_tree.qza \
--output-path exported_rooted_tree

qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted_tree.qza \
--i-table pet-table.qza \
--p-sampling-depth 10000 \
--m-metadata-file metadata.tsv \
--output-dir core_metrics_results
```
Now, once the different distances are obtained, from the above commands, the beta diversity box-plots are ploted.<br/>
```console
qiime diversity beta-group-significance \
--i-distance-matrix /your_path/core_metrics_results/distance_matrix.qza \
--m-metadata-file /your_path/metadata.tsv \
--m-metadata-column IBD_type \
--o-visualization /your_path/beta-group-significance.qzv
```
## Taxonomic analysis
The [SILVA database](https://www.arb-silva.de/) needs to be downloaded for the bacteria analysis. For the fungi analysis, the [UNITE databse](https://unite.ut.ee/repository.php) is the one requiered. <br/>
```console
wget -O "silva-138-99-seqs.qza" "https://data.qiime2.org/2024.2/common/silva-138-99-seqs.qza"

wget -O "silva-138-99-tax.qza" https://data.qiime2.org/2024.2/common/silva-138-99-tax.qza
```
Then, a model using the downloaded database needs to be trained. This should be done using a computer with a great processor and be aware that it is a long process so do not worry if it is taking a while.<br/>
```console
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads silva-138-99-seqs.qza \
--i-reference-taxonomy silva-138-99-tax.qza \
--o-classifier silva-classifier.qza
```
Finally, the taxonomic assignation is done. Basically, each otu is assigned to a taxonomic unit.<br/>
```console
qiime feature-classifier classify-sklearn \
--i-classifier /your_path/silva-classifier.qza \
--i-reads /your_path/rep_seqs.qza \
--o-classification /your_path/taxonomy.qza

qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv
```
Using [qiime2 view](https://view.qiime2.org/) the taxonomy.qzv file is imported so that it can be downloaded as a .tsv file, this will be the preprocessed tax_mat that will then be used for our R analysis.
