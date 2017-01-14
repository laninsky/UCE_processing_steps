# UCE_processing_steps
My general workflow for processing UCE data with Phyluce v1.5. The code also makes extensive use of phyluce (https://github.com/faircloth-lab/phyluce) and cloudforest (https://github.com/ngcrawford/CloudForest), written by Brant Faircloth and Nick Crawford respectively (and the excellent instructions for Phyluce at http://phyluce.readthedocs.io/en/latest/assembly.html). How I dealt with data using previous Phyluce versions is in prev_phyluce_versions.md in this current directory.

#Trimming and removing adaptor contamination
The first steps involve cleaning the reads and extracting the UCEs from the overall assembled contigs. One issue I have noticed with Illumiprocessor is that it expects samples to be named like the following {name}_L001_R1_001.fastq.gz and {name}_L001_R2_001.fastq.gz. Use the sed function to rename the files if need be. Check your cleaned reads through FastQC following the illumiprocessor step, because in same cases cutadapt may be needed to remove the adaptor sequence if trimmomatic doesn't get it all.
```
for i in *R1_001.fastq.gz; do basename=`echo $i | sed 's/R1_001.fastq.gz//g'`; cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o ${basename}adapttrimmed_R1_001.fastq.gz -p ${basename}adapttrimmed_R2_001.fastq.gz $i ${basename}R2_001.fastq.gz -q 5,15 -m 25 >> cutadapt.log; done
```

If there are issues with low percentages of trinity contigs mapping to UCEs, you can improve these numbers by removing over-represented sequences as identified through FastQC (filtering_overrepresented.R). If the FastQC comes back after this with identifiable spikes in kmers at the 5' end of each read, running filtering_kmers.R might be a good idea too. It might also be a good idea to do this before quality filtering so that these over-represented sequences are not shifted relative to the start of the read, however, this is a pretty time consuming step so I'd see how you go without it first.

In addition, if you have a lot of overlap between your F and R reads, it is a good idea to merge thenm using PEAR if you are planning on using Trinity as your assembler. ABYSS will not be able to utilize these merged reads, so you can skip this if you are going with ABYSS:
```
for i in *R1_001.fastq.gz; do /public/pear-0.9.10-bin-64/pear-0.9.10-bin-64 -f $i -r ${i/_R1_001/_R2_001} -o ${i/_adapttrimmed_R1_001.fastq.gz/} -n 33 >> pear.log 2>&1; done

#Moving the PEAR assemblies (after deleting the discarded files) into separate folders
for i in *.assembled.fastq; do basename=`echo $i | sed 's/.assembled.fastq//g'`; mkdir $basename; mv $i $basename/${basename}-READ-singleton.fastq; mv $basename.unassembled.forward.fastq $basename/${basename}-READ1.fastq; mv $basename.unassembled.reverse.fastq $basename/${basename}-READ2.fastq; done

for i in *; do gzip $i/*; done
```

After this, the first thing I do is rename all the files to lower case for the sample names, as this causes issues when converting files to phylip. I conduct the intialsteps on our local linux computers (e.g. the complabs). Things that need to be in paths: Phyluce needs to be in python path, raxml needs to be in path.

#Assembly
Running trinity assemblies using Phyluce 1.5 (http://phyluce.readthedocs.io/en/latest/assembly.html)
```
phyluce_assembly_assemblo_trinity --config trinity.conf --output trinity-assemblies/ --clean --cores 12 --log-path logs
```
Alternately, ABYSS:
```
#Navigate into the cutadaptrimmed folder
for i in *R1_001.fastq.gz; do basename=`echo $i | sed 's/_adapttrimmed_R1_001.fastq.gz//g'`; mkdir $basename; mv $i $basename/${basename}-READ1.fastq.gz; mv ${basename}_adapttrimmed_R2_001.fastq.gz $basename/${basename}-READ2.fastq.gz; done

cd ..

cp trinity.conf abyss.conf

sed -i 's/pearmergeddone/cutadapttrimmed/g' abyss.conf

#ABYSS needs unzipped data
for i in *R1_001.fastq.gz; do basename=`echo $i | sed 's/_adapttrimmed_R1_001.fastq.gz//g'`; mkdir $basename; mv $i $basename/${basename}-READ1.fastq.gz; mv ${basename}_adapttrimmed_R2_001.fastq.gz $basename/${basename}-READ2.fastq.gz; gunzip $basename/${basename}-READ1.fastq.gz; gunzip $basename/${basename}-READ2.fastq.gz; done

phyluce_assembly_assemblo_abyss --config abyss.conf --output abyss-assemblies --kmer 65 --cores 12 --clean --log-path logs
```

Before matching probes to contigs, I've found I get greater success if I filter the output of trinity to just the longest isoforms (so multiple isoforms of the same gene are not present in the output)
```
cd trinity-assemblies/contigs
mkdir longest_isoform
cd longest_isoform
for i in ../*.fasta; do newname=`echo $i | sed 's/..\///g'`; /public/trinityrnaseq-2.2.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl $i > $newname; done
```

#Matching probes to contigs
```
phyluce_assembly_match_contigs_to_probes --contigs trinity-assemblies/contigs/longest_isoform/ --probes coleoptera-v1-master-probe-list-DUPE-SCREENED.fasta --output match_contig_to_probes_longest_iso --log-path logs
```

#Generate data matrix
Set up a conf file with the names for all your taxa:
```
[dataset1]
sle0040
sle0047
sle0058
...
```
For incomplete matrix:
```
phyluce_assembly_get_match_counts --locus-db match_contig_to_probes_longest_iso/probe.matches.sqlite --taxon-list-config datasets.conf --taxon-group 'dataset1' --output dataset1.conf --log-path logs --incomplete-matrix
```
Note, if you get the following error, make sure your dataset.conf filename is correct:
```
Traceback (most recent call last):
  File "/public/uce/phyluce/bin/assembly/get_match_counts.py", line 351, in <module>
    main()
  File "/public/uce/phyluce/bin/assembly/get_match_counts.py", line 329, in main
    len(organisms),
TypeError: object of type 'NoneType' has no len()
```

#Extracting FASTA data
```
phyluce_assembly_get_fastas_from_match_counts --contigs trinity-assemblies/contigs/longest_isoform/ --locus-db match_contig_to_probes_longest_iso/probe.matches.sqlite --match-count-output dataset1.conf --incomplete-matrix dataset1.incomplete --output incomplete_fasta --log-path logs
```

#Aligning FASTA data
```
phyluce_align_seqcap_align --fasta incomplete_fasta --output incomplete_mafft_nexus --taxa 64 --aligner mafft --cores 12 --incomplete-matrix --log-path logs

phyluce_align_get_align_summary_data --alignments incomplete_mafft_nexus --cores 10 --log-path logs
```

#Locus name removal
```
phyluce_align_remove_locus_name_from_nexus_lines --alignments incomplete_mafft_nexus --output incomplete_mafft_fasta_no_locus_names --cores 12 --output-format fasta --log-path logs

phyluce_align_get_align_summary_data --alignments incomplete_mafft_fasta_no_locus_names --input-format fasta --cores 10 --log-path logs
```

#Gblocks on aligned data
(script will not work on nexus file input, even if you specify input-format nexus. You want to run this before adding missing data designators or it will attempt to implement it across the entire dataset, even across taxa the locus is not found in)
```
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed --alignments incomplete_mafft_fasta_no_locus_names --output incomplete_mafft_gblocks --output-format nexus --b2 0.5 --log-path logs --cores 10

phyluce_align_get_align_summary_data --alignments incomplete_mafft_gblocks --cores 10 --log-path logs
```

#Getting 50% complete data matrix
```
phyluce_align_get_only_loci_with_min_taxa --alignments incomplete_mafft_gblocks --taxa 64 --percent 0.5 --output 50perc_nexus --cores 12 --log-path logs
```

#Adding missing data designators 
```
phyluce_align_add_missing_data_designators --alignments 50perc_nexus --output 50perc_w_missing --match-count-output dataset1.conf --incomplete-matrix dataset1.incomplete --log-path logs --cores 10 --min-taxa 32 --verbatim --no-check-missing --log-path logs

phyluce_align_get_align_summary_data --alignments 50perc_w_missing --cores 10 --log-path logs
```

#Running RAxML on the concatenated dataset (without partitioning):
```
phyluce_align_format_nexus_files_for_raxml --alignments 50perc_w_missing --output concat_phylip

cd concat_phylip

raxmlHPC-PTHREADS-SSE3 -s 50perc_w_missing.phylip -n run1 -m GTRCATI -f a -N 100 -x $RANDOM -p $RANDOM -T 4

raxmlHPC-PTHREADS-SSE3 -s 50perc_w_missing.phylip -n run2 -m GTRCATI -f a -N 100 -x $RANDOM -p $RANDOM -T 4

```

#Making phylip alignments so we can use RAxML to estimate each of the gene trees
```
phyluce_align_convert_one_align_to_another --alignments 50perc_w_missing --output 50perc_w_missing_phylip --input-format nexus --output-format phylip --cores 4 --log-path logs

cd 50perc_w_missing_phylip

#Copy removing_missing.R from https://github.com/laninsky/Phase_hybrid_from_next_gen/blob/master/post-processing/removing_missing.R into the phylip directory
for i in `ls *.phylip`;
do echo $i > name;
Rscript removing_missing.R;
done;

cd ..
```

#Running RAxML to get the genetrees
I've had problems with phyluce's method for doing this, so this is free-hand code for doing this. After you have made your inputgenetree.tre, follow the code for "Getting the ASTRID and ASTRAL species trees based off gene trees" to generate species trees off your gene trees.
```
mkdir raxml_genetrees
wd=`pwd`
cd 50perc_w_missing_phylip
unset i
for i in `ls *.phylip`;
do mkdir ../raxml_genetrees/$i;
toraxml="raxmlHPC-SSE3 -m GTRGAMMA -n best -s $wd/50perc_w_missing_phylip/$i -p $RANDOM -w $wd/raxml_genetrees/$i";
$toraxml;
done;
 
#Summarizing the best trees 
cd ../raxml_genetrees;
touch inputgenetrees.tre;
for i in uce*; 
do cat $i/RAxML_bestTree.best >> inputgenetrees.tre
done;
```

#Running Cloudforest to get genetrees and models of substitution per locus
```
#cloudforest needs alignments to be ~9 bp longer than the number of bp in the alignment minus the number taxa in the alignment
cp -r 50perc_w_missing_phylip cloudforest_phylip
cd cloudforest_phylip
rm -rf *.reduced
rm removing_missing.R
rm name

for i in *.phylip;
do nochars=`head -n 1 $i | awk '{print $2}'`;
notaxa=`head -n 1 $i | awk '{print $1}'`;
nochars=$((nochars-notaxa))
nochars=$((nochars-9))
if [ $notaxa -gt $nochars ]
then rm -rf $i;
fi;
done;

cd ..

mkdir cloudforest_genetrees
python2 /usr/local/lib/python2.7/dist-packages/cloudforest/cloudforest_mpi.py cloudforest_phylip cloudforest_genetrees genetrees /public/PhyML-3.1/PhyML-3.1_linux64 --cores 5 --parallelism multiprocessing >> logs/cloudforest.log
```

#Partitioning loci by substitution model for concatenated RAxML runs, using cloudforest partitions
```
#cd into the genetrees folder. The following Phyluce code gets output for each model
phyluce_genetrees_split_models_from_genetrees --genetrees genetrees.tre --output output_models.txt

#R scripts for concatenating together loci with the same substitution model. Make sure you have the two cloudforestconcat Rscripts in the clouldforest_phylip directory 
cd ../cloudforest_phylip
cp ../cloudforest_genetrees/output_models.txt output_models.txt
Rscript cloudforestconcat1.R
bash concat_by_model.sh
Rscript cloudforestconcat2.R 

#Running RAxML on the partitioned dataset
raxmlHPC-PTHREADS-SSE3 -s concatphylip.phylip -q partitions.txt -n run1 -m GTRCATI -f a -N 100 -x $RANDOM -p $RANDOM -T 4
raxmlHPC-PTHREADS-SSE3 -s concatphylip.phylip -q partitions.txt -n run2 -m GTRCATI -f a -N 100 -x $RANDOM -p $RANDOM -T 4
```

#Extracting genetrees from cloudforest_genetrees output
```
#Inside the cloudforest_genetrees folder, getting just the genetrees so that the species tree gene tree methods can be run below
awk '{print $5}' genetrees.tre > inputgenetrees.tre

```

#Getting the ASTRID and ASTRAL species trees based off gene trees

Do this for RAxML and cloudforest gene trees - navigate into each of the gene tree folders, and then run the following code
```
/public/ASTRID-linux -i inputgenetrees.tre -o astrid.tre -m bionj

java -jar /public/Astral/astral.4.10.12.jar -i inputgenetrees.tre -o astral.tre

#Getting local support values for ASTRID tree through Astral
java -jar /public/Astral/astral.4.10.12.jar -q astrid.tre -i inputgenetrees.tre -o astrid_posterior.tre
```

#Running SVDquartets
https://github.com/laninsky/running_SVDquartets

#Submitting to Genbank
I've made some code that will take a tab-delimited file (called 'key') with the sample names as you have specified them in your UCE fasta files (1st column), as well as the R paste function for what you would like the genbank specifiers to be (2nd column). If you shove it in a folder full of fasta files labelled for each uce locus, it will use the name of the uce-loci in place of 'ucelocus' in the tab-delimited file to rename and print out the fasta sequences for each individual that has data (assuming missing is coded as "_").

```
>kaloula_baleata_jam3573	paste(">kaloula_baleata_jam3573","_",ucelocus," [organism=Kaloula baleata] [molecule=DNA] [mol_type=genomic DNA] [specimen_voucher=JAM3573] [note=Sampling location: Sulawesi] Kaloula baleata isolate JAM3573 ","ultra conserved element locus ",ucelocus," genomic sequence",sep="")
>kaloula_baleata_lsuhc5712	paste(">kaloula_baleata_lsuhc5712","_",ucelocus," [organism=Kaloula latidisca] [molecule=DNA] [mol_type=genomic DNA] [specimen_voucher=LSUHC5712] [note=Sampling location: Malaysia] Kaloula latidisca isolate LSUHC5712 ","ultra conserved element locus ",ucelocus," genomic sequence",sep="")
>kaloula_baleata_rmb2401	paste(">kaloula_baleata_rmb2401","_",ucelocus," [organism=Kaloula baleata] [molecule=DNA] [mol_type=genomic DNA] [specimen_voucher=TNHC67086] [note=Sampling location: Java] Kaloula baleata isolate TNHC67086 ","ultra conserved element locus ",ucelocus," genomic sequence",sep="")
```

After making sure the key file, your fasta files and the R-code are in the same folder. For the linux64.tbl2asn commands, you'll need ot have created a template file through sequin:
```
for i in *.fasta;
do echo $i > name;
Rscript genbankrename.R;
outputname=`echo $i | sed 's/.fasta//g'`;
/public/linux64.tbl2asn -i temp -t up_to_affiliation_template_7Nov2016 -o $outputname.sqn -s;
rm name;
rm temp;
done;
```

#Extracting mitogenomes - option 1: Pulling-out-mitogenomes-from-UCE-data 
Uses contig information from assembly step to generate mitogenomes: https://github.com/laninsky/Pulling-out-mitogenomes-from-UCE-data/blob/master/The_mito_pipeline.md

#Extracting mitogenomes - option 2: MITObim
Uses unassembled reads to generate mitogenomes. Install MIRA and MITObim: https://github.com/chrishah/MITObim

Interleave F and R reads (using trimmed reads from after illumiprocessor/cut-adapt stage)
```
/public/MITObim/misc_scripts/interleave-fastqgz-MITOBIM.py sle1004-READ1.fastq sle1004-READ2.fastq > sle1004_interleaved.fastq
```
Generate reference file based on mitogenome information available for same or related species on GenBank

Run MITObim:
```
/public/MITObim/MITObim_1.8.pl -end 100 -sample sle56 -ref hydrophiloidea --quick /home/a499a400/beetles/mitogenome_18_28S/sequence.fasta -readpool sle56_interleaved.fastq --pair --clean --denovo &> log
```
In the last iteration folder, there will be an assembly subfolder, and within that subfolder, a folder called *.info. Inside this folder will be a file called *_info_contigstats.txt. I pull this into my favorite spreadsheet program and sort on length. After identifying the contigs > 1,000bp, I copy them into a text file (each contig name on a new line) called "extract_contigs.txt". After putting extract_contigs.R in the same directory as "extract_contigs.txt", I set the file name for the contig file and then:
```
echo `pwd`/iteration26/sle638-hydrophiloidea_assembly/sle638-hydrophiloidea_d_results/sle638-hydrophiloidea_LargeContigs_out_sle638.unpadded.fasta > contig_file_name

Rscript extract_contigs.R
```
If you only have a single contig you can check for circularity of this sequence by:
```
/public/MITObim/misc_scripts/circules.py -f sle117_putative_mito.fasta -k 10-31
```
If no strong signal of circularity is found/we have multiple contigs, it is likely we have just recovered a partial mitogenome. I next use blast/R (make sure you have https://www.bioconductor.org/install/#update-bioconductor-packages installed - Biostrings in particular) to combine contigs that match to another contig within 50 bp of their end.
```
echo sle638-hydrophiloidea.fasta > contig_file_name 

old_line_number=`wc -l sle638-hydrophiloidea.fasta`
makeblastdb -in  sle638-hydrophiloidea.fasta -dbtype nucl

blastn -db  sle638-hydrophiloidea.fasta -query  sle638-hydrophiloidea.fasta -evalue 0.001 -outfmt "6 qacc sacc pident length qlen qstart qend slen sstart send evalue bitscore" | awk '$1!=$2 {print $0}' > blast_output.txt

Rscript blast_combine.R

new_line_number=`wc -l sle638-hydrophiloidea.fasta`
echo $old_line_number
echo $new_line_number
```

Keep looping through the above blast and R-code until your file doesn't change in length/size. Following this, I check check that the contigs seem to be mitochondrial in origin through the BLAST web-server, and that they don't overlap in my favorite assembler and by doing a command-line self-blast:
```

```
You are looking for matches between the beginning/end and end/beginning of the molecules (some may be reverse-complemented in the match). If you find them, delete the overlaps and combine the contigs (and then run them through the circularity script if you do end up with just one contig). I then run through the MITObim steps again, using these refined contigs as the starting reference. I do this to see if we can extend the contigs/obtain an entire mitogenome by just giving it just a limited number of references to work with (rather than splitting our reads among the many baits in the GenBank reference file):
```
mkdir original_run
mv iteration*/*assembly/*info/*contigreadlist.txt original_run
mv iteration*/*assembly/*info/*contigstats.txt original_run
#replace sle117 with your sample name
mv iteration*/*assembly/*results/sle117-hydrophiloidea_LargeContigs_out_sle117.unpadded.fasta original_run
mv log original_run
rm -rf iteration*

/public/MITObim/MITObim_1.8.pl -end 100 -sample sle117 -ref sle117 --quick sle117_putative_mito.fasta -readpool sle117_interleaved.fastq --pair --clean --denovo &> log
```
After this run, I follow the same steps as above (pulling out the largest contigs, checking whether they are circular/there is overlap between the beginning and the end, pulling out the next largest contig if not, confirming they are mtDNA through BLAST). 
```
head -n 718 /home/a499a400/beetles/mitogenome/sle117/iteration16/sle117-hydrophiloidea_assembly/sle117-hydrophiloidea_d_results/sle117-hydrophiloidea_LargeContigs_out_sle117.unpadded.fasta | tail -n 266 > sle117_final_mitobim.fasta

/public/MITObim/misc_scripts/circules.py -f sle117_final_mitobim.fasta -k 10-31

makeblastdb -in sle1004_final_mitobim.fasta -dbtype nucl
blastn -db sle1004_final_mitobim.fasta -query sle1004_final_mitobim.fasta -evalue 0.001 -outfmt 6 | awk '$7!=$9 {print $0}'

mkdir sle_specific_run
mv iteration*/*assembly/*info/*contigreadlist.txt sle_specific_run
mv iteration*/*assembly/*info/*contigstats.txt sle_specific_run
#replace sle117 with your sample name
mv iteration*/*assembly/*results/sle117-sle117_LargeContigs_out_sle117.unpadded.fasta sle_specific_run
mv log sle_specific_run
rm -rf iteration*
```

I then map our total paired end reads to these final MITObim contigs and generate a final consensus using bwa, samtools, gatk and picard (you'll need to have these installed). Change all instances of sle117 in the following code to your own sample name. In the  FindCoveredIntervals I've chosen to only output sequence where the coverage was >= 4 reads (contigs will be split into multiple contigs if they dip below this). You can tweak this if you like.
``` 
bwa index -a is sle117_final_mitobim.fasta 
samtools faidx sle117_final_mitobim.fasta 
/public/jdk1.8.0_112/bin/java -jar /public/picard.jar CreateSequenceDictionary R=sle117_final_mitobim.fasta O=sle117_final_mitobim.dict
bwa mem sle117_final_mitobim.fasta sle117-READ1.fastq sle117-READ2.fastq > temp.sam
/public/jdk1.8.0_112/bin/java -jar /public/picard.jar AddOrReplaceReadGroups I=temp.sam O=tempsort.sam SORT_ORDER=coordinate LB=rglib PL=illumina PU=phase SM=everyone
/public/jdk1.8.0_112/bin/java -jar /public/picard.jar MarkDuplicates MAX_FILE_HANDLES=1000 I=tempsort.sam O=tempsortmarked.sam M=temp.metrics AS=TRUE
/public/jdk1.8.0_112/bin/java -jar /public/picard.jar SamFormatConverter I=tempsortmarked.sam O=tempsortmarked.bam
samtools index tempsortmarked.bam
/public/jdk1.8.0_112/bin/java -jar /public/GenomeAnalysisTK.jar -T RealignerTargetCreator -R sle117_final_mitobim.fasta -I tempsortmarked.bam -o tempintervals.list
/public/jdk1.8.0_112/bin/java -jar /public/GenomeAnalysisTK.jar -T IndelRealigner -R sle117_final_mitobim.fasta -I tempsortmarked.bam -targetIntervals tempintervals.list -o temp_realigned_reads.bam
/public/jdk1.8.0_112/bin/java -jar /public/GenomeAnalysisTK.jar -T DepthOfCoverage -R sle117_final_mitobim.fasta -I temp_realigned_reads.bam -o temp.coverage

/public/jdk1.8.0_112/bin/java -jar /public/GenomeAnalysisTK.jar -T HaplotypeCaller -R sle117_final_mitobim.fasta -I temp_realigned_reads.bam --genotyping_mode DISCOVERY -stand_call_conf 30 -o temp_raw_variants.vcf --maxNumHaplotypesInPopulation 1
/public/jdk1.8.0_112/bin/java -jar /public/GenomeAnalysisTK.jar -T ReadBackedPhasing -R sle117_final_mitobim.fasta -I temp_realigned_reads.bam  --variant temp_raw_variants.vcf -o temp_phased_SNPs.vcf
/public/jdk1.8.0_112/bin/java -jar /public/GenomeAnalysisTK.jar -T FindCoveredIntervals -R sle117_final_mitobim.fasta -I temp_realigned_reads.bam -cov 4 -o temp_covered.list
/public/jdk1.8.0_112/bin/java -jar /public/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -V temp_phased_SNPs.vcf -R sle117_final_mitobim.fasta -L temp_covered.list -o sle117_final_mapped.fasta
mv temp.coverage coverage.txt
rm temp*
```
After verifying that the fragments are mitochondrial through BLAST, I annotate them via MITOS: http://mitos.bioinf.uni-leipzig.de/

You can then extract homologous regions across samples for downstream analyses based on the annotations.
