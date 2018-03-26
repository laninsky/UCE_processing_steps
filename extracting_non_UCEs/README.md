## Extracting mitogenomes - option 1: Pulling-out-mitogenomes-from-UCE-data 
Uses contig information from assembly step to generate mitogenomes: https://github.com/laninsky/Pulling-out-mitogenomes-from-UCE-data/blob/master/The_mito_pipeline.md

## Extracting mitogenomes, 18S/28S (or any other constant fragment) - option 2: MITObim
Uses unassembled reads to generate mitogenomes. Install MIRA and MITObim: https://github.com/chrishah/MITObim

Interleave F and R reads (using trimmed reads from after illumiprocessor/cut-adapt stage)
```
/public/MITObim/misc_scripts/interleave-fastqgz-MITOBIM.py sle1004-READ1.fastq sle1004-READ2.fastq > sle1004_interleaved.fastq
```
Generate reference file based on information available for same or related species on GenBank. You can include multiple fragments (e.g. mitogenomes and 18S/28S) in this reference file.

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
i=sle638-hydrophiloidea.fasta
echo $i > contig_file_name 

old_line_number=`wc -l $i`
makeblastdb -in $i -dbtype nucl

blastn -db $i -query $i -evalue 0.001 -outfmt "6 qacc sacc pident length qlen qstart qend slen sstart send evalue bitscore" | awk '$1!=$2 {print $0}' > blast_output.txt

Rscript blast_combine.R

new_line_number=`wc -l $i`
echo $old_line_number
echo $new_line_number
```

Keep looping through the above blast and R-code until your file doesn't change in length/size. Following this, I check check that the contigs seem to be mitochondrial (or your target region) in origin through the BLAST web-server. I then run through the MITObim steps again, using these refined contigs as the starting reference. I do this to see if we can extend the contigs/obtain an entire mitogenome by just giving it just a limited number of references to work with (rather than splitting our reads among the many baits in the GenBank reference file):
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
After this run, I follow the same steps as above (pulling out the largest contigs, checking whether they are circular/there is overlap between the beginning and the end, pulling out the next largest contig if not, confirming they are mtDNA/target fragments through BLAST). 
```
/public/MITObim/misc_scripts/circules.py -f sle117_final_mitobim.fasta -k 10-31

makeblastdb -in sle1004_final_mitobim.fasta -dbtype nucl
blastn -db sle1004_final_mitobim.fasta -query sle1004_final_mitobim.fasta -evalue 0.001 -outfmt 6 | awk '$1!=$2 {print $0}'

mkdir sle_specific_run
mv iteration*/*assembly/*info/*contigreadlist.txt sle_specific_run
mv iteration*/*assembly/*info/*contigstats.txt sle_specific_run
#replace sle117 with your sample name
mv iteration*/*assembly/*results/sle117-sle117_LargeContigs_out_sle117.unpadded.fasta sle_specific_run
mv log sle_specific_run
rm -rf iteration*
```

I then map our total paired end reads to these final MITObim contigs and generate a final consensus using bwa, samtools, gatk and picard (you'll need to have these installed). I do this separately for each locus I am trying to fish out, but if they had the same ploidy (i.e. all were nuclear) you could do them at the same time. Change all instances of sle117 in the following code to your own sample name. In the  FindCoveredIntervals I've chosen to only output sequence where the coverage was >= 4 reads (contigs will be split into multiple contigs if they dip below this). You can tweak this if you like. In addition, if you are fishing out nuclear regions (e.g. 18S etc), you probably want to change the --maxNumHaplotypesInPopulation flag to 2.

``` 
bwa index -a is sle117_final_mitobim.fasta 
samtools faidx sle117_final_mitobim.fasta 
/public/jdk1.8.0_112/bin/java -jar /public/picard.jar CreateSequenceDictionary R=sle117_final_mitobim.fasta O=sle117_final_mitobim.dict
bwa mem sle117_final_mitobim.fasta sle117-READ1.fastq sle117-READ2.fastq > temp.sam
/public/jdk1.8.0_112/bin/java -jar /public/picard.jar AddOrReplaceReadGroups I=temp.sam O=tempsort.sam SORT_ORDER=coordinate LB=rglib PL=illumina PU=phase SM=everyone
/public/jdk1.8.0_112/bin/java -jar /public/picard.jar MarkDuplicates MAX_FILE_HANDLES=1000 I=tempsort.sam O=tempsortmarked.sam M=temp.metrics AS=TRUE
/public/jdk1.8.0_112/bin/java -jar /public/picard.jar SamFormatConverter I=tempsortmarked.sam O=tempsortmarked.bam
samtools index tempsortmarked.bam
gatk DepthOfCoverage -R sle117_final_mitobim.fasta -I temp_realigned_reads.bam -o temp.coverage

gatk HaplotypeCaller -R sle117_final_mitobim.fasta -I temp_realigned_reads.bam -stand-call-conf 30 -o temp_raw_variants.vcf --maxNumHaplotypesInPopulation 1
gatk ReadBackedPhasing -R sle117_final_mitobim.fasta -I temp_realigned_reads.bam  --variant temp_raw_variants.vcf -o temp_phased_SNPs.vcf
gatk FindCoveredIntervals -R sle117_final_mitobim.fasta -I temp_realigned_reads.bam -cov 4 -o temp_covered.list
gatk FastaAlternateReferenceMaker -V temp_phased_SNPs.vcf -R sle117_final_mitobim.fasta -L temp_covered.list -o sle117_final_mapped.fasta
mv temp.coverage coverage.txt
rm temp*
```
This will do the same thing except using pre-version 4 of gatk
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

After verifying that the fragments were on target through BLAST, and replacing any spaces in the fasta header with undescores:
-- I annotated mitogenomes via MITOS: http://mitos.bioinf.uni-leipzig.de/
-- and 18S-28S via ITSX:
```
for i in *.fasta; basename=`echo $i | sed 's/.fasta//g'`; do ITSx --allow_single_domain --preserve --save_regions  all -i $i -o $basename; done 
```
\\\sed -ir 's|\(>[0-9]\)\(_sle\)\(.*\)|\2|g' test.fas
