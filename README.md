# UCE_processing_steps
The general workflow for processing UCE data that I use. One big acknowledgment: almost all the code I use comes from Carl Oliveros (https://github.com/carloliveros/uce-scripts). The code also makes extensive use of phyluce (https://github.com/faircloth-lab/phyluce) and cloudforest (https://github.com/ngcrawford/CloudForest), written by Brant Faircloth and Nick Crawford respectively. Any steps referred to are the steps from Carl's pipeline described at his github (https://github.com/carloliveros/uce-scripts/blob/master/UCE%20pipeline.md)

The first steps involve cleaning the reads and extracting the UCEs from the overall assembled contigs. One issue I have noticed with Illumiprocessor is that it expects samples to be named like the following {name}_L001_R1_001.fastq.gz and {name}_L001_R2_001.fastq.gz. Use the sed function to rename the files if need be. Check your cleaned reads through FastQC following the illumiprocessor step, because in same cases cutadapt may be needed to remove the adaptor sequence if trimmomatic doesn't get it all.
```
for i in *R1_001.fastq.gz; do basename=`echo $i | sed 's/R1_001.fastq.gz//g'`; cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o ${basename}adapttrimmed_R1_001.fastq.gz -p ${basename}adapttrimmed_R2_001.fastq.gz $i ${basename}R2_001.fastq.gz; done
```

If there are issues with low percentages of trinity contigs mapping to UCEs, it may also be a good idea to remove over-represented sequences as identified through FastQC (filtering_overrepresented.R). If the FastQC comes back after this with identifiable spikes in kmers at the 5' end of each read, running filtering_kmers.R might be a good idea too. I do this before quality filtering so that these over-represented sequences are not shifted relative to the start of the read.

After this, the first thing I do is rename all the files to lower case for the sample names, as this causes issues when converting files to phylip. I conduct the intialsteps on our local linux computers (e.g. the complabs). Things that need to be in paths: Phyluce needs to be in python path, raxml needs to be in path.

#STEP 3A
get the probes for matching the contigs to
```
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta
```
Making a log directory
```
mkdir logs
```
#STEP 3B
matching the trinity contigs to the probes
```
python /public/uce/phyluce/bin/assembly/match_contigs_to_probes.py --contigs /home/a499a400/gekko/trinity-assemblies/contigs/ --probes uce-5k-probes.fasta --output lastz --log-path logs
```
#STEP 4A
used Scott and Carl's previous dataset.conf file. Dataset name is 'All'. 
```
python /public/uce/phyluce/bin/assembly/get_match_counts.py --locus-db /home/a499a400/gekko/lastz/probe.matches.sqlite --taxon-list-config dataset.conf --taxon-group 'gekko' --output /home/a499a400/gekko/complete.conf --log-path logs

python /public/uce/phyluce/bin/assembly/get_match_counts.py --locus-db /home/a499a400/gekko/lastz/probe.matches.sqlite --taxon-list-config dataset.conf --taxon-group 'gekko' --output /home/a499a400/gekko/incomplete.conf --log-path logs --incomplete-matrix
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
#STEP 4B
```
python /public/uce/phyluce/bin/assembly/get_fastas_from_match_counts.py --contigs /home/a499a400/gekko/trinity-assemblies/contigs/ --locus-db /home/a499a400/gekko/lastz/probe.matches.sqlite --match-count-output /home/a499a400/gekko/complete.conf --output complete_fasta --log-path logs

python /public/uce/phyluce/bin/assembly/get_fastas_from_match_counts.py --contigs /home/a499a400/gekko/trinity-assemblies/contigs/ --locus-db /home/a499a400/gekko/lastz/probe.matches.sqlite --match-count-output /home/a499a400/gekko/incomplete.conf --output incomplete_fasta --incomplete-matrix incomplete.nostrict --log-path logs
```
#STEP 7A
```
python /public/uce/phyluce/bin/align/seqcap_align_2.py --fasta complete_fasta --output complete_aligned --output-format fasta --taxa 33 --aligner mafft --cores 12 --log-path logs

python /public/uce/phyluce/bin/align/seqcap_align_2.py --fasta incomplete_fasta --output incomplete_aligned --output-format fasta --taxa 33 --aligner mafft --cores 12 --log-path logs  --incomplete-matrix

python /public/uce/phyluce/bin/align/get_only_loci_with_min_taxa.py --alignments incomplete_aligned --taxa 33 --input-format fasta --output incomplete_75perc_aligned --percent 0.75 --cores 12 --log-path logs

python /public/uce/phyluce/bin/align/add_missing_data_designators.py --alignments incomplete_75perc_aligned --output incomplete_75perc_aligned_w_missing --input-format fasta --output-format fasta --match-count-output /home/a499a400/gekko/incomplete.conf --incomplete-matrix incomplete.nostrict --cores 9 --log-path logs
```
#STEP 7B
```
python /public/uce/phyluce/bin/align/get_gblocks_trimmed_alignments_from_untrimmed.py --alignments complete_aligned --output complete_gbtrimmed --b2 0.65 --cores 9 --log-path logs

python /public/uce/phyluce/bin/align/get_gblocks_trimmed_alignments_from_untrimmed.py --alignments incomplete_75perc_aligned_w_missing --output incomplete_75perc_gbtrimmed --b2 0.65 --cores 9 --log-path logs 
```
#STEP 7C
```
python /public/uce/phyluce/bin/align/remove_locus_name_from_nexus_lines.py --taxa 33 --alignment complete_gbtrimmed --output complete_renamed --cores 9 --log-path logs

python /public/uce/phyluce/bin/align/remove_locus_name_from_nexus_lines.py --taxa 33 --alignment incomplete_75perc_gbtrimmed --output incomplete_75perc_renamed --cores 9 --log-path logs
```
#STEP 7D
```
python /public/uce/phyluce/bin/align/get_align_summary_data.py --alignments complete_renamed --input-format nexus --cores 9 --log-path logs

python /public/uce/phyluce/bin/align/get_align_summary_data.py --alignments incomplete_75perc_renamed --input-format nexus --cores 9 --log-path logs
```
#STEP 8A
```
python /public/uce/phyluce/bin/align/format_nexus_files_for_raxml.py --alignments complete_renamed/ --output complete_raxml --log-path logs

python /public/uce/phyluce/bin/align/format_nexus_files_for_raxml.py --alignments incomplete_75perc_renamed/ --output incomplete_75perc_raxml --log-path logs
```
#CONCATENATED RAXML ANALYSIS
```
raxmlHPC-AVX -s outgroup_final_bylocus_nexus.phylip -n run1 -m GTRCATI -f a -N 100 -x $RANDOM -p $RANDOM
```

#SPECIES TREE STUFF
(COMPLETE DATASET ONLY)
```
python /public/uce/phyluce/bin/align/convert_one_align_to_another.py --alignments complete_renamed/ --output complete_phylip --input-format nexus --output-format phylip --shorten-names --name-conf  short.conf --cores 8

#short.conf file used for gekko data
[taxa]
gekko_kikuchi_hofh89053101:ki_9053101
gekko_kikuchi_hofh89053102:ki_9053102
gekko_mindorensis_acd1235:mi_acd1235
gekko_mindorensis_acd1454:mi_acd1454
gekko_mindorensis_acd1661:mi_acd1661
gekko_mindorensis_acd3392:mi_acd3392
gekko_mindorensis_acd3519:mi_acd3519
gekko_mindorensis_cds063:mi__cds063
gekko_mindorensis_cds099:mi__cds099
gekko_mindorensis_cds102:mi__cds102
gekko_mindorensis_cds1461:mi_cds1461
gekko_mindorensis_cds1462:mi_cds1462
gekko_mindorensis_cds1463:mi_cds1463
gekko_mindorensis_cds278:mi__cds278
gekko_mindorensis_cds382:mi__cds382
gekko_mindorensis_cds591:mi__cds591
gekko_mindorensis_cds607:mi__cds607
gekko_mindorensis_cds757:mi__cds757
gekko_mindorensis_cdsgs40:mi_cdsgs40
gekko_mindorensis_elr373:mi__elr373
gekko_mindorensis_elr395:mi__elr395
gekko_mindorensis_rmb13575:mi_mb13575
gekko_mindorensis_rmb2870:mi_rmb2870
gekko_mindorensis_rmb3969:mi_rmb3969
gekko_mindorensis_rmb4093:mi_rmb4093
gekko_mindorensis_rmb4981:mi_rmb4981
gekko_mindorensis_rmb5005:mi_rmb5005
gekko_mindorensis_rmb515:mi__rmb515
gekko_mindorensis_rmb5342:mi_rmb5342
gekko_mindorensis_rmb5568:mi_rmb5568
gekko_mindorensis_rmb6429:mi_rmb6429
gekko_monarchus_rmb2991:mo_rmb2991
gekko_palawanensis_rmb7588:pa_rmb7588

python /public/uce/phyluce/bin/genetrees/phyluce_genetrees_run_raxml_genetrees.py --input complete_phylip/ --output complete_genetrees --outgroup pa_rmb7588 --cores 6 --quiet 2>&1 | tee raxml_genetrees.txt

rm complete_genetrees/*/RAxML_log*
rm complete_genetrees/*/RAxML_parsimony*
rm complete_genetrees/*/RAxML_result*
```
If you have problems with the Phyluce script for making genetrees by locus, you can use the more manual shell script given here:
https://github.com/laninsky/Phase_hybrid_from_next_gen/tree/master/post-processing

For these next steps, I uploaded files to the KU cluster and ran summary species tree analyses up there. First, created 'gekko' folder on cluster, then ran following commands from complab2
```
scp -r complete_phylip a499a400@transfer.acf.ku.edu:/scratch/a499a400/gekko
scp -r complete_genetrees a499a400@transfer.acf.ku.edu:/scratch/a499a400/gekko
```
The following steps are completed on the cluster. Make sure phyluce is in your pythonpath e.g. export 
```
PYTHONPATH=/scratch/oliveros/phyluce:$PYTHONPATH
```
qlogin in to do this next step:
```
python /scratch/oliveros/phyluce/bin/genetrees/phyluce_genetrees_multilocus_bootstrap_count_replicates.py --input complete_phylip/ --bootstrap_replicates complete.bootstrap.replicates --bootstrap_counts complete.bootstrap.counts.csv --bootreps 500
```
Create PBS script to run genetrees for bootstrapped data
-- Edit number of loci in -t parameter (0 - [numloci-1])
-- Edit outgroup in python command

Make sure to create the output directory before trying to run the script. Check for the log files that will appear to see if anything is going wrong. You do not have to specify an output group with the --outgroup option if your downstream summary species methods don't need a fixed outgroup (and if you have an incomplete matrix you won't be able to use this option anyway).
```

#PBS -N complete.boot
#PBS -l nodes=1:ppn=1:avx,mem=2000m,walltime=5:00:00
#PBS -M a499a400@ku.edu
#PBS -t 0-1993
#PBS -r n
#PBS -m n
#PBS -j oe
#PBS -o /dev/null
#PBS -d /scratch/a499a400/gekko

LINENUM=$(echo "(${PBS_ARRAYID}+1)" | bc)
PARAMETERS=$(awk -v line=${LINENUM} '{if (NR == line) { print $0; };}' ${PBS_O_WORKDIR}/complete.bootstrap.counts.csv)
unbuffer python /scratch/oliveros/phyluce/bin/genetrees/phyluce_genetrees_run_raxml_per_locus_bootstraps.py --input /scratch/a499a400/gekko/complete_phylip/ --output /scratch/a499a400/gekko/complete_bootstraps/ --best-trees /scratch/a499a400/gekko/complete_genetrees/ --outgroup pa_rmb7588 ${PARAMETERS} > ${PBS_O_WORKDIR}/$PBS_JOBNAME.log
```
If you get the following error, make sure you have raxml installed in your path. Don't forget to delete anything inside your bootstrap folder before trying again, as it will not overwrite and will die with an error:
```
Traceback (most recent call last):
  File "/scratch/oliveros/phyluce/bin/genetrees/phyluce_genetrees_run_raxml_per_locus_bootstraps.py", line 130, in <module>
    main()
  File "/scratch/oliveros/phyluce/bin/genetrees/phyluce_genetrees_run_raxml_per_locus_bootstraps.py", line 126, in main
    run_raxml(args.locus, alignment, args.outgroup, args.bootreps, seed, args.output, time, patterns)
  File "/scratch/oliveros/phyluce/bin/genetrees/phyluce_genetrees_run_raxml_per_locus_bootstraps.py", line 105, in run_raxml
    seconds = time.search(stdout).groups()[0]
AttributeError: 'NoneType' object has no attribute 'groups'
```

While the bootstrapping genetrees script was running, I went ahead and ran the scripts on original data. First you need to pull out the trees into a single file (back on complab2):
```
cd outgroup_genetrees;
touch inputgenetrees.tre;
for i in uce*; 
do cat $i/RAxML_bestTree.best >> inputgenetrees.tre
done;
```

Then to run ASTRID:
```
 /public/ASTRID-linux -i inputgenetrees.tre -o astrid.tre -m bionj
 ```
To run ASTRAL-II:
```
java -jar ~/bin/ASTRAL/astral.4.10.5.jar -i input genetrees.tre -o astral.tre
```

Now, going back to the bootstrapped data (make the complete_boottree directory before executing the command):
```
python /scratch/oliveros/phyluce/bin/genetrees/phyluce_genetrees_sort_multilocus_bootstraps.py --input complete_bootstraps/ --bootstrap_replicates complete.bootstrap.replicates --output complete_boottree
```
Removing files we don't need
```
rm -rf complete_phylip; rm -f complete.boot*.log
```

Pull down bootstrapped data and run 

#Summarizing bootstrapped trees
Navigate to where the bootstrapped trees are. If you get any error messages, make sure you have DendroPy-4.0.3 installed. Make sure to change the number for the head command to whatever the number of lines you need is

```
/scratch/a499a400/bin/bin/sumtrees.py -o fours.astral.con.tre boot*.astral.tre
/scratch/a499a400/bin/bin/sumtrees.py -o fours.star.con.tre boot*.star.tre
/scratch/a499a400/bin/bin/sumtrees.py -o fours.steac.con.tre boot*.steac.tre
cat boot*.tre | grep "tree mpest" > summary
head -n 37 boot0.tre > mpest.500.tre
cat summary >> mpest.500.tre
echo "End;" >> mpest.500.tre
/scratch/a499a400/bin/bin/sumtrees.py -o fours.mpest.con.tre mpest.500.tre
```

Old species tree methods for posterity: running on original data

```
cd complete_genetrees
cp all-best-trees.tre genetrees.rooted

R
library("phybase")
source("/scratch/oliveros/NJst.R") 
outgrouptaxon <- "pa_rmb7588"

mytrees<-read.tree("genetrees.rooted")
taxaname<-mytrees[[1]]$tip.label
speciesname<-taxaname
ntaxa<-length(taxaname)
ngene<-length(mytrees)
print(paste("Outgroup taxon", outgrouptaxon))

treestringrooted<-read.tree.string("genetrees.rooted",format="phylip")
treesrooted<-treestringrooted$tree

# STAR
species.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
diag(species.structure)<-1
print("Estimating STAR tree")
star<-star.sptree(treesrooted, speciesname, taxaname, species.structure, outgroup=outgrouptaxon,method="nj")
write.table(star,"genetrees.star.tre",row.names=F,col.names=F,quote=F,append=TRUE)

# STEAC
species.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
diag(species.structure)<-1
print("Estimating STEAC tree")
steac<-steac.sptree(treesrooted, speciesname, taxaname, species.structure, outgroup=outgrouptaxon,method="nj")
write.table(steac,"genetrees.steac.tre",row.names=F,col.names=F,quote=F,append=TRUE)

# NJst
treestringphy<-read.tree.string("genetrees.rooted",format="phylip")
treesphy<-treestringphy$tree
species.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
diag(species.structure)<-1
print("Estimating NJst tree")
njsttree<-NJst(treesphy, speciesname, taxaname, species.structure)
write.table(njsttree,"genetrees.njst.tre",row.names=F,col.names=F,quote=F,append=TRUE)
```
For astral and mpest, these trees have to be run directly on the command line
```
java -jar /scratch/a499a400/bin/ASTRAL/Astral/astral.4.7.7.jar -i genetrees.rooted -o genetrees.astral.tre 2>&1 | tee genetrees.phy.astral.out
```
I find it easier to wait to construct the mpest tree until created the bootstrapped mpest control files (can just copy these over). So, after doing the steps on the boostrapped data summarized:
Changing directory to boots to create our STEAC and STAR bootstrap R scripts
```
R
numboot <- 500
outgroup <- "pa_rmb7588"
wd <- "/scratch/a499a400/gekko/four_pis_genetrees/boots"

numbers<-seq(from=0,to=numboot-1)  
phybase<-'library("phybase")'
workdir<-paste('setwd("',wd,'")',sep="")
outg1<-paste('outgrouptaxon<-"',outgroup,'"',sep="")
outg2<-'print(paste("Outgroup taxon", outgrouptaxon))'

for(i in numbers)
{
    fname<-paste("starsteac",i,".R",sep="")
    genetree<-paste('genetreefname<-paste("boot",',i,', sep="")')

    details<-'mytrees<-read.tree(genetreefname)
    taxaname<-mytrees[[1]]$tip.label
    speciesname<-taxaname
    ntaxa<-length(taxaname)
    ngene<-length(mytrees)
    
    treestring<-read.tree.string(genetreefname,format="phylip")
    trees<-treestring$tree

    ptm<-proc.time()
    species.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
    diag(species.structure)<-1
    star<-star.sptree(trees, speciesname, taxaname, species.structure, outgroup=outgrouptaxon,method="nj")
    starfname<-paste(genetreefname,".star.tre", sep="")
    write.table(star,starfname,row.names=F,col.names=F,quote=F,append=FALSE)
    diff<-proc.time() - ptm
    print(paste("STAR tree estimation for", genetreefname, "completed in", diff[3], "seconds"))

    ptm<-proc.time()
    species.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
    diag(species.structure)<-1
    steac<-steac.sptree(trees, speciesname, taxaname, species.structure, outgroup=outgrouptaxon,method="nj")
    steacfname<-paste(genetreefname,".steac.tre", sep="")
    write.table(steac,steacfname,row.names=F,col.names=F,quote=F,append=FALSE)
    diff<-proc.time() - ptm
    print(paste("STEAC tree estimation for", genetreefname, "completed in", diff[3], "seconds"))'

    a<-paste(phybase,workdir,genetree,outg1,outg2,details,sep="\n")
    write.table(a,fname, row.names=F,col.names=F,quote=F)
}

# CREATE CONTROL FILES FOR MPEST

nboot<-numboot-1 
ntaxa<-33  # number of taxa
ngenes<- 1254  # number of genes

#species-allele table below
c<-"ki_9053101 1 ki_9053101
ki_9053102 1 ki_9053102
mi_acd1235 1 mi_acd1235
mi_acd1454 1 mi_acd1454
mi_acd1661 1 mi_acd1661
mi_acd3392 1 mi_acd3392
mi_acd3519 1 mi_acd3519
mi__cds063 1 mi__cds063
mi__cds099 1 mi__cds099
mi__cds102 1 mi__cds102
mi_cds1461 1 mi_cds1461
mi_cds1462 1 mi_cds1462
mi_cds1463 1 mi_cds1463
mi__cds278 1 mi__cds278
mi__cds382 1 mi__cds382
mi__cds591 1 mi__cds591
mi__cds607 1 mi__cds607
mi__cds757 1 mi__cds757
mi_cdsgs40 1 mi_cdsgs40
mi__elr373 1 mi__elr373
mi__elr395 1 mi__elr395
mi_mb13575 1 mi_mb13575
mi_rmb2870 1 mi_rmb2870
mi_rmb3969 1 mi_rmb3969
mi_rmb4093 1 mi_rmb4093
mi_rmb4981 1 mi_rmb4981
mi_rmb5005 1 mi_rmb5005
mi__rmb515 1 mi__rmb515
mi_rmb5342 1 mi_rmb5342
mi_rmb5568 1 mi_rmb5568
mi_rmb6429 1 mi_rmb6429
mo_rmb2991 1 mo_rmb2991
pa_rmb7588 1 pa_rmb7588"

for(i in 0:nboot)
{
    file<-paste("control",i,sep="")   #filename of control file
    treefile<-paste("boot",i,sep="")  #filename of input rooted tree
    b<-floor(runif(1)*799736+1111)  # creates random seed number
    a<-paste(treefile,"0",b,paste(ngenes, ntaxa),c ,"0",sep="\n")  # contents of control file including num genetrees and num species
    write.table(a, file,row.names=F,col.names=F,quote=F)
}
```

After creating all of the *.R files and control files, execute them as job arrays through PBS script as below

# JOB ARRAY STAR, STEAC
```
#PBS -N fours.ss
#PBS -l nodes=1:ppn=1:avx,mem=5000m,walltime=2:00:00
#PBS -M a499a400@ku.edu
#PBS -t 0-499
#PBS -r n
#PBS -m n
#PBS -j oe
#PBS -o /dev/null
#PBS -d /scratch/a499a400/gekko/four_pis_genetrees/boots

R --vanilla < starsteac${PBS_ARRAYID}.R > ${PBS_O_WORKDIR}/$PBS_JOBNAME.log
```

# JOB ARRAY MPEST
```
#PBS -N fours.mp
#PBS -l nodes=1:ppn=1:avx,mem=5000m,walltime=4:00:00
#PBS -M a499a400@ku.edu
#PBS -t 0-499
#PBS -r n
#PBS -m n
#PBS -j oe
#PBS -o /dev/null
#PBS -d /scratch/a499a400/gekko/four_pis_genetrees/boots

mpest control${PBS_ARRAYID} > ${PBS_O_WORKDIR}/$PBS_JOBNAME.log
````

# JOB ARRAY ASTRAL
```
#PBS -N fours.ast
#PBS -l nodes=1:ppn=1:avx,mem=25000m,walltime=4:00:00
#PBS -M a499a400@ku.edu
#PBS -t 0-499
#PBS -r n
#PBS -m n
#PBS -j oe
#PBS -o /dev/null
#PBS -d /scratch/a499a400/gekko/four_pis_genetrees/boots

unbuffer java -jar /scratch/a499a400/bin/ASTRAL/Astral/astral.4.7.7.jar -i boot${PBS_ARRAYID} -o boot${PBS_ARRAYID}.astral.tre > ${PBS_O_WORKDIR}/$PBS_JOBNAME.log
```

After doing this, copied one of the MPEST control files and modified it so the input tree was genetrees.rooted and then copied it to the complete_genetrees folder. Ran MPEST on the original data by:
```
mpest control
```
