#!/bin/bash


#SBATCH --time=2-02:30:00 # ask that the job be allowed to run for 2 days, 2 hours, 30 minutes.
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   
#SBATCH --mem=8G # Maximum amount of memory this job will be given
#SBATCH --job-name="MysSpon_phy"
#SBATCH --open-mode=append
#SBATCH --output=MysSpon.out
#SBATCH --mail-user=makhaki@iastate.edu
#SBATCH --mail-type=FAIL,END


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load mafft/7.505-3bx4zsi
module load gblocks
module load cdhit/4.8.1-7sf2emj
module load raxml-ng/1.1.0-q3rg2p6


#rename sequences:

#for i in *.fa;
#do sed 's/>[A-Za-z0-9]\+|\([A-Za-z0-9_]\+\)|.\+ OS=\([A-Za-z\-]\+\) \([a-zA-Z0-9\.]\+\) .\+/>\1_\2_\3/' $i > $(basename $i .fasta).ren ;
#do sed 's/>\(\w\+\)\..* \[\(\w\+\) \(\w\+\) .*\].*/>\1_\2_\3/' $i >$(basename $1 .fa).ren;
#done

#cd-hit 
#for i in *.fa
 #       do cd-hit -i $i -o $(basename $i .fa).cdhit -c 0.95
#done

#alignment:
for i in *.cdhit;
      do mafft --anysymbol --reorder --auto $i > $(basename $i .ren).aln ;
done

#filter alignments:
for i in *.aln;
        do Gblocks $i -b5=a -p=n  ;# -b3=15 -b4=5  -b2=0.5 ;
done

#RAxML tree building:
for i in *.aln ;
        do raxml-ng --all --msa $i --model GTR+G ;
done
