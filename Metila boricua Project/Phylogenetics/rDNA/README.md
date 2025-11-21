# README

This folder contain the alignments and trees generated for phylogenetic analysis using the nuclear 18S and 28S ribosomal gene. 

The contents are as follows:

- README: This file
- 18S_Tree.aln: mutliple sequence alignment of 18S genes
- 18S_Tree.aln-gb: mutliple sequence alignment of 18S genes, with non-informative sites removed using Gblocks
- 18S_support.tre: Maximum likelihood tree using the above alignment, with bootstrap support values
- 28S_Tree.aln:  mutliple sequence alignment of 28S genes
- 28S_Tree.aln-gb: mutliple sequence alignment of 18S genes, with non-informative sites removed using Gblocks
- 28S_support.tre: Maximum likelihood tree using the above alignment, with bootstrap support values

Sequences were aligned using MAFFT (v7.508, Katoh and Standley 2013) with the “--auto” option, and highly variable sites were removed using GBlocks (Castresana 2000, -b5 = a, other settings default). The final alignments contained 880 positions for 18S and 2042 positions for 28S. Maximum likelihood (ML) trees were constructed using RAxML-NG (v1.1.0, Kozlov et al., 2019) using the GTR+G model, using 1000 bootstrap replicates. 