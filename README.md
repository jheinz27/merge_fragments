# merge_fragments

## Description 
An optional post-alignment modification tool that identifies gene fragments mapping next to each in the Liftoff generated target genome annotation and accordingly merges the fragments into a single annotation. Mosaic flags neighboring genes with identical gene names. It then determines whether the annotation is split in the middle of an exon, an intron or is correct in it being two distinct annotations. It then uses this information to merge the fragmented features into a single annotation. This results in a final output GFF file with all gene fragments that mapped next to each other annotated as one gene. 

The tool identifies neighboring  genes by checking if their gene name fields are identical. If the tool identifies a set of sequential genes on the same strand and chromosome, it begins to evaluate this set. The majority of the time, this set will only be two fragments, as it is unlikely that a gene was split across three or more contigs. The first evaluation checks whether any of the neighboring annotations are copies of each other using pairwise2 from the biopython project (Cock et al., 2009). The threshold proportion of identical sequence for genes to be considered copies of each other can be set by the user as a command line input, but the default is set to 1.00 or 100% identical sequence between the genes. Mosaic does not modify the copies and for further steps only considers the genes fragments that are not copies, if there are any. 

Next, the tool evaluates the sequence length between the fragments identified in the previous step. The length is calculated by finding the distance between the end location of the previous gene and the starting location of the current gene. Based on this distance, the tool determines if the annotation was split in the middle of an exon, intron, or neither if the fragments are too far apart to feasibly be part of the same gene model. Shorter distances between the mapped gene fragments is evidence of an exon split. The user can input the maximum distance for a split to qualify as an exon split, however, the default is set to 80 base pairs. In these cases, the tool merges the fragmented exon into a single exon annotation and subsequently merges the parent transcript/mRNA features and the parent gene. If the exon has an annotated CDS, it is also merged. An intron split is defined as a gene split that occurs in between two exons. This is characterized by longer distances in between sequential genes. The maximum length of an intron split can again be set by the user, as long as it is greater than the maximum exon split length, however, the default is set to 500,000 bp. In these cases, only transcript/mRNA features and gene features need to be merged. The tool's final output is a GTF annotation file identical to the input, but with all fragmented features  having been merged. 

## Usage
liftoffFragmentsV2.py [-h] [-i I] -fi <input.gff> -fa <target.fasta> [-min <min_intron_length_bp>] [-max <max_intron_length_bp>] -o <output.gff>

Check for copies in Liftoff fragments

optional arguments:

  -h, --help            show this help message and exit
  
  -i I                  Proportion threshold for copy identification
  
  -fi <input.gff>       desired input file
  
  -fa <target.fasta>    Fasta of target
  
  -min <min_intron_length_bp> min intron bp length
  
  -max <max_intron_length_bp> max intron bp length
  
  -o <output.gff>       desired output file
