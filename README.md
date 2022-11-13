# merge_fragments

usage: liftoffFragmentsV2.py [-h] [-i I] -fi <input.gff> -fa <target.fasta> [-min <min_intron_length_bp>] [-max <max_intron_length_bp>] -o <output.gff>

Check for copies in Liftoff fragments

optional arguments:

  -h, --help            show this help message and exit
  
  -i I                  Proportion threshold for copy identification
  
  -fi <input.gff>       desired input file
  
  -fa <target.fasta>    Fasta of target
  
  -min <min_intron_length_bp> min intron bp length
  
  -max <max_intron_length_bp> max intron bp length
  
  -o <output.gff>       desired output file
