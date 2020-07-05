# Python_scripts

Calculate_N50_gc_N_of_genome.py


```
usage: python3 Calculate_N50_gc_N_of_genome.py -i <genome.fa/genome.fa.gz> -o <output.txt> -m <N50|len>

This script is used to calculate the scaffold length, N50, GC% and N% of the genome.

required arguments:
  -i FASTA FILE, --in FASTA FILE
                        Input sequence file in FASTA format
  -o OUTPUT, --out OUTPUT
                        Output folders and files will be labelled with this name. WARNING: do not provide a path
  -m MODE, --mode MODE  Specify which analysis mode to run.
                        There are two valid modes:
                        -m len for caculating the length, GC%, and N% of each scaffold or genome
                        -m N50 for calcuting the scaffold N50 of the genome

optional arguments:
  -h, --help            Show this help message and exit
  ```
