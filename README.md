#	chimera


The purpose of these scripts is to locate the insertion points of one reference in another.
This is done by aligning a collection of reads to a reference, usually viral, and selecting those reads that soft aligned at and past the ends. 
We call this type of a read "chimeric" and can be at the beginning (pre) or ending (post) of the reference.
The part of the read that aligned is removed and new fasta files are created containing these chimeric reads.
These fasta files are the aligned to the second reference, usually human, and the insertion points are extracted.
The pre and post aligned insertion points are then compared and any within 10bp are considered an overlapper and are suggestive of a novel insertion of the first reference within the second.





##	Requirements

These scripts were written for use in a Unix/Linux/MacOSX environment
running bash-4.4.0, gawk-4.1.4, bowtie2-2.2.7 and samtools 1.2 (using htslib 1.2.1).




##	Installation

These scripts can be run from whereever they are. 
They don't need to be installed.
However, should you wish to install them in ~/.local/bin/ or elsewhere ...

Link or copy Makefile.example and edit as appropriate.

```BASH
ln -s Makefile.example Makefile
```

Then install. Installation will attempt to remove previous installs before installing.

```BASH
make install
```


##	Preparation

The alignments are done using bowtie2, so a bowtie2 index must be prepared for the human and viral references. The script will expect these indexes to either be in a directory pointed to by the environment variable $BOWTIE2_INDEXES or in $HOME/BOWTIE2_INDEXES.

For example ...
```BASH
bowtie2-build hg19.fasta hg19
bowtie2-build herv_k113.fasta herv_k113
bowtie2-build HCMV.fasta HCMV
```

##	Usage

If these scripts have not been installed to a directory in your path, prefix the command with the path to the location of this script.

```BASH
cd SOME/SAMPLE
chimera_insertions_and_overlappers.bash --human hg19 --viral HCMV Sample_1.fastq Sample_2.fastq
```

