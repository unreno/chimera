#	chimera (aka HERVnGoSeq)


The purpose of these scripts is to locate the insertion points of one reference in another.
This is done by aligning a collection of reads to a reference, usually viral, and selecting those reads that soft aligned at and past the ends. 
We call this type of a read "chimeric" and can be at the beginning (pre) or ending (post) of the reference.
The part of the read that aligned is removed and new fasta files are created containing these chimeric reads.
These fasta files are the aligned to the second reference, usually human, and the insertion points, the aligned position at which the initial read was snipped, are extracted.
The pre and post aligned insertion points are then compared and any within 10bp are considered an overlapper and are suggestive of a novel insertion of the first reference within the second.





##	Requirements

These scripts were written for use in a Unix/Linux/MacOSX environment
running bash-4.4.0, gawk-4.1.4,
[bowtie2](https://github.com/BenLangmead/bowtie2) 2.2.7 and
[samtools](https://github.com/samtools/samtools) 1.2 (using htslib 1.2.1).




##	Installation

There is only 1 script at the moment, so this is a bit excessive.
Simply copy it to wherever you like.

Eventually, however, ...

These scripts can be run from wherever they are. 
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
chimera_insertions_and_overlappers.bash --viral HCMV --human hg19 Sample_1.fastq Sample_2.fastq
```

will generate

```BASH
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.bam
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.pre.fasta
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.post.fasta
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.pre.bowtie2.hg19.bam
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.post.bowtie2.hg19.bam
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.pre.bowtie2.hg19.Q20.insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.post.bowtie2.hg19.Q20.insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.pre.bowtie2.hg19.Q20.rc_insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.post.bowtie2.hg19.Q20.rc_insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.both.bowtie2.hg19.Q20.insertion_points.overlappers
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.both.bowtie2.hg19.Q20.rc_insertion_points.rc_overlappers
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.pre.bowtie2.hg19.Q10.insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.post.bowtie2.hg19.Q10.insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.pre.bowtie2.hg19.Q10.rc_insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.post.bowtie2.hg19.Q10.rc_insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.both.bowtie2.hg19.Q10.insertion_points.overlappers
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.both.bowtie2.hg19.Q10.rc_insertion_points.rc_overlappers
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.pre.bowtie2.hg19.Q00.insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.post.bowtie2.hg19.Q00.insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.pre.bowtie2.hg19.Q00.rc_insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.post.bowtie2.hg19.Q00.rc_insertion_points
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.both.bowtie2.hg19.Q00.insertion_points.overlappers
SAMPLE.bowtie2.HCMV.__very_sensitive_local.aligned.both.bowtie2.hg19.Q00.rc_insertion_points.rc_overlappers
```




##	Docker

The following "Dockerfile" will create a functioning chimera install.

```BASH
FROM ubuntu

ENV DEBIAN_FRONTEND noninteractive
 
RUN apt-get update && apt-get install -y apt-utils dialog bzip2 gcc gawk zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev make libssl-dev libncurses5-dev zip g++ git libtbb-dev wget && apt-get clean

RUN cd / && wget https://github.com/samtools/htslib/releases/download/1.5/htslib-1.5.tar.bz2 && tar xvfj htslib-1.5.tar.bz2 && cd htslib-1.5 && ./configure && make && make install && cd ~ && /bin/rm -rf /htslib-1.5*
 
RUN cd / && wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 && tar xvfj samtools-1.5.tar.bz2 && cd samtools-1.5 && ./configure && make && make install && cd ~ && /bin/rm -rf /samtools-1.5*

RUN cd / && wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-source.zip/download -O bowtie2-2.3.2-source.zip && unzip bowtie2-2.3.2-source.zip && cd bowtie2-2.3.2 && make && make install && cd ~ && /bin/rm -rf /bowtie2-2.3.2*

RUN cd ~ && git clone http://github.com/unreno/chimera && cd chimera && ln -s Makefile.example Makefile && make BASE_DIR="/usr/local" install && cd ~ && /bin/rm -rf chimera
```




## Notes

More descriptive named script chimera_unpaired_local.bash is the same as chimera_insertions_and_overlappers.bash.


Deving chimera_paired_local.bash



