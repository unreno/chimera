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
chimera.bash --viral HCMV --human hg19 Sample_1.fastq Sample_2.fastq
```



##	Docker

Under development...


The following "Dockerfile" will create a functioning chimera install.

Install Docker, make a new empty temp directory and create a file "Dockerfile" with the following contents.

```BASH
FROM ubuntu

ENV DEBIAN_FRONTEND noninteractive
ENV SAMTOOLS_VERSION="1.5"
ENV BOWTIE2_VERSION="2.3.2"
WORKDIR=$HOME
 
RUN apt-get update && apt-get install -y apt-utils dialog bzip2 gcc gawk zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev make libssl-dev libncurses5-dev zip g++ git libtbb-dev wget && apt-get clean

RUN cd / && wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 && tar xvfj htslib-${SAMTOOLS_VERSION}.tar.bz2 && cd htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && cd ~ && /bin/rm -rf /htslib-${SAMTOOLS_VERSION}*
 
RUN cd / && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && tar xvfj samtools-${SAMTOOLS_VERSION}.tar.bz2 && cd samtools-${SAMTOOLS_VERSION} && ./configure && make && make install && cd ~ && /bin/rm -rf /samtools-${SAMTOOLS_VERSION}*

RUN cd / && wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-source.zip/download -O bowtie2-${BOWTIE2_VERSION}-source.zip && unzip bowtie2-${BOWTIE2_VERSION}-source.zip && cd bowtie2-${BOWTIE2_VERSION} && make && make install && cd ~ && /bin/rm -rf /bowtie2-${BOWTIE2_VERSION}*

RUN cd ~ && git clone http://github.com/unreno/chimera && cd chimera && ln -s Makefile.example Makefile && make BASE_DIR="/usr/local" install && cd ~ && /bin/rm -rf chimera
```

Afterward, from within this temp directory, run `docker build -t chimera .` to create the image.
It will take 5-10 minutes.
Then run `docker run -ti chimera` to start a new instance or container.

The resulting docker image will be nearly 500MB.


