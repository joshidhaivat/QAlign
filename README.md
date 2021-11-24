# QAlign 

### Contents <a id='contents'></a>

* <a href='#intro'>Introduction</a>
* <a href='#pub'>Publication</a>
* <a href='#setup'>Setup</a>
* <a href='#use'>Usage</a>

---

### Introduction <a id='intro'></a>

QAlign is a tool to apply quantization for nanopore long reads.

---

### Publication <a id='pub'></a>

If you find QAlign is helpful, we appreciate your citation of its publication [version](https://doi.org/10.1093/bioinformatics/btaa875):

Dhaivat Joshi, Shunfu Mao, Sreeram Kannan, Suhas Diggavi, QAlign: aligning nanopore reads accurately using current-level modeling, Bioinformatics, Volume 37, Issue 5, 1 March 2021, Pages 625–633.

---

### Setup <a id='setup'></a>

QAlign has been tested under Ubuntu 16.04. Please follow the below three steps to setup:

###### Step 1
Download QAlign.
```
$ git clone https://github.com/joshidhaivat/QAlign.git
```

###### Step 2
At the root directory of QAlign, run setup.py
```
$ python setup.py install
```

### Usage <a id='use'></a>

###### Convert
Convert a fasta file to its quantized version. At the root folder of qalign, run:

```
python qalign/main.py convert --input_fasta [/path/to/input/fasta/read]
                              --outdir [/path/to/output/folder]
                              --qlevel [quantization_level e.g. 2 or 3, default = 2]
                              --rc [enable_rev_complementary e.g. 1 or 0, default = 1]
                              --kmerpath [/path/to/kmermodel]
```

The quantized sequence (and its rev complementary) will be stored in outdir.


An example for conversion to quantized sequence using the test samples:

```
python qalign/main.py convert --input_fasta qalign/test_samples/reads.fasta --outdir qalign/test_samples/output/ --qlevel 2 --rc 1 --kmerpath qalign/kmermap_models/r9.4_6mer_nucleotide_template_model.txt
```

###### Alignment of quantized sequences
Align the quantized reads to quantized genome using minimap2

Alignment of forward strand Q2 reads to Q2 reference:
```
minimap2 -c -k 23 ref_q2.fasta res.fasta
```

Alignment of reverse strands Q2 reads to Q2 reference:
```
minimap2 -c -k 23 ref_q2.fasta rc_res.fasta
```
