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

If you find QAlign is helpful, we appreciate your citation of its pre-print [version](https://www.biorxiv.org/content/10.1101/862813v1):

Dhaivat Joshi, Shunfu Mao, Sreeram Kannan, Suhas Diggavi. QAlign: Aligning nanopore reads accurately using current-level modeling.

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
                              --qlevel [quantization_level e.g. 2 or 3]
                              --rc [enable_rev_complementary e.g. 1 or 0]
                              --kmerpath [/path/to/kmermodel]
```

The quantized sequence (and its rev complementary) will be stored in outdir.

