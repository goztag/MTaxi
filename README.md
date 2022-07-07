# MTaxi

A comparative tool for taxon identification of low coverage ancient genomes.

MTaxi is a method to classify pairs of closely related species using ancient DNA data, by utilizing mitochondrial tranversion-type substitutions between the two species.

It is developed and tested using both Linux operating system (Ubuntu 18.04 LTS) and WSL2 (Windows Subsystem for Linux 2), should work in other Unix-like platforms.

MTaxi is implemented with sheep and goat genomes in this repository, but can be used with arbitrary pairs if mtDNA references are available.


<br />

### **Requirements**


→ **Working environment**

- Python 3+ (tested with 3.6)

→ Bedtools ([https://bedtools.readthedocs.io/en/latest/](https://bedtools.readthedocs.io/en/latest/))
	- make sure it's on the path
    

→ **Required libraries**

 - Available in ```requirements.txt```. Install with ```pip3 install -r requirements.txt```.


<br />

### File Structure


**run_MTaxi.py :**
- run_MTaxi implementer file, call this to run MTaxi.


**settings.py :**

File locations are set here !

→ Edit the paths to sp1_ref_file, sp2_ref_file, transv_poly_sp1, transv_poly_sp2 files.
 - sp1_ref_file: species 1 mitochondrial reference file (in FASTA format)
 - transv_poly_sp1: species 1 substition positions file (in BED format)

<br />

### **How to install the program via GitHub**

```bash
git clone https://github.com/goztag/MTaxi.git
```

<br />

### **How to run** ###
**Command line arguments:**

→ Necessary arguments:
 - "-sp1": To specify the sample aligned to species 1 with .bam extension
 - "-sp2": To specify the .bam file aligned to species 2
 - "-out": To specify the output file name (format is set as .txt file)
	

→ Optional arguments:

"-shared" : To restrict the analysis to to those that aligned to both species’ references. As default, all reads (the reads that could be mapped to both and the ones that could map only to one of the species’ references) are included in the analysis.

"-d" : To use debug mode

**The order of arguments does not matter.**

→ Example run on bash:

```bash
python3 run_MTaxi.py -sp1 sample_sheep.bam -sp2 sample_goat.bam -out outf

```

```bash
python3 run_MTaxi.py -sp1 sample_sheep.bam -sp2 sample_goat.bam -shared -out outf

```

- sample_sheep.bam : Sample aligned to sheep mitochondrial reference
- sample_goat.bam : Sample aligned to goat mitochondrial reference
- outf: Name of the output file, which will be created as "outf.txt"


<br />

### How to run with debug mode

Use debug mode to see outputs of each step on the console.

In addition to this, debug mode enables you to have each intermediate file be created on data_out folder. Otherwise, only the required files are generated.

Debug mode is off by default. To activate debug mode run the following command :

```bash
python3 run_MTaxi.py -sp1 sample_sheep.bam -sp2 sample_goat.bam -out outf -d

```

### License

This project is licensed under the terms of the Creative Commons - Attribution 4.0 International (CC BY 4.0) license.
