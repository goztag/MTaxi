# MTaxi

A comparative tool for taxon identification of low coverage ancient genomes.

MTaxi is a method to classify pairs of closely related species using ancient DNA data, by utilizing mitochondrial tranversion-type substitions between the two species.

It is developed and tested using both Linux operating system (Ubuntu 18.04 LTS) and WSL2 (Windows Subsystem for Linux 2), should work in other Unix-like platforms.

MTaxi is implemented with sheep and goat genomes in this repository, but can be used with arbitrary pairs if mtDNA references are available.

For further details, please refer to .....

<br />

### **Requirements**


→ **Working environment**

- Python 3+ (tested with 3.6)

→ Bedtools ([https://bedtools.readthedocs.io/en/latest/](https://bedtools.readthedocs.io/en/latest/))
    

→ **Required libraries**

 - Available in ```requirements.txt```. Install with ```pip install -r requirements.txt```.


<br />

### File Structure


**choar_analysis.py :**
- choar_analysis implementer file. Call this to run MTaxi.


**settings.py :**

→ File locations are set here !

→ Change sp1_ref_file, sp2_ref_file, transv_poly_sp1, transv_poly_sp2 files.
 - sp1_ref_file: species 1 reference file
 - transv_poly_sp1: species 1 reference file

<br />

### **How to install the program via GitHub**

```bash
git clone https://github.com/...
```

<br />

### **How to run** ###
**command line arguments:**

→ Necessary arguments:
 - "-sp1": To specify the sample aligned to species 1 with .bam extension.
 - "-sp2": To specify the .bam file aligned to species 2
 - "-out": To specify the output file name.
	
Output format is set as .txt file.File name does not need to contain the extension, .txt. However, if it is to contain, it should only contain ".txt" extension.

→ Optional arguments:

"-d" : To use debug mode

**The order of arguments does not matter.**

→ example run on bash:

```bash
python3 choar_analysis.py -sp1 sample_sheep.bam -sp2 sample_goat.bam -out outf

```
- sample_sheep.bam : sample aligned to sheep mitochondrial reference
- sample_goat.bam : sample aligned to goat mitochondrial reference
- outf: Name of the output file. The output file will be created as "outf.txt" under the same directory as that of choar_analysis.py

<br />

### How to run with debug mode

Use debug mode to see outputs of each step on the console.

In addition to this, debug mode enables you to have each intermediate file be created on data_out folder. Otherwise, only needed files are generated.

Debug mode is off by default. To activate debug mode run the following command :

```bash
python3 choar_analysis.py -sp1 sample_sheep.bam -sp2 sample_goat.bam -out outf -d

```
