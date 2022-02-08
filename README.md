# TAXIDAG 

A comparative tool for taxon identification of low coverage ancient genomes.

TAXIDAG is an elaborative method  to distinguish between two closely related ancient genomes by utilizing mtDNA and finding different identical sites between these species.

TAXIDAG is developed and tested using WSL2 (Windows Subsystem for Linux 2). 

TAXIDAG is implemented with two similar species' ( Sheep and Goat ) genomes in this repository.
But can also be used with arbitrary two species.

For further details, please look at the poster in docs folder.

### **Requirements**

→ **Working environment**
- Python 3+ (tested with 3.6)

-Sam Tools
    Ubuntu/WSL installation :
    
    ``` sudo apt update
    sudo apt install python-pysam ```

- Bedtools ([https://bedtools.readthedocs.io/en/latest/](https://bedtools.readthedocs.io/en/latest/))
    Ubuntu/WSL installation : `sudo apt-get install bedtools` 

→ **Required libraries**

Available in ```requirements.txt```. Install with ```pip install -r requirements.txt```.

### File Structure

**docs/ :**

- Project documents available.

**choar_analysis.py :** 
- choar_analysis implementer file. Call this to run TAXIDAG.

**utils.py :** utility functions 

**settings.py :**
- File locations are set here !
- Change current_sample_dir to the directory name of sample data in tmp folder.
- Change sp1_ref_file, sp2_ref_file, transv_poly_sp1, transv_poly_sp2 files. 
	(sp1_ref_file: species 1 reference file)

### **How to install the program via GitHub**

```bash
git clone https://github.com/SevcanDogramaci/TAXIDAG.git
```

### **How to run** ###
command line arguments:

Necessary arguments:
"-sp1": To specify the aligned sample to species 1 with .bam extension.
"-sp2": To specify the .bam file of species 2
"-out": To specify the output file name.
	Output format is set as .txt file. 
	File name does not need to contain the extension, .txt. However, if it is to contain, it should only contain ".txt" extension.

Optional arguments:
"-d" : To use debug mode
"keep_trans": if given as command line argument, transitions will be kept. (By default, transitions are not kept).

**The order of arguments do not matter.

example run on bash:

```bash
python3 choar_analysis.py -sp1 o_dp_tps083_MT.bam -sp2 c_dp_tps083_MT.bam -out outf 
```
o_dp_tps083_MT.bam : bam file for species 1. Should be in the same path as choar_analysis.py. Or full path needs to be given as argument.
c_dp_tps083_MT.bam : bam file for species 2. Should be in the same path as choar_analysis.py. Or full path needs to be given as argument.
outf: Name of the output file. The output file will be created as "outf.txt" under the same directory as that of choar_analysis.py 


### How to run with debug mode

Use debug mode to see outputs of each step on the console. 

In addition to this, debug mode enables you to have each intermediate file be created on data_out folder. Otherwise, only needed files are generated.

Debug mode is off by default. To activate debug mode run the following command :

```bash
python3 choar_analysis.py -sp1 o_dp_tps083_MT.bam -sp2 c_dp_tps083_MT.bam -out outf -d
```

