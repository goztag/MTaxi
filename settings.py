import pathlib


# input file locations
# change file names according to the valid files to work with
 
current_dir = pathlib.Path(__file__).parent.absolute()

sp1_ref_file = f"{current_dir}/data/asm_mt.fasta"
sp2_ref_file = f"{current_dir}/data/equc_mt.fasta"
transv_poly_sp1_file = f"{current_dir}/data/tmp/asi_cab_bed/ttrvposasi.bed"
transv_poly_sp2_file = f"{current_dir}/data/tmp/asi_cab_bed/ttrvposcab.bed"