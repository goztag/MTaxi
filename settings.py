import pathlib


# input file locations
# change file names according to the valid files to work with
 
current_dir = pathlib.Path(__file__).parent.absolute()

sp1_ref_file = f"{current_dir}/sp1_mt.fasta"
sp2_ref_file = f"{current_dir}/sp2_mt.fasta"
transv_poly_sp1_file = f"{current_dir}/trvpos_sp1.bed"
transv_poly_sp2_file = f"{current_dir}/trvpos_sp2.bed"
