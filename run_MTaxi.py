import sys
import pysam
import pybedtools
import os
import shutil
import pathlib
import pandas as pd

import settings

from os import path, mkdir


import warnings
warnings.filterwarnings("ignore")

# outputs file locations
current_dir = pathlib.Path(__file__).parent.absolute()
data_out_dir = f"{current_dir}/tmp"
if not path.exists(data_out_dir):
    mkdir(data_out_dir)

sorted_aligned_sample_to_sp1_file = f"{data_out_dir}/sorted_sp1_dp_MT.bam"
sorted_aligned_sample_to_sp2_file = f"{data_out_dir}/sorted_sp2_dp_MT.bam"
transv_sample_sp1_file = f"{data_out_dir}/transv_sp1_samp.bed"
transv_sample_sp2_file = f"{data_out_dir}/transv_sp2_samp.bed"
sample_sp1_bamtobed_file = f"{data_out_dir}/samp_sp1_btb.bed"
sample_sp2_bamtobed_file = f"{data_out_dir}/samp_sp2_btb.bed"
shared_sp1_file = f"{data_out_dir}/shared_sp1.bed"
shared_sp2_file = f"{data_out_dir}/shared_sp2.bed"
uniq_intersections_sample_sp1_file = f"{data_out_dir}/uniq_sp1.bed"
uniq_intersections_sample_sp2_file = f"{data_out_dir}/uniq_sp2.bed"


def write_to_file(data, path):
    f = open(path, "w")
    
    for item in data:
        f.write(str(item))
    f.close()


def filter_indels(variant_calls):
    import re

    variant_call_ids_to_filter = []
    alt_base_column = []
    for i, read_result in variant_calls['READ_RESULTS'].iteritems():

        # search indels
        ins_search = "\+[0-9]+[ACGTNacgtn]+"
        del_search = "-[0-9]+[ACGTNacgtn]+"
        ins = re.search(ins_search, read_result)
        dels = re.search(del_search, read_result)

        if ins == None and dels == None:
            bases = "[ACGTNacgtn]"

            # get alt base by searching
            alt = re.search(bases, read_result)
            if alt != None:
                # get first alt base 
                alt = alt.group(0)
            else:
                alt = '.'
            alt_base_column.append(alt)
            if DEBUG:
                print(i, alt, variant_calls["POS"][i])
        else:
            row_id_to_delete = variant_calls.index[i]
            variant_call_ids_to_filter.append(row_id_to_delete)
    
    if DEBUG:
        print("\nLength of variant calls before: ", len(variant_calls), \
        " Length of variant calls to delete: ", len(variant_call_ids_to_filter))

    variant_calls = variant_calls.drop(variant_call_ids_to_filter,inplace=False)

    if DEBUG:
        print("Length of variant calls after: ", len(variant_calls), "\n")

    variant_calls.insert(len(variant_calls.columns), 'ALT', alt_base_column)
    return variant_calls
    

def filter_columns(variant_calls, columns_to_filter):
    return variant_calls.drop(columns_to_filter, axis=1, inplace=False)


def filter_variant_calls_from_pileup_format(variant_calls):
    from io import StringIO
    pileup_format_columns = ['CHROM', 'POS', 'REF', 'READ_COUNT', 'READ_RESULTS', 'READ_QUALITY']

    # convert variant calls to dataframe for easier data manipulation
    variant_calls = StringIO(variant_calls)
    variant_calls_dataframe = pd.read_table(variant_calls, names=pileup_format_columns)

    if DEBUG:
        print("\n\n>>> FILTER INSERTION & DELETIONS <<<")
    variant_calls_dataframe = filter_indels(variant_calls_dataframe)
    if DEBUG:
        print(variant_calls_dataframe.head(20))
    
    if DEBUG:
        print("\n\n>>> FILTER COLUMNS <<<")
    columns_to_delete = ['READ_COUNT', 'READ_QUALITY', 'READ_RESULTS']
    variant_calls_dataframe = filter_columns(variant_calls_dataframe, columns_to_delete)
    if DEBUG:
        print(variant_calls_dataframe.head(20))

    return variant_calls_dataframe



def get_variants_info(data, transv_sample_file):
    if DEBUG:
        print("\n\n>>> CREATE BED FILE <<<")

    new_data = [data["CHROM"], data["POS"].apply(lambda item: (item-1)), 
                data["POS"], data["ALT"]]
    new_data_headers = ["CHROM", "POS-1", "POS", "ALT"]
    dataframe = pd.concat(new_data, axis=1, keys=new_data_headers)


    if DEBUG:
        print(dataframe.head(20))
        dataframe.to_csv(transv_sample_file, header=False, index=False, sep="\t", mode="w")

    return dataframe


def create_shared_data_for_species(shared_data, columns, species_symbol):
    species_data = {}

    for column in columns:
        if column == "name":
            shared_column_name = column
        else:
            shared_column_name = f'{column}_{species_symbol}'
        species_data[column] = shared_data[shared_column_name]

    return pd.DataFrame(species_data)


def find_shared_reads(chi_sample, oar_sample):

    chi_btb  = chi_sample.to_dataframe()
    oar_btb  = oar_sample.to_dataframe()
    
    
    if SHARED :
        
        oar_btb["id"] = range(1, len(oar_btb)+1)
        
        shared2=chi_btb.merge(oar_btb, on="name")
        shared2.sort_values(by=["id"])

        shared_chi = create_shared_data_for_species(shared2, chi_btb.columns, "x")
        shared_oar = create_shared_data_for_species(shared2, chi_btb.columns, "y")
        
    if not SHARED:
    
        shared_chi = chi_btb
        shared_oar = oar_btb
        

    
    if DEBUG:
        shared_chi.to_csv(shared_sp1_file, header=False, index=False, sep="\t", mode="w")
        shared_oar.to_csv(shared_sp2_file, header=False, index=False, sep="\t", mode="w")



    return (shared_chi, shared_oar)


def find_alt_freqs(all_intersects_chi, all_intersects_oar):

    
    all_intersects = pd.merge(all_intersects_chi, all_intersects_oar, on="name")
    
    if not SHARED:
        sp1_unique = pd.concat([all_intersects_chi, all_intersects_oar,all_intersects_chi]).drop_duplicates(subset= 'name',keep=False)


        sp2_unique = pd.concat([all_intersects_oar, all_intersects_chi,all_intersects_oar]).drop_duplicates(subset= 'name',keep=False)


    dupl_indices = dict()
    df = pd.DataFrame(all_intersects)

    #print("dataframee \n", df)
    
    dups = df[all_intersects.duplicated('name', keep=False)]

  

    for i, name in zip(dups.index, dups['name']):
        if name not in dupl_indices.keys():
            dupl_indices[name] = [i]
        else: dupl_indices[name].append(i)

    #print("DUP INDICES TO BE DROPPED ARE::: \n", dupl_indices.values())
    #print()



    for indices in dupl_indices.values():
        i = indices[0]   

        for k in range(1,len(indices)):
            j = indices[k]

            all_intersects.at[i,'Ref_x'] += all_intersects.at[j,'Ref_x']
            all_intersects.at[i,'Alt_x'] += all_intersects.at[j,'Alt_x']
            all_intersects.at[i,'Total_x'] += all_intersects.at[j,'Total_x']
            all_intersects.at[i,'Ref_y'] += all_intersects.at[j,'Ref_y']
            all_intersects.at[i,'Alt_y'] += all_intersects.at[j,'Alt_y']
            all_intersects.at[i,'Total_y'] += all_intersects.at[j,'Total_y']

            all_intersects.at[j,'Alt_x'] = -1
            all_intersects.at[j,'Alt_y'] = -1


    
    if DEBUG:
        print("\n\n>>> ALL INTERSECTIONS <<<")
        print(all_intersects)

    all_intersects = all_intersects.query('(Alt_y<=Ref_x | Alt_x<=Ref_y) & (Alt_y>=Ref_x | Alt_x>=Ref_y)')

    all_intersects.insert(    
        len(all_intersects.columns),
        "Sp1_Alt_Freq", (all_intersects["Alt_y"] / all_intersects["Total_y"])
    )
    all_intersects.insert(    
        len(all_intersects.columns),
        "Sp2_Alt_Freq", (all_intersects["Alt_x"] / all_intersects["Total_x"])
    )

    if not SHARED:
        sp1_unique.insert(
            len(sp1_unique.columns),
            "Sp1_Alt_Freq", (sp1_unique["Alt"] / sp1_unique["Total"])
        )
        
        sp1_unique.insert(
            len(sp1_unique.columns),
            "Sp2_Alt_Freq", (1 - (sp1_unique["Alt"] / sp1_unique["Total"]))
        )
        
        sp2_unique.insert(
            len(sp2_unique.columns),
            "Sp1_Alt_Freq", (1 - (sp2_unique["Alt"] / sp2_unique["Total"]))
        )
        
        sp2_unique.insert(
            len(sp2_unique.columns),
            "Sp2_Alt_Freq", (sp2_unique["Alt"] / sp2_unique["Total"])
        )
    
    
    shared_alt_freqs = pd.concat([all_intersects["Sp1_Alt_Freq"], all_intersects["Sp2_Alt_Freq"]], axis=1)
    
    if SHARED:
        all_alt_freqs = shared_alt_freqs
    
    if not SHARED:
        unique_alt_freqs = pd.concat([sp2_unique[["Sp1_Alt_Freq","Sp2_Alt_Freq"]], sp1_unique[["Sp1_Alt_Freq","Sp2_Alt_Freq"]]])
    
        all_alt_freqs = pd.concat([shared_alt_freqs,unique_alt_freqs])
    
    
    
    return all_alt_freqs


def sort_and_index_aligned_file(sorted_aligned_sample_file, aligned_sample_file):
    pysam.sort("-o", sorted_aligned_sample_file, aligned_sample_file)
    pysam.index(sorted_aligned_sample_file)


def call_variants(ref_file, transv_poly_file, sorted_aligned_sample_file):
    # bcftools mpileup -B -f fastafile -R bedfile bamfile | 
    # bcftools call -mV indels -A --ploidy 1 -o tr_samp_oar.vcf
    if DEBUG:
        print("Snp calling is starting ...")
        
    variant_calls = pysam.mpileup("-f", ref_file, "-B", "-l", transv_poly_file, sorted_aligned_sample_file)
    if DEBUG:
        print("Snp calling finished ...")
    return variant_calls


def convert_bam_to_bed(bam_file, bam_to_bed_file):
    # bedtools bamtobed -i bamfile > samp_oarbtb.bed 
    sample_in_bam = pybedtools.example_bedtool(bam_file)
    sample_in_bed = sample_in_bam.bam_to_bed()
    if DEBUG:
        write_to_file(sample_in_bed, bam_to_bed_file) 

    return sample_in_bed


def find_intersections(shared, transv, uniq_intersections_sample_file):
    from re import sub
    
    # bedtools intersect -a sha.oar.bed -b transvoar_samp.bed -wb 
    shared = pybedtools.BedTool.from_dataframe(shared)
    transv_poly = pybedtools.BedTool.from_dataframe(transv)
    intersects = shared.intersect(transv_poly, wb=True)
    intersects = intersects.to_dataframe()

    intersects = pd.concat([intersects["name"], intersects["blockCount"]], axis=1)

    # |awk '{gsub(/A|T|G|C/,"N",$2)}1' | awk '{gsub(/N,N/,"N",$2)}1'
    intersects["blockCount"] = intersects["blockCount"].apply(
                                        lambda item: sub("A|a|T|t|G|g|C|c|(N,N)|(N,N,N)", "N", item))

    # sort | uniq -c |awk '{ print $1,'\t',$2,'\t',$3 }'
    intersects = intersects.value_counts().reset_index(name='counts')
    intersects = intersects.sort_values(by=["name"])

    if DEBUG:
        print("\n\n >>> INTERSECTIONS <<<")
        print(intersects.head(10))

    if DEBUG:
        intersects.to_csv(uniq_intersections_sample_file, header=False, sep="\t", mode="w")
        
    return intersects


def insert_ref_and_alt_allele_numbers(intersects):
    intersects["id"] = range(1, len(intersects)+1)

    intersects_with_ref_allele = intersects[intersects["blockCount"] == '.']
    intersects_with_alt_allele = intersects[intersects["blockCount"] == 'N']

    intersects_with_ref_allele.insert( len(intersects_with_ref_allele.columns), "Ref", intersects_with_ref_allele["counts"])
    intersects_with_ref_allele.insert( len(intersects_with_ref_allele.columns), "Alt", [0]*len(intersects_with_ref_allele))

    intersects_with_alt_allele.insert( len(intersects_with_alt_allele.columns),"Alt", intersects_with_alt_allele["counts"])
    intersects_with_alt_allele.insert(len(intersects_with_alt_allele.columns), "Ref", [0]*len(intersects_with_alt_allele))

    all_intersects = pd.concat([intersects_with_ref_allele, intersects_with_alt_allele])
    all_intersects.insert(len(all_intersects.columns), "Total", all_intersects["Ref"] + all_intersects["Alt"])
    all_intersects = filter_columns(all_intersects, ["counts", "blockCount"])
    all_intersects = all_intersects.sort_values(by=["id"])
    
    if DEBUG:
        print("\n\n>>> ALL INTERSECTIONS <<<")
        print(all_intersects.head(5))
    return all_intersects




SHARED = 0
DEBUG = 0
if __name__ == "__main__":    

    flag_sp1 = 0
    flag_sp2 = 0
    flag_out = 0
    for arg in sys.argv:
        if flag_sp1 == 1:
            aligned_sample_to_sp1_file = arg
            flag_sp1 = 0

        elif flag_sp2 == 1:
            aligned_sample_to_sp2_file = arg
            flag_sp2 = 0  

        elif arg == "-d":
            DEBUG = True
            
        elif arg == "-shared":
            SHARED = True

        elif arg == "-sp1":
            flag_sp1 = 1

        elif arg == "-sp2":
            flag_sp2 = 1


        elif arg == "-out":
            flag_out = 1

        elif flag_out == 1:  
            flag_out = 0  
            out_file = str(arg)
            if not out_file.endswith(".txt"):
                out_file = out_file + ".txt"


            


    print("Debug Mode On\n" if DEBUG else "Debug Mode Off\n")


# samtools faidx sequence1.fasta
pysam.faidx(settings.sp1_ref_file)
pysam.faidx(settings.sp2_ref_file)


# === STEP 1 ===
# sort and index the sample's aligned sequence to the sheep's or goat's mtDNA reference
sort_and_index_aligned_file(sorted_aligned_sample_to_sp1_file, aligned_sample_to_sp1_file)
sort_and_index_aligned_file(sorted_aligned_sample_to_sp2_file, aligned_sample_to_sp2_file)
# ==============


# === STEP 2 ===
variant_calls_sp1 = call_variants(settings.sp1_ref_file, settings.transv_poly_sp1_file, sorted_aligned_sample_to_sp1_file)
variant_calls_sp2 = call_variants(settings.sp2_ref_file, settings.transv_poly_sp2_file, sorted_aligned_sample_to_sp2_file)
# === ====== ===


# === STEP 3 ===
variant_calls_sp1 = filter_variant_calls_from_pileup_format(variant_calls_sp1)
variant_calls_sp2 = filter_variant_calls_from_pileup_format(variant_calls_sp2)


# awk '{print $1,$2-1,$2,$4}'


variants_info_sp1 = get_variants_info(variant_calls_sp1, transv_sample_sp1_file)
if variants_info_sp1.size == 0:

    print("no variants present in sp1")
    if path.exists(data_out_dir) and not DEBUG:
        shutil.rmtree(data_out_dir)
    if path.exists(data_out_dir) and not DEBUG:
        shutil.rmtree(data_out_dir)
    sys.exit()
	


    
  
	
variants_info_sp2 = get_variants_info(variant_calls_sp2, transv_sample_sp2_file)
if variants_info_sp2.size == 0:
    print("no variants present in sp2")
    if path.exists(data_out_dir) and not DEBUG:
        shutil.rmtree(data_out_dir)
    if path.exists(data_out_dir) and not DEBUG:
        shutil.rmtree(data_out_dir)    
    sys.exit()
	
    
        
    

# === ====== ===


# === STEP 4 ===
sample_sp1_in_bed = convert_bam_to_bed(sorted_aligned_sample_to_sp1_file, sample_sp1_bamtobed_file)   
sample_sp2_in_bed = convert_bam_to_bed(sorted_aligned_sample_to_sp2_file, sample_sp2_bamtobed_file)
# === ====== ===


# === STEP 5 ===
# sharead.R
shared_sp2, shared_sp1 = find_shared_reads(sample_sp2_in_bed, sample_sp1_in_bed)
if shared_sp1.size == 0 or shared_sp2.size == 0:
    print("no shared reads present")
    if path.exists(data_out_dir) and not DEBUG:
        shutil.rmtree(data_out_dir)
    if path.exists(data_out_dir) and not DEBUG:
        shutil.rmtree(data_out_dir)
    sys.exit()
    # === ====== ===



# === STEP 6 ===
intersect_sp1 = find_intersections(shared_sp1, variants_info_sp1, uniq_intersections_sample_sp1_file)
if intersect_sp1.size == 0:
	
    if path.exists(data_out_dir) and not DEBUG:
        shutil.rmtree(data_out_dir)
    if path.exists(data_out_dir) and not DEBUG:
        shutil.rmtree(data_out_dir)     
    print("no intersections between target positions and shared reads")
    sys.exit()
   
	
intersect_sp2 = find_intersections(shared_sp2, variants_info_sp2, uniq_intersections_sample_sp2_file)
if intersect_sp2.size == 0:
    if path.exists(data_out_dir) and not DEBUG:
        shutil.rmtree(data_out_dir)
    if path.exists(data_out_dir) and not DEBUG:
        shutil.rmtree(data_out_dir)     
    print("no intersections between target positions and shared reads")
    sys.exit()


# === ====== ===


# === STEP 7 ===
# create table of alternative, reference and total allele numbers
all_intersects_sp1 = insert_ref_and_alt_allele_numbers(intersect_sp1)
all_intersects_sp2 = insert_ref_and_alt_allele_numbers(intersect_sp2)
# === ====== ===


# === STEP 8 === 
all_alt_freqs = find_alt_freqs(all_intersects_sp2, all_intersects_sp1)
sp2_alt_freq_with_0_or_1 = all_alt_freqs.query('Sp2_Alt_Freq == 0 | Sp1_Alt_Freq == 0')



total_read_numbers = len(sp2_alt_freq_with_0_or_1)
sp2_read_numbers = len(sp2_alt_freq_with_0_or_1.query('Sp2_Alt_Freq == 0'))
sp1_read_numbers = len(sp2_alt_freq_with_0_or_1.query('Sp1_Alt_Freq == 0'))

print()
print("# of total reads: ", total_read_numbers)
print("# of Sp1 reads: ", sp1_read_numbers)
print("# of Sp2 reads: ", sp2_read_numbers)

from scipy.stats import binom_test
binom_test_result = binom_test(sp1_read_numbers, total_read_numbers, p=0.5, alternative='two-sided')

read_numbers_with_binom_result = pd.DataFrame({
                                        "Total_reads": [total_read_numbers],
                                        "Sp1_reads": [sp1_read_numbers],
                                        "Sp2_reads": [sp2_read_numbers],
                                        "p_value": [binom_test_result]
                                    })
print("\b\n>>> RESULTS <<<")
print(read_numbers_with_binom_result)
# === ====== ===




if path.exists(data_out_dir) and not DEBUG:
    shutil.rmtree(data_out_dir)


#read_numbers_with_binom_result.to_excel(out_file)
read_numbers_with_binom_result.to_csv(out_file, sep="\t", index = False, header=True)

