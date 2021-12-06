# FIRST FILTER IS OTHERINFO 10 = PASS AND FOR CHILDREN OTHERINFO 11 = DP > 15 and AF > .15 AND STRING FIRST 5 FIELDS
import pandas as pd
import os


def generate_file_names(sample_identification):
    file_names = []
    chromosome_names = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', "X", "Y"]
    for x in range(23):
        file_names.append("anno1.fixed.chr" + chromosome_names[x] + ".genes." + sample_identification + ".hg38_multianno.txt")
    return file_names


def pass_filter(file_df):
    rows = {}
    for x in range(len(file_df['Otherinfo10'].copy().index)):
        if file_df.loc[x, 'Otherinfo10'] == "PASS":
            rows[str(x + 2)] = "PASS"
        else:
            rows[str(x + 2)] = "REMOVED"
    return rows


def dp_af_filter(file_df, first_pass):
    # takes file_df, filters, then crosses with first pass from pass_filter
    rows = {}
    # for x in range(len(file_df['Otherinfo11'].copy().index)):
    # DP > 15 AF > .15
    # for each row of otherinfo11, separate by semicolon to find DP and AF
    col_len = len(file_df['Otherinfo11'].copy().index)
    # print(col_len)
    for x in range(len(file_df['Otherinfo11'].copy().index)):
        dp = False
        af = False
        for i in file_df.loc[x, 'Otherinfo11'].split(";"):
            if i[:3] == "DP=" and int(i[3:]) > 15:
                dp = True
            if i[:3] == "AF=" and float(i[3:]) > 0.15:
                af = True
        if dp and af:
            rows[str(x + 2)] = "PASS"
        else:
            rows[str(x + 2)] = "REMOVED"
    # cross with first pass from pass_filter
    # for x in range(len(first_pass)): if first_pass value(x) == pass and second_pass value(x) == pass: pass else: fail
    # second pass value is rows because rows is var of this function (dp_af_filter) not first_pass function
    for x in range(col_len):
        if first_pass[str(x+2)] == "PASS" and rows[str(x+2)] == "PASS":
            pass
        else:
            rows[str(x+2)] = "REMOVED"
    return rows


def combine(df, first_pass):
    # first_pass is the PASS REMOVED results to cross with UD, only need UDs that PASS
    start = df["Start"].tolist()
    start = [str(x) for x in start]
    end = df["End"].tolist()
    end = [str(x) for x in end]
    ref = df["Ref"].tolist()
    alt = df["Alt"].tolist()
    generated_ud = [j + k + l + m for j, k, l, m in zip(start, end, ref, alt)]
    # list version of first_pass
    pass_list = list(first_pass.values())
    return dict(zip(generated_ud, pass_list))

"""
def disjoint(combined_dict):
    disjoints = []
    for x in range(len(combined_dict):
        if x not in list(combined_dict.keys()):
            disjoints.append(x)
    return disjoints"""


"""def paternal_test():
    # will determine whether two sample ids and chromosome are cross-able (same chromosome and one is respective parent)
    return 0"""


# key number value is ud
# parent variant arguments are just placeholders not supposed to actually be parent and variant
def compare(parent, variant, parent_filter_output, variant_filter_output):

    disjoint = {}
    # disjoint parent - variant and variant - parent ONLY if they PASS
    # return list of every variant disjoint - how many variant uds not in parent
    for x in range(len(parent.keys())):
        if list(variant.values())[x] in list(parent.values()):
            # prints out what row is in variant
            pass
        else:
            disjoint[x] = list(variant.keys())[x]
    return disjoint


parent_ids = ['VE_22241', 'VE_22246', 'VE_22251']
folder_directory = r'C:\Users\kazij\Documents\Cisplatin_Resistance_Variant_Analysis\annovar-test'
"""
FOR EACH FILE IN ANNOVAR FOLDER
1) IF PARENT SAMPLE ID ONLY RUN PASS_FILTER NOT DP_FILTER AF_FILTER
2) IF VARIANT SAMPLE ID RUN PASS_FILTER AND DP_FILTER AND AF_FILTER
3) COMBINE THE FIRST FIVE FIELD STRINGS FOR REMAINING ROWS
"""
# MAIN FUNCTION
# dictionary of list type values, the keys will be crossed if they are parent and child, UD = Unique Identifier
# pre_filter is another way of saying the first pass
# file info is list of dictionaries, each dictionary says UD and whether that row PASS or REMOVED
# file_info = []
"""for filename in os.listdir(folder_directory):
    file_directory = os.path.join(folder_directory, filename)
    chromosome = filename.split(".")[0]
    sample_id = filename.split(".")[1]
    # if current file is a parent sample 1), else it's a variant sample 2)
    # print(chromosome, sample_id)

    # [:8] because cannot include "_test" part of filename
    if sample_id[:8] in parent_ids:
        pre_filter = pass_filter(pd.read_excel(r'' + file_directory))

    else:
        other_info_10 = pass_filter(pd.read_excel(r'' + file_directory))
        # print(other_info_10)
        pre_filter = dp_af_filter(pd.read_excel(r'' + file_directory), other_info_10)"""

# print(pre_filter)
# COMBINE(FILE) --> CREATE A UNIQUE VARIANT IDENTIFIER FOR EACH ROW AND RETURN IT AS A LIST
# file_info.append(combine(pd.read_excel(file_directory), pre_filter))

# compare tells you how many values in parent are in variant
# print(compare(file_info[0], file_info[1]))
# diff = [compare(file_info[1], file_info[0])]
# for each variant file, run it against the parent, append to diff
variant_list = generate_file_names("VE_22252")
parent_list = generate_file_names("VE_22251")
for x in range(len(variant_list)):
    folder_directory = "C:\\Users\\kazij\\Documents\\Cisplatin_Resistance_Variant_Analysis\\hn31_test\\source_files\\"

    file_name = folder_directory + variant_list[x]
    parent_file_name = folder_directory + parent_list[x]

    file = pd.read_csv(r"" + file_name, delimiter="\t")
    parent_file = pd.read_csv(r"" + parent_file_name, delimiter="\t")
    # run other 11 and 10 to variant file
    variant_filtered = dp_af_filter(file, pass_filter(file))
    parent_filtered = pass_filter(parent_file)
    # print(variant_filtered, parent_filtered)
    print()
"""    for y in variant_filtered:
        if variant_filtered[y] == "REMOVED":
            file.drop(file.index[int(y)], axis=0, inplace=True)
    for y in parent_filtered:
        if parent_filtered[y] == "REMOVED":
            parent_file.drop(parent_file.index[int(y)], axis=0, inplace=True)"""
    # for each row in parent and variant filtered if that rows key is removed then remove the row number index (key)
# diff will be list of all the uds that are disjoint of variant from parent
"""for sample_ud in ud:
    for row in sample_ud:
        # if that row passed and is in the other list, add the row number to parent - variant
        print(row)"""
# for every row in file combine(): combine returns list of unique variant identifier for each row in list
# variant_df = pd.read_excel(r'C:\Users\kazij\Documents\Cisplatin_Resistance_Variant_Analysis\annovar-test\chrY.VE_22334_test.xlsx')
# print(pass_filter(variant_df))
# UNIQUE IDENTIFIER IS START + END + REF + ALT str, ie. IDs ['27870242787024GT', '28416422841642TC', '28418792841879TA']
# Some REFs and ALTs are not single letter ('-', 'CGGC', etc.)
