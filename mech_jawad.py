import pandas as pd
import numpy as np
import os

usecols = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Otherinfo10", "Otherinfo11"]
dtype = {"Chr":"category", "Start":np.int64, "End":np.int64, "Ref":"category", "Alt":"category", "Func.refGene":"category", "Otherinfo10":"category", "Otherinfo11":"string"}

def add_dp_af(df):
    df = pd.eval("DP = df.Otherinfo11.str.split(';',3).str[2].str.split('DP=').str[1].astype('int32')", target=df)
    df = pd.eval("AF = df.Otherinfo11.str.split(';',5).str[4].str.split('AF=').str[1].astype('float')", target=df)
    return df

# works on Chr, Ref, Alt, Func.refGene, Otherinfo10
# for Func.refGene, enter 'FuncrefGene' instead
def category_filter(df, inplace=0, export=False, **kwargs):
    for key, value in kwargs.items():
        if key == "FuncrefGene":
            df = df[df["Func.refGene"] == value]
        else:
            df = df[df[key] == value]
    if inplace==1:
        globals()['dfa'] = df
    elif inplace==2:
        globals()['dfb'] = df
    if export:
        df.to_csv('cat_filter_dataframe.csv', index=False, header=True)
    return df

# works on DP,AF, assumes greater than i.e 'DP = 15' means keep DP vals > 15
def num_filter(df, inplace=0, export=False, **kwargs):
    for key, value in kwargs.items():
        df.loc[:,key].where(df.loc[:,key] > value,inplace=True)
    df = df.dropna()
    if inplace==1:
        globals()['dfa'] = df
    elif inplace==2:
        globals()['dfb'] = df
    if export:
        df.to_csv('cat_filter_dataframe.csv', index=False, header=True)
    return df

def subtraction(df1, df2, output_filename, export=True):
    uqid_df1 = df1[df1.columns[1:5]].apply(lambda x: ''.join(x.astype(str)),axis=1).squeeze()
    uqid_df2 = df2[df2.columns[1:5]].apply(lambda x: ''.join(x.astype(str)),axis=1).squeeze()  
    uqid = pd.concat([uqid_df1, uqid_df2], ignore_index=True)
    uqid_df1_len, uqid_df2_len= len(uqid_df1.index), len(uqid_df2.index)
    idx = np.append(np.arange(0, uqid_df1_len), np.arange(0, uqid_df2_len))
    uqid_df = pd.DataFrame({'idx':idx}, index=uqid)
    uqid_df.reset_index(inplace=True)
    uqid_df = uqid_df.sort_index().groupby('index').filter(lambda x: len(x) == 1)
    idx = uqid_df['idx'].to_numpy()
    mid_val = idx[len(idx)//5:].argmin() + len(idx)//5
    inv_arr = lambda max_val,idx_arr: np.array(sorted(set(range(0, max_val)).difference(idx_arr)))
    df1=df1.drop(df1.index[inv_arr(uqid_df1_len, idx[:mid_val])])
    df2=df2.drop(df2.index[inv_arr(uqid_df2_len, idx[mid_val:])])
    if export:
        df1.to_csv(output_filename + 'p.csv', index=False, header=True)
        df2.to_csv(output_filename + 'v.csv', index=False, header=True)
    return df1

# p - v, filter the variant and parent differently: rerun the algorithm v - p


# FOR BATCH PUT THESE FILES IN HERE
filenames = {'22241':['22334','22243','22244','22245'],"22241":['22334','22243','22244','22245']}
chrom_list = ['1','2', '3']
path = '../annovar-full/chr'

# '../annovar-full/chr(CHR#)_(SAMPLE ID).csv'
# '../annovar-full/chrY_22243.csv'

def main(filename_dict, parent, chromosomes):
    for chrom in chromosomes:
        for variant in filename_dict[parent]:
            parent_path
            # what is subtracting - (category criteria)
            # what is being subtracted from -(category criteria)
            parent_df = pd.read_csv(path + chrom + "_" + parent + '.csv', usecols=usecols, dtype=dtype, sep='\t')
            parent_df = category_filter(parent_df, FuncrefGene='exonic')
            
            subtraction(parent_df,
                        pd.read_csv(path + chrom + "_" + variant + '.csv', usecols=usecols, dtype=dtype, sep='\t'),
                        parent + "_" + variant + "_chr" + chrom)
            # subtraction(parent, variant)
            # 
            # subtraction(variant, parent)
    return None

# TODO: collate chromosomal outputs
# TODO: OR for inputs
# main(filenames, '22241', chrom_list)

