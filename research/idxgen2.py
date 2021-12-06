#!/usr/bin/env python
# coding: utf-8

# In[4]:


# All together for total time count
import pandas as pd
from bisect import bisect_left
import numpy as np
get_ipython().run_line_magic('pylab', 'inline')
from numba import jit
import pyarrow.feather as feather
import _pickle as cpickle
import pylzma as xz

chromosome_dict = dict(zip([x+1 for x in range(23)], [str(x+1) for x in range(23)]))
chromosome_dict[24], chromosome_dict[25] = "X", "Y"
sample_ids = [22241, 22334, 22243, 22244, 22245, 22246, 22247, 22248, 22249, 22250, 22251, 22252, 22253, 22254, 22255]
indexes = [0,1089538,1868287,2481623,2892120,3368329,3905767,4426173,4811381,5265718,5687521,6341258,6907171,7094109,7468963,7847065,8375551,9027612,9188821,9951403,10223830,10345825,10605338,10868333,10872903]
# inp_df = pd.read_csv("inp.txt", sep='\t', usecols=[0], nrows=1, dtype={"#Chr":'string'})

def BinarySearch(arr, low, high, x):
    if high >= low:
        mid = (high + low) // 2
        if arr[mid] == x:
            return mid
        elif arr[mid] > x:
            return BinarySearch(arr, low, mid - 1, x)
        else:
            return BinarySearch(arr, mid + 1, high, x)
    else:
        return -1

#     i = bisect_left(a, x)
#     if i != len(a) and a[i] == x:
#         return i
#     else:
#         return -1


# In[2]:


# input list of chr, output list of indexes
def get_idxs(db_chrs):
    # [1:4] because beginning of chromosome 2 index to beginning of chromosome Y, 1 is 0
    chrs = list(chromosome_dict.values())[1:]
    # db_chrs = ["Chr1", "Chr1", "Chr1", "Chr2", "Chr2", "Chr2", "Chr3","Chr3","Chr3", "Chr4", "Chr4"]
    db_chr_idx = [0]
    for x in chrs:
        db_chr_idx.append(np.searchsorted(db_chrs, "Chr" + x))
    return db_chr_idx


# In[4]:


# cp = new_df.to_pickle("exac03.pkl", compression=None) 63500521
# new_df = pd.read_csv("databases/gnomad30.txt", sep='\t', nrows=5,usecols=[1,3,5,6,7,8,9,10,11,12,13,14,15,16], dtype={"Start":'int',"Ref":"string","AF":"float","AF_raw":"float","AF_male":"float","AF_female":"float","AF_afr":"float","AF_ami":"float","AF_amr":"float","AF_asj":"float","AF_eas":"float","AF_fin":"float","AF_nfe":"float","AF_oth":"float","AF_sas":"float"})
# new_df
# cpdf = pd.read_pickle('exac03.pkl')
# oth_df = pd.read_csv("databases/gnomad30.txt", sep='\t', nrows=63500521,usecols=[1,3,6], dtype={"Start":'int',"Ref":"string","AF_raw":"string"})
# oth_df


# In[ ]:


# oth_df.to_pickle("gnomad30_chr1_af.pkl")
with open('gnomad30_chr1_af.pkl', 'rb') as gnomad1:
    data = cpickle.load(gnomad1)
with lzma.open()
cp = cpickle.load('gnomad30_chr1_af.pkl')
with lzma.open("lmza_test.xz", "wb") as f:
     pickle.dump(data, f)


# In[5]:


lzma.open()


# In[6]:


# list of database names as parameter, sidx eidx temp parameters
def get_db(sidx, eidx):
    return cpdf['Start'].values[sidx:eidx]

# FAKE INDEXES FOR TESTING PURPOSES RN
indexes = [0,1089538,1868287,2481623,2892120,3368329,3905767,4426173,4811381,5265718,5687521,6341258,6907171,7094109,7468963,7847065,8375551,9027612,9188821,9951403,10223830,10345825,10605338,10868333,10872903, 10872903]


# In[64]:


# list of indexes to get the db values ie AF_raw AF_amr etc
def get_db_info(indexes=None):
    return cp['AF_raw']


# In[39]:


def main(inp_start, out_start):    
    # locations = tuple(BinarySearch(out_start, 0, len(out_start) - 1, inp) for inp in inp_start)
    binsearch_results = []
    for x in range(len(inp_start)):
        binsearch_results.append((x, BinarySearch(out_start, 0, len(out_start) - 1, x)))
    # val_loc = zip(locations, [x for x in range(len(inp_start))])
    return binsearch_results


# In[48]:


# cp["Start"].values
sid = "22249"
chromosome_inp = "1"
inp_df = pd.read_csv("annovar-full/chr%s_%s.txt" % (chromosome_inp, sid), sep='\t', dtype={"Start":int})

with open("annovar-full/chr1_%s.txt" %sid) as inp:
    inp_chr = ''.join(inp.readlines()[1:2][0][3:5].split())

index_val = list(chromosome_dict.values()).index(inp_chr)
start_idx = indexes[index_val]
end_idx = indexes[index_val+1]

df_start = inp_df["Start"].values


# In[49]:


final_output = main(df_start, cp["Start"].values) # get_db(start_idx, end_idx))
# del final_output[-1]
final_output


# In[65]:


# df_start[final_output.get(1)]
# len 15798
# GET THE ROWS FROM THE TUPLE OUTPUT
# creates column af_raw and fills with .

# VECTORIZE THIS FOR LOOP ITS TOO SLOW BOOTLENECK IN CODE
AF = get_db_info()
inp_df["AF_raw"] = ['.'] * len(inp_df)
for final in range(len(final_output)):
    if final_output[final][1] > -1:
        inp_df["AF_raw"].iloc[final] = AF[final_output[final][0]]
# numpy.where(a > -1, a, get value from index)
# numpy.where(replace all == -1 with .)


# In[ ]:


inp_df.to_csv('test_output.csv')


# In[27]:


# from 0 to last line in the file for gnomad30genome
# ['0', '63500523'],  first chromosome already done
for x in range(len(gg_idx)):
    for i in range(2):
        gg_idx[x][i] = int(gg_idx[x][i])
gg_idx


# In[33]:


gg_idx = [[458633875, 496794134], [496794134, 534168131], [534168131, 571116550], [571116550, 597829579], [597829579, 623024981], [623024981, 646931121], [646931121, 673364087], [673364087, 697350022], [697350022, 718159091], [718159091, 737255085], [737255085, 754960502], [754960502, 766773369], [766773369, 779364375], [779364375, 811571117], [811571117, 813080256]]
for idx in range(len(gg_idx)):
    gg_df = pd.read_csv("databases/gnomad30.txt", sep='\t', skiprows=gg_idx[idx][0], nrows=gg_idx[idx][1]-gg_idx[idx][0],usecols=[1,3,6], dtype={"Start":'int',"Ref":"string","AF_raw":"string"})
    gg_df.to_pickle("gnomad30_chr%s_af.pkl" % str(idx+10))
    print("Dont with pickle gnomad30_chr%s_af" % str(idx+10))
    # 1000 17 files run through the person

