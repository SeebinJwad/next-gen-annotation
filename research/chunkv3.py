import pandas as pd
import numpy as np
from numba import njit
import os
%pylab inline

test_inp = pd.read_csv("first50/gnomad30_genome_first_50.txt", sep=',')
test_inp_start = test_inp["Start"].to_numpy()


test_inp = pd.read_csv("first50/gnomad30_genome_first_50.txt", sep=',')
test_inp_start = test_inp["Start"].to_numpy()

test_out = [-1] * len(test_inp)

chunksize = 10 ** 6
dtypes={
    "Chr#": "string",
    "Start": "int",
    "End": "int",
    "Ref": "string",
    "Alt": "string",
    "AF": "float"
}
# 63_500_523
curr_inp_idx = 0
for chunk in pd.read_csv("databases/gnomad/gnomad30.txt", 
                         dtype=dtypes, index_col=0, sep='\t', 
                         usecols=[1,2,3,4,5], na_filter=True,
                         na_values=".", chunksize=chunksize):
    # test_inp end value not start value ideally
    for idx, val in enumerate(test_inp_start):
        try:
            value = chunk.loc[val]
        except (KeyError) as error:
            continue
        else:
            if type(value) == pd.core.frame.DataFrame:
                # test_out[idx] = 
                print(checks := value.values[:-1]) # value.iloc[0][-1]
                # print(value)
                for check in checks:
                    if check[:-1] == test_inp.iloc[idx][1:4]:
                        test_out[idx] = checks[-1]
            else:
                # test_out[idx] = value[-1]
                print(checks := value.tolist())
                # gnomad line without the AF
                if checks[:-1] == test_inp.iloc[idx][1:4]:
                    test_out[idx] = checks[-1]
    break