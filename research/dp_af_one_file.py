# FIRST FILTER IS OTHERINFO 10 = PASS AND FOR CHILDREN OTHERINFO 11 = DP > 15 and AF > .15 AND STRING FIRST 5 FIELDS
import pandas as pd
import os


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
            # REMOVE THE ROW FROM EXCEL
    # cross with first pass from pass_filter
    # for x in range(len(first_pass)): if first_pass value(x) == pass and second_pass value(x) == pass: pass else: fail
    # second pass value is rows because rows is var of this function (dp_af_filter) not first_pass function
    for x in range(col_len):
        if first_pass[str(x+2)] == "PASS" and rows[str(x+2)] == "PASS":
            pass
        else:
            rows[str(x+2)] = "REMOVED"
    return rows


chr_1 = pass_filter(pd.read_excel(r"C:\Users\kazij\Documents\Cisplatin_Resistance_Variant_Analysis\annovar-full\anno1.fixed.chr1.genes.VE_22249.hg38_multianno.txt"))
