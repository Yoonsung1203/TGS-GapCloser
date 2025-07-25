#%%
import gzip,re,json,os,sys,argparse

import pandas as pd


def run(path_paf, path_save, idy, align, hang_wrong, align_ratio):
    with open(path_paf, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            list_values = line.strip('\n').split('\t')

            list_columns = [ 
                "qryname",
                "qrylen",
                "qry_alignstart",
                "qry_alignend",
                "strand",
                "refname",
                "reflen",
                "ref_alignstart",
                "ref_alignend",
                "match",
                "align",
                "quality"
            ]

            list_columns_int = [
                "qrylen",
                "qry_alignstart",
                "qry_alignend",
                "reflen",
                "ref_alignstart",
                "ref_alignend",
                "match",
                "align",
                "quality"
            ]

            dict_row = dict(zip(
                list_columns, list_values[:12]
            ))

            for col_int in list_columns_int:
                dict_row[col_int] = int(dict_row[col_int])

            hang_size = get_wrong_hang_length(dict_row["strand"], dict_row["reflen"], dict_row["qrylen"], dict_row["ref_alignstart"], dict_row["ref_alignend"], dict_row["qry_alignstart"], dict_row["qry_alignend"])
            identity = dict_row["match"] / dict_row["align"]

            is_identity = identity > idy
            is_wronghang = hang_size < hang_wrong
            is_align = dict_row["align"] > align
            is_alignratio = (dict_row["align"] / (dict_row["align"] + hang_size)) > align_ratio

            if is_identity and is_wronghang and is_align and is_alignratio:
                fw.write(line)
    

def read_paf(path_paf):
    table_paf = pd.read_csv(path_paf, sep = '\t', names = list(range(1, 13)), usecols = list(range(1, 13)))
    table_paf = table_paf.rename(
        columns = {
            1:"qryname",
            2:"qrylen",
            3:"qry_alignstart",
            4:"qry_alignend",
            5:"strand",
            6:"refname",
            7:"reflen",
            8:"ref_alignstart",
            9:"ref_alignend",
            10:"match",
            11:"align",
            12:"quality"
        }
    )
    return table_paf

def get_wrong_hang_length(strand, reflen, qrylen, ref_alignstart, ref_alignend, qry_alignstart, qry_alignend):
    is_reverse_complement = strand == '-'

    ref_left_hang = ref_alignstart
    ref_right_hang = reflen - ref_alignend

    if is_reverse_complement:
        qry_left_hang = qrylen - qry_alignend
        qry_right_hang = qry_alignstart
    else:
        qry_left_hang = qry_alignstart
        qry_right_hang = qrylen - qry_alignend
    
    left_hang = min(ref_left_hang, qry_left_hang)
    right_hang = min(ref_right_hang, qry_right_hang)

    return left_hang + right_hang


#%%
if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--paf")
    argument_parser.add_argument("--out")
    argument_parser.add_argument("--idy", type = float)
    argument_parser.add_argument("--align", type = float)
    argument_parser.add_argument("--hang_wrong", type = float)
    argument_parser.add_argument("--align_ratio", type = float)
    
    args = argument_parser.parse_args()
    
    path_paf = args.paf
    path_save = args.out
    idy = args.idy
    align = args.align
    hang_wrong = args.hang_wrong
    align_ratio = args.align_ratio
    
    run(path_paf, path_save, idy, align, hang_wrong, align_ratio)