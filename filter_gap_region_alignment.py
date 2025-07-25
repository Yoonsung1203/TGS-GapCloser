#%%
import gzip,re,json,os,sys,argparse

import pandas as pd


def run(path_paf, path_scaff_info, path_save):
    dict_scaff_info = read_scaff_info(path_scaff_info)
    list_contigs_between_gap = get_contigs_bewteen_gap(dict_scaff_info)
    dict_qry_to_gapcontigs = get_gap_associated_aligned_contignames_per_readname(path_paf, list_contigs_between_gap)
    
    list_qrynames_assigned_multiple_contigs = list(filter(lambda name: len(dict_qry_to_gapcontigs[name]) > 1, dict_qry_to_gapcontigs.keys()))
    write_paf_with_desired_qrynames(path_paf, path_save, list_qrynames_assigned_multiple_contigs, list_contigs_between_gap)

def read_scaff_info(path_scaff_info):
    dict_scaff_info = dict()
    
    with open(path_scaff_info, 'r') as fr:
        for line in fr:
            if line.startswith('>'): continue
            values = line.strip('\n').split('\t')
            
            contigname = values[0]
            scaffname = values[6]
            contiglen = int(values[3])
            contigind = int(values[5])

            if dict_scaff_info.get(scaffname) == None:
                dict_scaff_info[scaffname] = dict()
            assert dict_scaff_info[scaffname].get(contigind) == None
            dict_scaff_info[scaffname][contigind] = {
                "name" : contigname,
                "length" : contiglen
            }
    return dict_scaff_info

def get_contigs_bewteen_gap(dict_scaff_info):
    list_contigs = list()

    for scaffname, dict_contigs in dict_scaff_info.items():
        list_contig_order = list(dict_contigs.keys())
        if len(list_contig_order) > 1:
            list_contigs.extend(
                list(map(lambda dict_info: dict_info["name"], dict_contigs.values()))
            )
    return list_contigs

def get_gap_associated_aligned_contignames_per_readname(path_paf, list_contigs_between_gap):
    contignames_search = set(list_contigs_between_gap)
    
    dict_qry_to_contigs = dict()
    
    with open(path_paf, 'r') as fr:
        for line in fr:
            list_values = line.strip('\n').split('\t')
            qryname = list_values[0]
            refname = list_values[5]
            
            if refname in contignames_search:
                if dict_qry_to_contigs.get(qryname) == None:
                    dict_qry_to_contigs[qryname] = set()
                dict_qry_to_contigs[qryname].add(refname)
    return dict_qry_to_contigs

def write_paf_with_desired_qrynames(path_paf, path_save, list_qrynames, list_contigs_between_gap):
    qrynames_search = set(list_qrynames)
    contignames_search = set(list_contigs_between_gap)
    
    with open(path_paf, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            list_values = line.strip('\n').split('\t')
            qryname = list_values[0]
            refname = list_values[5]
            
            if qryname in qrynames_search and refname in contignames_search:
                fw.write(line)
                
            
            
            
# %%
if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--paf")
    argument_parser.add_argument("--out")
    argument_parser.add_argument("--scaff_info")
    
    args = argument_parser.parse_args()
    
    path_paf = args.paf
    path_save = args.out
    path_scaff_info = args.scaff_info
    
    run(path_paf, path_scaff_info, path_save)
