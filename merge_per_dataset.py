import re
import pandas as pd
import glob
from statistics import NormalDist
import numpy as np
import itertools

class AlignedDiffResultsCollector():
    def __init__(self, dir_with_aligned_diffresults):
        self._dir_with_aligned_diffresults = dir_with_aligned_diffresults
        self._list_of_files = None
        self.gene2info = {}

        self._define_list_of_files()
        self._iterate_through_files()

    
    def _define_list_of_files(self):
        self._list_of_files = glob.glob(f"{self._dir_with_aligned_diffresults}/*")

    def _iterate_through_files(self):
        for file in self._list_of_files:
            collector = ResultsFileCollector(file)
            collector.update_gene2info_w_results_df(self.gene2info)


class ResultsFileCollector():
    def __init__(self, results_file):
        self._results_file = results_file
        self._nd = NormalDist()

        self._condpair = None
        self._results_df = None

        self._define_condpair()
        self._define_results_df()
        
    def update_gene2info_w_results_df(self, gene2info):
        for _, row in self._results_df.iterrows():
            protein = row["protein"]
            pval = row["p-val"]
            log2fc = row['log2FC']
            z_val = self._convert_p_fc_to_zval(pval, log2fc, self._nd)
            self._update_geneinfo(protein, z_val, gene2info)

    def _define_condpair(self):
        condpairname = self._results_file.split("/")[-1].replace(".tsv", "")
        split_condpairname = condpairname.split("_VS_")
        self._condpair = {split_condpairname[0], split_condpairname[1]}

    def _define_results_df(self):
        self._results_df = pd.read_csv(self._results_file, sep = "\t")

    
    @staticmethod
    def _convert_p_fc_to_zval(p_val, log2fc, nd):
        p_val = max(p_val, 1e-9)
        z_abs = abs(nd.inv_cdf(p_val/2))
        direction = np.sign(log2fc)
        return z_abs*direction
    
    def _update_geneinfo(self, gene, z_val, gene2info):
        geneinfo = gene2info.get(gene, GeneInfo(gene))
        geneinfo.zvals.append(z_val)
        geneinfo.condpairs.append(self._condpair)
        gene2info[gene] = geneinfo


class GeneNameUpdater():
    def __init__(self, protid2gene_file):
        self._protid2gene_file = protid2gene_file
        self._protid2gene = None

        self._define_protid2gene()

    def _define_protid2gene(self):
        protid2gene_df = pd.read_csv(self._protid2gene_file, sep = "\t", header=None)
        self._protid2gene = dict(zip(protid2gene_df[0], protid2gene_df[1]))

    def update_collection_of_geneinfos(self, collection_of_geneinfos):
        for geneinfo in collection_of_geneinfos:
            self._change_protname_to_reference_name(geneinfo)

    def _change_protname_to_reference_name(self, geneinfo):
        if "sp|" in geneinfo.name:
            geneinfo.name = geneinfo.name.split("|")[1]
        geneinfo.name = self._protid2gene.get(geneinfo.name)


class GeneInfo():
    def __init__(self, name):
        self.name = name
        self.condpairs = []
        self.zvals = []


class GeneToZvalMerger():
    def __init__(self, gene2info):
        self._gene2info = gene2info
        self.gene2mergedz_df = None
        self._define_gene2mergedz_df()
    
    def _define_gene2mergedz_df(self):
        genes = [x.name for x in self._gene2info.values()]
        zvals = [self._merge_zvals(x) for x in self._gene2info.values()]
        self.gene2mergedz_df = pd.DataFrame({'protein' : genes, 'z-value':zvals})
    
    # def _merge_zvals(self, geneinfo):

    #     num_zs = len(zvals)
    #     if num_zs ==1:
    #         return zvals[0]
    #     z_summed = sum(zvals)
    #     num_elems = 0.5*(np.sqrt(8*num_zs +1)+1)
    #     sigma = np.sqrt(0.5*num_elems*((num_elems-1)**2))
    #     z_merged = z_summed/sigma
    #     return z_merged

    def _merge_zvals(self, geneinfo):
        z_summed = sum(geneinfo.zvals)
        variance = self.calc_summed_covariances(geneinfo.condpairs)
        merged_z = z_summed/np.sqrt(variance)
        merged_cut_z = np.sign(merged_z)*min(7.3487, abs(merged_z))
        return merged_cut_z

    @staticmethod
    def calc_summed_covariances(list_of_pairs):
        num_overlaps = 0
        for pair_of_pairs in itertools.combinations(list_of_pairs, 2):
            elems1 = pair_of_pairs[0]
            elems2 = pair_of_pairs[1]
            num_intersects = len(elems1.intersection(elems2))
            if num_intersects ==1:
                num_overlaps+=1
            elif num_intersects == 0:
                continue
            else:
                raise Exception("two overlaps should not exist!")
        result = num_overlaps+len(list_of_pairs)
        return result