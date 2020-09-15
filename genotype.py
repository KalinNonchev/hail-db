import hail as hl
import glob
import pandas as pd
from tqdm import tqdm
import numpy as np
import sample_anno as anno

class GenotypeHAIL():
    
    def __init__(self):
        hl.init(quiet=True)
    
    def get_genotype(self, var_df, samples):
        return self._get(var_df, samples, 'GT')
    
    def get_gtquality(self, var_df, samples):
        return self._get(var_df, samples, 'GQ')
    
    def get_sequencingdepth(self, var_df, samples):
        return self._get(var_df, samples, 'DP')

    def _sum(self, obj):
        if obj is None:
            return -1
        else:
            return sum(obj)

    def _get_Tables_paths(self, samples):
        sample_anno = anno.SampleAnnotation()
        paths = sample_anno.get_vcfs_for_samples(samples)
        mt_paths = [f"db/Tables/{(path.split('/')[-1]).split('.')[0]}.ht" for path in paths]
        return mt_paths
    
    def _get(self, var_df, samples, field):
        # array with samples in HAIL
        if type(samples) is str:
            samples = [samples]

        # create table with vars in HAIL
        ht = hl.Table.from_pandas(var_df)
                    
        # create table with samples
        df = pd.DataFrame({
            's': samples
        })
        ht_samples = hl.Table.from_pandas(df)
                
        ht = ht.join(ht_samples)
        ht = ht.annotate(pos=hl.int32(ht.pos))
        ht = ht.add_index()
        ht = ht.key_by(locus=hl.struct(contig=ht.chrom, position=ht.pos), alleles=hl.array([ht.ref, ht.alt]), s=ht.s)
        
        # all variants per sample
        res_table = None
        
        ht_paths = self._get_Tables_paths(samples)
        # iterate through ht_vcfs with samples
        for ht_path in ht_paths:
            ht_vcf = hl.read_table(ht_path)
            
            ht_n = ht.join(ht_vcf, how='left')
            
            if res_table is None:
                res_table = ht_n
                res_table = res_table.checkpoint('db/checkpoint/ht1.ht', overwrite=True)
            else:
                res_table = res_table.union(ht_n)


        # all variants per sample
        res_table = res_table.annotate(GT = hl.coalesce(res_table.GT, 0))
        res_table = res_table.annotate(DP = hl.coalesce(res_table.DP, 0))
        res_table = res_table.annotate(GQ = hl.coalesce(res_table.GQ, 0))
        res_table = res_table.order_by(res_table.idx)
        res_table = res_table.checkpoint('db/checkpoint/ht2.ht', overwrite=True)
        return np.column_stack([np.array(res_table.filter(res_table.s==sample)[field].collect()).reshape(-1, 1) for sample in samples])
"""
    def _get(self, var_df, samples, field):
        # create table with vars in HAIL
        ht = hl.Table.from_pandas(var_df)
        #ht_filter = ht.key_by(**hl.parse_variant(ht.chrom + hl.literal(':') + hl.str(ht.pos) + hl.literal(':') + ht.ref + hl.literal(':') + ht.alt))
        #ht_filter = hl.MatrixTable.from_rows_table(ht_filter)
        #ht_filter = ht_filter.select_rows()
        #ht_filter = ht_filter.select_cols()
        # ht = ht.key_by(locus=hl.struct(contig=ht.chrom, position=ht.pos), alleles=hl.array([ht.ref, ht.alt]))

        # array with samples in HAIL
        if type(samples) is str:
            samples = [samples]
            
        hl_samples = hl.array(samples)
        
        # create table with samples
        df = pd.DataFrame({
            's': samples
        })
        ht_samples = hl.Table.from_pandas(df)
        ht_paths = self._get_Tables_paths(samples)
        
        ht = ht.join(ht_samples)
        ht = ht.annotate(pos=hl.int32(ht.pos))
        ht = ht.add_index()
        ht = ht.key_by(locus=hl.struct(contig=ht.chrom, position=ht.pos), alleles=hl.array([ht.ref, ht.alt]), s=ht.s)
        #ht = ht.persist()
        # checkpoints to speed up
        #ht.write('db/variants_checkpoint.ht', overwrite=True)
        #ht = hl.read_table('db/variants_checkpoint.ht')
        #ht = ht.checkpoint('db/checkpoint/variants_checkpoint.ht', overwrite=True)
        #ht_filter = ht_filter.checkpoint('db/checkpoint/variants_checkpoint.mt', overwrite=True)
        
        # all variants per sample
        res_table = None
        
        # iterate through mt
        for ht_path in ht_paths:
            mt = hl.read_table(ht_path)
            # get only lines with variants of interest
            #mt = mt.filter(~hl.is_missing(ht_filter.index_rows(hl.locus(mt.locus.contig, mt.locus.position), mt.alleles)))

            # filter samples
            #mt = mt.filter_cols(hl_samples.contains(mt.s))
            
            #mt = mt.annotate_entries(GT = mt.GT[0]+mt.GT[1])
            #mt = mt.annotate_entries(GT = hl.coalesce(mt.GT, -1))
            
            ht_n = ht.join(mt, how='left')
            
            if res_table is None:
                res_table = ht_n
            else:
                #res_table = res_table.checkpoint('db/checkpoint/res_table_in.ht', overwrite=True)
                res_table = res_table.union(ht_n)
                #res_table = res_table.checkpoint('db/checkpoint/res_table_out.ht', overwrite=True)

            #df = mt.entries().to_pandas()

            # format df
            #df[['ref','alt']] = pd.DataFrame(df.alleles.tolist(), index=df.index)
            #df['GT'] = df['GT.alleles'].apply(self._sum)  #mt = mt.annotate_entries(GT = mt.GT[0]+mt.GT[1])
            #                                            #mt = mt.annotate_entries(GT = hl.coalesce(mt.GT, 0))
            #                                            #mt = mt.annotate_entries(DP = hl.coalesce(mt.DP, 0))
            #                                            #mt = mt.annotate_entries(GQ = hl.coalesce(mt.GQ, 0))
            #df = df.drop(['alleles', 'GT.alleles', 'GT.phased'], axis=1)
            #df.columns = ['chrom', 'pos', 'sample', 'DP', 'GQ', 'ref', 'alt', 'GT']
            #df = df[['sample', 'chrom', 'pos', 'ref', 'alt', 'GT', 'DP', 'GQ']]

            # append to final df with all variants per sample
            #res_df = res_df.append(df, ignore_index = True)
        res_table = res_table.checkpoint('db/checkpoint/ht4.ht', overwrite=True)
        res_table = res_table.annotate(GT = hl.coalesce(res_table.GT, 0))
        res_table = res_table.annotate(DP = hl.coalesce(res_table.DP, 0))
        res_table = res_table.annotate(GQ = hl.coalesce(res_table.GQ, 0))
        res_table = res_table.order_by(res_table.idx)
        res_table = res_table
        #res_table.write('db/checkpoint.ht', overwrite=True)
        #res_table = hl.read_table('db/checkpoint.ht')
        #res_df = res_table.to_pandas()
        #res_df = res_df[['chrom', 'pos', 'ref', 'alt', 's', 'GT', 'DP', 'GQ']]

        res_np = list()
        for sample in samples:
            # only variant per sample
            #res_sample = res_df.groupby('s').get_group(sample)
            # get info
            #arr_sample = pd.merge(res_sample, var_df, on=['chrom', 'pos', 'ref', 'alt'], how='right')[field].astype('int32').to_numpy().reshape(-1, 1)
            #res_np.append(arr_sample)
            res_np.append(np.array(res_table.filter(res_table.s==sample)[field].collect()).reshape(-1, 1))

        return np.column_stack(res_np)

"""