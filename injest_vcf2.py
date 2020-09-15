import hail as hl
import glob
import pandas as pd
from tqdm import tqdm

def injest_vcf1(input_vcfs, output_folder='db/MatrixTables'):
    hl.init(quiet=True)

    # sanitize chr
    recode = {f"chr{i}":f"{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
    
    if type(input_vcfs) == str:
        input_vcfs = [input_vcfs]
    
    for input_vcf in tqdm(input_vcfs):
        vcf_name = (input_vcf.split('/')[-1]).split('.')[0]
        try:
            # load vcf
            mt = hl.methods.import_vcf(input_vcf, contig_recoding=recode, force_bgz=True, reference_genome=None)
            
            # clean information
            mt = (mt.select_entries(mt.GT, mt.DP, mt.GQ))
            mt = mt.select_rows()
            mt = mt.select_cols()
            
            valid = hl.is_valid_locus(mt.locus.contig, mt.locus.position, reference_genome='GRCh37')
            
            mt_wrong = mt.filter_rows(~valid)
            mt_correct = mt.filter_rows(valid)
            
            # store correct variants in MatrixTables
            if mt_correct.rows().count() > 0:
                mt_correct.write(f'{output_folder}/{vcf_name}.mt', overwrite=True)
            
            # store incorrect variants in tsv table
            if mt_wrong.rows().count() > 0:
                mt_wrong.rows().export(f'db/errors/{vcf_name}.tsv')
            
        except Exception as e:
            with open('invalid_vcf.txt', 'a') as f:
                f.write(f'{input_vcf}\n')
            with open('error_log.txt', 'a') as log:
                log.write(f'{input_vcf}\n{e}\n')
    
    return 'done'


def injest_vcf2(input_vcfs, output_folder='db/Tables'):

    # sanitize chr
    recode = {f"chr{i}":f"{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
    
    if type(input_vcfs) == str:
        input_vcfs = [input_vcfs]
    table = None
    for input_vcf in tqdm(input_vcfs):
        vcf_name = (input_vcf.split('/')[-1]).split('.')[0]
        try:
            # load vcf
            mt = hl.methods.import_vcf(input_vcf, contig_recoding=recode, force_bgz=True, reference_genome=None)
            
            # clean information
            mt = (mt.select_entries(mt.GT, mt.DP, mt.GQ))
            mt = mt.select_rows()
            mt = mt.select_cols()
            
            valid = hl.is_valid_locus(mt.locus.contig, mt.locus.position, reference_genome='GRCh37')
            
            mt_wrong = mt.filter_rows(~valid)
            mt_correct = mt.filter_rows(valid)
            
            mt_correct = mt_correct.annotate_entries(GT = mt_correct.GT[0]+mt_correct.GT[1])
            mt_correct = mt_correct.annotate_entries(GT = hl.coalesce(mt_correct.GT, -1))
            mt_correct = mt_correct.annotate_entries(DP = hl.coalesce(mt_correct.DP, 0))
            mt_correct = mt_correct.annotate_entries(GQ = hl.coalesce(mt_correct.GQ, 0))
            
            # store correct variants in MatrixTables
            if mt_correct.rows().count() > 0:
                if table is None:
                    table = mt_correct.entries()
                else:
                    table = table.join(mt_correct.entries())
                #mt_correct.entries().write(f'{output_folder}/{vcf_name}.ht', overwrite=True)
            
            # store incorrect variants in tsv table
            if mt_wrong.rows().count() > 0:
                mt_wrong.rows().export(f'db/errors/{vcf_name}.tsv')
            
        except Exception as e:
            with open('invalid_vcf.txt', 'a') as f:
                f.write(f'{input_vcf}\n')
            with open('error_log.txt', 'a') as log:
                log.write(f'{input_vcf}\n{e}\n')
    table.write('db/Table/table.ht', overwrite=True)
    
    return 'done'