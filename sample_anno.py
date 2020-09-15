import os
import pandas as pd
import sqlite3
from cyvcf2 import VCF
import numpy as np
from tqdm import tqdm

class SampleAnnotation():
    
    def __init__(self, sample_path='db/'):
        samples_file = os.path.join(sample_path, 'samples.db')

        if not os.path.exists(samples_file):
            if not os.path.exists(sample_path):
                os.mkdir(sample_path)
            self.conn = sqlite3.connect(samples_file)
            self.create_table()
        else:
            self.conn = sqlite3.connect(samples_file)
            self.load_anno_to_pd()
    
    def create_table(self):
        sql_create = """
        CREATE TABLE sample (
            sample_idx INTEGER,
            file_path VARCHAR,
            file_name VARCHAR,
            sample_id VARCHAR,
            dataset VARCHAR,
            read_status INTEGER,
            PRIMARY KEY (sample_idx));
        """
        c = self.conn.cursor()
        c.executescript(sql_create)
    
    def load_anno_to_pd(self):
        self.anno = pd.read_sql_query("SELECT * from sample", self.conn)
    
    
    def register_sample(self, file_paths, dataset=None):
        files_in_db = list(pd.read_sql_query("SELECT file_path from sample", self.conn)['file_path'])
        if isinstance(file_paths, str) and file_paths not in files_in_db:
            self._register_sample(file_paths, dataset=dataset)
        else:
            for file_path in tqdm(file_paths):
                if file_path not in files_in_db:
                    self._register_sample(file_path, dataset=dataset)

        self.load_anno_to_pd()
    
    def _register_sample(self, file_path, dataset):
        file_name = file_path.split('/')[-1].split('.')[0]
        sample_id = VCF(file_path, lazy=True).samples
        df = pd.DataFrame({'file_path': file_path,
                           'file_name': file_name,
                           'sample_id': sample_id,
                           'read_status': 0,
                           'dataset': dataset})
        df.to_sql('sample', self.conn, if_exists='append', index=False)
    
    
    def get_vcfs_for_samples(self, samples):
        if type(samples) == str:
            samples = [samples]
        
        sample_req = ','.join([f"'{sample}'" for sample in samples])
        sql_request = f"""
        SELECT DISTINCT file_path from sample
        WHERE sample_id in ({sample_req});
        """
        return pd.read_sql_query(sql_request, self.conn).file_path
        
    
    def set_file_path_as_read(self, file_path: str):
        sql_statement = f"""
        UPDATE sample
        SET read_status = 1
        WHERE file_path = '{file_path}';
        """
        c = self.conn.cursor()
        c.executescript(sql_statement)
    
    
    def list_unread_files(self):
        sql_request = """
        SELECT DISTINCT file_path from sample
        WHERE read_status = 0;
        """
        return pd.read_sql_query(sql_request, self.conn).file_path