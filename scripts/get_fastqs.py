import pandas as pd
import sys
import os 
import numpy as np
import logging 
import time
from requests.api import get 

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.sra import *

class Query(object) :
    def __init__(self):
        # self.species = 'Delphinapterus leucas'
        self.resourcedir = '/data/johlee/cross_mammal_xci/resource'
        self.metadir = '/data/johlee/cross_mammal_xci/metadata'
        self.tempdir = '/data/johlee/cross_mammal_xci/temp'
        self.runs_df = f'{self.resourcedir}/all_mammal_runs.tsv'
        self.xid_batchsize = 20
        self.log = logging.getLogger('sra')

        self.sra_efetch= 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra'

        pass

    def execute(self, species = 'Delphinapterus leucas'):
        pdf = get_available_runs(file_path = self.runs_df ,sciname = species)
        explist = list(set(pdf.Experiment))
        proj_rows = []
        samp_rows = []
        exp_rows = []
        run_rows = []

        nexps = len(explist)
        nbatches = int(np.ceil(nexps / self.xid_batchsize ))
        for i in range(nbatches): 
            exps = explist[ i*self.xid_batchsize: (i+1)*self.xid_batchsize  ]
            exps = ",".join(exps)
            # if exp in self.expids:
            exd = self.query_experiment_package_set(exps)
            (projrows, samprows, exprows,
                runs) = self.parse_experiment_package_set(exd)
            proj_rows = itertools.chain(proj_rows, projrows)
            samp_rows = itertools.chain(samp_rows, samprows)
            exp_rows = itertools.chain(exp_rows, exprows)
            run_rows = itertools.chain(run_rows, runs)
            # else:
            #     self.log.debug(f'expid {exp} not in set of search expids.  ')
        proj_rows = list(proj_rows)
        samp_rows = list(samp_rows)
        exp_rows = list(exp_rows)
        run_rows = list(run_rows)

        pass


    def query_experiment_package_set(self, xid):
        """
        Query XML data for this experiment ID. 

        """
        xmldata = None
        try:
            url = f"{self.sra_efetch}&id={xid}"
            self.log.debug(f"fetch url={url}")

            # XXX could run infinitely???
            while True:
                r = requests.post(url)
                if r.status_code == 200:
                    xmldata = r.content.decode()
                    self.log.debug(f'good HTTP response for {xid}')
                    break
                else:
                    self.log.warn(
                        f'bad HTTP response for id {xid}. retry in 10s')
                    time.sleep(10)

        except Exception as ex:
            self.log.error(f'problem with NCBI id {xid}')
            logging.error(traceback.format_exc(None))

        finally:
            self.log.debug(
                f"sleeping {1} secs between fetch calls...")
            time.sleep(1)
        return xmldata


    def parse_experiment_package_set(self, xmlstr):
        """
        package sets should have one package per uid pulled via efetch, e.g.

        https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=12277089,12277091
        https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=13333495 

        """
        root = et.fromstring(xmlstr)
        self.log.debug(f"root={root}")
        proj_rows = []
        samp_rows = []
        exp_rows = []
        run_rows = []

        n_processed = 0
        for exp in root.iter("EXPERIMENT_PACKAGE"):
            (projrow, samprow, exprow, newruns) = self.parse_experiment_package(exp)
            proj_rows.append(projrow)
            samp_rows.append(samprow)
            exp_rows.append(exprow)
            # run_rows.append(newruns)
            run_rows = itertools.chain(run_rows, newruns)
            n_processed += 1
        self.log.debug(f"processed {n_processed} experiment package(s).")
        run_rows = list(run_rows)
        self.log.debug(
            f'returning\n    proj_rows: {proj_rows}\n    exp_rows: {exp_rows} \n    run_rows: {run_rows}')
        return (proj_rows, samp_rows, exp_rows, run_rows)

    def parse_experiment_package(self, root):
        """
        NCBI provides no XSD, so we shouldn't rely on order

        """
        self.log.debug('parsing experiment package...')
        exp = root.find('EXPERIMENT')
        sub = root.find('SUBMISSION')
        proj = root.find('STUDY')
        samp = root.find('SAMPLE')
        runs = root.find('RUN_SET')

        # get experiment properties.

        # get submission properties
        sra_id = sub.get('accession')

        # get study/project properties title, abstract
        projrow = self.parse_proj(proj)
        projrow.append(sra_id)
        proj_id = projrow[0]

        # get sample properties - append project and sra ids
        samprow = self.parse_sample(samp)
        samprow.append(proj_id)
        samprow.append(sra_id)
        samp_id = samprow[0]

        # get experiment properties
        exprow = self.parse_exp(exp)
        exprow.append(sra_id)
        exp_id = exprow[0]

        # get run properties - list of lists
        runrows = self.parse_run_set(runs, proj_id, sra_id)

        self.log.debug(
            f'exprow: proj_id={proj_id} exp_id={exp_id} sra_id={sra_id} samp_id={samp_id}')
        # projrow = [proj_id, title, pubdate, abstract]
        # exprow = [proj_id, exp_id, sra_id, gsm, gse, lcp, sample_attributes]
        self.log.debug(
            f'\n  projrow: {projrow}\n   exprow: {exprow} \n  runrows: {runrows}')
        return(projrow, samprow, exprow, runrows)

    def parse_run_set(self, runs, proj_id, sra_id):
        """

        """
        runrows = []
        for run in runs.findall('RUN'):
            runrow = self.parse_run(run)
            runrow.append(proj_id)
            runrow.append(sra_id)
            runrows.append(runrow)
        return runrows

    # need to append project id
    def parse_run(self, run):

        run_ext_ids = {}
        ids = run.find('IDENTIFIERS')
        run_id = ids.find('PRIMARY_ID').text
        for elem in ids.findall('EXTERNAL_ID'):
            tag = elem.get('namespace')
            val = elem.text
            run_ext_ids[tag] = val
        run_ext_ids = str(run_ext_ids)
        avail_status = run.get('unavailable')
        if avail_status == 'true':
            raise RunUnavailableException(f'run data unavailable for {run_id}')

        total_spots = run.get('total_spots')
        total_bases = run.get('total_bases')
        run_size = run.get('size')
        pdate = run.get('published')

        expid = run.find('EXPERIMENT_REF').get('accession')
        pool = run.find('Pool')
        sampleid = pool.find('Member').get('accession')
        taxon = pool.find('Member').get('tax_id')
        organism = pool.find('Member').get('organism')

        nreads = run.find('Statistics').get('nreads')

        srafiles = run.find('SRAFiles')

        (file_url, file_size) = self.parse_srafiles(srafiles, run_id)

        bases = run.find('Bases')
        basecounts = {}
        for base in bases:
            tag = base.get('value')
            val = base.get('count')
            basecounts[tag] = val

        basecounts = str(basecounts)

        runrow = [run_id, run_ext_ids, total_spots, total_bases, run_size, pdate,
                  taxon, organism, nreads,  basecounts, file_url, file_size , expid, sampleid ]

        return runrow

    def parse_srafiles(self, srafiles, runid):
        """
        Get file info. Choose Amazon URL if available. 
        
        """
        file_url = None
        file_size = None
        for srafile in srafiles.findall('SRAFile'):
            if srafile.get('filename') == runid:
                file_size = srafile.get('size')
                file_url = srafile.get('url')
                for altern in srafile.findall('Alternatives'):
                    url = altern.get('url')
                    if 'amazonaws.com' in url:
                        self.log.debug(f'found AWS alternate: {url}  Using...')
                        file_url = url
        self.log.info(f'got info for {runid}: {file_size} {file_url}')
        return (file_url, file_size)

    def parse_exp(self, exp):

        exp_ext_ids = {}
        ids = exp.find('IDENTIFIERS')
        exp_id = ids.find('PRIMARY_ID').text
        self.log.debug(f'parsing exp_id: {exp_id}')
        for elem in ids.findall('EXTERNAL_ID'):
            tag = elem.get('namespace')
            val = elem.text
            exp_ext_ids[tag] = val
        exp_ext_ids = str(exp_ext_ids)
        projid = exp.find('STUDY_REF').get('accession')
        des = exp.find('DESIGN')
        sampid = des.find('SAMPLE_DESCRIPTOR').get('accession')
        ldes = des.find('LIBRARY_DESCRIPTOR')

        lcp = ""
        strat = ""
        source = ""
        try:
            lcp = ldes.find('LIBRARY_CONSTRUCTION_PROTOCOL').text
            lcp = lcp.strip()
            strat = ldes.find('LIBRARY_STRATEGY').text
            source = ldes.find('LIBRARY_SOURCE').text
        except AttributeError as ae:
            self.log.warn(f'attribute error parsing LIBRARY_DESCRIPTOR children.')

        exprow = [exp_id, exp_ext_ids,  strat, source, lcp,  sampid, projid]
        return exprow

    # need to append project id - done in execute
    def parse_sample(self, samp):
        samp_ext_ids = {}
        ids = samp.find('IDENTIFIERS')
        samp_id = ids.find('PRIMARY_ID').text
        for elem in ids.findall('EXTERNAL_ID'):
            tag = elem.get('namespace')
            val = elem.text
            samp_ext_ids[tag] = val
        samp_ext_ids = str(samp_ext_ids)
        
        try:
            samptitle = samp.find('TITLE').text
        except:
            samptitle = samp_id

        sample_attributes = {}
        for elem in samp.find('SAMPLE_ATTRIBUTES').findall('SAMPLE_ATTRIBUTE'):
            tag = elem.find('TAG').text
            val = elem.find('VALUE').text
            sample_attributes[tag] = val
        sample_attributes = str(sample_attributes)

        taxid = samp.find('SAMPLE_NAME').find('TAXON_ID').text
        sciname = samp.find('SAMPLE_NAME').find('SCIENTIFIC_NAME').text
        samprow = [samp_id, samp_ext_ids,  taxid,
                   sciname, samptitle, sample_attributes]

        return samprow

    def parse_proj(self, proj):

        proj_ext_ids = {}
        ids = proj.find('IDENTIFIERS')
        proj_id = ids.find('PRIMARY_ID').text
        for elem in ids.findall('EXTERNAL_ID'):
            tag = elem.get('namespace')
            val = elem.text
            proj_ext_ids[tag] = val
        proj_ext_ids = str(proj_ext_ids)  # convert to strings to store in df
        d_elem = proj.find('DESCRIPTOR')
        title = d_elem.find('STUDY_TITLE').text
        abstract = d_elem.find('STUDY_ABSTRACT').text

        projrow = [proj_id, proj_ext_ids, title, abstract]
        return projrow

    pass


class GetFastq(object):

    def __init__(self) :
        self.resourcedir = '/data/johlee/cross_mammal_xci/resource'
        self.tempdir = '/data/johlee/cross_mammal_xci/temp'
        self.runs_df = f'{self.resourcedir}/all_mammal_runs.tsv'
        # species in ncbi refseq have both fna and gtf files

        self.spec = 'Delphinapterus leucas'

    def execute(self):
        # get all runs associatated with the species. 
        rdf = get_available_runs(file_path = self.runs_df ,sciname = self.spec)
        runs = rdf.Run
        projs = rdf.SRAStudy.unique()
        # 'SRP097694', 'SRP120137', 'SRP106871'
        

        for run_id in runs :
            # prefetch

            # fasterq dump

            pass




def get_available_runs(file_path = '/data/johlee/cross_mammal_xci/resource/all_mammal_runs.csv', sciname = None):
    df = pd.read_csv(file_path, sep="\t",index_col=0)
    df = df.loc[~ df.duplicated(),:].reset_index(drop=True)
    print(sciname)
    if sciname :
        print(sciname)

        df = df.loc[df.ScientificName == sciname, : ].reset_index(drop=True)

    return(df)