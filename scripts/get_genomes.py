import subprocess 
import requests
import sys
import glob 
import re
import os

gitpath = os.path.expanduser("~/git/cross_mammal_xci")
sys.path.append(gitpath)
from scripts.utils import *


class GetGenomes(object):
    def __init__(self):
            self.tempdir = '/home/johlee/cross_mammal_xci/temp/'
            
            self.ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/'
            self.ens_url = 'http://ftp.ensembl.org/pub/current_fasta/'
            self.ens_url_gtf  = 'http://ftp.ensembl.org/pub/current_gtf/'



    def get_available_genomes(self):
            
        ens_cmd = ['wget','-r','--directory-prefix',f'{self.tempdir}', f'{self.ens_url}' ]
        ens_stderr, ens_stdout, ens_returncode = run_command( ens_cmd )     # from utils.py

        ens_htmlpath = self.ens_url.replace('http://',f'{self.tempdir}') 
        ens_htmlpath = f'{ens_htmlpath}index.html'

        cmd_ncbi = ['wget','-r', '--directory-prefix',f'{self.tempdir}', f'{self.ncbi_url}' ]
        ncbi_stderr, ncbi_stdout,ncbi_returncode = run_command( cmd_ncbi )    # from utils.py
        
        ncbi_htmlpath = self.ncbi_url.replace('https://',f'{self.tempdir}') 
        ncbi_htmlpath = f'{ncbi_htmlpath}index.html'


        ncbi_specs = get_links_from_html( ncbi_htmlpath )[1:-3]
        ens_specs = get_links_from_html( ens_htmlpath )[1:]

        all_species = set([spec.upper() for spec in ens_specs])  | set([spec.upper() for spec in ncbi_specs])

        rdf = self.get_available_runs()
        avail_specs = { sciname.replace(' ','_').upper() +'/' for sciname in rdf.ScientificName.unique() }
        specs_to_download = list(avail_specs & all_species)

        # look for links to the fa files.
        # prioritize fasta from ncbi
        for spec in specs_to_download :
            specens = spec.lower()
            specncbi = spec[0].upper() + spec[1:].lower()

            # does the species have a fasta file in ncbi?
            if specncbi in ncbi_specs :

                srcurl = f'{self.ncbi_url}{specncbi}latest_assembly_versions/'
                cmd = ['wget','-r','--directory-prefix',f'{self.tempdir}',f'{srcurl}']
                stderr, stdout, rc = run_command(cmd)

                if str(rc) =='0' :
                    htmlpath = srcurl.replace('https://',f'{self.tempdir}/') + 'index.html'
                    if os.path.isfile(htmlpath) :
                        # get the name of the latest version
                        print(f'downloading {htmlpath}')
                        version_name = get_links_from_html(htmlpath)[1][:-1]

                        fa_cmd = ['wget','-nc','-nv', '--directory-prefix', f'{self.tempdir}/{specncbi}/',f'{srcurl}{version_name}/{version_name}_genomic.fna.gz']
                        gtf_cmd= ['wget','-nc','-nv', '--directory-prefix', f'{self.tempdir}/{specncbi}/',f'{srcurl}{version_name}/{version_name}_genomic.gtf.gz']

                        fastderr, fastdout, farc = run_command(fa_cmd)
                        gtfstderr, gtfstdout, gtfrc = run_command(gtf_cmd)

                        print(fastdout)
                        print(gtfstdout)

            # ... if not, does it have a fasta file in ensembl?
            if specens in ens_specs :
                # srcurl = f'{self.ens_url}{specens}/latest_assembly_versions/'
                # cmd = ['wget','-r','--directory-prefix',f'{self.tempdir}',f'{srcurl}']
                # stderr, stdout, rc = run_command(cmd)


                pass
            


    def get_species_from_html(self,htmlpath) :
        with open(htmlpath) as f :
            html_lines=f.readlines()


        species = [ line.split("/</a>")[0] for line in html_lines ]
        species = [ spec.split('/">')[-1]  for spec in species]
        species = [ spec for spec in species if (not '<' in spec) and (not spec =='..') ]
        return(species)

    def get_available_runs(self,file_path = '/data/johlee/cross_mammal_xci/resource/all_mammal_runs.csv', sciname = None):
        df = pd.read_csv(file_path, sep=",",index_col=0)
        if sciname :
            df.loc[df.ScientificName == sciname, : ]

        return(df)


    def get_runs(self):
        '''
        reads the sra metadata for genus specific runs, concatonates and filters 
        only needs to be done once.
        '''

        paths = glob.glob('/data/johlee/cross_mammal_xci/resource/mammaldata/*')
        alldfs = []
        print('reading genus specific dataframes.')
        for fpath in paths :
            alldfs.append(pd.read_csv(fpath,sep=",", dtype = str))
    
        print('joining dataframes and filtering')
        df = pd.concat(alldfs)
        df = df.sort_values('ScientificName')
        # remove rows where the species is not provided/ failed to 
        df = df.loc[~ df.ScientificName.isna(),:]
        # remove rows where the gender is confirmed male
        df = df.loc[~ df.Sex.isin(['Male','male','male (gelding)','maie','gelding','M','male and neuter' ]) ,: ].reset_index(drop=True)
        # remove monotremes
        df.loc[ ~ df.ScientificName.str.contains('Tachyglossus|Ornithorhynchus') ,:] 
        # save to disk
        print('saving to disk')
        df.to_csv('/data/johlee/cross_mammal_xci/resource/all_mammal_metadata.csv',sep=",")
        df[['Run','Experiment','Sample','SRAStudy','TaxID','ScientificName','Sex']].to_csv('/data/johlee/cross_mammal_xci/resource/all_mammal_runs.csv',sep=",")

        return(df)


def get_links_from_html(htmlpath):
    with open(htmlpath) as f:
        htmltxt = f.read()

    html_tag_regex = re.compile(r'<a[^<>]+?href=([\'\"])(.*?)\1', re.IGNORECASE)
    return [match[1] for match in html_tag_regex.findall(htmltxt)]




