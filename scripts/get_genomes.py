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
        # ens downloads  ~ 1MB / sec
        # ncbi downloads ~ 40MB / sec
         
        self.tempdir = '/home/johlee/cross_mammal_xci/temp/'
        # species in ncbi refseq have both fna and gtf files
        self.ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/'
        # need to intersect which species have both fasta and gtf from ens (different from ncbi)
        self.ens_url_fa = 'http://ftp.ensembl.org/pub/current_fasta/'
        self.ens_url_gtf  = 'http://ftp.ensembl.org/pub/current_gtf/'



    def get_available_genomes(self):
            
        ens_cmd_fa = ['wget','-r','--directory-prefix',f'{self.tempdir}', f'{self.ens_url_fa}' ]
        ens_stderr, ens_stdout, ens_returncode = run_command( ens_cmd_fa )     # from utils.py

        ens_htmlpath_fa = self.ens_url_fa.replace('http://',f'{self.tempdir}') 
        ens_htmlpath_fa = f'{ens_htmlpath_fa}index.html'

        ens_cmd_gtf = ['wget','-r','--directory-prefix',f'{self.tempdir}', f'{self.ens_url_gtf}' ]
        ens_stderr, ens_stdout, ens_returncode = run_command( ens_cmd_gtf )     # from utils.py

        ens_htmlpath_gtf = self.ens_url_gtf.replace('http://',f'{self.tempdir}') 
        ens_htmlpath_gtf = f'{ens_htmlpath_gtf}index.html'

        cmd_ncbi = ['wget','-r', '--directory-prefix',f'{self.tempdir}', f'{self.ncbi_url}' ]
        ncbi_stderr, ncbi_stdout,ncbi_returncode = run_command( cmd_ncbi )    # from utils.py
        
        ncbi_htmlpath = self.ncbi_url.replace('https://',f'{self.tempdir}') 
        ncbi_htmlpath = f'{ncbi_htmlpath}index.html'


        ncbi_specs = get_links_from_html( ncbi_htmlpath )[1:-3]
        ens_specs_fa = get_links_from_html( ens_htmlpath_fa )[1:]
        ens_specs_gtf = get_links_from_html( ens_htmlpath_gtf )[1:]

        ens_specs = set(ens_specs_gtf) &set(ens_specs_fa)
        all_species = set([spec.upper() for spec in ens_specs])  | set([spec.upper() for spec in ncbi_specs])

        rdf = get_available_runs()
        avail_specs = { sciname.replace(' ','_').upper() +'/' for sciname in rdf.ScientificName.unique() }
        specs_to_download = list(avail_specs & all_species)

        # look for links to the fa files.
        # prioritize fasta from ncbi
        for spec in specs_to_download :
            specens = spec.lower()
            specncbi = spec[0].upper() + spec[1:].lower()
                        # does the species have a fasta file in ncbi?
            if specncbi in ncbi_specs :
                try : 
                    srcurl = f'{self.ncbi_url}{specncbi}latest_assembly_versions/'
                    cmd = ['wget','-r','--directory-prefix',f'{self.tempdir}',f'{srcurl}']
                    stderr, stdout, rc = run_command(cmd)

                    if str(rc) =='0' :
                        htmlpath = srcurl.replace('https://',f'{self.tempdir}') + 'index.html'
                        if os.path.isfile(htmlpath) :
                            # get the name of the latest version
                            print(f'downloading {htmlpath}')
                            version_name = get_links_from_html(htmlpath)[1][:-1]

                            fa_cmd = ['wget','-nc','-nv', '--directory-prefix', f'{self.tempdir}{specncbi}/',f'{srcurl}{version_name}/{version_name}_genomic.fna.gz']
                            gtf_cmd= ['wget','-nc','-nv', '--directory-prefix', f'{self.tempdir}{specncbi}/',f'{srcurl}{version_name}/{version_name}_genomic.gtf.gz']

                            fastderr, fastdout, farc = run_command(fa_cmd)
                            gtfstderr, gtfstdout, gtfrc = run_command(gtf_cmd)

                            print(fastdout)
                            print(gtfstdout)
                except :
                    print(f'failed for {srcurl}')
            # ... if not, does it have a fasta file in ensembl?
            elif specens in ens_specs :
                print(f'downloading {specens} from ensembl')
                srcurl = f'{self.ens_url_fa}{specens}dna/'
                cmd = ['wget','-r','--directory-prefix',f'{self.tempdir}',f'{srcurl}']
                stderr, stdout, rc = run_command(cmd)
                htmlpath = srcurl.replace('http://',f'{self.tempdir}') + 'index.html'
                fa_files = get_links_from_html(htmlpath)
                # Cercocebus_atys.Caty_1.0.dna.toplevel.fa.gz 
                fa_file = [ fa for fa in fa_files if '.dna.toplevel.fa.gz' in fa  ][0]
                fa_cmd = ['wget','-nc','-nv', '--directory-prefix', f'{self.tempdir}{specncbi}/', f'{srcurl}{fa_file}']

                srcurl = f'{self.ens_url_gtf}{specens}'
                cmd = ['wget','-r','--directory-prefix',f'{self.tempdir}',f'{srcurl}']
                stderr, stdout, rc = run_command(cmd)
                htmlpath = srcurl.replace('http://',f'{self.tempdir}') + 'index.html'
                gtf_files = get_links_from_html(htmlpath)
                # Cercocebus_atys.Caty_1.0.dna.toplevel.fa.gz 
                gtf_files = [ gtf for gtf in gtf_files if '.gtf.gz' in gtf  ]
                gtf_file = [ gtf for gtf in gtf_files if 'abinitio' not in gtf ][0]
                gtf_cmd = ['wget','-nc','-nv', '--directory-prefix', f'{self.tempdir}{specncbi}/', f'{srcurl}{gtf_file}']
                stderr, stdout, rc = run_command(gtf_cmd)

                
                
            


    def get_species_from_html(self,htmlpath) :
        with open(htmlpath) as f :
            html_lines=f.readlines()


        species = [ line.split("/</a>")[0] for line in html_lines ]
        species = [ spec.split('/">')[-1]  for spec in species]
        species = [ spec for spec in species if (not '<' in spec) and (not spec =='..') ]
        return(species)






def get_links_from_html(htmlpath):
    with open(htmlpath) as f:
        htmltxt = f.read()

    html_tag_regex = re.compile(r'<a[^<>]+?href=([\'\"])(.*?)\1', re.IGNORECASE)
    return [match[1] for match in html_tag_regex.findall(htmltxt)]


# not used here
def get_available_runs(file_path = '/data/johlee/cross_mammal_xci/resource/all_mammal_runs.csv', sciname = None):
    df = pd.read_csv(file_path, sep="\t",index_col=0)
    df = df.loc[~ df.duplicated(),:].reset_index(drop=True)
    print(sciname)
    if sciname :
        print(sciname)

        df = df.loc[df.ScientificName == sciname, : ].reset_index(drop=True)

    return(df)

def merge_queries():
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
    print('removing monotremes')
    df.loc[ ~ df.ScientificName.str.contains('Tachyglossus|Ornithorhynchus') ,:] 

    df = df.loc[~ df.duplicated() ,:].reset_index(drop=True)

    # save to disk
    print('saving to disk')
    df.to_csv('/data/johlee/cross_mammal_xci/resource/all_mammal_metadata.tsv',sep="\t")
    df[['Run','Experiment','Sample','SRAStudy','TaxID','ScientificName','Sex']].to_csv('/data/johlee/cross_mammal_xci/resource/all_mammal_runs.tsv',sep="\t")
    # SRR15209859,SRR15209849,SRR15209874,SRR15209875,SRR15206524
    return(df)


