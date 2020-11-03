import os
import requests
from ruamel.yaml import YAML

script_dir = os.path.dirname(os.path.realpath(__file__))

annotation_dir = os.path.join(script_dir, '..', 'snakemake_wrapper', 'annotation')

CGI_LOCATION_FILE_TEMPLATE    = os.path.normpath(os.path.join(annotation_dir, 'cgi-locations-{}.csv.gz'))
REPEAT_MASKER_FILE_TEMPLATE   = os.path.normpath(os.path.join(annotation_dir, 'repeat-masker-{}.csv.gz'))
GTF_ANNOTATION_FILE_TEMPLATE  = os.path.normpath(os.path.join(annotation_dir, 'genes.{}.gtf.gz'))

ANNOTATION_DOWNLOAD_LINKS_FILE = os.path.join(script_dir, 'annotation_download_links.yaml')
DEFAULT_OPTIONALS_FILE = os.path.join(script_dir, 'optionals.yaml')

yaml = YAML(typ='safe')
yaml.default_flow_style = False


def download_if_necessary(file_name):

    if not os.path.isfile(file_name):

        print("[WARNING] " + file_name + " not found. Attempting download from public repositories.")

        with open(ANNOTATION_DOWNLOAD_LINKS_FILE) as f:

            annotation_links = yaml.load(f)

            base_file = os.path.basename(file_name)

            target_url = annotation_links['download_links'][base_file]

            print("[INFO] Downloading file " + base_file + " from " + target_url)

            request = requests.get(target_url, allow_redirects = True)
            open(file_name, 'wb').write(request.content)


def get_reference_annotation_files(genome_build):

    if genome_build == 'None':

        print("[WARNING] Segmentation with MethylSeekR requires CGI annotation.")
        print("[WARNING] Set 'cgi_annotation_file' accordingly or remove 'segmentation/umr-lmr-all.csv' from 'target_files' in configuration file to prevent errors.")

        return {
            'cgi_annotation_file': None,
            'gtf_annotation_file': None,
            'repeat_masker_annotation_file': None
        }

    cgi_location_file = CGI_LOCATION_FILE_TEMPLATE.format(genome_build)
    gtf_annotation_file = GTF_ANNOTATION_FILE_TEMPLATE.format(genome_build)
    repeat_masker_annotation_file = REPEAT_MASKER_FILE_TEMPLATE.format(genome_build)

    download_if_necessary(gtf_annotation_file)
    download_if_necessary(repeat_masker_annotation_file)

    return {
        'cgi_annotation_file': cgi_location_file,
        'gtf_annotation_file': gtf_annotation_file,
        'repeat_masker_annotation_file': repeat_masker_annotation_file
    }


def get_default_optional_parameters():

    with open(DEFAULT_OPTIONALS_FILE) as f:

        default_optionals = yaml.load(f)

        return default_optionals


def get_methylseekr_calibration_chromosome(genome_build):

    if genome_build is None:
        print("[WARNING] Genome build not set. Make sure to set chromosome for MethylSeekR calibration accordingly with parameter 'methylseekr_pmd_chromosome' in config.")

    options = {
        'hg19': 'chr22',
        'hg38': 'chr22',
        'None': 'chr22',
        'mmul10': 'chr20',
    }

    return options[genome_build]


def get_default_config(fastq_dir, fasta_ref, group1, group2, genome_build, cores_per_job, output_dir):

    mandatory_parameters = {
        'rawdir': os.path.abspath(fastq_dir),
        'output_dir': os.path.abspath(output_dir),
        'group1': group1,
        'group2': group2,
        'samples': group1 + group2,
        'ref': os.path.abspath(fasta_ref),
        'computing_threads': cores_per_job,
        'io_threads': cores_per_job
    }

    optional_parameters = get_default_optional_parameters()

    reference_annotation_files = get_reference_annotation_files(genome_build)

    methylseekr_calibration_chr = get_methylseekr_calibration_chromosome(genome_build)

    return {
        **mandatory_parameters,
        **optional_parameters,
        **reference_annotation_files,
        'methylseekr_pmd_chromosome': methylseekr_calibration_chr
    }


def dump_config(config_dict, target_file):

    with (open(target_file, 'w')) as f:

        yaml.dump(config_dict, f)


def read_samples_from_file(sample_file):

    with open(sample_file) as f:

        return f.read().splitlines()


def create_config(use_sample_files, genome_build, cores_per_job, fastq_dir, reference_fasta, group1, group2, output_dir, target_yaml):

    if (use_sample_files):

        samples_in_group1 = read_samples_from_file(group1)
        samples_in_group2 = read_samples_from_file(group2)

    else:

        samples_in_group1 = group1.split(',')
        samples_in_group2 = group2.split(',')

    config_yaml = get_default_config(fastq_dir, reference_fasta, samples_in_group1, samples_in_group2, genome_build, cores_per_job, output_dir)

    dump_config(config_yaml, target_yaml)
