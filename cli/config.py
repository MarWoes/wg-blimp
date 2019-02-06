from ruamel.yaml import YAML

CGI_LOCATION_FILE    = 'annotation/cgi-locations-{}.csv.gz'
GENE_ANNOTATION_FILE = 'annotation/gene-locations-{}.csv.gz'
REPEAT_MASKER_FILE   = 'annotation/repeat-masker-{}.csv.gz'
TSS_FILE             = 'annotation/transcription-start-sites-{}.csv.gz'

DEFAULT_OPTIONALS_FILE = 'cli/optionals.yaml'

yaml = YAML(typ='safe')
yaml.default_flow_style = False

def get_reference_annotation_files(genome_build):

    return {
        'cgi_annotation_file':
          CGI_LOCATION_FILE.format(genome_build),

        'gene_annotation_file':
          GENE_ANNOTATION_FILE.format(genome_build),

        'repeat_masker_annotation_file':
          REPEAT_MASKER_FILE.format(genome_build),

        'transcript_start_site_file':
          TSS_FILE.format(genome_build)
    }

def get_default_optional_parameters():

    with open(DEFAULT_OPTIONALS_FILE) as f:

        default_optionals = yaml.load(f)

        return default_optionals

def get_default_config(fastq_dir, fasta_ref, group1, group2, genome_build, output_dir):

    mandatory_parameters = {
        'rawdir': fastq_dir,
        'output_dir': output_dir,
        'group1': group1,
        'group2': group2,
        'samples': group1 + group2,
        'ref': fasta_ref
    }

    optional_parameters = get_default_optional_parameters()

    reference_annotation_files = get_reference_annotation_files(genome_build)

    return {
        **mandatory_parameters,
        **optional_parameters,
        **reference_annotation_files
    }

def dump_config(config_dict, target_file):

    with (open(target_file, 'w')) as f:

        yaml.dump(config_dict, f)
