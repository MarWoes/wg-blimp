import os
import snakemake
import cli.config

# get right directory, see https://stackoverflow.com/questions/4934806/how-can-i-find-scripts-directory-with-python
script_dir = os.path.dirname(os.path.realpath(__file__))
snakefile_location = os.path.join(script_dir,'..','snakemake', 'Snakefile')


def run_snakemake_from_config(dry_run, config_yaml, delete_all_output=False, cores=1):

    config_dir= os.path.dirname(os.path.realpath(config_yaml))

    snakemake.snakemake(
        snakefile=snakefile_location,
        configfile=config_yaml,
        workdir=config_dir,
        dryrun=dry_run,
        cores=cores,
        printshellcmds=True,
        delete_all_output=delete_all_output
    )


def run_snakemake(dry_run, use_sample_files, cores, genome_build, fastq_dir, reference_fasta, group1, group2, output_dir):

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    yaml_config_file = os.path.join(output_dir, 'config.yaml')

    cli.config.create_config(use_sample_files, genome_build, cores, fastq_dir, reference_fasta, group1, group2, output_dir, yaml_config_file)

    run_snakemake_from_config(dry_run, yaml_config_file, cores=1)


def delete_all_output(dry_run, config_yaml):

    run_snakemake_from_config(dry_run, config_yaml, True)
