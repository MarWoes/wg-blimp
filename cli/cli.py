import click
import cli.config

def read_samples_from_file(sample_file):

    with open(sample_file) as f:

        return f.readlines()

@click.group()
def main():
    pass

@main.command(help='Create a config YAML file for running the Snakemake pipeline.')
@click.option('--use-sample-files', is_flag=True, default=False, help='Load sample names from text files instead of passing them as a comma-seperated list.')
@click.option('--genome_build', type=click.Choice(['hg19','hg38']), default='hg38', help='Build of the reference used for annotation.')
@click.argument('fastq_dir')
@click.argument('reference_fasta')
@click.argument('group1')
@click.argument('group2')
@click.argument('output_dir')
@click.argument('target_yaml')
def create_config(use_sample_files, genome_build, fastq_dir, reference_fasta, group1, group2, output_dir, target_yaml):

    if (use_sample_files):

        samples_in_group1 = read_samples_from_file(group1)
        samples_in_group2 = read_samples_from_file(group2)

    else:

        samples_in_group1 = group1.split(',')
        samples_in_group2 = group2.split(',')

    config_yaml = cli.config.get_default_config(fastq_dir, reference_fasta, samples_in_group1, samples_in_group2, genome_build, output_dir)

    cli.config.dump_config(config_yaml, target_yaml)
