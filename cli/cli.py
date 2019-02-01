import click
import cli.config

@click.group()
def main():
    pass

@main.command()
def create_config():

    click.echo('welcome to wg-blimp')
    test_config = cli.config.get_default_config('/scripts/test/fastq', '/scripts/test/lambda-phage.fa', '/output', ['simulated1', 'simulated2'], ['simulated3', 'simulated4'], 'hg38')

    cli.config.dump_config(test_config, 'test-config.yaml')
