import click

@click.group()
def main():
    pass

@main.command()
def create_config():

    click.echo('welcome to wg-blimp')
