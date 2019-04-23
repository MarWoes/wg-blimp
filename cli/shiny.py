import subprocess
import tempfile
import shutil
import os
import atexit


def dump_config_paths_to_file(config_files, target_file):

    with (open(target_file, 'w')) as f:

        for config_file in config_files:

            f.write('{}\n'.format(os.path.abspath(config_file)))


def get_shiny_starting_command_string(host, port):

    script_dir = os.path.dirname(os.path.realpath(__file__))
    shiny_code_dir = os.path.join(script_dir,'..','shiny')

    return 'shiny::runApp(appDir = "{}", host = "{}", port = {})'.format(shiny_code_dir, host, port)


def start_shiny(config_files, host, port):

    shiny_workdir = tempfile.mkdtemp()
    atexit.register(shutil.rmtree, shiny_workdir)

    projects_file = os.path.join(shiny_workdir, 'projects.txt')
    multiqc_dir = os.path.join(shiny_workdir, 'multiqc')

    dump_config_paths_to_file(config_files, projects_file)

    shiny_running_command = [
        'Rscript',
        '-e', 'projectsFileArgument <- "{}"'.format(projects_file),
        '-e', 'multiqcDirArgument <- "{}"'.format(multiqc_dir),
        '-e', get_shiny_starting_command_string(host, port)
    ]

    subprocess.run(shiny_running_command)
