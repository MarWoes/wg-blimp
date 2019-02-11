import subprocess
import tempfile
import shutil
import os
import atexit


def dump_dirs_to_file(project_dirs, target_file):

    with (open(target_file, 'w')) as f:

        f.write('{}\n'.format(project_dirs))


def get_shiny_starting_command_string(host, port):

    script_dir = os.path.dirname(os.path.realpath(__file__))
    shiny_code_dir = os.path.join(script_dir,'..','shiny')

    return 'shiny::runApp(appDir = "{}", host = "{}", port = {})'.format(shiny_code_dir, host, port)


def start_shiny(project_dirs, host, port):

    shiny_workdir = tempfile.mkdtemp()
    atexit.register(shutil.rmtree, shiny_workdir)

    projects_file = os.path.join(shiny_workdir, 'projects.txt')
    multiqc_dir = os.path.join(shiny_workdir, 'multiqc')

    dump_dirs_to_file(project_dirs, projects_file)

    shiny_running_command = [
        'Rscript',
        '-e', 'projectsFileArgument <- "{}"'.format(projects_file),
        '-e', 'multiqcDirArgument <- "{}"'.format(multiqc_dir),
        '-e', get_shiny_starting_command_string(host, port)
    ]

    subprocess.run(shiny_running_command)
