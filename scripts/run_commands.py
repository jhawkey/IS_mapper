from subprocess import call
from subprocess import check_output, CalledProcessError, STDOUT
import logging

# Exception classes
class CommandError(Exception):
    pass


class BedtoolsError(Exception):
    pass

def run_command(command, **kwargs):
    """
    Execute a shell command and check the exit status and any O/S exceptions.
    """

    command_str = ' '.join(command)
    logging.info('Running: {}'.format(command_str))
    try:
        exit_status = call(command_str, **kwargs)
    except OSError as e:
        message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
        raise CommandError({"message": message})
    if exit_status == 139 and command[0] == 'closestBed':
        raise BedtoolsError({'message':'One or more bed files are empty. Writing out empty results table.'})
    if exit_status != 0:
        message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
        raise CommandError({"message": message})

def check_command(command_call, command_name):
    '''
    Check that the dependency is installed.
    Exits the program if it can't be found.
        - command_list is the command to run to determine the version.
        - command_name is the name of the command to show in the error message.
    '''
    try:
        command_stdout = check_output(command_call, stderr=STDOUT)
        logging.info('Found dependency %s', command_name)
    except OSError as e:
        logging.error("Failed command: %s", command_call)
        logging.error(str(e))
        logging.error("Do you have %s installed in your PATH?", command_name)
        raise CommandError
    except CalledProcessError as e:
        # some programs such as samtools return a non zero exit status
        # when you ask for the version. We ignore it here.
        command_stdout = e.output
        logging.info('Found dependency %s', command_name)

def make_directories(dir_list):

    """
    Take a list of folders and make each directory.
    """
    for directory in dir_list:
        run_command(['mkdir', '-p', directory], shell=True)