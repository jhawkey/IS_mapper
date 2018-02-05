from subprocess import call
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

def make_directories(dir_list):

    """
    Take a list of folders and make each directory.
    """
    for directory in dir_list:
        run_command(['mkdir', '-p', directory], shell=True)