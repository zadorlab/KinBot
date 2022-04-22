import os
import logging


def load_launchpad(lpad_file, **kwargs):  # TODO build a LaunchPad from kinbot parameters
    """Load a launchpad from a .yaml file and remove all its workflows."""
    from datetime import date
    from fireworks import LaunchPad

    if os.path.isfile(lpad_file):
        pass
    elif os.path.isfile("my_launchpad.yaml"):
        lpad_file = "my_launchpad.yaml"
    else:
        err_msg = "When using fireworks, either 'my_launchpad.yaml' should " \
                  "be present in the currrent working directory or the path " \
                  "to the launchpad .yaml file must be specified with the " \
                  "'lpad_file' keyword.\n" \
                  "(eg. \"lpad_file\": \"/home/my_user/my_new_launchpad.yaml\")"
        logging.error(err_msg)
        raise FileNotFoundError(err_msg)

    if not os.path.isfile('my_qadapter.yaml'):
        err_msg = "When using fireworks, 'my_qadapter.yaml' should " \
                  "be present in the currrent working directory."
        #           " or a path " \
        #           "to the launchpad .yaml file must be provided as value of " \
        #           "the fireworks keyword in the .json input file. " \
        #           "(eg. \"fireworks\" : \"/home/my_user/my_new_launchpad.yaml)"
        logging.error(err_msg)
        raise FileNotFoundError(err_msg)

    try:
        lpad = LaunchPad.from_file(lpad_file)
    except ValueError:
        err_msg = f"Incorrect format for fireworks launchpad file: {lpad_file}."
        logging.error(err_msg)
        raise ValueError(err_msg)

    if kwargs['reset']:
        lpad.reset(str(date.today()))   # TODO Implement a pid-specific reset.

    return lpad


def load_qadapter(qadapter_file, **kwargss):  # TODO build a QueueAdapter from kinbot parameters
    from fireworks.queue.queue_adapter import QueueAdapterBase
    return QueueAdapterBase.from_file(qadapter_file)


def load_fworker(fworker_file, **kwargs):  # TODO build a FWorker from kinbot parameters
    from fireworks import FWorker
    if os.path.isfile(fworker_file):
        pass
    elif os.path.isfile("my_fworker.yaml"):
        fworker_file = "my_fworker.yaml"
    try:
        return FWorker.from_file(fworker_file)
    except (FileNotFoundError, ValueError) as e:
        return FWorker()


def place_job(lpad, cmd, name="KinBot workflow"):
    from fireworks import Workflow, Firework, ScriptTask
    firetask = ScriptTask.from_str(cmd)
    firework = Firework(firetask)
    firework.spec['_launch_dir'] = os.getcwd()
    firework.name = os.getpid()
    workflow = Workflow([firework], name=name)
    lpad.add_wf(workflow)


def setup_fireworks(lpad_file='my_launchpad.yaml',
                    fworker_file='my_fworker.yaml',
                    qadapter_file='my_qadapter.yaml',
                    **kwargs):
    from subprocess import Popen
    import time

    if 'num_jobs' not in kwargs:
        num_jobs = 1
    elif 'num_jobs' == -1:
        num_jobs = 1
    else:
        num_jobs = int(kwargs['num_jobs'])

    lpad = load_launchpad(lpad_file, **kwargs)
    fworker = load_fworker(fworker_file)  # Not used now
    # qadapter = load_qadapter(qadapter_file)

    # TODO use python API
    cmd = f"qlaunch -fm singleshot"
    # cmd = f"qlaunch --launch_dir {os.getcwd()} --block_dir block_fake -fm rapidfire -m 3 --nlaunches infinite"

    for i in range(num_jobs):
        Popen(cmd, shell=True)
        time.sleep(1)

    return lpad
