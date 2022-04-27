"""Module to manage the execution of FireWorks.

"""
import os
import logging


def reset_lpad(lpad):
    """Eliminates the jobs of former kinbot executions from the launchapd.

    @param lpad: LaunchPad object to manage Workflows, FireWorks and tasks.

    This is mainly used to prevent fireworks tasks to remain in the RUN status
    after a job reaches its walltime.
    """
    import json
    pids = []
    if os.path.isfile('pids'):
        with open('pids', 'r') as pids_fh:
            pids = [int(pid) for pid in pids_fh]
    elif os.path.isfile('FW.json'):
        with open('FW.json', 'r') as json_fh:
            fw_dict = json.load(json_fh)
            pids = [fw_dict['name']]

    for fw_id in lpad.get_fw_ids():
        fw = lpad.get_fw_by_id(fw_id)
        if fw.name in pids:
            lpad.delete_fws([fw_id])


def load_launchpad(lpad_file, **kwargs):  # TODO build a LaunchPad from kinbot parameters
    """Load a launchpad from a .yaml file and remove reset the present fireworks.

    @param lpad_file: Path to the yaml file containing the settings for Fireworks
        to connect to the MongoDB FireServer storing all the tasks to carry out.
    """
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
        reset_lpad(lpad)
        # lpad.reset(str(date.today()))  # Old way

    return lpad


def load_qadapter(qadapter_file, **kwargss):  # TODO build a QueueAdapter from kinbot parameters
    """Load a QAdapter from a yaml file or build one from KinBot input parameters.

    @param qadapter_file: Path to the yaml file containing the settings for Fireworks
        to manage the scheduler operations.
    """
    from fireworks.queue.queue_adapter import QueueAdapterBase
    return QueueAdapterBase.from_file(qadapter_file)


def load_fworker(fworker_file, **kwargs):  # TODO build a FWorker from kinbot parameters
    """Load a FWorker from a yaml file or build one from KinBot input parameters

    @param fworker_file: Path to the yaml file containing the settings for Fireworks
        to manage the worker carrying out the calculations.
    """
    from fireworks import FWorker
    if os.path.isfile(fworker_file):
        pass
    elif os.path.isfile("my_fworker.yaml"):
        fworker_file = "my_fworker.yaml"
    try:
        return FWorker.from_file(fworker_file)
    except (FileNotFoundError, ValueError) as e:
        return FWorker()


def place_job(lpad, cmd, wf_name="KinBot workflow"):
    """Place a jobn in the LaunchPad server for a FireWorker to execute.

    @param lpad: LaunchPad object to manage Workflows, FireWorks and tasks.
    @param cmd: Command to execute in the FireWorks task.
    @param wf_name: Name of the Task to be placed.

    """
    from fireworks import Workflow, Firework, ScriptTask
    firetask = ScriptTask.from_str(cmd)
    firework = Firework(firetask)
    firework.spec['_launch_dir'] = os.getcwd()
    firework.name = os.getpid()
    workflow = Workflow([firework], name=wf_name)
    lpad.add_wf(workflow)


def setup_fireworks(lpad_file='my_launchpad.yaml',
                    fworker_file='my_fworker.yaml',
                    qadapter_file='my_qadapter.yaml',
                    **kwargs):
    """

    """
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
