"""Run an unchanged ROTD_py input with KinBot's shared Slurm admission.

Only the controller's submission command and local job cap are adapted.
Native calculation settings, sampling input, and existing job IDs stay intact.
"""
import runpy
import shlex
import sys

from kinbot import run_policy


def cooperative_multi(base):
    class SharedPoolMulti(base):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            if self.sub_cmd != 'sbatch':
                raise ValueError('Shared ROTD_py admission requires Slurm.')
            self.sub_cmd = shlex.join([sys.executable, '-m', 'kinbot.run_policy', 'submit'])
            self._refresh_limit()

        def _refresh_limit(self):
            # Zero is enforced by the shared admission gate. Keep ROTD's
            # local polling loop able to reconcile completed jobs meanwhile.
            self.calculator['max_jobs'] = max(1, run_policy.limit(self.calculator['max_jobs']))

        def submit_work(self, *args, **kwargs):
            self._refresh_limit()
            return super().submit_work(*args, **kwargs)

        def check_running_jobs(self, *args, **kwargs):
            self._refresh_limit()
            return super().check_running_jobs(*args, **kwargs)

    return SharedPoolMulti


def main():
    from rotd_py import new_multi
    new_multi.Multi = cooperative_multi(new_multi.Multi)
    filename = sys.argv[1]
    sys.argv = sys.argv[1:]
    runpy.run_path(filename, run_name='__main__')


if __name__ == '__main__':
    main()
