"""Task graph framework."""
import os
import math
import logging
import multiprocessing
import traceback
import sys

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('ipbes-cv')


def worker(input_queue):
    """A worker function for distributed processing.

    Parameters:
        input_queue (multiprocessing.JoinableQueue): work queue containing
            either a (func, args) pair, or the token 'STOP'.  If the former
            then `func(*args)` is invoked, a try/except catches any errors
            and logs them, and `input_queue.task_complete` is invoked.
            Otherwise if 'STOP', `input_queue.task_complete` is invoked and
            the function returns.

    Returns:
        None."""
    for func, args in iter(input_queue.get, 'STOP'):
        try:
            func(*args)
        except Exception as exception:
            LOGGER.error(
                "".join(traceback.format_exception(*sys.exc_info())))
            LOGGER.error(exception)
        input_queue.task_done()
    input_queue.task_done()


class Task(object):
    """Encapsulates work/task state for multiprocessing."""

    def __init__(
            self, func, args, expected_output_path_list, dependant_task_list):
        """Make a task.

            Parameters:
                func (function): a function that takes the argument list
                    `args`
                args (tuple): a list of arguments to pass to `func`.
                expected_output_path_list (list): a list of strings
                    representing expected file path outputs.
                dependant_task_list (list of function): a list of other
                    functions class by `Task` that contain an .execute()
                    method: These are all invoked before `func(args)` is
                    invoked.
        """
        self.func = func
        self.args = args
        self.expected_output_path_list = expected_output_path_list
        self.dependant_task_list = dependant_task_list

    def __call__(self):
        """Invoke this method when ready to execute task."""
        if all(os.path.exists(p) for p in self.expected_output_path_list):
            LOGGER.info(
                "All expected files exist for %s so not executing",
                self.func.__name__)
            return

        for task in self.dependant_task_list:
            task()
        self.func(*self.args)
