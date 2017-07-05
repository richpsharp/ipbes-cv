"""Task graph framework."""
import traceback
import os
import math
import logging
import multiprocessing
import traceback
import sys
import threading

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('ipbes-cv')


class Task(object):
    """Encapsulates work/task state for multiprocessing."""

    def __init__(
            self, func, args, expected_output_path_list, dependent_task_list,
            global_lock_dict, global_lock):
        """Make a task.

        Parameters:
            func (function): a function that takes the argument list
                `args`
            args (tuple): a list of arguments to pass to `func`.
            expected_output_path_list (list): a list of strings
                representing expected file path outputs.
            dependent_task_list (list of function): a list of other
                functions class by `Task` that contain an .execute()
                method: These are all invoked before `func(args)` is
                invoked.
            global_lock_dict (dict): maps task id to a lock.  If task is
                scheduled it exists in this dict and can be synchronized
                with this lock.
        """
        self.func = func
        self.args = args
        self.expected_output_path_list = expected_output_path_list
        self.dependent_task_list = dependent_task_list
        self.global_lock_dict = global_lock_dict
        self.global_lock = global_lock
        with self.global_lock:
            if not hasattr(Task, 'next_task_id'):
                Task.next_task_id = 0
            self.task_id = Task.next_task_id
            Task.next_task_id += 1

    def __call__(self):
        """Invoke this method when ready to execute task."""
        # if this Task is currently running somewhere, then wait for it.
        LOGGER.debug("Entering task %s", self.task_id)
        with self.global_lock:
            if self.task_id in self.global_lock_dict:
                LOGGER.debug("Waiting for task %s", self.task_id)
                # it's already scheduled, just block until it's done
                task_lock = self.global_lock_dict[self.task_id]
                with task_lock:
                    return
            else:
                # otherwise this is the worker thread, register a lock and
                # acquire it
                task_lock = threading.Lock()
                task_lock.acquire()
                self.global_lock_dict[self.task_id] = task_lock

        LOGGER.debug("Starting task %s", self.task_id)
        # finally block releases `task_lock` and de-registers it
        try:
            if all(os.path.exists(p) for p in self.expected_output_path_list):
                LOGGER.info(
                    "All expected files exist for %s so not executing",
                    self.func.__name__)
                return

            # Otherwise execute dependencies
            thread_list = []
            dependant_lock_list = []
            for task in self.dependent_task_list:
                LOGGER.debug("Checking dependent task %d", task.task_id)
                with self.global_lock:
                    if task.task_id in self.global_lock_dict:
                        # If the task is already got a thread assigned to it,
                        # then wait for that thread to terminate
                        LOGGER.debug(
                            "%d is a task in the global lock dict",
                            task.task_id)
                        dependant_lock_list.append(
                            self.global_lock_dict[task.task_id])
                    else:
                        # No thread exists yet, regster to start outside lock
                        thread = threading.Thread(
                            target=task, name=task.task_id)
                        thread_list.append(thread)

            # start the threads
            LOGGER.debug("Starting dependent threads %d", self.task_id)
            for thread in thread_list:
                thread.start()

            # Block on all the dependent locks
            LOGGER.debug("Waiting for dependent locks %d", self.task_id)
            for lock in dependant_lock_list:
                lock.acquire()
                lock.release()

            # wait for all the dependent threads to quit
            LOGGER.debug("Joining dependant threads %d", self.task_id)
            for thread in thread_list:
                thread.join()

            # Do this task's work
            self.func(*self.args)
        finally:
            with self.global_lock:
                del self.global_lock_dict[self.task_id]
                task_lock.release()
