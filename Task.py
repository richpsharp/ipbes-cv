"""Task graph framework."""
import datetime
import hashlib
import json
import pickle
import os
import logging
import multiprocessing
import threading
import errno

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('ipbes-cv')


class TaskGraph(object):
    """Encapsulates the woker and tasks states for parallel processing."""

    def __init__(self, token_storage_path, n_workers):
        """Create a task graph.

        Creates an object for building task graphs, executing them,
        parallelizing independent work notes, and avoiding repeated calls.

        Parameters:
            token_storage_path (string): path to a directory where work tokens
                (files) can be stored.  Task graph checks this directory to
                see if a task has already been completed.
            n_workers (int): number of parallel workers to allow during
                task graph execution.
        """
        # https://stackoverflow.com/questions/273192/how-can-i-create-a-directory-if-it-does-not-exist
        try:
            os.makedirs(token_storage_path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        self.token_storage_path = token_storage_path
        self.worker_pool = multiprocessing.Pool(n_workers)

        # used to lock global resources
        self.global_lock = threading.Lock()
        # if a Task is in here, it's currently executing
        self.global_working_task_dict = {}

        # list of tasks for the graph
        self.task_list = []

    def add_task(
            self, target=None, args=None, kwargs=None,
            dependent_task_list=None):
        """Add a task to the task graph.

        Parameters:
            target (function): target function
            args (list): argument list for `target`
            kwargs (dict): keyword arguments for `target`
            dependent_task_list (list): list of `Task`s that this task is
                dependent on.

        Returns:
            Task which was just added to the graph.
        """
        if args is None:
            args = []
        if kwargs is None:
            kwargs = {}
        if dependent_task_list is None:
            dependent_task_list = []
        if target is None:
            target = lambda: None
        self.task_list.append(
            Task(target, args, kwargs, dependent_task_list,
                 self.token_storage_path))
        return self.task_list[-1]

    def execute(self):
        """Execute the task graph.

        At the top level, function iterates through all tasks and triggers a
        `__call__` on each Task with an empty dependency graph.

        Returns:
            None
        """
        thread_list = []
        for task in self.task_list:
            thread = threading.Thread(
                target=task, args=(
                    self.global_lock,
                    self.global_working_task_dict,
                    self.worker_pool))
            thread_list.append((thread, task))
            thread.start()
        for thread, task in thread_list:
            thread.join()
            if not task.is_complete():
                raise RuntimeError(
                    "Task %s didn't complete, discontinuing task graph." % (
                        task.task_id))


class Task(object):
    """Encapsulates work/task state for multiprocessing."""

    def __init__(
            self, target, args, kwargs, dependent_task_list,
            token_storage_path):
        """Make a task.

        Parameters:
            target (function): a function that takes the argument list
                `args`
            args (tuple): a list of arguments to pass to `target`.  Can be
                None.
            kwargs (dict): keyword arguments to pass to `target`.  Can be
                None.
            dependent_task_list (list of Task): a list of other
                `Task`s that are to be invoked before `target(args)` is
                invoked.
            token_storage_path (string): path to a directory that exists
                where task can store a file to indicate completion of task.
        """
        self.target = target
        self.args = args
        self.kwargs = kwargs
        self.dependent_task_list = dependent_task_list

        # Make a unique hash of the input parameters of the function call
        task_string = '%s:%s:%s' % (
            target.__name__, pickle.dumps(args),
            json.dumps(kwargs, sort_keys=True))
        self.task_id = '%s_%s' % (
            target.__name__, hashlib.sha1(task_string).hexdigest())

        # https://stackoverflow.com/questions/273192/how-can-i-create-a-directory-if-it-does-not-exist
        try:
            os.makedirs(token_storage_path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        # The following file will be written when work is complete
        self.token_path = os.path.join(token_storage_path, self.task_id)

        # used to lock execution of the current Task
        self.task_lock = threading.Lock()

    def __call__(
            self, global_lock, global_working_task_dict,
            global_worker_pool):
        """Invoke this method when ready to execute task.

        Parameters:
            global_lock (threading.Lock): use this to lock global
                the global resources to the task graph.
            global_working_task_dict (dict): contains a dictionary of task_ids
                to Tasks that are currently executing.  Global resource and
                should acquire lock before modifying it.
            global_worker_pool (multiprocessing.Pool): a process pool used to
                execute subprocesses.

        Returns:
            None
        """
        # if this Task is currently running somewhere, then wait for it.
        LOGGER.debug("Entering task %s", self.task_id)
        wait_for_task = False
        with global_lock:
            if self.task_id not in global_working_task_dict:
                # acquire our own lock so nobody else runs and register in
                # the global task_dict
                self.task_lock.acquire()
                global_working_task_dict[self.task_id] = self
            else:
                wait_for_task = True
        if wait_for_task:
            # Task is running somewhere else, block on this call in case
            # something else is waiting for it and return when lock is
            # released
            LOGGER.debug(
                "Task is running somewhere else, waiting for it to finish %s",
                self.task_id)
            with self.task_lock:
                return

        LOGGER.debug("Starting task %s", self.task_id)
        # finally block below releases `self.task_lock` and deregisters it
        try:
            if self.is_complete():
                LOGGER.info(
                    "Completion token exists for %s so not executing",
                    self.task_id)
                return

            # Otherwise execute dependencies
            pending_dependent_thread_list = []
            pending_dependent_tasks = []
            with global_lock:
                for task in self.dependent_task_list:
                    LOGGER.debug("Checking dependent task %s", task.task_id)
                    if task.task_id in global_working_task_dict:
                        # task is already executing, wait for it to terminate
                        LOGGER.debug(
                            "%s is a task in the global_working_task_dict",
                            task.task_id)
                        pending_dependent_tasks.append(task)
                    else:
                        # No thread exists yet, register to start outside lock
                        thread = threading.Thread(
                            target=task,
                            args=(global_lock, global_working_task_dict))
                        global_working_task_dict[self.task_id] = self
                        pending_dependent_thread_list.append((thread, task))

            # start the threads
            if len(pending_dependent_thread_list) > 0:
                LOGGER.debug("Starting dependent threads %s", self.task_id)
                for thread, _ in pending_dependent_thread_list:
                    thread.start()

            # Block on all the dependent locks
            if len(pending_dependent_tasks) > 0:
                LOGGER.debug("Waiting for dependent locks %s", self.task_id)
                for task in pending_dependent_tasks:
                    with task.task_lock:
                        pass
                    if not task.is_complete():
                        raise RuntimeError(
                            "Task %s didn't complete, discontinuing "
                            "execution of %s" % (task.task_id, self.task_id))

            # wait for all the dependent threads
            if len(pending_dependent_thread_list) > 0:
                LOGGER.debug("Joining dependent threads %s", self.task_id)
                for thread, task in pending_dependent_thread_list:
                    thread.join()
                    if not task.is_complete():
                        raise RuntimeError(
                            "Task %s didn't complete, discontinuing "
                            "execution of %s" % (task.task_id, self.task_id))

            # Do this task's work
            LOGGER.debug("Starting process for %s", self.task_id)
            result = global_worker_pool.apply_async(
                func=self.target, args=self.args, kwds=self.kwargs)
            result.get()
            with open(self.token_path, 'w') as token_file:
                token_file.write(str(datetime.datetime.now()))
        finally:
            with global_lock:
                del global_working_task_dict[self.task_id]
            self.task_lock.release()

    def is_complete(self):
        """Return true if complete token exists."""
        return os.path.exists(self.token_path)
