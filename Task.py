"""Task graph framework."""
import types
import collections
import traceback
import time
import datetime
import hashlib
import json
import pickle
import os
import logging
import multiprocessing
import threading
import errno
import Queue
import inspect

import psutil

LOGGER = logging.getLogger('Task')


def _worker(work_queue):
    """Thread worker.  `work_queue` has func/args tuple or 'STOP'."""
    for func, args in iter(work_queue.get, 'STOP'):
        func(*args)


class TaskGraph(object):
    """Encapsulates the worker and tasks states for parallel processing."""

    def __init__(self, token_storage_path, n_workers):
        """Create a task graph.

        Creates an object for building task graphs, executing them,
        parallelizing independent work notes, and avoiding repeated calls.

        Parameters:
            token_storage_path (string): path to a directory where work tokens
                (files) can be stored.  Task graph checks this directory to
                see if a task has already been completed.
            n_workers (int): number of parallel workers to allow during
                task graph execution.  If set to 0, use current process.
        """
        # https://stackoverflow.com/questions/273192/how-can-i-create-a-directory-if-it-does-not-exist
        try:
            os.makedirs(token_storage_path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        self.token_storage_path = token_storage_path
        self.work_queue = Queue.Queue()
        self.n_workers = n_workers
        for thread_id in xrange(n_workers):
            threading.Thread(
                target=TaskGraph.worker, args=(self.work_queue,),
                name=thread_id).start()

        if n_workers > 0:
            self.worker_pool = multiprocessing.Pool(n_workers)
            parent = psutil.Process()
            parent.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            for child in parent.children():
                child.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        else:
            self.worker_pool = None

        # used to lock global resources
        self.global_lock = threading.Lock()

        # if a Task is in here, it's been previously created
        self.global_working_task_set = set()

        self.closed = False

    def __del__(self):
        """Clean up task graph by injecting STOP sentinels."""
        self.close()

    @staticmethod
    def worker(work_queue):
        """Worker taking (func, args, kwargs) tuple from `work_queue`."""
        for func, args, kwargs in iter(work_queue.get, 'STOP'):
            try:
                func(*args, **kwargs)
            except Exception as e:
                LOGGER.error(traceback.format_exc())

    def close(self):
        """Prevent future tasks from being added to the work queue."""
        self.closed = True
        for thread_id in xrange(self.n_workers):
            self.work_queue.put('STOP')

    def add_task(
            self, target=None, args=None, kwargs=None,
            target_path_list=None,
            ignore_path_list=None,
            dependent_task_list=None,
            ignore_directories=True):
        """Add a task to the task graph.

        See the docstring for Task.__call__ to determine how it determines
        if it should execute.

        Parameters:
            target (function): target function
            args (list): argument list for `target`
            kwargs (dict): keyword arguments for `target`
            target_path_list (list): if not None, a list of file paths that
                are expected to be output by `target`.  If any of these paths
                don't exist, or their timestamp is earlier than an input
                arg or work token, target will be executed.  If None, not
                considered when scheduling task.
            ignore_path_list (list): list of file paths that could be in
                args/kwargs that should be ignored when considering timestamp
                hashes.
            dependent_task_list (list): list of `Task`s that this task must
                `join` before executing.
            ignore_directories (boolean): if the existence/timestamp of any
                directories discovered in args or kwargs is used as part
                of the work token hash.

        Returns:
            Task which was just added to the graph.
        """
        if self.closed:
            raise ValueError(
                "The task graph is closed and cannot accept more tasks.")
        if args is None:
            args = []
        if kwargs is None:
            kwargs = {}
        if dependent_task_list is None:
            dependent_task_list = []
        if target_path_list is None:
            target_path_list = []
        if ignore_path_list is None:
            ignore_path_list = []
        if target is None:
            target = lambda: None
        task = Task(
            target, args, kwargs, target_path_list, ignore_path_list,
            dependent_task_list, ignore_directories, self.token_storage_path)

        with self.global_lock:
            self.global_working_task_set.add(task)
        if self.n_workers > 0:
            self.work_queue.put(
                (task,
                 (self.global_lock,
                  self.global_working_task_set,
                  self.worker_pool),
                 {}))
        else:
            task(
                self.global_lock, self.global_working_task_set,
                self.worker_pool)
        return task

    def join(self):
        """Join all threads in the graph."""
        for task in self.global_working_task_set:
            task.join()


class Task(object):
    """Encapsulates work/task state for multiprocessing."""

    def __init__(
            self, target, args, kwargs, target_path_list, ignore_path_list,
            dependent_task_list, ignore_directories, token_storage_path):
        """Make a Task.

        Parameters:
            target (function): a function that takes the argument list
                `args`
            args (tuple): a list of arguments to pass to `target`.  Can be
                None.
            kwargs (dict): keyword arguments to pass to `target`.  Can be
                None.
            target_path_list (list): a list of filepaths that this task
                should generate.
            dependent_task_list (list of Task): a list of other
                `Task`s that are to be invoked before `target(args)` is
                invoked.
            target_path_list (list): list of filepaths expected
                to be generated by this target and args/kwargs.
            ignore_path_list (list): list of file paths that could be in
                args/kwargs that should be ignored when considering timestamp
                hashes.
            ignore_directories (boolean): if the existence/timestamp of any
                directories discovered in args or kwargs is used as part
                of the work token hash.
            token_storage_path (string): path to a directory that exists
                where task can store a file to indicate completion of task.
        """
        self.target = target
        self.args = args
        self.kwargs = kwargs
        self.target_path_list = target_path_list
        self.dependent_task_list = dependent_task_list
        self.target_path_list = target_path_list
        self.ignore_path_list = ignore_path_list
        self.ignore_directories = ignore_directories
        self.token_storage_path = token_storage_path
        self.token_path = None  # not set until dependencies are blocked

        # Used to ensure only one attempt at executing and also a mechanism
        # to see when Task is complete
        self.lock = threading.Lock()
        self.lock.acquire()  # the only release is at the end of __call__

        # https://stackoverflow.com/questions/273192/how-can-i-create-a-directory-if-it-does-not-exist
        try:
            os.makedirs(token_storage_path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    def _calculate_token(self):
        """Make a unique hash of the call. Standalone so it can be threaded."""
        try:
            if not hasattr(Task, 'target_source_map'):
                Task.target_source_map = {}
            # memoize target source code because it's likely we'll import
            # the same target many times and reflection is slow
            if self.target not in Task.target_source_map:
                Task.target_source_map[self.target] = (
                    inspect.getsource(self.target))
            source_code = Task.target_source_map[self.target]
        except IOError:
            # we might be in a frozen binary, so just leave blank
            source_code = ''

        file_stat_list = list(_get_file_stats(
            [self.args, self.kwargs],
            self.target_path_list+self.ignore_path_list,
            self.ignore_directories))

        task_string = '%s:%s:%s:%s:%s:%s' % (
            self.target.__name__, pickle.dumps(self.args),
            json.dumps(self.kwargs, sort_keys=True), source_code,
            self.target_path_list, str(file_stat_list))

        return task_string

    def __call__(
            self, global_lock, global_working_task_dict,
            global_worker_pool):
        """Invoke this method when ready to execute task.

        This function will execute `target` on `args`/`kwargs` under the
        following circumstances:
            * if no work token exists (a work token is a hashed combination
                of the source code, arguments, target files, and associated
                time stamps).
            * if any input filepath arguments have a newer timestamp than the
              work token or any path in `target_path_list`.
            * AND all the tasks in `dependant_task_list` have been joined.


        Parameters:
            global_lock (threading.Lock): use this to lock global
                the global resources to the task graph.
            global_working_task_dict (dict): contains a dictionary of task_ids
                to Tasks that are currently executing.  Global resource and
                should acquire lock before modifying it.
            global_worker_pool (multiprocessing.Pool): a process pool used to
                execute subprocesses.  If None, use current process.

        Returns:
            None
        """
        try:
            if len(self.dependent_task_list) > 0:
                for task in self.dependent_task_list:
                    task.join()
                    if not task.is_complete():
                        raise RuntimeError(
                            "Task %s didn't complete, discontinuing "
                            "execution of %s" % (
                                task.task_id, self.target.__name__))

            token = self._calculate_token()
            task_id = '%s_%s' % (
                self.target.__name__, hashlib.sha1(token).hexdigest())
            self.token_path = os.path.join(self.token_storage_path, task_id)
            LOGGER.debug("Starting task %s", task_id)

            if self.is_complete():
                LOGGER.debug(
                    "Completion token exists for %s so not executing",
                    task_id)
                return

            LOGGER.debug("Starting process for %s", task_id)
            if global_worker_pool is not None:
                result = global_worker_pool.apply_async(
                    func=self.target, args=self.args, kwds=self.kwargs)
                result.get()
            else:
                self.target(*self.args, **self.kwargs)
            LOGGER.debug("Complete process for %s", task_id)
            with open(self.token_path, 'w') as token_file:
                token_file.write(str(datetime.datetime.now()))
        finally:
            self.lock.release()

    def is_complete(self):
        """Return true if complete token and expected files exist."""
        return all([
            os.path.exists(path)
            for path in [self.token_path] + self.target_path_list])

    def join(self):
        """Block until task is complete."""
        with self.lock:
            pass


def _get_file_stats(base_value, ignore_list, ignore_directories):
    """Iterate over any values that are filepaths by getting filestats.

    Parameters:
        base_value: any python value.
        ignore_list (list): any paths found in this list are not included
            as part of the file stats
        ignore_directories (boolean): If True directories are not
            considered for filestats.

    Return:
        list of (timestamp, filesize) tuples for any filepaths found in
            base_value or nested in base value that are not otherwise
            ignored by the input parameters.
        """
    if isinstance(base_value, types.StringType):
        try:
            if base_value not in ignore_list:
                if (not isinstance(base_value, types.StringType) or
                        not os.path.isdir(base_value) or
                        not ignore_directories):
                    yield (
                        os.path.getmtime(base_value), os.path.getsize(
                            base_value), base_value)
        except OSError:
            pass
    elif isinstance(base_value, collections.Mapping):
        for key in sorted(base_value.iterkeys()):
            value = base_value[key]
            for x in _get_file_stats(value, ignore_list, ignore_directories):
                yield x
    elif isinstance(base_value, collections.Iterable):
        for value in base_value:
            for x in _get_file_stats(value, ignore_list, ignore_directories):
                yield x
