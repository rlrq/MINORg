import os
import logging

from datetime import datetime
# from pathlib import Path

## notes after experimentation:
## - multiple loggers can write to the same file
##  - take advantage of this to format different types of things to log
##   - we don't want everything to be prefixed w/ logging level, but the option would be nice esp for warnings
##   - e.g. logger_args can handle logging arguments, logger_map can log index-query mapping, etc.

## custom exceptions (in case we need to catch and ignore specific problems)
class InvalidArgs(Exception): pass
class InvalidName(Exception): pass
class InvalidType(Exception): pass
class DuplicateObject(Exception): pass
class UnknownChild(Exception): pass
class UnknownInitialiser(Exception): pass


## handle filename changes
class DynamicFileHandler(logging.FileHandler):
    
    def __init__(self, filename, check_path = True):
        super().__init__(filename)
        if check_path: self._check_path()
    
    ## if file already exists, append unix time (seconds) to new logfile name
    def _check_path(self):
        if os.path.exists(self.baseFilename) and os.stat(self.baseFilename).st_size > 0:
            ## close existing file
            self.close()
            ## update new logfile name
            p = ('.'.join(self.baseFilename.split('.')[:-1]) + '_' + \
                 str(int(datetime.utcnow().timestamp())) + ".log")
            self.baseFilename = p
            ## open logfile in updated location
            self._open()
        return
    
    ## moves log file that has been written to new location
    def update_filename(self, filename, check_path = True):
        p_old = self.baseFilename
        ## close current logfile
        self.close()
        ## update logfile location and check availability
        self.baseFilename = os.path.abspath(filename)
        if check_path and p_old != filename: self._check_path()
        # self.close() ## close to avoid potential problems when moving logfile to this location in next step
        ## move written logfile
        if os.path.exists(p_old): os.rename(p_old, self.baseFilename)
        ## reopen moved logfile for further logging
        self._open()
        return


class MultiChild():
    
    def __init__(self, child_adj = "child", child_obj = "Child", child_class = None,
                 make_child = None, unique = False):
        self._unique = unique
        self._children = {}
        self._child_adj = child_adj
        self._child_obj = child_obj
        self._child_class = child_class
        self._make_child = make_child
    
    def _add_child(self, child, child_name, replace = False):
        if self.is_child(child_name):
            if replace:
                print (f"{self._child_obj} name '{child_name}' is already in use."
                       f" Updating with new {self._child_obj} object.")
            else:
                raise InvalidName(f"{self._child_obj} name '{child_name}' is already in use.")
        elif child_name is not None and child_name in dir(self):
            raise InvalidName(f"'{child_name}' is a reserved name for class {self.__class__.__name__}")
        elif self._unique and child in self.children:
            raise DuplicateObject( (f"Repeated {self._child_obj} objects not allowed."
                                    f" This {self._child_obj} object has already been indexed as"
                                    f" '{self.get_child_name(child)[0]}'") )
        self._children[child_name] = child
        setattr(self, child_name, child)
        return
    
    def _remove_child(self, child_name):
        del self._children[child_name]
        delattr(self, child_name)
        return
    
    @property
    def children(self): return list(self._children.values())
    @property
    def children_names(self): return list(self._children.keys())
    @property
    def children_map(self): return self._children
    
    def is_child_class(self, other):
        return isinstance(other, self._child_class)
    
    def get_child(self, child):
        if child in self.children_names: return self.children_map[child]
        elif child in self.children: return child
        else: raise UnknownChild( (f"'{child}' is not a known {self._child_obj} name or object"
                                   " Use DynamicFileLogger.add_<file/stream>_handler"
                                   " to add it to known Handlers.") )
    
    def get_child_name(self, child):
        names = [k for k, v in self.children_map.items() if v is child]
        if self._unique: return names[0]
        return names
    
    def is_child(self, other):
        try: tmp = self.get_child(other)
        except UnknownChild: return False
        return True
    
    def add_child(self, child = None, child_name: str = None, replace = False, **kwargs):
        """
        accepts Child obj or child_name (str)
        """
        ## check that at least child or child_name is provided
        if child is None and child_name is None:
            raise InvalidArgs( (f"{self._child_obj} object or name of {self._child_obj} object"
                                " to be created required") )
        ## if only Child object provided, use Child.name to index
        elif child is not None and child_name is None:
            try:
                child_name = child.name
            except AttributeError:
                raise InvalidArgs( (f"{self._child_adj}_name is required if {self._child_obj} object"
                                    " does not have attribute 'name'.") )
            if child_name in self._children:
                raise InvlaidName( (f"{self._child_obj} name '{child_name}' is already in use."
                                    f" Please rename the {self._child_obj} object.") )
            else:
                self._add_child(child, child_name, replace = replace)
                return
        ## if only child_name is provided, create Child object and return it
        elif child is None and child_name is not None:
            if self._make_child is None:
                raise UnknownInitialiser( (f"Function to make child object is not defined.") )
            child = self._make_child(child_name, **kwargs)
            self._add_child(child, child_name, replace = replace)
            return child
        ## if both Child object and child_name are provided, use child_name to index regardless of Child.name
        else: self._add_child(child, child_name, replace = replace)
        return
    
    def remove_child(self, child, class_name = None):
        """
        accepts both child_name (str) and Child obj
        """
        if class_name is None: class_name = self.__class__.__name__
        ## if child_name provided
        if isinstance(child, str):
            if child in self._children: self._remove_child(child)
            else: raise UnknownChild( (f"Cannot remove {self._child_obj}: '{child}' is not a"
                                       f" known {self._child_obj} name") )
        ## if Child object provided
        else:
            ## if Child obj is in self.children, retrieve key
            inv_children = dict((v, k) for k, v in self._children.items())
            if child in inv_children: self._remove_child(inv_children[child])
            else: raise UnknownChild( (f"Cannot remove {self._child_obj}: The {self._child_obj} object"
                                       f" provided is not known to this {class_name}") )
        return
    

class PublicLogger(logging.Logger, MultiChild):
    """
    Logger object, but it adds Handler objects as attributes for easy access
    """
    
    def __init__(self, name):
        MultiChild.__init__(self, child_adj = "handler", child_obj = "Handler")
        logging.Logger.__init__(self, name)
    
    # @property
    # def handlers(self): return self.children
    @property
    def handler_names(self): return self.children_names
    @property
    def handler_map(self): return self.children_map
    
    def addHandler(self, handler, handler_name, replace = False):
        """
        Accepts a Handler object (e.g. StreamHandler, FileHandler, DynamicFileHandler)
        """
        super()._add_child(handler, handler_name, replace = replace)
        super().addHandler(handler)
        return
    
    def removeHandler(self, handler):
        ## get Handler object if handler name (str) provided
        if isinstance(handler, str): handler = self.get_child(handler)
        elif not self.is_child_class(handler): raise InvalidType("Unknown object class")
        super().remove_child(handler)
        super().removeHandler(handler)
        return

class MultiLogger(MultiChild):
    
    def __init__(self):
        super().__init__(child_adj = "logger", child_obj = "Logger", child_class = PublicLogger,
                         make_child = lambda name: PublicLogger(name))
    
    @property
    def loggers(self): return self.children
    @property
    def logger_names(self): return self.children_names
    @property
    def logger_map(self): return self.children_map
    
    def is_logger(self, other):
        return self.is_child_class(other)
    
    def add_logger(self, logger, name: str = None, replace = False):
        if self.is_logger(logger):
            return super().add_child(child = logger, child_name = name, replace = replace)
        elif isinstance(logger, str) and name is None:
            return super().add_child(child_name = logger, replace = replace)
        elif isinstance(logger, str) and name is not None:
            print("Both arguments are strings. Ignoring 2nd argument.")
            return super().add_child(child_name = logger, replace = replace)
        else:
            raise InvalidArgs("Unexpected 1st argument type. Must be Logger object or str.")
        
    def remove_logger(self, logger):
        return super().remove_child(logger, class_name = self.__class__.__name__)
    
    def get_logger(self, logger):
        return self.get_child(logger)


class DynamicFileLogger(MultiLogger):
    """
    User specifies a group when passing message
    (E.g. DynamicFileLogger.cmdname('minimumset') and the relevant handlers will do their job)
    (Note that this works by creating a Logger called 'cmdname' and handlers linked to that Logger 
       will do their job. Maybe a StreamHandler will not print it cuz it's not relevant, 
       but the FileHandler may wish to do so for thoroughness. There may be different formatting.)
    User may also change log file location whenever desired
    (E.g. DynamicFileLogger.update_filename('/mnt/chaelab/rachelle/logs/loggertest.log'))
    
    Groups are effectively Logger objects
    """
    
    def __init__(self, filename, **config_defaults):
        super().__init__()
        self._base_handler = DynamicFileHandler(filename, check_path = True)
        ## stores preset formats
        self.format = MultiChild(child_adj = "formatter", child_obj = "Formatter",
                                 child_class = logging.Formatter,
                                 make_child = lambda name, format: logging.Formatter(format))
        ## stores resuable handlers
        self.handler = MultiChild(child_adj = "handler", child_obj = "Handler", unique = True)
        return
    
    ## some getters for users
    @property
    def filename(self): return self._base_handler.baseFilename
    @property
    def groups(self): return self.logger_names
    @property
    def handlers(self): return self.handler.children
    @property
    def handler_names(self): return self.handler.children_names
    @property
    def handler_map(self): return self.handler.children_map
    @property
    def formats(self): return self.format.children
    @property
    def format_names(self): return self.format.children_names
    @property
    def format_map(self): return self.format.children_map
    def active_handlers(self, name = True):
        import itertools
        return list(set(itertools.chain(*[(logger.handler_names if name else logger.handlers)
                                          for logger in self.loggers])))
    def orphan_handlers(self, name = True):
        return list(set(self.handler.children_names if name else self.handler.children) -
                    set(self.active_handlers(name = name)))
    def group_handlers(self, name = True, group_name = False, handler_name = False):
        if group_name or handler_name: name = False
        if name: group_name = handler_name = True
        return {(logger_name if group_name else logger):
                (logger.handler_names if handler_name else logger.handlers)
                for logger_name, logger in self.logger_map.items()}
    
    def _add_handler(self, name, handler, format = None, level = None, replace = False):
        self.handler._add_child(child = handler, child_name = name, replace = replace)
        if isinstance(format, logging.Formatter):
            handler.setFormatter(format)
        elif isinstance(format, str) and self.format.is_child(format):
            handler.setFormatter(self.format.get_child(format))
        elif isinstance(format, str):
            handler.setFormatter(logging.Formatter(format))
        else:
            raise InvalidArgs(f"'{format}' is an unknown format")
        if level is not None:
            handler.setLevel(level)
        return
    
    def _add_handler_to_group(self, group, handler):
        group_logger = self.get_logger(group)
        group_name = group_logger.name
        handler_obj = self.get_handler(handler)
        handler_name = self.handler.get_child_name(handler_obj)
        group_logger.addHandler(handler_obj, handler_name)
        return
    
    def _remove_handler_from_group(self, group, handler):
        group_logger = self.get_logger(group)
        handler_obj = self.get_handler(handler)
        if handler_obj not in group_logger.handlers:
            print(f"'{name}' is not in group '{group_logger.name}'. Ignoring command.")
        group_logger.removeHandler(handler_obj)
        return
    
    def add_format(self, name, format, replace = False):
        self.format.add_child(child_name = name, format = format, replace = replace)
        return
        
    def add_file_handler(self, name, format = None, level = None, replace = False):
        new_handler = DynamicFileHandler(self.filename, check_path = False)
        self._add_handler(name, new_handler, format = format, level = level, replace = replace)
        return
    
    def add_stream_handler(self, name, format = None, level = None, replace = False):
        new_handler = logging.StreamHandler()
        self._add_handler(name, new_handler, format = format, level = level, replace = replace)
        return
    
    def get_handler(self, handler):
        return self.handler.get_child(handler)
        
    def add_group(self, name, *handlers):
        """
        Positional 'handlers' for Handlers already added to DynamicFileLogger object.
        """
        if self.is_child(name): raise InvalidName(f"'{name}' is an existing group.")
        new_logger = self.add_logger(name)
        ## parse existing handlers
        for i, handler in enumerate(handlers):
            self._add_handler_to_group(name, handler)
        return
    
    def remove_group(self, name):
        if not self.is_child(name): print(f"'{name}' is not an existing group. Ignoring command.")
        self.remove_logger(name)
        return
    
    def update_group(self, group, add = None, remove = None):
        """
        Add or remove handlers from group
        """
        logger = self.get_child(group)
        group = logger.name
        if add is not None:
            if isinstance(add, str) or isinstance(add, logging.Handler):
                add = [add]
            for handler in add:
                self._add_handler_to_group(group, handler)
        if remove is not None:
            if instance(remove, str) or isinstance(remove, logging.Handler):
                remove = [remove]
            for handler in remove:
                self._remove_handler_from_group(group, handler)
        return
    
    def update_filename(self, filename):
        """
        Closes connections from DynamicFileHandler to current logfile,
        replaces DynamicFileHandler.baseFilename, copies closed logfile to new location,
        and opens connection between all DynamicFileHandler objects to new location
        """
        ## update path for DynamicFileHandler, from which we will retrieve the new path for other Handlers
        self._base_handler.update_filename(filename, check_path = True)
        file_handlers = [h for h in self.handler.children if isinstance(h, DynamicFileHandler)]
        ## get new log file location (raise check_path only for first Handler processed)
        for file_handler in file_handlers:
            file_handler.update_filename(self._base_handler.baseFilename, check_path = False)
        return

class DynamicMultiFileLogger(DynamicFileLogger):
    """
    Allows writing to multiple files simultaneously, while sharing group names and Formatters
    Possibly only really useful for reusing headers such as "--Initiating <subcmd>--"
    
    dmflogger = DynamicMultiFileLogger() ## instantiate w/ dummy log location, close immediately & remove file
    dmflogger.add_file("log1", "/mnt/chaelab/rachelle/tmp/log1.log") ## check & store filename in dict
    dmflogger.add_file("log2", "/mnt/chaelab/rachelle/tmp/log2.log") ## check & store filename in dict
    dmflogger.add_format("cmd", "command: %(message)s")
    dmflogger.add_format("full", "%(asctime)s - %(name)s - %(levelname)s:%(message)s")
    dmflogger.add_stream_handler("cheader", level = logging.INFO, format = "cmd")
    dmflogger.add_file_handler("fheader", level = logging.INFO, format = "cmd") ## mk template DynamicFileHandler
    dmflogger.add_stream_handler("cfull", level = logging.INFO, format = "full")
    dmflogger.add_file_handler("ffull", level = logging.INFO, format = "full") ## mk template DynamicFileHandler
    dmflogger.add_group("header", "cheader", "fheader")
    dmflogger.add_group("full", "cfull", "ffull")
    dmflogger.add_group("all", *dmflogger.handlers) ## perhaps have a built-in 'all' group?
    dmflogger.log1.add_handler("header") ## duplicate fheader's DynamicFileHandler 
    dmflogger.log2.add_handler("header", "full") ## duplicate ffull AND fheader's DynamicFileHandlers
    dmflogger.header.info("executing minimumset")
    dmflogger.full.debug("we're doing something")
    dmflogger.log2.update_filename("/mnt/chaelab/rachelle/log3.log")
    
    IMPT: ensure no duplicate calls of StreamHandlers for groups w/ StreamHandler that are
      associated with multiple files
    - each 'file' is represented by a DynamicFileLogger (allows file-specific logging)
      - ensure all changes to DynamicMultiFileLogger groups are propagated to each file's DynamicFileLogger
    - if logging directly from DynamicMultiFileLogger, DON"T call each file's DynamicFileLogger.error...
    - file handlers (DynamicFileHandlers) to be handled distinctly from StreamHandlers
      - when new file handlers are added, create template file handler w/ dummy log location, close, then rm
      - when file handlers are assigned to a file, duplicate handler template, update log location w/ filename, then add duplicated handler to group's MultiLogger and file's DynamicFileLogger
      - ensure that groups are also propagated to each file's DynamicFileLogger
    - when removing handlers from groups, ensure that this is propagated to group's MultiLogger AND file's DynamicFileLogger
    """
    pass
    
## test
fmt_cmd = "command: %(message)s"
fmt_timed = "%(asctime)s %(message)s"
fmt_full = "%(asctime)s - %(name)s - %(levelname)s:%(message)s"
logger = DynamicFileLogger("/mnt/chaelab/rachelle/tmp/testlog1.log")
logger.add_format(name = "cmd", format = fmt_cmd)
logger.add_format(name = "timed", format = fmt_timed)
logger.add_format(name = "full", format = fmt_full)
logger.format.full = fmt_full ## problematic assignment
logger.add_format(name = "full", format = fmt_full, replace = True)
logger.format.full._fmt == fmt_full ## True
logger.add_stream_handler("sinfo", level = logging.INFO, format = "timed")
logger.add_file_handler("finfo", level = logging.INFO, format = logger.format.full)
logger.add_stream_handler("sheader", level = logging.ERROR, format = "cmd")
logger.add_file_handler("fheader", level = logging.ERROR, format = "cmd")
logger.add_stream_handler("random", level = logging.DEBUG,
                          format = "SPECIAL %(asctime)s - %(name)s - %(levelname)s:%(message)s")
# ## raise error here
# logger.add_stream_handler("random", level = logging.DEBUG,
#                           format = "SPECIAL %(asctime)s - %(name)s - %(levelname)s:%(message)s")
logger.add_stream_handler("random", level = logging.DEBUG,
                          format = "SPECIAL %(asctime)s - %(name)s - %(levelname)s:%(message)s", replace = True)
# ## raise error here (reserved name)
# logger.add_stream_handler("get_child", "sinfo")

logger.add_group("header", "fheader")
logger.update_group("header", add = ["sheader"])
# ## raise error here (unknown Handler)
# logger.update_group("header", add = ["ohno"])

logger.add_group("info", logger.handler.sinfo, "finfo")
logger.add_group("all", *logger.handler.children)
# ## raise error here (reserved name)
# logger.add_group("get_child", "sinfo")

## abstraction breaking. (Maybe find a way to stop this? Or ensure that even if this happens
##   we can propagate it up to DynamicFileLogger so it's accessibly to other groups as well AND
##   shows up when logger.handlers/handler_names/handler_map is called (alt, add new function so user
##   can retrieve handlers that were added in this abstraction-breaking way)?)
logger.all.addHandler(logging.StreamHandler(), "breakabstraction")
logger.all.breakabstraction.setLevel(logging.ERROR)
logger.all.breakabstraction.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(message)s"))

logger.header.info("executing 'minimumset'")
logger.info.info("hello we're executing something")
logger.all.warning("hey stop and pay attention")
logger.all.debug("it's a mess here is what it is")

logger.handlers
logger.handler_names
logger.handler_map
logger.groups
logger.formats
logger.format_names
{k: fmt._fmt for k, fmt in logger.format_map.items()}
logger.active_handlers(name = True)
logger.active_handlers(name = False)
logger.orphan_handlers(name = True)
logger.group_handlers(name = True)
logger.header.handler_names

## moved from previous location to new location
logger.update_filename("/mnt/chaelab/rachelle/tmp/testlog2.log")
logger.filename

logger2 = DynamicFileLogger("/mnt/chaelab/rachelle/tmp/testlog2.log")
assert(logger2.filename != "/mnt/chaelab/rachelle/tmp/testlog2.log", "Successfully moved to avoid conflict.")
