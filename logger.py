import os
import logging

from datetime import datetime
# from pathlib import Path

## TODO: use more specific Exceptions


## notes after experimentation:
## - multiple loggers can write to the same file
##  - take advantage of this to format different types of things to log
##   - we don't want everything to be prefixed w/ logging level, but the option would be nice esp for warnings
##   - e.g. logger_args can handle logging arguments, logger_map can log index-query mapping, etc.

# ## how to use CustomFileHandler w/ MultiLogger:
# mlogger = MultiLogger()
# mlogger.add_logger("subcmd1")
# mlogger.subcmd1.addHandler(DynamicFileHandler("/mnt/chaelab/rachelle/test/loggertrial.log"))
# mlogger.subcmd1.basicConfig(level = logging.DEBUG, format = '%(levelname)s:%(message)s')

# mlogger.add_logger("subcmd2")

# mlogger.subcmd1.error("OH NO")

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
            p = '.'.join(self.baseFilename.split('.')[:-1]) + str(int(datetime.utcnow().timestamp())) + ".log"
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
        if check_path: self._check_path()
        # self.close() ## close to avoid potential problems when moving logfile to this location in next step
        ## move written logfile
        if os.path.exists(p_old): os.rename(p_old, self.baseFilename)
        ## reopen moved logfile for further logging
        self._open()
        return


class MultiChild():
    
    # def __init__(self, child_adj = "child", child_obj = "Child", child_class = object,
    #              make_child = lambda child_name: child_name):
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
                raise Exception(f"{self._child_obj} name '{child_name}' is already in use.")
        elif child_name is not None and child_name in dir(self):
            raise Exception(f"'{child_name}' is a reserved name")
        elif self._unique and child in self.children():
            raise Exception( (f"Repeated {self._child_obj} objects not allowed."
                              f" This {self._child_obj} object has already been indexed as"
                              f" '{self.get_child_name()[0]}'") )
        self._children[child_name] = child
        setattr(self, child_name, child)
        return
    
    def _remove_child(self, child_name):
        del self._children[child_name]
        delattr(self, child_name)
        return
    
    def children(self):
        return list(self._children.values())
    
    def children_names(self):
        return list(self._children.keys())
    
    def children_map(self):
        return self._children
    
    def is_child_class(self, other):
        return isinstance(other, self._child_class)
    
    def get_child(self, child):
        if child in self._children: return self._children[child]
        elif child in self._children.values(): return child
        else: raise Exception( (f"'{child}' is not a known {self._child_obj} name or object"
                                " Use DynamicFileLogger.add_<file/stream>_handler"
                                " to add it to known Handlers.") )
    
    def get_child_name(self, child):
        names = [k for k, v in self._children.items() if v is child]
        if self._unique: return names[0]
        return names
    
    def is_child(self, other):
        try: tmp = self.get_child(other)
        except Exception: return False
        return True
    
    def add_child(self, child = None, child_name: str = None, replace = False, **kwargs):
        """
        accepts Child obj or child_name (str)
        """
        ## check that at least child or child_name is provided
        if child is None and child_name is None:
            raise Exception( (f"{self._child_obj} object or name of {self._child_obj} object"
                              " to be created required") )
        # ## check that child_name is not already in use or reserved
        # elif child_name is not None and child_name in self._children:
        #     raise Exception(f"{self._child_obj} name '{child_name}' is already in use")
        # elif child_name is not None and child_name in dir(self):
        #     raise Exception(f"'{children}' is a reserved name")
        ## if only Child object provided, use Child.name to index
        elif child is not None and child_name is None:
            try:
                child_name = child.name
                if child_name in self._children:
                    raise Exception( (f"{self._child_obj} name '{child_name}' is already in use."
                                      f" Please rename the {self._child_obj} object.") )
                else:
                    self._add_child(child, child_name, replace = replace)
                    return
            except AttributeError:
                raise Exception( (f"{self._child_adj}_name is required if {self._child_obj} object"
                                  " does not have attribute 'name'.") )
        ## if only child_name is provided, create Child object and return it
        elif child is None and child_name is not None:
            if self._make_child is None:
                raise Exception( (f"Function to make child object is not defined.") )
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
            else: raise Exception( (f"Cannot remove {self._child_obj}: '{child}' is not a"
                                    f" known {self._child_obj} name") )
        ## if Child object provided
        else:
            ## if Child obj is in self.children, retrieve key
            inv_children = dict((v, k) for k, v in self._children.items())
            if child in inv_children: self._remove_child(inv_children[child])
            else: raise Exception( (f"Cannot remove {self._child_obj}: The {self._child_obj} object"
                                    f" provided is not known to this {class_name}") )
        return
    

class PublicLogger(logging.Logger, MultiChild):
    """
    Logger object, but it adds Handler objects as attributes for easy access
    """
    
    def __init__(self, name):
        MultiChild.__init__(self, child_adj = "handler", child_obj = "Handler")
        logging.Logger.__init__(self, name)
    
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
        elif not self.is_child_class(handler): raise Exception("Unknown object class")
        super().remove_child(handler)
        super().removeHandler(handler)
        return

class MultiLogger(MultiChild):
    
    def __init__(self):
        # super().__init__(child_adj = "logger", child_obj = "Logger", child_class = logging.Logger,
        #                  make_child = lambda name: logging.getLogger(name)):
        super().__init__(child_adj = "logger", child_obj = "Logger", child_class = PublicLogger,
                         make_child = lambda name: PublicLogger(name))
    
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
            raise Exception("Unexpected 1st argument type. Must be Logger object or str.")
    
    # def add_logger(self, untyped_logger = None, logger: logging.Logger = None, logger_name: str = None):
    #     if untyped_logger is not None and (logger is not None or logger_name is not None):
    #         raise Exception( (f":(") )
    #     elif isinstance(untyped_logger, str): logger_name = untyped_logger
    #     elif untyped_logger is not None: logger = untyped_logger
    #     return super().add_child(child = logger, child_name = logger_name)
    
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
    """
    
    def __init__(self, filename, **config_defaults):
        super().__init__()
        # self._groups = {} ## {<logger object name>: [list of Handler obj associated w/ logger object]}
        ## stores preset formats
        self.format = MultiChild(child_adj = "formatter", child_obj = "Formatter",
                                 child_class = logging.Formatter,
                                 make_child = lambda name, format: logging.Formatter(format))
        ## stores resuable handlers
        self.handler = MultiChild(child_adj = "handler", child_obj = "Handler", unique = True)
        self._private_handler = DynamicFileHandler(filename, check_path = True)
        return
    
    def filename(self):
        return self._private_handler.baseFilename
    
    def _add_handler(self, name, handler, format = None, level = None, replace = False):
        self.handler._add_child(child = handler, child_name = name, replace = replace)
        if isinstance(format, logging.Formatter):
            handler.setFormatter(format)
        elif isinstance(format, str) and self.format.is_child(format):
            handler.setFormatter(self.format.get_child(format))
        elif isinstance(format, str):
            handler.setFormatter(logging.Formatter(format))
        else:
            raise Exception(f"'{format}' is an unknown format")
        if level is not None:
            handler.setLevel(level)
        return
    
    def _add_handler_to_group(self, group, handler):
        group_logger = self.get_logger(group)
        group_name = group_logger.name
        handler_obj = self.get_handler(handler)
        handler_name = self.handler.get_child_name(handler_obj)
        # self._groups[group_name].append(handler_obj)
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
        new_handler = DynamicFileHandler(self.filename(), check_path = False)
        self._add_handler(name, new_handler, format = format, level = level, replace = replace)
        return
    
    def add_stream_handler(self, name, format = None, level = None, replace = False):
        new_handler = logging.StreamHandler()
        self._add_handler(name, new_handler, format = format, level = level, replace = replace)
        return
    
    def get_handler(self, handler):
        try: return(self.handler.get_child(handler))
        except Exception as e:
            # print( ( f"Handler ({handler}) is not known."
            #          " Use DynamicFileLogger.add_<file/stream>_handler"
            #          " to add it to known Handlers." ) )
            raise e
        return
        
    def add_group(self, name, *handlers):
        """
        Positional 'handlers' for Handlers already added to DynamicFileLogger object.
        """
        if self.is_child(name): raise Exception(f"'{name}' is an existing group.")
        # self._groups[name] = []
        new_logger = self.add_logger(name)
        ## parse existing handlers
        for i, handler in enumerate(handlers):
            self._add_handler_to_group(name, handler)
        return
    
    def remove_group(self, name):
        if not self.is_child(name): print(f"'{name}' is not an existing group. Ignoring command.")
        self.remove_logger(name)
        # del self._groups[name]
        return
    
    def update_group(self, group, add = [], remove = []):
        """
        Add or remove handlers from group
        """
        logger = self.get_child(group)
        group = logger.name
        for handler in add:
            self._add_handler_to_group(group, handler)
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
        self._private_handler.update_filename(filename, check_path = True)
        file_handlers = [h for h in self.handler.children() if isinstance(h, DynamicFileHandler)]
        ## get new log file location (raise check_path only for first Handler processed)
        for file_handler in file_handlers:
            file_handler.update_filename(self._private_handler.baseFilename, check_path = False)
        return

## test
logger = DynamicFileLogger("/mnt/chaelab/rachelle/tmp/testlog1.log")
logger.add_format(name = "cmd", format = "command: %(message)s")
logger.add_format(name = "timed", format = "%(asctime)s %(message)s")
logger.add_format(name = "full", format = "%(asctime)s - %(name)s - %(levelname)s:%(message)s")
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

logger.add_group("header", "fheader")
logger.update_group("header", add = ["sheader"])
# ## raise error here
# logger.update_group("header", add = ["ohno"])

logger.add_group("info", logger.handler.sinfo, "finfo")
logger.add_group("all", *logger.handler.children())

logger.header.info("executing 'minimumset'")
logger.info.info("hello we're executing something")
logger.all.warning("hey stop and pay attention")
logger.all.debug("it's a mess here is what it is")

logger.update_filename("/mnt/chaelab/rachelle/tmp/testlog2.log")

# class MultiFormatLogger(MultiLogger):
    
#     def __init__(self, **config_defaults):
#         super().__init__()
#         self._template = self._make_child("template")
#         self._template.config
#         self._config_defaults = config_defaults
    
#     def _add_format(self, name, format, **kwargs):
#         new_logger = self._make_child(name)
        
    
#     def add_format(self, *args, **kwargs):
        

    
# class GroupLogger:
    
#     def __init__(self, name):
#         self._name = name
#         self._loggers = []

# class MultiGroupLogger(MultiChild):

# ## something like multilogger? idk if we need to handle the, uh, handlers separately if we're already handling the loggers separately. we can just keep making new loggers
# class CustomLogger(logging.Logger):
    
#     def __init__(self, name):
#         super().__init__(name)
#         self._handlers = {}
#     def __call__(self, handler_name):
#         return self._handlers[handler_name]
    

# class MultiLogger():
    
#     def __init__(self):
#         self.loggers = {}
    
#     def _add_logger(self, logger, logger_name):
#         self.loggers[logger_name] = logger
#         setattr(self, logger_name, logger)
#         return
    
#     def _remove_logger(self, logger_name):
#         del self.loggers[logger_name]
#         delattr(self, logger_name)
#         return
    
#     def add_logger(self, logger: logging.Logger = None, logger_name: str = None):
#         """
#         accepts Logger obj or logger_name (str)
#         """
#         ## check that at least logger or logger_name is provided
#         if logger is None and logger_name is None:
#             raise Exception("Logger object or name of Logger object to be created required")
#         ## check that logger_name is not already in use or reserved
#         elif logger_name is not None:
#             if logger_name in self.loggers: raise Exception(f"'{logger_name}' is already in use")
#             if logger_name == "loggers": raise Exception(f"'{loggers}' is a reserved name")
#         ## if only Logger object provided, use Logger.name to index
#         elif logger is not None and logger_name is None:
#             logger_name = logger.name
#             if logger_name is in self.loggers:
#                 raise Exception(f"'{logger_name}' is already in use. Please rename the Logger object.")
#             else:
#                 self._add_logger(logger, logger_name)
#         ## if only logger_name is provided, create Logger object from 
#         elif logger is None and logger_name is not None:
#             logger = logging.getLogger(logger_name)
#             self._add_logger(logger, logger_name)
#         ## if both Logger object and logger_name are provided, use logger_name to index regardless of Logger.name
#         else: self._add_logger(logger, logger_name)
#         return
    
#     def remove_logger(self, logger):
#         """
#         accepts both logger_name (str) and logger obj (Logger)
#         """
#         ## if logger_name provided
#         if isinstance(logger, str):
#             if logger in self.loggers: self._remove_logger(logger)
#             else: raise Exception(f"Cannot remove logger: '{logger}' is not a known Logger name")
#         ## if Logger object provided
#         elif isinstance(logger, logging.Logger):
#             ## if logger obj is in self.loggers, retrieve key
#             inv_loggers = dict((v, k) for k, v in self.loggers.items())
#             if logger in inv_loggers: self._remove_logger(inv_loggers[logger])
#             else: raise Exception( ("Cannot remove logger:  The Logger object provided is not"
#                                     " known to this MultiLogger") )
#         return
        
