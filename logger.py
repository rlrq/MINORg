import os
import logging

from datetime import datetime
from pathlib import Path


## notes after experimentation:
## - multiple loggers can write to the same file
##  - take advantage of this to format different types of things to log
##   - we don't want everything to be prefixed w/ logging level, but the option would be nice esp for warnings
##   - e.g. logger_args can handle logging arguments, logger_map can log index-query mapping, etc.

## handle filename changes
class CustomFileHandler(logging.FileHandler):
    
    def __init__(self, filename):
        super().__init__(filename)
        self._check_path()
    
    ## if file already exists, append unix time (seconds) to new logfile name
    def _check_path(self):
        if os.path.exists(self.baseFilename):
            ## close existing file
            self.close()
            ## update new logfile name
            p = '.'.join(self.baseFilename.split('.')[:-1]) + str(int(datetime.time())) + ".log"
            self.baseFilename = p
            ## open logfile in updated location
            self._open()
        return
    
    ## moves log file that has been written to new location
    def update_filname(self, filename):
        p_old = self.baseFilename
        ## close current logfile
        self.close()
        ## update logfile location and check availability
        self.baseFilename = os.path.abspath(filename)
        self._check_path()
        self.close() ## close to avoid potential problems when moving logfile to this location in next step
        ## move written logfile
        os.rename(p_old, self.baseFilename)
        ## reopen moved logfile for further logging
        self._open()
    return


class MultiChild():
    
    def __init__(self, child_adj = "child", child_obj = "Child", child_class = object,
                 make_child = lambda child_name: child_name):
        self._children = {}
        self._child_adj = child_adj
        self._child_obj = child_obj
        self._child_class = child_class
        self._make_child = make_child
    
    def _add_child(self, child, child_name):
        self._children[child_name] = child
        setattr(self, child_name, child)
        return
    
    def _remove_child(self, child_name):
        del self._children[child_name]
        delattr(self, child_name)
        return
    
    def add_child(self, child: = None, child_name: str = None):
        """
        accepts Child obj or child_name (str)
        """
        ## check that at least child or child_name is provided
        if child is None and child_name is None:
            raise Exception( (f"{self._child_obj} object or name of {self._child_obj} object"
                              " to be created required") )
        ## check that child_name is not already in use or reserved
        elif child_name is not None:
            if child_name in self._children: raise Exception(f"'{child_name}' is already in use")
            if child_name in dir(self): raise Exception(f"'{children}' is a reserved name")
        ## if only Child object provided, use Child.name to index
        elif child is not None and child_name is None:
            try:
                child_name = child.name
                if child_name is in self._children:
                raise Exception( (f"'{child_name}' is already in use."
                                  f" Please rename the {self._child_obj} object.") )
                else:
                    self._add_child(child, child_name)
            except AttributeError:
                raise Exception( (f"{self._child_adj}_name is required if {self._child_obj} object"
                                  " does not have attribute 'name'.") )
        ## if only child_name is provided, create Child object and return it
        elif child is None and child_name is not None:
            child = self._child_class(child_name)
            self._add_child(child, child_name)
            return chld
        ## if both Child object and child_name are provided, use child_name to index regardless of Child.name
        else: self._add_child(child, child_name)
        return
    
    def remove_child(self, child, class_name = self.__class__.__name__):
        """
        accepts both child_name (str) and Child obj
        """
        ## if child_name provided
        if isinstance(child, str):
            if child in self.children: self._remove_child(child)
            else: raise Exception( (f"Cannot remove {self._child_obj}: '{child}' is not a"
                                    f" known {self._child_obj} name") )
        ## if Child object provided
        else:
            ## if Child obj is in self.children, retrieve key
            inv_children = dict((v, k) for k, v in self.children.items())
            if child in inv_children: self._remove_child(inv_children[child])
            else: raise Exception( (f"Cannot remove {self._child_obj}: The {self._child_obj} object"
                                    f" provided is not known to this {class_name}") )
        return
        

class MultiLogger(MultiChild):
    
    def __init__(self):
        super().__init__(child_adj = "logger", child_obj = "Logger", child_class = logging.Logger,
                         make_child = lambda name: logging.getLogger(name)):
    
    def add_logger(self, logger: logging.Logger = None, logger_name: str = None):
        super().add_child(child = logger, child_name = logger_name)
    
    def remove_logger(logger):
        super().remove_child(logger)


class MultiLogger(MultiChild):
    
    def __init__(self):
        super().__init__(child_adj = "logger", child_obj = "Logger", child_class = logging.Logger,
                         make_child = lambda name: logging.getLogger(name)):
    
    def add_logger(self, logger: logging.Logger = None, logger_name: str = None):
        return super().add_child(child = logger, child_name = logger_name)
    
    def remove_logger(logger):
        return super().remove_child(logger, class_name = self.__class__.__name__)


class CustomLogger(logging.Logger):
    
    def __init__(self, name):
        super().__init__(name)
        self._handlers = {}
    def __call__(self, handler_name):
        return self._handlers[handler_name]
    

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
        
