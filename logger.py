import os
import logging

from datetime import datetime
from pathlib import Path

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


class CustomLogger(logging.Logger):
    
    def __init__(self, name):
        super().__init__(name)
        self._handlers = {}
    def __call__(self, handler_name):
        return self._handlers[handler_name]

class MultiLogger():
    
    def __init__(self):
        self.loggers = {}
    
    def _add_logger(self, logger, logger_name):
        self.loggers[logger_name] = logger
        setattr(self, logger_name, logger)
        return
    
    def _remove_logger(self, logger_name):
        del self.loggers[logger_name]
        delattr(self, logger_name)
        return
    
    def add_logger(self, logger: logging.Logger = None, logger_name: str = None):
        """
        accepts Logger obj or logger_name (str)
        """
        ## check that at least logger or logger_name is provided
        if logger is None and logger_name is None:
            raise Exception("Logger object or name of Logger object to be created required")
        ## check that logger_name is not already in use or reserved
        elif logger_name is not None:
            if logger_name in self.loggers: raise Exception(f"'{logger_name}' is already in use")
            if logger_name == "loggers": raise Exception(f"'{loggers}' is a reserved name")
        ## if only Logger object provided, use Logger.name to index
        elif logger is not None and logger_name is None:
            logger_name = logger.name
            if logger_name is in self.loggers:
                raise Exception(f"'{logger_name}' is already in use. Please rename the Logger object.")
            else:
                self._add_logger(logger, logger_name)
        ## if only logger_name is provided, create Logger object from 
        elif logger is None and logger_name is not None:
            logger = logging.getLogger(logger_name)
            self._add_logger(logger, logger_name)
        ## if both Logger object and logger_name are provided, use logger_name to index regardless of Logger.name
        else: self._add_logger(logger, logger_name)
        return
    
    def remove_logger(self, logger):
        """
        accepts both logger_name (str) and logger obj (Logger)
        """
        ## if logger_name provided
        if isinstance(logger, str):
            if logger in self.loggers: self._remove_logger(logger)
            else: raise Exception(f"Cannot remove logger: '{logger}' is not a known Logger name")
        ## if Logger object provided
        elif isinstance(logger, logging.Logger):
            ## if logger obj is in self.loggers, retrieve key
            inv_loggers = dict((v, k) for k, v in self.loggers.items())
            if logger in inv_loggers: self._remove_logger(inv_loggers[logger])
            else: raise Exception( ("Cannot remove logger:  The Logger object provided is not"
                                    " known to this MultiLogger") )
        return
        
