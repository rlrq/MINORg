import os
import copy
import shutil
import logging
import tempfile
import warnings
import itertools

lvl = logging.INFO

# logging.basicConfig(level=slogger.debug)
slogger = logging.Logger("slogger")
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
slogger.addHandler(ch)
slogger.setLevel(lvl)


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

class ReservedName(InvalidName): pass

class InvalidType(Exception): pass

class DuplicateObject(Exception): pass

class DuplicateName(InvalidName, DuplicateObject): pass

class UnknownChild(Exception): pass

class UnknownInitialiser(Exception): pass

class UnknownReferenceWarning(UserWarning): pass

class HiddenAttribute(Exception):
    """
    Used to hide superclass attributes by throwing error upon attempted access to specific method names
    Used in conjuction with '@property' decorator, e.g.:
     @property
     def method_to_hide: raise HiddenMethod(super().method_to_hide)
    When the method to be hidden is already a property, use:
     @property
     def method_to_hide: raise HiddenMethod(super(self.__class__, self.__class__).method_to_hide.fget, self)
    """
    def __init__(self, attribute_type, method, method_bound_obj = None):
        import re
        ## get method_name
        if isinstance(method, str):
            method_name = re.search("(?<=\.)[^.]+(?= at )", method).group(0)
        else:
            method_name = method.__name__
        ## get method_bound_class
        if method_bound_obj is not None:
            method_bound_class = method_bound_obj.__class__.__name__
        elif re.search("^<bound method ", str(method)):
            method_bound_class = re.search("(?<=\.).+(?= object at)",
                                           str(method).split(f".{method_name} of <")[1]).group(0)
        else:
            method_bound_class = None
        ## get method_definition_class
        if re.search("^<bound method ", str(method)):
            method_definition_class = re.search(f"(?<=^<bound method ).+(?=.{method_name} of <)",
                                                str(method)).group(0)
        elif re.search("^<function ", str(method)):
            method_definition_class = re.search(f"(?<=^<function ).+(?=.{method_name} at )",
                                                str(method)).group(0)
        else:
            method_definition_class = None
        ## build message
        self.message = ( f"{attribute_type} '{method_name}'" +
                         ('' if method_definition_class is None else \
                          f" (defined in class '{method_definition_class}')") +
                         (" is hidden" if method_bound_class is None else \
                          f" is hidden from class '{method_bound_class}'"))
        super().__init__(self.message)

class HiddenMethod(HiddenAttribute):
    def __init__(self, method, method_bound_obj = None):
        super().__init__("Method", method, method_bound_obj = method_bound_obj)

class HiddenProperty(HiddenAttribute):
    def __init__(self, method, method_bound_obj = None):
        super().__init__("Property", method, method_bound_obj = method_bound_obj)


## handler functions
def file_in_use(path):
    return os.path.exists(path) and os.stat(path).st_size > 0

def append_time(path):
    directory, basename = os.path.split(path)
    has_suffix = '.' in basename
    new_basename = '.'.join(path.split('.')[:(-1 if has_suffix else 1)]) + '_' + \
                   str(int(datetime.utcnow().timestamp())) + \
                   (('.' + basename.split('.')[-1]) if has_suffix else '')
    return os.path.join(directory, new_basename)

def check_path_and_return(path):
    if file_in_use(path): return append_time(path)
    return path

## handle filename changes
class DynamicFileHandler(logging.FileHandler):
    
    def __init__(self, filename, check_path = True):
        slogger.debug(f"DynamicFileHandler super().__init__")
        logging.FileHandler.__init__(self, filename)
        if check_path: self._check_path()
    
    ## if file already exists, append unix time (seconds) to new logfile name
    def _check_path(self):
        if file_in_use(self.baseFilename):
            ## close existing file
            self.close()
            ## update new logfile name
            p = append_time(self.baseFilename)
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
        ## move written logfile
        if os.path.exists(p_old): shutil.move(p_old, self.baseFilename)
        ## reopen moved logfile for further logging
        self._open()
        return

class DynamicMultiFileHandler(logging.FileHandler):
    """
    Writes to multple files
    """
    def __init__(self, *filenames, check_path = True):
        slogger.debug(f"DynamicMultiFileHandler super().__init__")
        ## filenames are indexed by positional index (+1) if given as <path> instead of {<alias>: <path>} dict
        ##  we don't use the filenames alone because they may be modified by _check_paths and it'll become hard
        ##  to retrieve if the user doesn't know the new filename
        if len(filenames) == 0: self._filenames_original = {}
        elif isinstance(filenames[0], dict): self._filenames_original = filenames[0]
        else: self._filenames_original = {i+1: filename for i, filename in enumerate(filenames)}
        self._filenames = copy.deepcopy(self._filenames_original)
        ## check all filenames if 'check_path' raised before calling super().__init__ with the first of them
        if check_path: self._check_paths()
        ## call super().__init__ with first filename (or dummy filepath if none provided)
        tmp_filename = (tempfile.mkstemp()[1] if len(filenames) == 0 else tuple(self._filenames.values())[0])
        logging.FileHandler.__init__(self, tmp_filename)
        ## close and remove file (if empty) so all files start with nothing
        self.close()
        if os.path.exists(tmp_filename) and os.stat(tmp_filename).st_size == 0: os.remove(tmp_filename)
    
    def _check_paths(self):
        for k, filename in self._filenames.items():
            self._filenames[k] = check_path_and_return(filename)
        return
    
    def _next_index(self, name):
        if name is not None: return name
        if all(map(lambda x: not isinstance(x, int), self._filenames.keys())): return 1
        else: return max([x for x in self._filenames.keys() if isinstance(x, int)]) + 1
    
    def _get_indices(self, name):
        if name in self._filenames: indices = [name]
        elif name in self._filenames.values() or self._filenames_original.values():
            indices = [k for k, v in self._filenames.items() if v == name] + \
                      [k for k, v in self._filenames_original.items() if v == name]
        else: indices = []
        return indices
    
    def _get_index(self, name):
        indices = self._get_indices(name)
        if not indices: return None
        else: return indices[0]
    
    def add_file(self, filename, name = None, check_path = True, replace = False, quiet = False):
        def local_print(msg):
            if not quiet: print(msg)
        name = self._next_index(name)
        if (name in self._filenames
            and self._filenames[name] != filename
            and self._filenames_original[name] != filename):
            if not replace: raise DuplicateObject(f"Filename alias '{name}' is already in use.")
            elif self._filenames_original[name] != filename and self._filenames[name] != filename:
                local_print ( (f"Filename alias '{name}' is already in use."
                               f" Updating with new filename ({filename}).") )
        elif filename in self._filenames.values() or name in self._filenames_original.values():
            local_print(f"'{filename}' is already used for logging.")
            return
        self._filenames_original[name] = filename
        self._filenames[name] = filename if not check_path else check_path_and_return(filename)
        
    def remove_file(self, name):
        indices = self._get_indices(name)
        for index in indices:
            del self._filenames[index]
            del self._filenames_original[index]
        return
    
    def update_filename(self, old, new, check_path = True, original_new = None):
        indices = self._get_indices(old)
        if not indices: raise InvalidName(f"'{old}' is not known to '{self}'.")
        checked_new = new if not check_path else check_path_and_return(new)
        for index in indices:
            p_old = self._filenames[index]
            self._filenames_original[index] = new if original_new is None else original_new
            self._filenames[index] = checked_new
            if os.path.exists(p_old): shutil.move(p_old, checked_new)
        return
    
    def emit(self, record):
        """
        Call FileHandler's emit method for each filename
        """
        ## open each file, write, and close
        for filename in self._filenames.values():
            self.baseFilename = filename
            self._open()
            logging.FileHandler.emit(self, record)
            self.close()
        return

# ## dmfh test
# dmfh = DynamicMultiFileHandler("/mnt/chaelab/rachelle/tmp/log4.log", "/mnt/chaelab/rachelle/tmp/log5.log")
# dmfh.setLevel(slogger.debug)
# logger = logging.Logger("base")
# logger.addHandler(dmfh)
# logger.info("hello")
# logger.debug("oh no")

# dmfh.setFormatter(logging.Formatter(fmt_full))
# logger.debug("eheh")

# dmfh.remove_file(1)
# logger.debug("ho")

# dmfh.add_file("/mnt/chaelab/rachelle/tmp/log5.log", check_path = False)
# logger.debug("anyway")

# dmfh.add_file("/mnt/chaelab/rachelle/tmp/log4.log", check_path = True)
# logger.debug("hooboy")

# dmfh.add_file("/mnt/chaelab/rachelle/tmp/log4.log", check_path = False)
# logger.debug("anyhoo")


class MultiChild():
    
    def __init__(self, child_adj: str = "child", child_obj: str = "Child", child_alias: str = "Child",
                 child_class = None, make_child = None, unique: bool = False,
                 children: dict = {}, **args_for_instantiating_children_upon_creation):
        slogger.debug(f"MultiChild __init__")
        self._unique = unique
        self._children = {}
        self._child_adj = child_adj
        self._child_obj = child_obj
        self._child_alias = child_adj if child_alias is None else child_alias
        self._child_class = child_class
        self._make_child = make_child
        ## if children already provided, add them
        for name, child in children.items():
            self.add_child(child = child, child_name = name, replace = False,
                           **args_for_instantiating_children_upon_creation)
    
    def _add_child(self, child, child_name, replace = False, ignore_duplicate = False, quiet = False):
        def local_print(msg):
            if not quiet: print(msg)
        if child_name is not None:
            try:
                self.name_available(child_name, raise_exception = True)
            except DuplicateName as e:
                if replace:
                    local_print((f"{self._child_alias} name '{child_name}' is already in use."
                                f" Updating with new {self._child_obj} object."))
                elif child is self.get_child(child_name):
                    local_print((f"{self._child_alias} '{child_name}' has already been added."))
                else:
                    # raise InvalidName(f"{self._child_obj} name '{child_name}' is already in use.")
                    raise e
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
    
    ## some filtering functions
    def filter_children_map(self, f = lambda name, obj: True):
        return {name: child for name, child in self.children_map.items() if f(name, child)}
    def filter_children_names(self, f):
        return list(self.filter_children_map(f).keys())
    def filter_children(self, f):
        return list(self.filter_children_map(f).values())
    
    def is_child_class(self, other):
        return isinstance(other, self._child_class)
    
    def get_child(self, child):
        if child in self.children_names: return self.children_map[child]
        elif child in self.children: return child
        else: raise UnknownChild( (f"'{child}' is not a known {self._child_obj} name or object.") )
    
    def get_child_name(self, child):
        names = [k for k, v in self.children_map.items() if v is child]
        if self._unique: return names[0]
        return names
    
    def is_child(self, other):
        try: tmp = self.get_child(other)
        except UnknownChild: return False
        return True
    
    def name_available(self, name, raise_exception = False):
        if self.is_child(name):
            if raise_exception: raise DuplicateName(f"{self._child_alias} name '{name}' is already in use.")
            else: return False
        elif name in dir(self):
            if raise_exception: raise ReservedName((f"'{name}' is a reserved name for"
                                                    f" class {self.__class__.__name__}."))
            else: return False
        else: return True
    
    def add_child(self, child = None, child_name: str = None, replace = False, quiet = False, **kwargs):
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
                raise DuplicateName( (f"{self._child_alias} name '{child_name}' is already in use."
                                      f" Please rename the {self._child_obj} object.") )
            else:
                self._add_child(child, child_name, replace = replace)
                return
        ## if only child_name is provided, create Child object and return it
        elif child is None and child_name is not None:
            if self._make_child is None:
                raise UnknownInitialiser( (f"Function to make child object is not defined.") )
            child = self._make_child(child_name, **kwargs)
            self._add_child(child, child_name, replace = replace, quiet = quiet)
            return child
        ## if both Child object and child_name are provided, use child_name to index regardless of Child.name
        else: self._add_child(child, child_name, replace = replace, quiet = quiet)
        return
    
    def remove_child(self, child, class_name = None):
        """
        accepts both child_name (str) and Child obj
        """
        if class_name is None: class_name = self.__class__.__name__
        ## if child_name provided
        if isinstance(child, str):
            if child in self._children: self._remove_child(child)
            else: raise UnknownChild( (f"Cannot remove {self._child_adj}: '{child}' is not a"
                                       f" known {self._child_alias} name") )
        ## if Child object provided
        else:
            ## if Child obj is in self.children, retrieve key
            inv_children = dict((v, k) for k, v in self._children.items())
            if child in inv_children: self._remove_child(inv_children[child])
            else: raise UnknownChild( (f"Cannot remove {self._child_adj}: The {self._child_obj} object"
                                       f" provided is not known to this {class_name}") )
        return
    

class MultiHandlerBase(MultiChild):
    
    def __init__(self, **kwargs):
        slogger.debug("MultiHandlerBase super().__init__")
        MultiChild.__init__(self, child_adj = "handler", child_obj = "Handler",
                            child_class = logging.Handler,  **kwargs)
        self.__is_stream_handler = lambda n, h: (type(h) is logging.StreamHandler)
        self.__is_file_handler = lambda n, h: isinstance(h, logging.FileHandler)
    
    @property
    def handler_map(self): return self.children_map
    @property
    def stream_handler_map(self): return self.filter_children_map(self.__is_stream_handler)
    @property
    def file_handler_map(self): return self.filter_children_map(self.__is_file_handler)
    @property
    def all_handlers(self): return self.children
    @property
    def stream_handlers(self): return self.filter_children(self.__is_stream_handler)
    @property
    def file_handlers(self): return self.filter_children(self.__is_file_handler)
    @property
    def handler_names(self): return self.children_names
    @property
    def stream_handler_names(self): return self.filter_children_names(self.__is_stream_handler)
    @property
    def file_handler_names(self): return self.filter_children_names(self.__is_file_handler)
    
    def add_handler(self, handler, handler_name, replace = False, ignore_duplicate = False, quiet = False):
        """
        Accepts a Handler object (e.g. StreamHandler, FileHandler, DynamicFileHandler)
        """
        super()._add_child(handler, handler_name, replace = replace, ignore_duplicate = ignore_duplicate,
                           quiet = quiet)
        return
    
    def remove_handler(self, handler):
        ## get Handler object if handler name (str) provided
        if isinstance(handler, str): handler = self.get_child(handler)
        elif not self.is_child_class(handler): raise InvalidType("Unknown object class")
        super().remove_child(handler)
        return

class MultiHandler(MultiHandlerBase):
    def __init__(self, unique = True):
        MultiHandlerBase.__init__(self, unique = unique)
    
    @property
    def handlers(self): return super().all_handlers

class PublicLogger(logging.Logger, MultiHandlerBase):
    """
    Logger object, but it adds Handler objects as attributes for easy access
    """
    
    def __init__(self, name, unique = True, default_log_level = logging.INFO):
        slogger.debug(f"PublicLogger MultiHandlerBase.__init__")
        MultiHandlerBase.__init__(self, unique = unique)
        slogger.debug(f"PublicLogger logging.Logger.__init__")
        logging.Logger.__init__(self, name)
        self._default_log_level = default_log_level
    
    def addHandler(self, handler, handler_name, replace = False, ignore_duplicate = False, quiet = False):
        super().add_handler(handler, handler_name, replace = replace, ignore_duplicate = ignore_duplicate,
                            quiet = quiet)
        super().addHandler(handler)
        return
    
    def removeHandler(self, handler):
        super().remove_handler(handler)
        super().removeHandler(handler)
        return
    
    def __call__(self, msg, *args, **kwargs):
        if self.isEnabledFor(self._default_log_level):
            self._log(self._default_log_level, msg, args, **kwargs)

class MultiLogger(MultiChild):
    
    def __init__(self, logger_class = PublicLogger, make_logger = lambda name: PublicLogger(name),
                 default_log_level = logging.INFO):
        slogger.debug(f"MultiLogger super().__init__")
        MultiChild.__init__(self, child_adj = "logger", child_obj = "Logger", child_class = logger_class,
                            make_child = make_logger)
        self._default_log_level = default_log_level
    
    @property
    def loggers(self): return self.children
    @property
    def logger_names(self): return self.children_names
    @property
    def logger_map(self): return self.children_map
    
    def is_logger(self, other):
        return self.is_child_class(other)
    
    def add_logger(self, logger, name: str = None, replace = False, quiet = False, level = lvl,
                   default_log_level = None):
        if self.is_logger(logger):
            logger = super().add_child(child = logger, child_name = name, replace = replace, quiet = quiet)
        elif isinstance(logger, str) and name is None:
            logger = super().add_child(child_name = logger, replace = replace, quiet = quiet)
        elif isinstance(logger, str) and name is not None:
            warnings.warn("Both arguments are strings. Ignoring 2nd argument.")
            logger = super().add_child(child_name = logger, replace = replace, quiet = quiet)
        else:
            raise InvalidArgs("Unexpected 1st argument type. Must be Logger object or str.")
        logger.setLevel(level)
        logger._default_log_level = self._default_log_level if default_log_level is None else default_log_level
        return logger
        
    def remove_logger(self, logger):
        return super().remove_child(logger, class_name = self.__class__.__name__)
    
    def get_logger(self, logger):
        return self.get_child(logger)

class GroupLogger(MultiLogger):
    
    def __init__(self, **kwargs):
        slogger.debug(f"GroupLogger super().__init__")
        MultiLogger.__init__(self, **kwargs)
        ## stores preset formats
        slogger.debug(f"GroupLogger._format MultiChild.__init__")
        self._format = MultiChild(child_adj = "formatter", child_obj = "Formatter",
                                  child_class = logging.Formatter,
                                  make_child = lambda name, format: logging.Formatter(format))
        ## stores resuable handlers
        slogger.debug(f"GroupLogger._handler MultiChild.__init__")
        self._handler = MultiHandler(unique = True)
        return
    
    @property
    def group_names(self): return self.logger_names
    @property
    def handler(self): return self._handler
    # ## commented out to avoid conflict w/ Logger when DynamicFileLoggerBase is
    # ##  inherited w/ it in DynamicFileSingleLogger
    # @property
    # def handlers(self): return self.handler.children
    @property ## placeholder for def handlers
    def handler_children(self): return self.handler.all_handlers
    @property
    def handler_names(self): return self.handler.handler_names
    @property
    def handler_map(self): return self.handler.handler_map
    @property
    def format(self): return self._format
    @property ## left in despite inconsistency/nonuniformity w/ handlers.
    def formats(self): return self.format.children
    @property ## for consistency w/ handler_children
    def format_children(self): return self.format.children
    @property
    def format_names(self): return self.format.children_names
    @property
    def format_map(self): return self.format.children_map
    @property
    def all_handlers(self): return self.handler_children
    @property
    def stream_handlers(self): return self.handler.stream_handlers
    ## some properties from class GroupLogger
    @property
    def file_handlers(self): return self.handler.file_handlers
    @property
    def stream_handler_names(self): return self.handler.stream_handler_names
    @property
    def file_handler_names(self): return self.handler.file_handler_names
    @property
    def active_handlers(self): return list(set(itertools.chain(*[logger.all_handlers
                                                                 for logger in self.loggers])))
    @property
    def active_handler_names(self): return list(set(itertools.chain(*[logger.handler_names
                                                                      for logger in self.loggers])))
    @property
    def orphan_handlers(self): return list(set(self.all_handlers) - set(self.active_handlers))
    @property
    def orphan_handler_names(self): return list(set(self.handler_names) - set(self.active_handler_names))
    def groups(self, name = True, group_name = False, handler_name = False):
        if group_name or handler_name: name = False
        if name: group_name = handler_name = True
        return {(logger_name if group_name else logger):
                (logger.handler_names if handler_name else logger.all_handlers)
                for logger_name, logger in self.logger_map.items()}
    
    def _add_handler(self, name, handler, format = "%(asctime)s %(message)s",
                     level = None, replace = False, quiet = False):
        self.handler._add_child(child = handler, child_name = name, replace = replace, quiet = quiet)
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
    
    def _add_handler_to_group(self, group, handler, ignore_duplicate = False):
        group_logger = self.get_logger(group)
        group_name = group_logger.name
        slogger.debug(f"{self}, {group}, {handler}")
        handler_obj = self.get_handler(handler)
        handler_name = self.handler.get_child_name(handler_obj)
        if not ignore_duplicate and handler_obj in group_logger.handlers:
            warnings.warn("Handler '{handler_name}' has already been added.", DuplicateWarning)
        group_logger.addHandler(handler_obj, handler_name)
        return
    
    def _remove_handler_from_group(self, group, handler, name = None):
        group_logger = self.get_logger(group)
        handler_obj = self.get_handler(handler)
        if handler_obj not in group_logger.all_handlers:
            handler_name = self.handler.get_child_name(handler_obj)
            warnings.warn(f"'{handler_name}' is not in group '{group_logger.name}'. Ignoring command.",
                          UnknownReferenceWarning)
        group_logger.removeHandler(handler_obj)
        return
    
    def add_format(self, name, format, replace = False, quiet = False):
        self.format.add_child(child_name = name, format = format, replace = replace, quiet = quiet)
        return
        
    def add_file_handler(self, *names, format = "%(asctime)s %(message)s",
                         level = None, replace = False, quiet = False):
        for name in names:
            new_handler = DynamicFileHandler(self.filename, check_path = False)
            self._add_handler(name, new_handler, format = format, level = level,
                              replace = replace, quiet = quiet)
        return
    
    def add_stream_handler(self, *names, format = "%(asctime)s %(message)s",
                           level = None, replace = False, quiet = False):
        for name in names:
            new_handler = logging.StreamHandler()
            self._add_handler(name, new_handler, format = format, level = level,
                              replace = replace, quiet = quiet)
        return
    
    def get_handler(self, handler):
        return self.handler.get_child(handler)
        
    def add_group(self, name, *handlers, level = lvl, default_log_level = None):
        """
        Positional 'handlers' for Handlers already added to DynamicFileLogger object.
        """
        try: self.name_available(name, raise_exception = True)
        except DuplicateName: raise DuplicateName(f"'{name}' is an existing group.")
        new_logger = self.add_logger(name, level = level, default_log_level = default_log_level)
        ## parse existing handlers
        for i, handler in enumerate(handlers):
            self._add_handler_to_group(name, handler)
        return
    
    def remove_group(self, name):
        if not self.is_child(name): warning.warn(f"'{name}' is not an existing group. Ignoring command.",
                                                 UnknownReferenceWarning)
        self.remove_logger(name)
        return
    
    def update_group(self, group, add = None, remove = None, level = None):
        """
        Add or remove handlers from group
        """
        logger = self.get_child(group)
        group = logger.name
        if add is not None:
            slogger.debug(f"add: {add}")
            if isinstance(add, str) or isinstance(add, logging.Handler):
                add = [add]
            slogger.debug(f"add2: {add}")
            for handler in add:
                self._add_handler_to_group(group, handler, ignore_duplicate = True)
        if remove is not None:
            if isinstance(remove, str) or isinstance(remove, logging.Handler):
                remove = [remove]
            for handler in remove:
                self._remove_handler_from_group(group, handler)
        if level is not None:
            logger.setLevel(level)
        return    


class PublicAwareLogger(PublicLogger):
    """
    Handles crosstalk and logging level
    """
    
    def __init__(self, name, parent, crosstalk = True):
        self._parent = parent
        self._crosstalk = crosstalk
        super().__init__(name)
    
    @property
    def crosstalk(self): return self._crosstalk
    @crosstalk.setter
    def crosstalk(self, v: bool): self._crosstalk = v
    
    def _log(self, level, *args, **kwargs):
        """
        Modify logging.Logger._log so that files not linked to self's group won't be written to if
         self.crosstalk is False
        """
        if isinstance(self._parent, DynamicFileLogger):
            if self._parent.is_enabled_for(level):
                super()._log(level, *args, **kwargs)
        elif isinstance(self._parent, DynamicMultiFileLogger):
            if self.crosstalk:
                linked_files = {self._parent._files.get_child_name(f)[0]: f.filename
                                for f in self._parent._files.children
                                if f.is_enabled_for(level)}
            else:
                linked_files = {self._parent._files.get_child_name(f)[0]: f.filename
                                for f in self._parent._files.children
                                if f.is_enabled_for(level) and self in f._active_groups}
            store = {}
            ## settle each handler individually
            for handler in self.handlers:
                if isinstance(handler, DynamicMultiFileHandler):
                    store[handler] = handler._filenames
                    filename_intersection = (set(handler._filenames.keys()) & set(linked_files.keys()))
                    handler._filenames = {fname: linked_files[fname] for fname in filename_intersection}
                elif isinstance(handler, logging.FileHandler): ## FileHandler or DynamicFileHandler
                    store[handler] = handler
                    ## remove handler so it doesn't get called if wrong level or file not linked to group
                    if handler.baseFilename not in linked_files.values():
                        self.removeHandler(handler)
                else: continue ## don't do anything if is StreamHandler
            ## call logging.Logger._log
            super()._log(level, *args, **kwargs)
            ## restore previous filename(s) for all handlers and add handlers back to self
            for handler, stored in store.items():
                if isinstance(handler, DynamicMultiFileHandler):
                    handler._filenames = store[handler]
                elif isinstance(handler, logging.FileHandler):
                    self.addHandler(handler)
        return


class DynamicFileLoggerBase(GroupLogger):
    """
    User specifies a group when passing message
    (E.g. DynamicFileLoggerBase.cmdname('minimumset') and the relevant handlers will do their job)
    (Note that this works by creating a Logger called 'cmdname' and handlers linked to that Logger 
       will do their job. Maybe a StreamHandler will not print it cuz it's not relevant, 
       but the FileHandler may wish to do so for thoroughness. There may be different formatting.)
    User may also change log file location whenever desired
    (E.g. DynamicFileLoggerBase.update_filename('/mnt/chaelab/rachelle/logs/loggertest.log'))
    
    Groups are effectively Logger objects
    """
    
    def __init__(self, filename, level = lvl, default_log_level = logging.INFO,
                 **config_defaults): ## config_defaults currently not used
        slogger.debug(f"DynamicFileLoggerBase super().__init__")
        GroupLogger.__init__(self, logger_class = PublicAwareLogger,
                             make_logger = lambda name: PublicAwareLogger(name, self),
                             default_log_level = default_log_level)
        slogger.debug(f"DynamicFileLoggerBase._base_handler DynamicFileHandler.__init__")
        self._base_handler = DynamicFileHandler(filename, check_path = True)
        self._level = level
        return
    
    @property
    def filename(self): return self._base_handler.baseFilename
    @property
    def level(self): return self._level
    @level.setter
    def level(self, level): self._level = level
    
    def set_level(self, level): self.level = level
    def is_enabled_for(self, level): return level >= self._level
    
    def update_filename(self, filename, check_path = True):
        """
        Closes connections from DynamicFileHandler to current logfile,
        replaces DynamicFileHandler.baseFilename, copies closed logfile to new location,
        and opens connection between all DynamicFileHandler objects to new location
        """
        # checked_filename = check_path_and_return(filename) if check_path else filename
        ## update path for DynamicFileHandler, from which we will retrieve the new path for other Handlers
        self._base_handler.update_filename(filename, check_path = check_path)
        file_handlers = [h for h in self.handler.children if isinstance(h, DynamicFileHandler)]
        ## get new log file location (raise check_path only for first Handler processed)
        for file_handler in file_handlers:
            file_handler.update_filename(self._base_handler.baseFilename, check_path = False)
        return


class DynamicFileLogger(DynamicFileLoggerBase):
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
    
    def __init__(self, filename, level = lvl, **config_defaults):
        slogger.debug(f"DynamicFileLogger super().__init__")
        DynamicFileLoggerBase.__init__(self, filename, level = level, **config_defaults)
    
    @property
    def handlers(self): return self.handler_children

# class DynamicFileParasiticSingleLogger(PublicLogger, DynamicFileLoggerBase):
class DynamicFileParasiticSingleLogger(DynamicFileLoggerBase):
    
    def __init__(self, name, filename, parent, level = lvl, **config_defaults):
        self.name = name
        self._parent = parent
        self._base_logger = PublicLogger(name, unique = True)
        self._base_logger.setLevel(level)
        slogger.debug("DynamicFileSingleLogger DynamicFileLoggerBase.__init__")
        DynamicFileLoggerBase.__init__(self, filename, **config_defaults)
        ## replace self._format and self._handler objs w/ parent DynamicMultiFileLogger's
        self._format = self._parent._format
        self._handler = self._parent._handler
        self._active_groups = set()
        ## add self's logger as (one and only) child Logger
        slogger.debug("DynamicFileSingleLogger adding self._base_logger to list of loggers now")
        super().add_child(self._base_logger, name = name)
    
    def _get_group(self, grp): return (None if not self._parent.is_child(grp) else
                                       self._parent.get_logger(grp))
    def _to_handler_group(self, group_or_handler, file_handler_only = True):
        """
        If group_or_handler is a group, returns ([<handlers in group>], [<group>])
        Else, returns ([<handler>], [])
        """
        if self._handler.is_child(group_or_handler):
            handlers = [self._handler.get_child(group_or_handler)]
            groups = []
        elif self._parent.is_child(group_or_handler):
            handlers = list(self._parent.get_logger(group_or_handler).all_handlers)
            groups = [self._parent.get_logger(group_or_handler)]
        else:
            handlers = []
            groups = []
        if len([handler for handler in handlers if not isinstance(handler, logging.FileHandler)]) != 0:
            print("Handlers that are not FileHandler objects will be ignored.")
        if file_handler_only: handlers = [h for h in handlers if isinstance(h, logging.FileHandler)]
        slogger.debug(handlers)
        return handlers, groups
    def _to_handler(self, group_or_handler, file_handler_only = True):
        return self._to_handler_group(group_or_handler, file_handler_only = file_handler_only)[0]
    def _to_group(self, group_or_handler): return self._to_handler_group(group_or_handler)[1]
    
    @property
    def level(self): return self._base_logger.level
    @level.setter
    def level(self, level): self._base_logger.setLevel(level)
    def set_level(self, level): self.level = level
    def is_enabled_for(self, level): return self._base_logger.isEnabledFor(level)
    
    ## hide group-related properties/methods
    @property
    def group_names(self): raise HiddenProperty(super(self.__class__, self.__class__).group_names.fget, self)
    # @property
    # def groups(self): raise HiddenMethod(super(self.__class__, self.__class__).groups, self)
    @property
    def add_group(self): raise HiddenMethod(super(self.__class__, self.__class__).add_group, self)
    @property
    def add_logger(self): raise HiddenMethod(super(self.__class__, self.__class__).add_logger, self)
    @property
    def remove_logger(self): raise HiddenMethod(super(self.__class__, self.__class__).remove_logger, self)
    @property
    def update_group(self): raise HiddenMethod(super(self.__class__, self.__class__).update_group, self)
    
    ## hide handler-related methods so user can't accidentally modify handlers when using DynamicMultiFileLogger
    @property
    def add_stream_handler(self):
        raise HiddenMethod(super(self.__class__, self.__class__).add_stream_handler, self)
    @property
    def add_file_handler(self):
        raise HiddenMethod(super(self.__class__, self.__class__).add_file_handler, self)
        
    ## returns groups that are explicitly included in self.
    @property
    def groups(self): return self._active_groups
    @property
    def group_names(self): return [self._parent.get_child_name(group) for group in self.groups]
    ## returns groups that are implicitly included in self. (I.e. group's file handlers is subset of self's)
    @property
    def handler_groups(self):
        self_handlers = set(handler for handler in self.active_handlers
                            if isinstance(handler, logging.FileHandler))
        parent_groups = self._parent.groups(name = False)
        parent_groups_file_handlers = {group: set(handler for handler in handlers
                                                  if isinstance(handler, logging.FileHandler))
                                       for group, handlers in parent_groups.items()}
        self_groups = [group for group, handlers in parent_groups_file_handlers.items()\
                       if handlers.issubset(self_handlers)]
        return self_groups
    @property
    def handler_group_names(self): return [self._parent.get_child_name(group) for group in self.handler_groups]
        
    ## 1 master function and...4 other functions that are basically the master function but more specific lol
    ## note that these 'add' and 'remove' functions only modify what's active for this logger
    ## - they do not modify the shared pool of handlers/groups. to do that, use the the _parent's functions
    def update_level(self, level): self._base_logger.setLevel(level)
    def addHandler(self, handler, handler_name): ## for use by _add_handler_to_group
        if isinstance(handler, DynamicMultiFileHandler): handler.add_file(self.filename, self.name)
        super().addHandler(handler, handler_name)
    def add_handler(self, *handlers):
        """
        Accepts both group names/obj and handler names/obj
         - if group name/obj is received, handlers will be added to file but group will NOT be
        """
        handlers = list(set(itertools.chain(*[self._to_handler(h) for h in handlers])))
        for h in handlers:
            h.add_file(self.filename, self.name)
        self._update(add = list(itertools.chain(*[self._to_handler(h) for h in handlers])))
    def removeHandler(self, handler): ## for use by _remove_handler_from_group
        if isinstance(handler, DynamicMultiFileHandler): handler.remove_file(self.name)
        super().removeHandler(handler)
    def remove_handler(self, *handlers):
        """
        Accepts both group names/obj and handler names/obj
         - if group name/obj is received, handlers will be removed from file but group will NOT be
          - that is, if a handler already added to file is then added to group, group.emit will cause that
            newly added handler to write to this file
        """
        handlers = list(set(itertools.chain(*[self._to_handler(h) for h in handlers])))
        for h in handlers:
            h.remove_file(self.name)
        self._update(remove = handlers)
    def update_handlers(self, add = None, remove = None):
        self._update(add = add, remove = remove)
    
    ## active group modifiers
    def add_group(self, *groups):
        groups_to_add = []
        handlers_to_add = []
        for group in groups:
            if self._get_group(group) is not None:
                groups_to_add.append(group)
                self._active_groups |= {self._get_group(group)}
                handlers_to_add.extend(list(self._get_group(group).file_handlers))
            else: warnings.warn(f"'{group}' is not a known group. Will be ignored.", UnknownReferenceWarning)
        self.add_handler(*handlers_to_add)
    def remove_group(self, *groups):
        groups_to_remove = []
        handlers_to_remove = []
        for group in groups:
            if self._get_group(group) is not None:
                if self._get_group(group) in self._active_groups:
                    groups_to_remove.append(group)
                    handlers_to_remove.extend(list(self._get_group(group).file_handlers))
                    self._active_groups -= {self._get_group(group)}
                else: warnings.warn(f"'{group}' is not active for this file '{self._base_logger.name}'.")
            else: warnings.warn(f"'{group}' is not a known group. Will be ignored.", UnknownReferenceWarning)
        self.remove_handler(*handlers_to_remove)
    
    def update_filename(self, filename, check_path = True):
        """
        Replaces DynamicMultiFileHandler records for current file with new location
        """
        ## get new log file location
        checked_filename = filename if not check_path else check_path_and_return(filename)
        ## move file
        shutil.move(self.filename, checked_filename)
        ## separately update each handler's stored filename(s)
        file_handlers = [h for h in self.handler.children if isinstance(h, logging.FileHandler)]
        for file_handler in file_handlers:
            if isinstance(file_handler, DynamicMultiFileHandler):
                if self.filename in file_handler._filenames.values():
                    file_handler.update_filename(self.filename, checked_filename,
                                                 original_new = filename, check_path = False)
            elif isinstance(file_handler, DynamicFileHandler):
                if file_handler.filename == self.filename:
                    file_handler.update_filename(checked_filename, check_path = False)
            else:
                if file_handler.baseFilename == self.filename:
                    file_handler.baseFilename = checked_filename
        ## update self's filename location
        self._base_handler.baseFilename = checked_filename
        return
    
    def _update(self, **kwargs):
        super().update_group(self._base_logger, **kwargs)
    # def update(self, add = None, remove = None, level = None, filename = None):
    def update(self, add = None, remove = None, level = None, filename = None):
        """
        Both 'add' and 'remove' accept groups, not just handlers? Maybe? TODO
        If a name (provided to 'add' or 'remove' refers to both a handler and a group, handler is prioritised
        """
        if filename is not None: self.update_filename(filename)
        if level is not None: self.update_level(level)
        if add is not None:
            add_g = list(set(itertools.chain(*[self._to_group(e) for e in add])))
            add_h = list(set(itertools.chain(*[self._to_handler(e) for e in add])) -
                         set(itertools.chain(*[self._to_handler(e) for e in add_g])))
            self.add_group(*add_g) ## adds groups + all (file) handlers associated w/ the groups
            self.add_handler(*add_h)
        if remove is not None:
            remove_g = list(set(itertools.chain(*[self._to_group(e) for e in remove])))
            remove_h = list(set(itertools.chain(*[self._to_handler(e) for e in remove])) -
                            set(itertools.chain(*[self._to_handler(g) for g in remove_g])))
            self.remove_group(*remove_g) ## remove groups + all (file) handlers associated w/ the groups
            self.remove_handler(*remove_h)
        return

## this needs to be reorganised to use DynamicMultiFileHandler
class DynamicMultiFileLogger(GroupLogger):
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
    # dmflogger.log1.add_handler("header") ## duplicate fheader's DynamicFileHandler 
    # dmflogger.log2.add_handler("header", "full") ## duplicate ffull AND fheader's DynamicFileHandlers
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
    def __init__(self, crosstalk = True):
        """
        group_crosstalk:
         - if True, groups that were not added to file can write to file via shared handlers.
           - Handlers can be added by both add_group and add_handler.
         - if False, groups can only write to file if group has been added to file
           - Groups can only be added to file using:
             - DynamicMultiFileLogger.file1.add_group(<group>)
             - DynamicMultiFileLogger.update_file(<file1>, add = <group>)
             - DynamicMultiFileLogger.update_group(<file1>, add = <group>)
           - DynamicMultiFileLogger.file1.add_handler(<group>) and 
             DynamicMultiFileLogger.update_file(<file1>, add = <group handlers>)
             DynamicMultiFileLogger.update_group(<file1>, add = <group handlers>)
             will link individual handlers in the group but not the group to file
        """
        self._crosstalk = crosstalk
        slogger.debug(f"DynamicMultiFileLogger GroupLogger.__init__")
        GroupLogger.__init__(self,
                             logger_class = PublicAwareLogger,
                             make_logger = (lambda name:
                                            PublicAwareLogger(name, self, crosstalk = crosstalk)))
        ## Logger obj will be created by self <DynamicMultiFileLogger> then added to self._files.
        self._files = MultiChild(child_adj = "file", child_obj = "DynamicFileParasiticSingleLogger",
                                 child_class = DynamicFileParasiticSingleLogger,
                                 make_child = (lambda name, filename, **kwargs: \
                                               DynamicFileParasiticSingleLogger(name, filename, self, **kwargs)))
    
    @property
    def handlers(self): return self.handler_children
    @property
    def crosstalk(self): return self._crosstalk
    @crosstalk.setter
    def crosstalk(self, crosstalk: bool):
        for logger in self.loggers:
            logger._crosstalk = crosstalk
        self._crosstalk = crosstalk
    
    def add_file_handler(self, name, format = "%(asctime)s %(message)s",
                         level = None, replace = False, quiet = False):
        new_handler = DynamicMultiFileHandler() ## differs from GroupLogger.add_file_handler only at this line
        self._add_handler(name, new_handler, format = format, level = level, replace = replace, quiet = quiet)
        return
    
    def add_file(self, name, filename, *handlers, replace = False, check_path = True):
        """
        Positional 'handlers' for Handlers already added to DynamicFileLogger object.
        """
        ## check if existing group/file
        if self.is_child(name):
            if name in self._files:
                if replace: print((f"File alias '{name}' is already in use."
                                   f" Updating with new filename ('{filename}')."))
                else: raise DuplicateName((f"'{name}' is an existing file at"
                                           f" '{self.get_child(name).filename}'"))
            else: raise DuplicateName(f"'{name}' is an existing group.")
        ## make new DynamicFileSingleLogger
        new_dynamicfilesinglelogger = self._files.add_child(child_name = name, filename = filename,
                                                            check_path = check_path)
        ## add new PublicLogger to self, so it can now be treated as a 'group'
        self.add_child(new_dynamicfilesinglelogger, name = name)
        ## parse handlers
        for i, handler in enumerate(handlers):
            self._add_handler_to_group(name, handler)
            handler.add_file(filename, name = name)
        return
    
    def update_group(self, group, add = None, remove = None, level = None):
        """
        Add or remove handlers from group
        """
        if self._files.is_child(group):
            Warning(f"'{group}' is a file.")
            self.update_file(group, add = add, remove = remove, level = level)
        else: super().update_group(group, add = add, remove = remove, level = level)
    
    def update_file(self, fname, add = None, remove = None, level = None, filename = None):
        """
        Both 'add' and 'remove' accept groups, not just handlers
        If a name (provided to 'add' or 'remove' refers to both a handler and a group, handler is prioritised
        """
        f = self._files.get_child(fname)
        # super().update_group(fname, add = add, remove = remove, level = level, filename = filename)
        f.update_group(fname, add = add, remove = remove, filename = filename)
    
    def update_filename(self, name, filename):
        f.update_file(name, filename = path)

## test
def test_logger(test: str = "dmfl", dir = None, test_getters: bool = True, logger = None):
    
    ## TODO: test file logging level restriction for DynamicFileLogger
    
    def print_file(fname, msg):
        with open(fname, 'r') as f:
            print(f"\n{msg} ({fname}):\n\t" + str(f.read()).replace('\n', "\n\t"))
    
    test = ("dmfl" if isinstance(logger, DynamicMultiFileLogger) \
            else "dfl" if isinstance(logger, DynamicFileLogger) \
            else test)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        if dir is not None: mkfname = lambda f: os.path.join(dir, f)
        else: mkfname = lambda f: os.path.join(tmpdir, f)
        
        ## instantiate formats
        fmt_cmd = "command: %(message)s"
        fmt_timed = "%(asctime)s - %(levelname)s - %(message)s"
        fmt_full = "%(asctime)s - %(name)s - %(levelname)s:%(message)s"
        fmt_rand = "SPECIAL %(asctime)s - %(name)s - %(levelname)s:%(message)s"
        
        ## create logger
        if logger is None:
            if test == "dmfl":
                logger = DynamicMultiFileLogger(crosstalk = False)
            elif test == "dfl":
                logger = DynamicFileLogger(mkfname("dfl_testlog1.log"))
            else:
                print("Unknown test schema '{test}'")
                return
        
        ## add formats
        logger.add_format(name = "cmd", format = fmt_cmd)
        logger.add_format(name = "timed", format = fmt_timed)
        logger.add_format(name = "full", format = fmt_full)
        logger.format.full = fmt_full ## problematic assignment
        logger.add_format(name = "full", format = fmt_full, replace = True)
        assert (isinstance(logger.format.full, logging.Formatter)
                and logger.format.full._fmt == fmt_full), "Problematic assignment not resolved"
        logger.add_stream_handler("sinfo", level = logging.INFO, format = "timed")
        logger.add_file_handler("finfo", level = logging.INFO, format = logger.format.full)
        logger.add_stream_handler("sheader", level = logging.DEBUG)
        logger.add_file_handler("fheader", level = logging.DEBUG, format = "cmd")
        logger.add_stream_handler("soops", level = logging.ERROR, format = "timed")
        logger.add_file_handler("foops", level = logging.ERROR, format = "timed")
        logger.add_stream_handler("random", level = logging.DEBUG, format = fmt_rand)
        
        ## raise error here
        try: logger.add_stream_handler("random", level = logging.DEBUG, format = fmt_rand)
        except DuplicateName: print("DuplicateName exception caught")
        logger.add_stream_handler("random", level = logging.DEBUG, format = fmt_rand, replace = True)
        
        ## raise error here (reserved name)
        try: logger.add_stream_handler("get_child", "sinfo")
        except ReservedName: print("ReservedName exception caught")
        
        logger.add_group("header", "fheader")
        logger.update_group("header", add = ["sheader"])
        ## raise error here (unknown Handler)
        try: logger.update_group("header", add = ["ohno"])
        except UnknownChild: print("UnknownChild exception caught")
        
        logger.add_group("errors", "soops", "foops")
        logger.add_group("info", logger.handler.sinfo, "finfo")
        logger.add_group("all", *logger.handlers)
        
        logger.update_group(logger.all, remove = ["soops", "foops"])
        logger.update_group("all", add = logger.handlers)
        
        ## abstraction breaking. (Maybe find a way to stop this? Or ensure that even if this happens
        ##   we can propagate it up to DynamicFileLogger so it's accessibly to other groups as well AND
        ##   shows up when logger.handlers/handler_names/handler_map is called (alt, add new function so user
        ##   can retrieve handlers that were added in this abstraction-breaking way)?)
        logger.all.addHandler(logging.StreamHandler(), "breakabstraction")
        logger.all.breakabstraction.setLevel(logging.ERROR)
        logger.all.breakabstraction.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(message)s"))
        
        if test == "dmfl":
            ## add first file
            logger.add_file("log1", mkfname("dmfl_testlog1.log"))
            logger.log1.active_handlers
            logger.log1.add_handler("finfo")
            logger.log1.active_handlers
            print("\tlogger.log1.groups:", logger.log1.groups)
            print("\tlogger.log1.handler_names:", logger.log1.handler_names)
            logger.log1.add_group("errors")
            logger.log1.active_handlers
            logger.log1.active_handlers[0]._filenames
            logger.log1.add_group("info")
            logger.log1.add_handler("info")
            logger.log1.update(add = "info")
            print("\tlogger.log1.group_names:", logger.log1.group_names)
            # ## add second file
            # logger.add_file("log2", mkfname("dmfl_testlog2.log"))
            # logger.log2.add_handler("errors")
            # logger.log2.add_group("info")
            # logger.log2.set_level(logging.WARNING) ## TODO: verify logger.info.info(xxx) does not write to log2
            
        
        ## test handler outputs
        ## print filenames for reference
        if test == "dmfl": print("logger.log1.filename:", logger.log1.filename)
        elif test == "dfl": print("logger.filename:", logger.filename)
        logger.errors.error("houston we have a problem")
        logger.header.info("executing 'minimumset'")
        logger.info.info("hello we're executing something")
        logger.all.warning("hey stop and pay attention")
        logger.all.debug("it's a mess here is what it is")
        logger.update_group("all", level = logging.DEBUG)
        ## 'all' should not show up in file
        
        if test == "dmfl":
            
            ## check no crosstalk between group 'all' and logger.log1
            with open(logger.log1.filename, 'r') as f:
                assert "- all -" not in f.read(), "Crosstalk inappropriately enabled"
            
            logger.crosstalk = True
            logger.errors.error("houston we have a problem")
            logger.header.info("executing 'minimumset'")
            logger.info.info("hello we're executing something")
            logger.all.warning("hey stop and pay attention")
            logger.all.debug("it's a mess here is what it is")
            ## 'all' should show up in file now
            
            ## check crosstalk between group 'all' and logger.log1
            print_file(logger.log1.filename, "Crosstalk enabled")
            with open(logger.log1.filename, 'r') as f:
                assert "- all -" in f.read(), "Crosstalk inexplicably disabled"
    
            ## check that updated filename doesn't break 
            logger.log1.update_filename(mkfname("dmfl_testlog3.log"))
            print("Updated log1 filename:", logger.log1.filename)
            print_file(logger.log1.filename, "File moved")
            
            logger.errors.error("houston we have a problem")
            logger.header.info("executing 'minimumset'")
            logger.info.info("hello we're executing something again")
            logger.all.warning("hey stop and pay attention again!")
            logger.all.debug("it's a mess here is what it is")
            # logger.update_group("all", level = logging.DEBUG)
            ## these should be written to the new file
            
            ## check crosstalk between group 'all' and logger.log1
            print_file(logger.log1.filename, "File moved and written")
            with open(logger.log1.filename, 'r') as f:
                assert f.read().count("- all -") == 2, "Filename not properly updated"
        
        if test == "dfl":
            ## move from previous location to new location (for DynamicFileLogger)
            ftestlog2 = mkfname("dfl_testlog2.log")
            logger.update_filename(ftestlog2)
            logger.filename
            logger2 = DynamicFileLogger(ftestlog2)
            assert logger2.filename!=ftestlog2, "Move unsuccessful. Name conflict between logger1 and logger2."
            
            ## test logging level restriction
            logger.level = logging.INFO
            logger.header.debug("log level lower than logger level")
            with open(logger.filename, 'r') as f:
                assert f.read().count("log level lower") == 0, "Logger level restriction fail"
            
            logger.header.error("log level higher than logger level")
            with open(logger.filename, 'r') as f:
                assert f.read().count("log level higher") == 1, "Logger level enable fail"
            
            logger.level = logging.DEBUG
            logger.header.debug("log level equal to logger level")
            with open(logger.filename, 'r') as f:
                assert f.read().count("log level equal to") == 1, "Logger equivalent level enable fail"
        
        if test_getters:
            print("Testing getters")
            print("handlers:", logger.handlers)
            print("handler_names:", logger.handler_names)
            print("handler_map:", logger.handler_map)
            print("groups:", logger.groups)
            print("formats:", logger.formats)
            print("format_names:", logger.format_names)
            print("<fmt name>: <fmt>", {k: fmt._fmt for k, fmt in logger.format_map.items()})
            print("active_handlers:", logger.active_handlers)
            print("active_handler_names:", logger.active_handler_names)
            print("orphan_handlers:", logger.orphan_handlers)
            print("orphan_handler_names:", logger.orphan_handler_names)
            print("groups(name = True):", logger.groups(name = True))
            print("group_names:", logger.group_names)
            print("header.handler_names (i.e. names of handlers in group 'header'):",
                  logger.header.handler_names)
            print("header.stream_handler_names:", logger.header.stream_handler_names)
            print("[g.name for g in logger.groups(name = False).keys()]:",
                  [g.name for g in logger.groups(name = False).keys()])
    
    return

# tmp_dir = "/mnt/chaelab/rachelle/tmp"
# tmp_f = tempfile.mkstemp(dir = "/mnt/chaelab/rachelle/tmp")[1]

# logger = DynamicMultiFileLogger(crosstalk = False)
# test_logger(dir = tmp_dir, test_getters = False, logger = logger)

# logger = DynamicFileLogger(tmp_f[1])
# test_logger(dir = tmp_dir, test_getters = False, logger = logger)

def minorg_logger(fname = None, level = logging.INFO, default_log_level = logging.INFO):
    if fname is None: fname = tempfile.mkstemp()[1]
    logger = DynamicFileLogger(fname, level = level, default_log_level = default_log_level)
    ## formats
    logger.add_format("header", "-- %(message)s --")
    logger.add_format("plain", "%(message)s")
    logger.add_format("full", "%(asctime)s - %(name)s - %(levelname)s:%(message)s")
    ## generic handlers (set level to most permissive)
    logger.add_stream_handler("splain", format = "plain", level = logging.DEBUG)
    logger.add_file_handler("fplain", format = "plain", level = logging.DEBUG)
    logger.add_stream_handler("sfull", format = "full", level = logging.DEBUG)
    logger.add_file_handler("ffull", format = "full", level = logging.DEBUG)
    ## generic groups
    logger.add_group("plain", "splain", "fplain")
    logger.add_group("splain", "splain")
    logger.add_group("fplain", "fplain")
    logger.add_group("full", "sfull", "ffull")
    logger.add_group("sfull", "sfull")
    logger.add_group("ffull", "ffull")
    ## wrap preset log calls (except exception, which seems to require special handling)
    logger.add_group("critical", "sfull", "ffull", default_log_level = logging.CRITICAL)
    logger.add_group("error", "sfull", "ffull", default_log_level = logging.ERROR, level = logging.ERROR)
    logger.add_group("warning", "sfull", "ffull", default_log_level = logging.WARNING)
    logger.add_group("info", "sfull", "ffull", default_log_level = logging.INFO)
    logger.add_group("debug", "sfull", "ffull", default_log_level = logging.DEBUG, level = logging.DEBUG)
    logger.add_group("sdebug", "sfull", level = logging.DEBUG, default_log_level = logging.DEBUG) # console
    logger.add_group("fdebug", "ffull", level = logging.DEBUG, default_log_level = logging.DEBUG) # file
    ## args
    # logger.add_file_handler("fargs", format = "plain", level = logging.DEBUG)
    logger.add_group("args", "fplain")
    ## header
    logger.add_stream_handler("sheader", format = "header", level = logging.DEBUG)
    logger.add_file_handler("fheader", format = "header", level = logging.DEBUG)
    logger.add_group("header", "sheader", "fheader")
    return logger

# logger = minorg_logger(level = logging.DEBUG)
# logger.filename
# logger.header.info("full") ## writes to both file and console
# logger.args.info("--gene AT1G01010 --acc ref") ## writes only to file
# logger.update_filename("/mnt/chaelab/rachelle/tmp/testminorglogger.log")
# logger.debug.debug("testing debug level")
# logger.args("testing default logging in args")
# logger.error("test error default log (ERROR)")
# logger.warning("test warning default log (WARNING)")
# logger.error.info("test error with info call") ## does nothing cuz we set error's level to ERROR
# logger.warning.info("test warning with info call") ## does something cuz warning's level is defualt (INFO)
# logger.error.info("test error with higher call level (INFO) than error group level (ERROR)") ## prints nothing

# logger2 = minorg_logger(level = logging.INFO)
# logger2.filename
# logger2.header.info("full") ## writes to both file and console
# logger2.args.info("--gene AT1G01010 --acc ref") ## writes only to file
# logger2.update_filename("/mnt/chaelab/rachelle/tmp/testminorglogger.log")
# logger2.debug("test debug with default log (DEBUG)")
# logger2.debug.info("test debug with info call") ## prints something
