from dynamic_logger import *

# def minorg_logger(fname = None, level = logging.INFO, default_log_level = logging.INFO):
#     if fname is None: fname = tempfile.mkstemp()[1]
#     logger = DynamicFileLogger(fname, level = level, default_log_level = default_log_level)
#     ## formats
#     logger.add_format("header", "-- %(message)s --")
#     logger.add_format("plain", "%(message)s")
#     logger.add_format("full", "%(asctime)s - %(name)s - %(levelname)s:%(message)s")
#     ## generic handlers (set level to most permissive)
#     logger.add_stream_handler("splain", format = "plain", level = logging.DEBUG)
#     logger.add_file_handler("fplain", format = "plain", level = logging.DEBUG)
#     logger.add_stream_handler("sfull", format = "full", level = logging.DEBUG)
#     logger.add_file_handler("ffull", format = "full", level = logging.DEBUG)
#     ## generic groups
#     logger.add_group("plain", "splain", "fplain")
#     logger.add_group("splain", "splain")
#     logger.add_group("fplain", "fplain")
#     logger.add_group("full", "sfull", "ffull")
#     logger.add_group("sfull", "sfull")
#     logger.add_group("ffull", "ffull")
#     ## wrap preset log calls (except exception, which seems to require special handling)
#     logger.add_group("critical", "sfull", "ffull", default_log_level = logging.CRITICAL)
#     logger.add_group("error", "sfull", "ffull", default_log_level = logging.ERROR, level = logging.ERROR)
#     logger.add_group("warning", "sfull", "ffull", default_log_level = logging.WARNING)
#     logger.add_group("info", "sfull", "ffull", default_log_level = logging.INFO)
#     logger.add_group("debug", "sfull", "ffull", default_log_level = logging.DEBUG, level = logging.DEBUG)
#     logger.add_group("sdebug", "sfull", level = logging.DEBUG, default_log_level = logging.DEBUG) # console
#     logger.add_group("fdebug", "ffull", level = logging.DEBUG, default_log_level = logging.DEBUG) # file
#     ## args
#     # logger.add_file_handler("fargs", format = "plain", level = logging.DEBUG)
#     logger.add_group("args", "fplain")
#     ## header
#     logger.add_stream_handler("sheader", format = "header", level = logging.DEBUG)
#     logger.add_file_handler("fheader", format = "header", level = logging.DEBUG)
#     logger.add_group("header", "sheader", "fheader")
#     return logger


class MINORgLogger(DynamicFileLogger):
    
    def __init__(self, filename = None, level = logging.INFO, default_log_level = logging.INFO):
        
        if filename is None: filename = tempfile.mkstemp()[1]
        DynamicFileLogger.__init__(self, filename, level = level, default_log_level = default_log_level)
        
        ## formats
        self.add_format("header", "-- %(message)s --")
        self.add_format("plain", "%(message)s")
        self.add_format("full", "%(asctime)s - %(name)s - %(levelname)s:%(message)s")
        ## generic handlers (set level to most permissive)
        self.add_stream_handler("splain", format = "plain", level = logging.DEBUG)
        self.add_file_handler("fplain", format = "plain", level = logging.DEBUG)
        self.add_stream_handler("sfull", format = "full", level = logging.DEBUG)
        self.add_file_handler("ffull", format = "full", level = logging.DEBUG)
        ## generic groups
        self.add_group("plain", "splain", "fplain")
        self.add_group("splain", "splain")
        self.add_group("fplain", "fplain")
        self.add_group("full", "sfull", "ffull")
        self.add_group("sfull", "sfull")
        self.add_group("ffull", "ffull")
        ## developer-only groups
        import getpass
        if getpass.getuser() in {"rachelle"}:
            self.add_group("devsplain", "splain")
        else:
            self.add_group("devsplain")
        ## wrap preset log calls (except exception, which seems to require special handling)
        self.add_group("critical", "sfull", "ffull", default_log_level = logging.CRITICAL)
        self.add_group("error", "sfull", "ffull", default_log_level = logging.ERROR, level = logging.ERROR)
        self.add_group("warning", "sfull", "ffull", default_log_level = logging.WARNING)
        self.add_group("info", "sfull", "ffull", default_log_level = logging.INFO)
        self.add_group("debug", "sfull", "ffull", default_log_level = logging.DEBUG, level = logging.DEBUG)
        self.add_group("sdebug", "sfull", level = logging.DEBUG, default_log_level = logging.DEBUG) # console
        self.add_group("fdebug", "ffull", level = logging.DEBUG, default_log_level = logging.DEBUG) # file
        ## args (keep either group args of fmt_args (rename if fmt_args is kept)
        # self.add_file_handler("fargs", format = "plain", level = logging.DEBUG)
        self.add_group("args", "fplain",
                       format_msg = lambda type_args: f"(({type_args[0]} args))\n{' '.join(type_args[1])}")
        args_to_ignore = {"genomes", "references", "clusters", "members"}
        mk_expanded_args_log = lambda args, params:'\n'.join([(k + ':\t' + str(v) +
                                                               ('' if v != vars(params)[k].default
                                                                else ' (default)'))
                                                              for k, v in vars(args).items()
                                                              if k not in args_to_ignore])
        self.add_group("args_expanded", "fplain",
                       format_msg = lambda args_params: ("\n((expanded args))\n" + \
                                                         mk_expanded_args_log(*args_params)))
        # self.add_format("format_args", "%(message)s",
        #                 format_msg = lambda type_args: f"(({type_args[0]} args))\n{' '.join(type_args[1])}")
        # self.add_stream_handler("format_args", "format_args")
        # self.add_group("args", "format_args")
        ## header
        self.add_stream_handler("sheader", format = "header", level = logging.DEBUG)
        self.add_file_handler("fheader", format = "header", level = logging.DEBUG)
        self.add_group("header", "sheader", "fheader")
    
    def move(self, config, args, check_path = True):
        fout = os.path.join(config.directory, f"{args.prefix}.log")
        self.update_filename(fout, check_path = check_path)
        

# logger = MINORgLogger(level = logging.DEBUG)
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

# logger2 = MINORgLogger(level = logging.INFO)
# logger2.filename
# logger2.header.info("full") ## writes to both file and console
# logger2.args.info("--gene AT1G01010 --acc ref") ## writes only to file
# logger2.update_filename("/mnt/chaelab/rachelle/tmp/testminorglogger.log")
# logger2.debug("test debug with default log (DEBUG)")
# logger2.debug.info("test debug with info call") ## prints something

# import sys, importlib
# sys.path.append("/mnt/chaelab/rachelle/scripts/minorgpy")
# import dynamic_logger

# importlib.reload(dynamic_logger)

# x = dynamic_logger.DynamicFileLogger("/mnt/chaelab/rachelle/tmp/testcustomformatter.log")
# x.add_format("plain", "%(message)s")
# x.add_format("format_args", "%(message)s", format_msg = lambda type_args: f"(({type_args[0]}))\n{' '.join(type_args[1])}", replace = True)
# x.add_stream_handler("splain", "plain")
# x.add_stream_handler("fmt_args", format = "format_args")
# x.add_group("args", "splain", format_msg = lambda type_args: f"(({type_args[0]}))\n{' '.join(type_args[1])}")
# x.add_group("fmtargs", "fmt_args")
# args = ["raw", ['--gene', 'AT1G12210']]
# ## for some reason all of the below defaulted to %(asctime)s %(message)s?
# x.args(args) ## incorrect formatting
# x.args.info(args) ## correct formatting
# x.fmtargs(args)
# x.fmtargs.info(args) ## same as above
# ## TODO: figure out which function is being called in Formatter
