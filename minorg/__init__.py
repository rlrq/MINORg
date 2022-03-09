"""
gRNA set generator for multiple genes in multiple genomes.

https://rlrq.github.io/MINORg
"""

import logging

_logging_level = logging.DEBUG

__version__ = "0.2.1.4alpha1"

class MINORgWarning(Warning):
    """
    MINORg warning.
    
    MINORg should use this warning (or subclasses of it), making it easy to
    silence all its warnings should you wish to.
    
    >>> import warnings
    >>> from minorg import MINORgWarning
    >>> warnings.simplefilter('ignore', MINORgWarning)
    
    Consult the warnings module documentation for more details.
    """

def _warning(message, category, filename, lineno, file=None, line=None):
    print(f"{filename}:{lineno}: {category.__name__}: {message}")

class MINORgError(Exception):
    """
    MINORg exception.
    
    MINORg should use this exception (or subclasses of it) to raise
    non-typer or non-click exceptions without printing traceback.
    (i.e. Traceback is reserved for errors of unknown origin.)
    """
    def __init__(self, message):
        self.message = message
        self.prefix = "MINORgError"
    def __repr__(self):
        return self.message
    def print_message(self):
        print(f"{self.prefix}: {self.message}")
        return

# class MessageError(MINORgError):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)

class InputFormatError(MINORgError):
    def __init__(self, error_src = None, error_src_type = None, hint = None,
                 message = None):
        super().__init__( ( f"The format of the input {error_src} is incorrect."
                            if message is None else message ) +
                            ('' if hint is None else '\nHint: ' + hint) )

class InvalidPath(MINORgError):
    def __init__(self, path):
        super().__init__( f"{path} is not a valid path." )

class InvalidFile(MINORgError):
    def __init__(self, path):
        super().__init__( f"{path} is not a file." )

class UnreadableFile(MINORgError):
    def __init__(self, path):
        super().__init__( f"You do not have permission to read this file: {path}" )

