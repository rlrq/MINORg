class MINORgError(Exception):
    def __init__(self, message):
        self.message = message
    def __repr__(self):
        return self.message
    def print_message(self):
        print(self.message)
        return


class MessageError(MINORgError):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class InputFormatError(MINORgError):
    def __init__(self, error_src = None, error_src_type = None, hint = None,
                 message = None):
        super().__init__( ( f"Error: The format of the input {error_src} is incorrect."
                            if message is None else message ) +
                            ('' if hint is None else '\nHint: ' + hint) )

class InvalidPath(MINORgError):
    def __init__(self, path):
        super().__init__( f"Error: {path} is not a valid path." )

class InvalidFile(MINORgError):
    def __init__(self, path):
        super().__init__( f"Error: {path} is not a file." )

class UnreadableFile(MINORgError):
    def __init__(self, path):
        super().__init__( f"Error: You do not have permission to read this file: {path}" )
