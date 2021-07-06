import typer

## functions for displaying messages

def print_indent(msg, lvl: int = 0, c: str = ' ', overwrite: bool = False):
    '''
    Print with indentation
    '''
    msg = (c * lvl) + str(msg)
    if overwrite:
        print(msg, end = '\r')
    else:
        print(msg)
    return

def print_overwrite_multi(msgs, lvl: int = 0, overwrite_last: bool = False, pause: float = 0):
    '''
    Print multiple lines with indentation
    '''
    from time import sleep
    for i, msg in enumerate(msgs):
        if i < len(msgs) - 1 or overwrite_last:
            print_indent(msg, lvl = lvl, overwrite = True)
        else:
            print_indent(msg, lvl = lvl, overwrite = False)
        sleep(pause)
    return

def make_print_preindent(initial_lvl: int = 0):
    '''
    Generate print_indent function w/ predefined initial indent level
    '''
    def print_preindent(msg, lvl: int = 0, **kwargs):
        return print_indent(msg, initial_lvl + lvl, **kwargs)
    return print_preindent

def make_print_overwrite_multi_preindent(initial_lvl: int = 0):
    '''
    Generate print_overwrite_multi function w/ predefined initial indent level
    '''
    def print_overwrite_multi_preindent(msgs, lvl: int = 0, **kwargs):
        return print_overwrite_multi(msgs, initial_lvl = lvl, **kwargs)
    return print_overwrite_multi_preindent

