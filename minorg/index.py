import copy
from Bio.Seq import Seq
from pyfaidx import Fasta

from typing import Union

###################
##  INDEX FILES  ##
##   (CLASSES)   ##
###################

class IndexedFile:
    """
    File indexed by line.
    
    Stores position of every ``chunk_lines`` th line.
    During retrieval, jumps to largest ``chunk_lines`` th line before or equal to desired line number
    and iterates from that line.
    """
    def __init__(self, filename, chunk_lines = 10000, skip = lambda x: False):
        """
        Create IndexedFile object.
        
        Arguments:
            filename (str): required, path to file to index
            chunk_lines (int): number of lines between each stored position.
                The larger the number the fewer number of positions that will have to be stored in memory
                BUT the slower the lookup.
            skip (func): function that accepts line (str) from file and outputs whether to skip that line. 
                Skipped lines are effectively hidden from self and do not contribute to line count 
                during indexing or retrieval. (default=lambda x: False)
        """
        self.filename = filename
        self.chunk_lines = chunk_lines
        self.indices = None
        ## if skip(line) returns True, the line will not be indexed
        ##  and will be ignored (and will not add to line count) when using get_line
        self.skip = skip
        self.index()
    def __iter__(self):
        """Iterates through lines and yields contents of line if ``not self.skip(line)``."""
        with open(self.filename, 'r') as f:
            for line in f:
                if self.skip(line): continue
                else: yield line
    def index(self, chunk_lines = None) -> None:
        """
        Index file.
        
        Arguments:
            chunk_lines (int): see __init__.
        """
        if chunk_lines is not None:
            self.chunk_lines = chunk_lines
        indices = [0]
        with open(self.filename, 'r') as f:
            i = 0
            while True:
                line = f.readline()
                if not line: break
                if self.skip(line): continue
                i += 1
                if i % self.chunk_lines == 0:
                    indices.append(f.tell())
        self.indices = indices
    def get_line(self, *indices, strip_newline = False, output_fmt = None) -> Union[str, list]:
        """
        Retrieve line content by line number(s).
        
        Arguments:
            *indices (int): line numbers of lines to retrieve
            strip_newline (bool): remove newline from returned lines
            output_fmt (type): output format. Valid value: ``list``.
                If ``output_fmt = list``, returns list even if ``len(indices) == 1``.
        
        Returns
        -------
        str
            If ``len(indices) == 1 and output_fmt != list``
        list
            If ``len(indices) > 1 or output_fmt == list``
        """
        single_i = len(indices) == 1
        ## remove invalid indices (<0)
        is_valid = lambda i: i >= 0
        invalid = [i for i in indices if not is_valid(i)]
        if invalid: print(f"Invalid indices (will be ignored): {','.join(map(str, sorted(invalid)))}")
        indices = [i for i in indices if is_valid(i)]
        ## sort indices and split into bins according to self._chunk_lines
        bins = {}
        for i in sorted(indices):
            i_bin = i // self.chunk_lines
            i_offset = i % self.chunk_lines
            bins[i_bin] = bins.get(i_bin, []) + [i_offset]
        ## jump to each bin and iterate within it
        output = []
        with open(self.filename, 'r') as f:
            for i_bin, i_offsets in bins.items():
                last_i_offset = i_offsets[-1]
                i_offsets = set(i_offsets)
                ## jump to first line of bin
                f.seek(self.indices[i_bin])
                i = 0
                for line in f:
                    if self.skip(line): continue
                    if i in i_offsets:
                        output.append(line)
                        ## exit bin if all desired lines in bin have been visited
                        if i == last_i_offset: break
                    i += 1
        if strip_newline: output = [line.replace('\n', '') for line in output]
        return output if ((not single_i) or (output_fmt is list)) else output[0]

## basically pyfaidx.Fasta class, but returns Bio.Seq.Seq item when sliced and filename for __repr__
class IndexedFasta(Fasta):
    """
    pyfaidx.Fasta class object wrapped to return Bio.Seq.Seq when sliced.
    """
    def __init__(self, fasta, *args, **kwargs):
        if isinstance(fasta, self.__class__):
            self.__dict__ = copy.copy(fasta.__dict__) ## no need for deep copy
        elif isinstance(fasta, Fasta):
            super().__init__(fasta.filename, *args, **kwargs)
        else:
            super().__init__(fasta, *args, **kwargs)
    def __repr__(self):
        return f"IndexedFasta({self.filename})"
    def __str__(self):
        return self.filename
    def get_seq(self, *args, **kwargs) -> Seq:
        pyfaidx_seq = super().get_seq(*args, **kwargs)
        return Seq(pyfaidx_seq.seq)
    def get_spliced_seq(self, *args, **kwargs) -> Seq:
        pyfaidx_seq = super().get_spliced_seq(*args, **kwargs)
        return Seq(pyfaidx_seq.seq)

