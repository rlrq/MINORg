## suppress SearchIO.parse warning
from Bio import SearchIO, BiopythonWarning
import warnings

## wrapper for SearchIO.parse that suppresses the specific Biopython warning: "/usr/local/lib/python3.6/dist-packages/Bio/SearchIO/_legacy/NCBIStandalone.py:45: BiopythonWarning: Parsing BLAST plain text output file is not a well supported functionality anymore. Consider generating your BLAST output for parsing as XML or tabular format instead." Depending on your Biopython version, this warning may be erroneously raised even if you are parsing an xml or tab file.
def searchio_parse(*args, **kwargs):
    def generator():
        with warnings.catch_warnings():
            if args[1] in {"blast-xml", "blast-tab"}:
                warnings.filterwarnings("ignore", "Parsing BLAST plain text output file is not a well supported functionality anymore. Consider generating your BLAST output for parsing as XML or tabular format instead.")
            for x in SearchIO.parse(*args, **kwargs):
                yield x
    return generator()
