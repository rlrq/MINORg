"""GFF3 class"""

import itertools
from typing import Union, Generator, Callable

from minorg.index import IndexedFile

from minorg.functions import (
    assign_alias, make_local_print,
    non_string_iter
)

from minorg.display import print_indent

#################
##  GFF_MANIP  ##
##  (CLASSES)  ##
#################

def get_recursively(d, default, *keys):
    def helper(d, keys):
        key = keys[0]
        if key in d:
            if len(keys) == 1: return d[key]
            else: return helper(d[key], keys[1:])
        else:
            return default
    return helper(d, keys)

class GFF:
    """
    Representation of GFF3 file.
    
    Large files can be indexed instead of read to memory by using ``memsave=True``.
    
    >>> my_gff = GFF('/path/to/large_gff.gff', memsave = True)
    
    Attributes:
        _fname (str): path to GFF3 file
        _fmt (str): GFF3 file format
        _data (list): stores annotation data as list of :class:`minorg.annotation.Annotation` objects.
            Not used if fname is not None but memsave=True.
        _string (str): stores raw string data if string!=None
        _attr_mod (dict): attribute modification mapping for non-standard attribute field names
        _attr_fields (dict): full attribute field name mapping
        _quiet (bool): print only essential messages
        _kwargs: stores additional arguments when parsing GFF3 entries 
            to :class:`minorg.annotation.Annotation` objects
        _seqids (dict): stores order of seqids/molecules/chromosomes
        _chunk_lines (int): number of lines between each stored position (used for indexing when memsave=True)
        _indexed_file (:class:`minorg.index.IndexedFile`): indexed GFF file
        _memsave (bool): index file instead of reading data to memory
    
    .. automethod:: __iter__
    """
    def __init__(self, fname = None, data = [], string = None, attr_mod = {}, genetic_code = 1, fmt = None,
                 quiet = False, memsave = False, chunk_lines = 1000, **kwargs):
        """
        Create a GFF object.
        
        Arguments:
            fname (str): optional, path to GFF3 file or BED file generated using gff2bed
            data (list): optional, list of :class:`minorg.annotation.Annotation` objects
            string (str): optional, string contents in the format of a GFF3 file or 
                BED file generated using gff2bed
            attr_mod (dict): optional, dictionary of mapping for non-standard attribute field names
            genetic_code (int or str): NCBI genetic code name or number
            fmt (str): optional, valid values: BED, GFF, GFF3.
                If not provided and fname != None, will be inferred from fname extension.
            quiet (bool): print only essential messages
            memsave (bool): index file instead of reading data to memory
            chunk_lines (int): number of lines between each stored line index (default=1000)
            **kwargs: additional arguments when parsing GFF3 entries 
                to :class:`minorg.annotation.Annotation` objects
        """
        self._fname = fname
        if fmt or not fname:
            self._fmt = (fmt if fmt else None)
        else:
            self._infer_fmt()
        self._data = data
        self._string = string
        self._attr_mod = attr_mod
        self._attr_fields = {"all": {"ID": "ID", "Name": "Name", "Alias": "Alias", "Parent": "Parent",
                                     "Target": "Target", "Gap": "Gap", "Derives_from": "Derives_from",
                                     "Note": "Note", "Dbxref": "Dbxref", "Ontology_term": "Ontology_term",
                                     "Is_circular": "Is_circular"}}
        self._attr_fields_inv = {}
        self._genetic_code = genetic_code
        self._quiet = quiet
        self._kwargs = kwargs
        # self._seqids = {}
        self._seqids = {}
        self._chunk_lines = chunk_lines
        self._indexed_file = None
        self._memsave = memsave
        self.update_attr_fields()
        if not self._memsave:
            self.parse()
        elif self.has_file():
            self.index()
    
    def __iter__(self):
        """
        Read from data stored in memory if self._data is not empty, else read directly from file.
        
        Yields
        ------
        :class:`minorg.annotation.Annotation`
        """
        if self._data:
            for entry in self._data:
                yield entry
        else:
            mk_annotation = self.make_annotation_from_str_gen(strip_newline = True)
            for entry in self.iter_raw():
                yield mk_annotation(entry)
        return
    
    def __len__(self):
        """
        Length of data stored in memory. Returns 0 if no data or file is only indexed.
        """
        return len(self._data)
    
    def _infer_fmt(self):
        """
        Infer file format from file name and store format at self._fmt.
        """
        self._fmt = ("BED" if self._fname.split('.')[-1].upper() == "BED" else "GFF3")
    
    def empty_copy(self, other = None) -> Union['GFF', None]:
        """
        Shallow copy self's attributes (BUT NOT DATA) to another :class:`minorg.annotation.GFF` object.
        
        If ``other=None``, create new :class:`minorg.annotation.GFF` object, 
        copy attributes to it, and return the new object.
        
        Arguments:
            other (:class:`minorg.annotation.GFF`): optional.
                If not provided, creates a new :class:`minorg.annotation.GFF` object with copied attributes.
        
        Returns
        -------
        :class:`minorg.annotation.GFF`
            If ``other=None``
        """
        if other:
            other._attr_mod = self._attr_mod
            other._genetic_code = self._genetic_code
            other._quiet = self._quiet
            other._chunk_lines = self._chunk_lines,
            other._memsave = self._memsave
            other._kwargs = self._kwargs
        else:
            return GFF(attr_mod = self._attr_mod, genetic_code = self._genetic_code,
                       quiet = self._quiet, chunk_lines = self._chunk_lines,
                       memsave = self._memsave, **self._kwargs)
    
    def is_bed(self) -> bool:
        """
        Whether input file format is BED.
        """
        if not self._fmt: return False
        return self._fmt.upper() == "BED"
    def is_gff(self) -> bool:
        """
        Whether input file format is GFF.
        """
        if not self._fmt: return False
        return self._fmt.upper() in {"GFF", "GFF3"}
    def has_file(self) -> bool:
        """
        Whether self is associated with a file.
        """
        return self._fname is not None
    def index(self, chunk_lines = None) -> None:
        """
        Index file.
        
        Arguments:
            chunk_lines (int): optional, number of lines between each indexed line.
                If not provided, defaults to self._chunk_lines.
        """
        if chunk_lines:
            self._chunk_lines = chunk_lines
        if self.has_file():
            self._indexed_file = IndexedFile(self._fname, chunk_lines = self._chunk_lines,
                                             skip = lambda line: (tuple(line[:1]) == ('#',) or line == '\n'))
    def is_indexed(self) -> bool:
        """
        Whether file has been indexed.
        """
        return bool(self._indexed_file)
    
    def make_annotation_from_str_gen(self, strip_newline = True) -> Callable[[str], 'Annotation']:
        """
        Create function to parse raw string entries into :class:`minorg.annotation.Annotation` objects 
        based on inferred data format (GFF3 or BED generated by gff2bed)
        
        Arguments:
            strip_newline (bool): whether to strip newline if it present in string entry when read from file
        
        Returns
        -------
        func
        """
        if strip_newline:
            split_line = lambda entry: entry.replace('\n', '').split('\t')
        else:
            split_line = lambda entry: entry.split('\t')
        if self.is_bed():
            def mk_annotation(entry):
                entry = split_line(entry)
                gff_fmt = [entry[0], entry[6], entry[7], str(int(entry[1]) + 1), entry[2],
                           entry[4], entry[5], entry[8], entry[9]]
                return Annotation(gff_fmt, self, **self._kwargs)
        else: ## GFF3
            def mk_annotation(entry):
                entry = split_line(entry)
                return Annotation(entry, self, **self._kwargs)
        return mk_annotation
    
    def iter_raw(self) -> Generator[str, None, None]:
        """
        Yields raw string read from file.
        """
        if self._indexed_file is not None:
            for entry in self._indexed_file:
                yield entry
        elif self._fname is not None:
            with open(self._fname, 'r') as f:
                for entry in f:
                    if tuple(entry[:1]) == ('#',) or entry == '\n': continue
                    yield entry
        elif self._string is not None:
            for entry in self._string.split('\n'):
                if tuple(entry[:1]) == ('#',) or entry == '\n': continue
                yield entry
        return
    
    def parse(self) -> None:
        """
        Read file stored at self._fname to memory.
        """
        if self._fname is not None or self._string is not None:
            ## start parsing
            self._data = []
            for entry in self:
                self.add_entry(entry)
            self.index_seqids()
    
    def read_file(self, fname) -> None:
        """
        Read file to memory.
        
        Arguments:
            fname (str): path to GFF3 file
        """
        self._fname = fname
        self._infer_fmt()
        self.parse()
    
    def index_seqids(self) -> None:
        """
        Store order of seqid.
        """
        self._seqids = {}
        index = 0
        curr = None
        for entry in self:
            if entry.seqid != curr and entry.seqid not in self._seqids:
                curr = entry.seqid
                index += 1
                self._seqids[entry.seqid] = index
        return
    
    def sort(self) -> None:
        """
        Sort data stored in memory. Does NOT sort indexed files.
        
        Sort key: seqid, start, -end, source, feature type, str of attributes with standardised field names
        """
        if self._seqids:
            seqid_key = lambda seqid: self._seqids[seqid]
        else:
            seqid_key = lambda seqid: seqid
        self._data.sort(key = lambda entry: (seqid_key(entry.seqid), entry.start, -entry.end,
                                             entry.source, entry.feature, entry.attributes.standardise_fields()))
    
    def add_entry(self, gff_entry, duplicate_check = False) -> None:
        """
        Add :class:`minorg.annotation.Annotation` object to self's data at self._data.
        
        Arguments:
            gff_entry (:class:`~minorg.annotation.Annotation`): required, Annotation object to add
            duplicate_check (bool): check for duplicates and only add ``gff_entry`` if not already in data
        """
        if duplicate_check:
            for entry in self:
                if gff_entry == entry: return
        self._data.append(gff_entry)
    
    def update_attr_fields(self, attr_mod = None) -> None:
        """
        Update attribute fields w/ user-provided attribute field modification dictionary.
        (Otherwise, use self._attr_fields)
        
        Arguments:
            attr_mod (dict): optional, dictionary of attribute modifications
        """
        if attr_mod is None: attr_mod = self._attr_mod
        for feature, mods in attr_mod.items():
            if feature in self._attr_fields:
                for canon, atypical in mods.items():
                    self._attr_fields[feature][canon] = atypical
            else:
                self._attr_fields[feature] = mods
        self._attr_fields_inv = self.invert_attr_fields()
        return
    
    def invert_attr_fields(self) -> dict:
        """
        Generate {<feature>: {<NONSTANDARD attribute field name>: <STANDARD attribute field name>}}
        mapping from self._attr_fields (which is in format
        {<feature>: {STANDARD attribute field name>: <NONSTANDARD attribute field name>}}
        
        Returns
        -------
        dict
            Of reversed attribute field name mapping
        """
        output = {feature: {v: k for k, v in attr_map.items()}
                 for feature, attr_map in self._attr_fields.items()}
        return output
    
    def get_features(self, *feature_types, index = False):
        """
        Get entries of specific feature types.
        
        Arguments:
            *feature_types (str): GFF3 feature types to retrieve
            index (bool): return line number (index) instead of Annotation objects
        
        Returns
        -------
        list
            Of :class:`minorg.annotation.Annotation` objects if index=False
        list
            Of int line number of entries if index=True
        """
        feature_types = set(feature_types)
        if index:
            return [i for i, ann in enumerate(self) if ann.type in feature_types]
        else:
            return [ann for ann in self if ann.type in feature_types]
    
    def get_subfeatures(self, *feature_ids, feature_types = [], index = False):
        """
        Get all features that are subfeatures of user-provided feature_ids.
        
        Arguments:
            *feature_ids (list): list of feature IDs
            feature_types (str): feature type(s) to retain
            index (bool): return line number of feature instead of :class:`minorg.annotation.Annotation` objects
        
        Returns
        -------
        list
            Of :class:`minorg.annotation.Annotation` objects if index=False
        list
            Of int line number of entries if index=True
        """
        feature_ids = set(feature_ids)
        if type(feature_types) is str:
            feature_types = {feature_types}
        else:
            feature_types = set(feature_types)
        indices = [i for i, entry in enumerate(self) if
                   ((not feature_types or entry.type in feature_types) and
                    entry.has_attr("Parent", feature_ids))]
        if index: return indices
        else: return self.get_i(*indices, output_list = True)
    
    def get_subfeatures_full(self, *feature_ids, feature_types = [], index = False, preserve_order = True):
        """
        Get all features that are subfeatures of user-provided feature_ids
        AND subfeatures of those subfeatures, until there are no sub-sub...sub-features left.
        
        Arguments:
            *feature_ids (str): parent feature IDs
            feature_types (list of str): GFF3 feature type(s) to retrieve
            index (bool): return line number of feature instead of :class:`minorg.annotation.Annotation` objects
            preserve_order (bool): sort output by line number (preserve original order)
        
        Returns
        -------
        list
            Of :class:`minorg.annotation.Annotation` objects if index=False
        list
            Of int line number of entries if index=True
        """
        printi = make_local_print(quiet = self._quiet, printf = print_indent)
        feature_types = set(feature_types) if non_string_iter(feature_types) else {feature_types}
        output_indices = []
        curr_indices = []
        iteration_n = 0
        while True:
            iteration_n += 1
            printi(f"Subsetting iteration: {iteration_n}", overwrite = True)
            curr_indices = self.get_subfeatures(*feature_ids, index = True)
            feature_ids = set(itertools.chain(*[entry.get_attr("ID")
                                                for entry in self.get_i(*curr_indices, output_list = True)]))
            # feature_ids = []
            # for i in curr_indices:
            #     try:
            #         feature_ids.append(self.get_i(i).get_attr("ID"))
            #     except Exception as e:
            #         print(i, self.get_i(i, verbose = True))
            #         raise e
            # feature_ids = set(itertools.chain(*[self.get_i(i).get_attr("ID") for i in curr_indices]))
            output_indices.extend(curr_indices)
            if not feature_ids: break
        printi("", overwrite = True) ## clear printed progress
        if feature_types:
            output_indices = [i for i in output_indices if self.get_i(i).type in feature_types]
        ## prepare output
        if index: return sorted(output_indices)
        else: return self.get_i(*output_indices, sort = preserve_order, output_list = True)
    
    def get_features_and_subfeatures(self, *feature_ids, index = False, full = True,
                                     preserve_order = True):
        """
        Gets features w/ feature_ids AND subfeatures of those features.
        If full = True, executes get_subfeatures_full for subfeature discovery, else get_subfeatures
        
        Arguments:
            feature_ids (str): parent feature IDs
            index (bool): return line number of feature instead of :class:`minorg.annotation.Annotation` objects
            full (bool): return feature(s) as well as its/their subfeatures
            preserve_order (bool): sort output by line number (i.e. preserve original order)
        
        Returns
        -------
        list
            Of :class:`minorg.annotation.Annotation` objects if index=False
        list
            Of int line number of entries if index=True
        """
        features_indices = self.get_id(*feature_ids, index = True, output_list = True)
        if full: subfeatures_indices = self.get_subfeatures_full(*feature_ids, index = True)
        else: subfeatures_indices = self.get_subfeatures(*feature_ids, index = True)
        ## prepare output
        final_features = sorted(features_indices + subfeatures_indices)
        if index: return final_features
        else: return self.get_i(*final_features,
                                output_list = True,
                                sort = preserve_order) ## it's already sorted by it won't hurt to raise it
    
    def get_i_raw(self, *indices, strip_newline = True, output_list = True,
                  sort = True) -> Union[list, str, None]:
        """
        Get raw string of entry/entries by line index.
        
        Arguments:
            indices (int): line numbers (indices) of entries to retrieve
            strip_newline (bool): remove newline from returned lines
            output_list (bool): return list even if ony one line number is provided
            sort (bool): sort output by line number
        
        Returns
        -------
        list
            Of str of entries if ``output_list=True`` or multiple lines were requested
        str
            Of entry if ``output_list=False`` and only one line was requested
        None
            If ``output_list=False`` and the specified line does not exist
        """
        lines = self._indexed_file.get_line(*indices, strip_newline = strip_newline, output_fmt = list)
        if len(indices) == 1 and not output_list:
            if lines: return lines[0]
            else: return None
        else: return lines
    
    def get_i(self, *indices, output_list = False, sort = True) -> Union[list, 'Annotation', None]:
        """
        Get :class:`~minorg.annotation.Annotation` of GFF entry/entries by line index.
        
        Arguments:
            indices (int): line number(s) (indices) of entries to retrieve
            output_list (bool): return list even if ony one line number is provided
            sort (bool): sort output by line number
        
        Returns
        -------
        list
            Of :class:`~minorg.annotation.Annotation` if ``output_list=True`` or multiple indices were requested
        :class:`~minorg.annotation.Annotation`
            If ``output_list=False`` and only one index was requested
        None
            If ``output_list=False`` and the specified line does not exist
        """
        if sort:
            indices = sorted(indices)
        if self._data:
            output = [self._data[i] for i in indices]
        elif self.is_indexed():
            raw_entries = self.get_i_raw(*indices, strip_newline = True)
            ## newline already stripped by get_i_raw. No need to strip them again.
            mk_annotation = self.make_annotation_from_str_gen(strip_newline = False)
            output = [mk_annotation(entry) for entry in raw_entries]
        else:
            output = []
        if len(indices) == 1 and not output_list:
            if output: return output[0]
            else: return None
        else: return output
    
    def get_id(self, *feature_ids, index = False, output_list = False, preserve_order = True):
        ## if even if output_list is only used when len(indices) == 1 or == 0 AND type(feature_ids) is str.
        ##  Always returns list otherwise.
        """
        Get :class:`~minorg.annotation.Annotation` of GFF entry/entries by feature ID.
        
        Arguments:
            *feature_ids (str): feature ID(s)
            index (bool): return line number(s) of feature(s) 
                instead of :class:`minorg.annotation.Annotation` object(s)
            output_list (bool): return list even if ony one feature ID is provided
            preserve_order (bool): sort output by line number (preserve original order)
        
        Returns
        -------
        list
            If ``output_list=True`` or more than one feature ID was provided.
            List of :class:`minorg.annotation.Annotation` objects if ``index=False``.
            List of int if ``index=True``.
        :class:`~minorg.annotation.Annotation`
            If ``output_list=False`` and only one feature ID was provided and ``index=False``
        int
            If ``output_list=False`` and only one feature ID was provided and ``index=True``
        None
            If ``output_list=False`` and no feature with the specified feature ID can be found
        """
        if type(feature_ids) is str:
            feature_ids = {feature_ids}
        else:
            feature_ids = set(feature_ids)
        indices = [i for i, entry in enumerate(self) if entry.has_attr("ID", feature_ids)]
        if not output_list and len(indices) <= 1:
            if not indices: return None
            elif index: return indices[0]
            else: return self.get_i(*indices, output_list = False)
        elif index: return indices
        else: return self.get_i(*indices, sort = preserve_order, output_list = True)
    
    def subset(self, feature_ids = None, feature_types = None, subfeatures = True, preserve_order = True):
        """
        Subset data by feature ID (feature_ids) and/or feature type (feature_types)
        and generate new GFF object from them.
        
        Arguments:
        
        """
        new_gff = self.empty_copy()
        ## start subsetting
        if feature_ids:
            if isinstance(feature_ids, str):
                feature_ids = [feature_ids]
            if subfeatures:
                indices = self.get_features_and_subfeatures(*feature_ids, index = True)
            else:
                indices = self.get_id(*feature_ids, index = True, output_list = True)
            new_gff._data = self.get_i(*indices, sort = preserve_order)
            if feature_types:
                return new_gff.subset(feature_types = feature_types)
            else:
                return new_gff
        elif feature_types:
            if isinstance(feature_types, str):
                feature_types = [feature_types]
            data = self.get_features(*feature_types, index = False)
            new_gff._data = data
            if sort:
                new_gff.sort()
            return new_gff
        else:
            return self
    
    def write(self, fout = None, entries = None, **kwargs):
        '''
        Writes entries to file
        '''
        if not fout:
            fout = self._fname
        if not entries:
            entries = self
        with open(fout, "w+") as f:
            f.write('\n'.join([entry.generate_str(**kwargs) for entry in entries]) + '\n')
    
    def write_i(self, fout, indices, **kwargs):
        '''
        Executes get_i, then writes output to file
        '''
        self.write(fout, self.get_i(*indices), **kwargs)
        return
    
    def write_id(self, fout, feature_ids, **kwargs):
        '''
        Executes get_id, then writes output to file
        '''
        self.write(fout, self.get_id(*feature_ids), **kwargs)

class Annotation:
    """
    Representation of GFF3 annotation/feature
    """
    def __init__(self, entry, gff, **kwargs):
        """
        Create an Annotation object.
        """
        self._gff = gff
        try:
            self.molecule = entry[0]
            self.seqid = entry[0]
            self.source = entry[1]
        except Exception as e:
            print("entry:", entry)
            print("entry[0]:", entry[0])
            raise e
        self.feature = entry[2]
        self.type = entry[2]
        self.start = int(entry[3]) - 1
        # self.start0 = self.start ## 0-indexed start
        self.end = int(entry[4])
        self.score = entry[5] if entry[5] == '.' else float(entry[5])
        self.strand = entry[6]
        self.phase = entry[7] if entry[7] == '.' else int(entry[7])
        self.attributes = Attributes(entry[8], gff, self, **kwargs)
        self._fields = {"molecule": self.molecule,
                        "seqid": self.seqid,
                        "source": self.source,
                        "feature": self.feature,
                        "type": self.type,
                        "start": self.start,
                        "end": self.end,
                        "score": self.score,
                        "strand": self.strand,
                        "phase": self.phase,
                        "attributes": self.attributes}
    def __eq__(self, other):
        return (self.generate_gff(standardise = True) == other.generate_gff(standardise = True))
    @property
    def plus(self):
        return self.strand == '+'
    def generate_attr(self, original = True, fields = None):
        if original: return self.attributes._raw
        else: return self.attributes.standardise_fields()
    def generate_str(self, fmt = "GFF"):
        if fmt.upper() in {"GFF", "GFF3"}:
            output = self.generate_gff()
        elif fmt.upper() in {"BED"}:
            output = self.generate_bed()
        return '\t'.join(map(str, output))
    def generate_gff(self, standardise = False):
        return list(map(str, [self.seqid, self.source, self.type, self.start + 1, self.end,
                              self.score, self.strand, self.phase,
                              self.generate_attr(original = (not standardise))]))
    def generate_bed(self, standardise = False):
        return list(map(str, [self.seqid, self.start, self.end, self.attributes.get("ID", fmt = str),
                              self.score, self.strand, self.source, self.type, self.phase,
                              self.generate_attr(original = (not standardise))]))
    def get(self, *fields):
        return [self.f_dict[field] for field in fields]
    
    ## wrappers for Attribute methods
    def get_attr(self, a, **kwargs): return self.attributes.get(a, **kwargs)
    def is_attr(self, a, val, **kwargs): return self.attributes.is_attr(a, val, **kwargs)
    def has_attr(self, a, vals, **kwargs): return self.attributes.has_attr(a, vals, **kwargs)

class Attributes:
    """
    Reprsentation of GFF3 feature attributes
    """
    def __init__(self, val, gff, entry, field_sep_inter = ';', field_sep_intra = ',',
                 **for_dummy_gff):
        self._entry = entry
        self._gff = gff if gff is not None else GFF(**for_dummy_gff)
        self._raw = val
        self._sep_inter = field_sep_inter
        self._sep_intra = field_sep_intra
        self._data = self._parse()
    def __eq__(self, other):
        return (self.standardise_fields(return_str = False) == other.standardise_fields(return_str = False))
    def __repr__(self):
        return self._raw
    def __str__(self):
        return self._raw
    def _parse(self):
        import re
        ## if multiple separate entries for same field (e.g. 'Parent=abc.1;Parent=abc.2'), parse properly
        attributes_l = [x.split('=') for x in re.findall(f"[^{self._sep_intra}{self._sep_inter}=]+=.+?(?=[^{self._sep_intra}{self._sep_inter}=]+=.+?|$)", self._raw)]
        attributes = {}
        for attribute in attributes_l:
            attributes[attribute[0]] = attributes.get(attribute[0], []) + \
                                       re.search(f'^.+?(?=[{self._sep_intra}{self._sep_inter}]?$)',
                                                 attribute[1]).group(0).split(self._sep_intra)
        return attributes
    def get(self, a, fmt = list):
        feature = self._entry.feature
        attr_fields = self._gff._attr_fields
        if feature in attr_fields and a in attr_fields[feature]:
            mapped_field = attr_fields[feature][a]
        elif a in attr_fields["all"]:
            mapped_field = attr_fields["all"][a]
        else:
            mapped_field = a
        iterable = self._data.get(mapped_field, [])
        ## format return
        if fmt is str and iterable: return self._sep_intra.join(iterable)
        elif fmt is str and not iterable: return '.'
        else: return fmt(iterable)
    def is_attr(self, a, val):
        '''
        Checks if 'val' is a value of attribute 'a'
        '''
        return val in self.get(a)        
    def has_attr(self, a, vals):
        '''
        Checks if at least 1 value in 'vals' is also a value of attribute 'a'
        '''
        return (set(self.get(a)) & set(vals)) ## checks intersection != empty set
    def standardise_fields(self, return_str = True):
        def get_mod(field):
            feature_mod = get_recursively(self._gff._attr_fields_inv, field, self._entry.feature, field)
            if feature_mod != field: return feature_mod
            else: return get_recursively(self._gff._attr_fields_inv, field, "all", field)
        if return_str:
            return ';'.join([f"{get_mod(k)}={'.'.join(v)}" for k, v in self._data.items()])
        else:
            return {get_mod(k): v for k, v in self._data.items()}


# #############################
# ##  EXTRACT_GFF_FEAUTURES  ##
# #############################

# def extract_features_and_subfeatures(gff_fname, feature_ids, fout, quiet = False, attr_mod = {},
#                                      fin_fmt = "GFF", fout_fmt = "GFF", memsave = True,
#                                      preserve_order = True):
#     printi = make_local_print(quiet = quiet)
#     printi("Reading files")
#     gff = GFF(fname = gff_fname, fmt = fin_fmt, quiet = quiet, memsave = memsave, attr_mod = attr_mod)
#     # feature_ids = splitlines(feature_id_fname)
#     printi("Retrieving features")
#     gff_subset = gff.subset(feature_ids = feature_ids, subfeatures = True, preserve_order = preserve_order)
#     # output_features = gff.get_features_and_subfeatures(feature_ids, index = True, full = True)
#     printi("Writing features")
#     gff_subset.write(fout)
#     return


# #################
# ##  ANN MANIP  ##
# #################

# def subset_ann(gff_beds, ids, fout = None, mk_tmpf_name = None, fout_fmt = "BED", attr_mods = {},
#                memsave = True, preserve_order = True):
#     if mk_tmpf_name is None:
#         import tempfile
#         mk_tmpf_name = lambda x: tempfile.mkstemp()[1]
#     gff_beds = assign_alias(gff_beds)
#     gff_subsets = {}
#     for alias, gff_bed in gff_beds.items():
#         gff_subset = mk_tmpf_name(alias)
#         gff_subsets[alias] = gff_subset
#         fmt = "BED" if gff_subset.split('.')[-1].upper() == "BED" else "GFF"
#         extract_features_and_subfeatures(gff_bed, ids, gff_subset, quiet = True, memsave = memsave,
#                                          attr_mod = attr_mods.get(alias), fin_fmt = fmt, fout_fmt = fout_fmt,
#                                          preserve_order = preserve_order)
#     if fout:
#         cat_files(tuple(gff_subsets.values()), fout, remove = True)
#     else:
#         return gff_subsets

# def reduce_ann(gff_beds, ids, fout = None, mk_tmpf_name = None, fout_fmt = "BED", attr_mods = {},
#                memsave = True, preserve_order = True):
#     return subset_ann(gff_beds, ids, fout = fout, mk_tmpf_name = mk_tmpf_name,
#                       fout_fmt = fout_fmt, memsave = memsave, attr_mods = attr_mods, sort = preserve_order)
