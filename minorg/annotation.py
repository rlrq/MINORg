import itertools

from minorg.index import IndexedFile

from minorg.functions import (
    assign_alias, make_local_print
)

#################
##  GFF_MANIP  ##
##  (CLASSES)  ##
#################
    
class GFF:
    
    def __init__(self, fname = None, data = [], attr_mod = {}, fmt = None, quiet = False,
                 memsave = False, chunk_lines = 1000, **kwargs):
        self._fmt = (fmt if fmt
                     else None if not fname
                     else ("BED" if fname.split('.')[-1].upper() == "BED" else "GFF3"))
        self._fname = fname
        self._data = data
        self._attr_mod = attr_mod
        self._attr_fields = {"all": {"ID": "ID", "Name": "Name", "Alias": "Alias", "Parent": "Parent",
                                     "Target": "Target", "Gap": "Gap", "Derives_from": "Derives_from",
                                     "Note": "Note", "Dbxref": "Dbxref", "Ontology_term": "Ontology_term",
                                     "Is_circular": "Is_circular"}}
        self._attr_fields_inv = {}
        self._quiet = quiet
        self._kwargs = kwargs
        self._molecules = {}
        self._chunk_lines = chunk_lines
        self._indexed_file = None
        self.update_attr_fields()
        if not memsave:
            self.parse()
        elif self.has_file():
            self.index()
    
    def __iter__(self):
        if self._data:
            for entry in self._data:
                yield entry
        else:
            mk_annotation = self.make_annotation_from_str_gen(strip_newline = True)
            for entry in self.iter_raw():
                yield mk_annotation(entry)
        return
    
    def __len__(self):
        return len(self._data)
    
    def is_bed(self):
        return self._fmt.upper() == "BED"
    def is_gff(self):
        return self._fmt.upper() in {"GFF", "GFF3"}
    def has_file(self):
        return self._fname is not None
    def index(self, chunk_lines = None):
        if chunk_lines:
            self._chunk_lines = chunk_lines
        if self.has_file():
            self._indexed_file = IndexedFile(self._fname, chunk_lines = self._chunk_lines,
                                             skip = lambda line: (tuple(line[:1]) == ('#',) or line == '\n'))
    def is_indexed(self):
        return bool(self._indexed_file)
    
    def make_annotation_from_str_gen(self, strip_newline = True):
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
    
    def iter_raw(self):
        if self._indexed_file is not None:
            for entry in self._indexed_file:
                yield entry
        elif self._fname is not None:
            with open(self._fname, 'r') as f:
                for entry in f:
                    if tuple(entry[:1]) == ('#',) or entry == '\n': continue
                    yield entry
        return
    
    def parse(self):
        if self._fname is not None:
            ## start parsing
            self._data = []
            for entry in self:
                self.add_entry(entry)
            self.index_molecules()
    
    def read_file(self, fname):
        self._fname = fname
        self.parse()
    
    def index_molecules(self, reset = True):
        if reset:
            self._molecules = {}
        index = 0
        curr = None
        for entry in self:
            if entry.molecule != curr and entry.molecule not in self._molecules:
                curr = entry.molecule
                index += 1
                self._molecules[entry.molecule] = index
        return
    
    def sort(self):
        self._data.sort(key = lambda entry: (self._molecules[entry.molecule], entry.start, -entry.end,
                                             entry.source, entry.feature, entry.attributes))
    
    def add_entry(self, gff_entry):
        self._data.append(gff_entry)
    
    def update_attr_fields(self, attr_mod = None):
        '''
        Update attribute fields w/ user-provided attribute field modification dictionary.
        (Otherwise, use self._attr_fields)
        '''
        if attr_mod is None: attr_mod = self._attr_mod
        for feature, mods in attr_mod.items():
            if feature in self._attr_fields:
                for canon, atypical in mods.items():
                    self._attr_fields[feature][canon] = atypical
            else:
                self._attr_fields[feature] = mods
        self._attr_fields_inv = self.invert_attr_fields()
        return
    
    def invert_attr_fields(self):
        output = {feature: {v: k for k, v in attr_map.items()}
                 for feature, attr_map in self._attr_fields.items()}
        return output
    
    def get_subfeatures(self, feature_ids, *features, index = False):
        '''
        Get subfeatures of features w/ user-provided feature_ids.
        '''
        features = set(features)
        if type(feature_ids) is str:
            feature_ids = {feature_ids}
        else:
            feature_ids = set(feature_ids)
        indices = [i for i, entry in enumerate(self) if
                   ((not features or entry.feature in features) and
                    entry.has_attr("Parent", feature_ids))]
        if index: return indices
        else: return self.get_i(indices)
    
    def get_subfeatures_full(self, feature_ids, *features, index = False):
        '''
        Get all features that are subfeatures of user-provided feature_ids
        AND subfeatures of those subfeatures, until there are no sub-sub...sub-features left.
        '''
        printi = make_local_print(quiet = self._quiet)
        features = set(features)
        output_indices = []
        curr_indices = []
        iteration_n = 0
        while True:
            iteration_n += 1
            printi(f"Executing iteration {iteration_n}")
            curr_indices = self.get_subfeatures(feature_ids, index = True)
            feature_ids = set(itertools.chain(*[self.get_i(i).get_attr("ID") for i in curr_indices]))
            output_indices.extend(curr_indices)
            if not feature_ids: break
        if features:
            output_indices = [i for i in output_indices if self.get_i(i).feature in features]
        ## prepare output
        output_indices = sorted(output_indices)
        if index: return output_indices
        else: return self.get_i(sorted(output_indices))
    
    def get_features_and_subfeatures(self, feature_ids, index = False, full = True):
        '''
        Gets features w/ feature_ids AND subfeatures of those features.
        If full = True, executes get_subfeatures_full for subfeature discovery, else get_subfeatures
        '''
        features = self.get_id(feature_ids, index = True, output_list = True)
        if full: subfeatures = self.get_subfeatures_full(feature_ids, index = True)
        else: subfeatures = self.get_subfeatures(feature_ids, index = True)
        ## prepare output
        final_features = sorted(features + subfeatures)
        if index: return final_features
        else: return self.get_i(final_features)

    def get_i_raw(self, indices, strip_newline = True):
        indices = indices if not isinstance(indices, int) else [indices]
        return self._indexed_file.get_line(*indices, strip_newline = strip_newline)
    
    def get_i(self, indices):
        '''
        Get GFF entry/entries by line index
        '''
        if self._data:
            if type(indices) is int: return self._data[indices]
            else: return [self._data[i] for i in indices]
        elif self.is_indexed():
            raw_entries = self.get_i_raw(indices, strip_newline = True)
            ## newline already stripped by get_i_raw. No need to strip them again.
            mk_annotation = self.make_annotation_from_str_gen(strip_newline = False)
            if type(indices) is int: return mk_annotation(raw_entries)
            else: return [mk_annotation(entry) for entry in raw_entries]
    
    def get_id(self, feature_ids, index = False, output_list = False):
        ## if even if output_list is only used when len(indices) == 1 or == 0 AND type(feature_ids) is str.
        ##  Always returns list otherwise.
        '''
        Get GFF entry by feature ID
        '''
        if type(feature_ids) is str:
            indices = [i for i, entry in enumerate(self) if entry.has_attr("ID", [feature_ids])]
        else:
            indices = [i for i, entry in enumerate(self) if entry.has_attr("ID", feature_ids)]
        if type(feature_ids) and not output_list and len(indices) <= 1:
            if not indices: return None
            elif index: return indices[0]
            else: return self.get_i(indices[0])
        elif index: return indices
        else: return list(self.get_i(indices))
    
    def write(self, fout, entries, **kwargs):
        '''
        Writes entries to file
        '''
        with open(fout, "w+") as f:
            f.write('\n'.join([entry.generate_str(**kwargs) for entry in entries]) + '\n')
    
    def write_i(self, fout, indices, **kwargs):
        '''
        Executes get_i, then writes output to file
        '''
        self.write(fout, self.get_i(indices), **kwargs)
        return
    
    def write_id(self, fout, feature_ids, **kwargs):
        '''
        Executes get_id, then writes output to file
        '''
        self.write(fout, self.get_id(feature_ids), **kwargs)
        
class Annotation:
    def __init__(self, entry, gff, **kwargs):
        self._gff = gff
        self.molecule = entry[0]
        self.source = entry[1]
        self.feature = entry[2]
        self.start = int(entry[3])
        self.end = int(entry[4])
        self.score = entry[5] if entry[5] == '.' else float(entry[5])
        self.strand = entry[6]
        self.phase = entry[7] if entry[7] == '.' else int(entry[7])
        self.attributes = Attributes(entry[8], gff, self, **kwargs)
        self._fields = {"molecule": self.molecule,
                        "source": self.source,
                        "feature": self.feature,
                        "start": self.start,
                        "end": self.end,
                        "score": self.score,
                        "strand": self.strand,
                        "phase": self.phase,
                        "attributes": self.attributes}
    def generate_attr(self, original = True, fields = None):
        if original: return self.attributes._raw
        else: return self.attributes.standardise_fields()
    def generate_str(self, fmt = "GFF"):
        if fmt.upper() in {"GFF", "GFF3"}:
            output = self.generate_gff()
        elif fmt.upper() in {"BED"}:
            output = self.generate_bed()
        return '\t'.join(map(str, output))
    def generate_gff(self):
        return list(map(str, [self.molecule, self.source, self.feature, self.start, self.end,
                              self.score, self.strand, self.phase, self.generate_attr()]))
    def generate_bed(self):
        return list(map(str, [self.molecule, self.start - 1, self.end, self.attributes.get("ID", fmt = str),
                              self.score, self.strand, self.source, self.feature, self.phase,
                              self.generate_attr()]))
    def get(self, *fields):
        return [self.f_dict[field] for field in fields]
    
    ## wrappers for Attribute methods
    def get_attr(self, a, **kwargs): return self.attributes.get(a, **kwargs)
    def is_attr(self, a, val, **kwargs): return self.attributes.is_attr(a, val, **kwargs)
    def has_attr(self, a, vals, **kwargs): return self.attributes.has_attr(a, vals, **kwargs)

class Attributes:
    def __init__(self, val, gff, entry, field_sep_inter = ';', field_sep_intra = ',',
                 **for_dummy_gff):
        self._entry = entry
        self._gff = gff if gff is not None else GFF(**for_dummy_gff)
        self._raw = val
        self._sep_inter = field_sep_inter
        self._sep_intra = field_sep_intra
        self._data = self._parse()
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
    def standardise_fields(self):
        def get_mod(field):
            feature_mod = get_recursively(self._gff._attr_fields_inv, field, self.feature, field)
            if feature_mod != field: return feature_mod
            else: return get_recursively(self._gff._attr_fields_inv, field, "all", field)
        return ';'.join([f"{get_mod(k)}={'.'.join(v)}" for k, v in self._data.items()])


#############################
##  EXTRACT_GFF_FEAUTURES  ##
#############################

def extract_features_and_subfeatures(gff_fname, feature_ids, fout, quiet = False,
                                     fin_fmt = "GFF", fout_fmt = "GFF", memsave = True):
    printi = make_local_print(quiet = quiet)
    printi("Reading files")
    gff = GFF(fname = gff_fname, fmt = fin_fmt, quiet = quiet, memsave = memsave)
    # feature_ids = splitlines(feature_id_fname)
    printi("Retrieving features")
    output_features = gff.get_features_and_subfeatures(feature_ids, index = True, full = True)
    printi("Writing features")
    gff.write_i(fout, output_features, fmt = fout_fmt)
    return


#################
##  ANN MANIP  ##
#################

def reduce_ann(gff_beds, ids, fout = None, mk_tmpf_name = None, fout_fmt = "BED", memsave = True):
    if mk_tmpf_name is None:
        import tempfile
        mk_tmpf_name = lambda x: tempfile.mkstemp()[1]
    gff_beds = assign_alias(gff_beds)
    bed_reds = {}
    for alias, gff_bed in gff_beds.items():
        bed_red = mk_tmpf_name(alias)
        bed_reds[alias] = bed_red
        fmt = "BED" if gff_bed.split('.')[-1].upper() == "BED" else "GFF"
        extract_features_and_subfeatures(gff_bed, ids, bed_red, quiet = True, memsave = memsave,
                                         fin_fmt = fmt, fout_fmt = fout_fmt)
    if fout:
        cat_files(tuple(bed_reds.values()), fout, remove = True)
    else:
        return bed_reds

