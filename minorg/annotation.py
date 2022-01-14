import itertools

from minorg.index import IndexedFile

from minorg.functions import (
    assign_alias, make_local_print
)

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
    
    def __init__(self, fname = None, data = [], attr_mod = {}, genetic_code = 1, fmt = None,
                 quiet = False, memsave = False, chunk_lines = 1000, **kwargs):
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
        self._genetic_code = genetic_code
        self._quiet = quiet
        self._kwargs = kwargs
        self._molecules = {}
        self._chunk_lines = chunk_lines
        self._indexed_file = None
        self._memsave = memsave
        self.update_attr_fields()
        if not self._memsave:
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
    
    def empty_copy(self, other = None):
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
    
    def is_bed(self):
        if not self._fmt: return False
        return self._fmt.upper() == "BED"
    def is_gff(self):
        if not self._fmt: return False
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
        if self._molecules:
            mol_key = lambda mol: self._molecules[mol]
        else:
            mol_key = lambda mol: mol
        self._data.sort(key = lambda entry: (mol_key(entry.molecule), entry.start, -entry.end,
                                             entry.source, entry.feature, entry.attributes.standardise_fields()))
    
    def add_entry(self, gff_entry, duplicate_check = False):
        if duplicate_check:
            for entry in self:
                if gff_entry == entry: return
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
    
    def get_features(self, *features, index = False):
        '''
        Get entries of specific feature types.
        '''
        features = set(features)
        if index:
            return [i for i, ann in enumerate(self) if ann.feature in features]
        else:
            return [ann for ann in self if ann.feature in features]
    
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
            # iteration_n += 1
            # printi(f"Executing iteration {iteration_n}")
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
        return self._indexed_file.get_line(*indices, strip_newline = strip_newline, output_fmt = list)
    
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
            output = [mk_annotation(entry) for entry in raw_entries]
            if type(indices) is int: return output[0]
            else: return output
    
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
    
    def subset(self, feature_ids = None, features = None, subfeatures = True, sort = True):
        '''
        Subset data by feature ID (feature_ids) and/or feature type (features)
          and generates new GFF object from them
        '''
        new_gff = self.empty_copy()
        ## start subsetting
        if feature_ids:
            if subfeatures:
                data = self.get_features_and_subfeatures(feature_ids, index = False)
            else:
                data = self.get_id(feature_ids, index = False, output_list = True)
            new_gff._data = data
            if sort:
                new_gff.sort()
            if features:
                return new_gff.subset(features = features)
            else:
                return new_gff
        else:
            if features:
                if isinstance(features, str):
                    features = [features]
                data = self.get_features(*features, index = False)
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
        try:
            self.molecule = entry[0]
            self.source = entry[1]
        except Exception as e:
            print("entry:", entry)
            print("entry[0]:", entry[0])
            raise e
        self.feature = entry[2]
        self.start = int(entry[3])
        self.start0 = self.start - 1 ## 0-indexed start
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
    def __eq__(self, other):
        return (self.generate_gff() == other.generate_gff())
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
    def __eq__(self, other):
        return (self.standardise_fields() == other.standardise_fields())
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
            feature_mod = get_recursively(self._gff._attr_fields_inv, field, self._entry.feature, field)
            if feature_mod != field: return feature_mod
            else: return get_recursively(self._gff._attr_fields_inv, field, "all", field)
        return ';'.join([f"{get_mod(k)}={'.'.join(v)}" for k, v in self._data.items()])


#############################
##  EXTRACT_GFF_FEAUTURES  ##
#############################

def extract_features_and_subfeatures(gff_fname, feature_ids, fout, quiet = False, attr_mod = {},
                                     fin_fmt = "GFF", fout_fmt = "GFF", memsave = True, sort = True):
    printi = make_local_print(quiet = quiet)
    printi("Reading files")
    gff = GFF(fname = gff_fname, fmt = fin_fmt, quiet = quiet, memsave = memsave, attr_mod = attr_mod)
    # feature_ids = splitlines(feature_id_fname)
    printi("Retrieving features")
    gff_subset = gff.subset(feature_ids, subfeatures = True, sort = sort)
    # output_features = gff.get_features_and_subfeatures(feature_ids, index = True, full = True)
    printi("Writing features")
    gff_subset.write(fout)
    return


#################
##  ANN MANIP  ##
#################

def subset_ann(gff_beds, ids, fout = None, mk_tmpf_name = None, fout_fmt = "BED", attr_mods = {},
               memsave = True, sort = True):
    if mk_tmpf_name is None:
        import tempfile
        mk_tmpf_name = lambda x: tempfile.mkstemp()[1]
    gff_beds = assign_alias(gff_beds)
    gff_subsets = {}
    for alias, gff_bed in gff_beds.items():
        gff_subset = mk_tmpf_name(alias)
        gff_subsets[alias] = gff_subset
        fmt = "BED" if gff_subset.split('.')[-1].upper() == "BED" else "GFF"
        extract_features_and_subfeatures(gff_bed, ids, gff_subset, quiet = True, memsave = memsave,
                                         attr_mod = attr_mods.get(alias), fin_fmt = fmt, fout_fmt = fout_fmt,
                                         sort = sort)
    if fout:
        cat_files(tuple(gff_subsets.values()), fout, remove = True)
    else:
        return gff_subsets

def reduce_ann(gff_beds, ids, fout = None, mk_tmpf_name = None, fout_fmt = "BED", attr_mods = {},
               memsave = True, sort = True):
    return subset_ann(gff_beds, ids, fout = fout, mk_tmpf_name = mk_tmpf_name,
                      fout_fmt = fout_fmt, memsave = memsave, attr_mods = attr_mods, sort = sort)
