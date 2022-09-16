"""
Implementation of the weight-based set cover algorithm presented in Ajami, Z. and Cohen, S. (2019) Enumerating Minimal Weight Set Covers. Proceedings - International Conference on Data Engineering, 518-529
"""

## based on Ajami and Cohen (2019) Enumerating Minimal Weight Set Covers. Proceedings - International Conference on Data Engineering, 518-529

import copy

def make_examples():
    ## make example
    a, b, c, d, e = 'a', 'b', 'c', 'd', 'e'
    global s1, s2, s3, s4, S, U
    s1 = Set('s1', 2, [a, b, c])
    s2 = Set('s2', 3, [a, b, d])
    s3 = Set('s3', 2, [c, d])
    s4 = Set('s4', 5, [d, e])
    S = SetOfSets(s1, s2, s3, s4)
    U = S.elements
    return

def make_examples_2():
    a, b, c, d, e, f, g, h, i, j = range(10)
    global s1, s2, s3, s4, S, U
    s1 = Set('s1', 9, [a, b, c, d, e, f, g, h, i])
    s2 = Set('s2', 7, [a, b, c, d, e, f, g])
    s3 = Set('s3', 3, [h, i, j])
    s4 = Set('s4', 6, [e, f, g, h, i, j])
    S = SetOfSets(s1, s2, s3, s4)
    U = S.elements
    return

class e:
    ## a, b, c are all SetOfSets objects
    def __init__(self, a, b, c, tiebreaker = None):
        self.a = a
        self.b = b
        self.c = c
        self.tiebreaker = tiebreaker
    def __repr__(self):
        list_names = lambda setofsets:sorted(s.name for s in setofsets)
        return (f"{self.__class__.__name__}({list_names(self.a)},"
                f" {list_names(self.b)},"
                f" {list_names(self.c)})")
    @property
    def C(self): return self.a
    @property
    def Sfc(self): return self.b
    @property
    def So(self): return self.c
    @property
    def w(self): return self.C.w
    def copy(self):
        return self.__class__(self.a.copy(),
                              self.b.copy(),
                              self.c.copy(),
                              tiebreaker = self.tiebreaker)
    def _make_sorting_key(self):
        if self.tiebreaker:
            return lambda s: (s.w, self.tiebreaker(s))
        return lambda s: s.w
    ## dedicated functions so this class can be adapted to work when So is a SetOfSets of CollapsedNamedSets
    ## (or something like it that tracks equivalent sets and only removes one set of the CollapsedNamedSets
    ## but not the entire CollapsedNamedSets)
    def So_min(self):
        sorting_key = self._make_sorting_key()
        return min(self.So - self.C, key = sorting_key)
    def remove(self, *args, **kwargs):
        super().remove(*args, **kwargs)

class e1(e):
    def __init__(self, a, b, c):
        super().__init__(a, b, c)
    @property
    def Sf(self): return self.Sfc

class e2(e):
    def __init__(self, a, b, c):
        super().__init__(a, b, c)
    @property
    def Sc(self): return self.Sfc

class Q(list):
    def __init__(self, set_of_sets_class = None):
        super().__init__()
        self.set_of_sets_class = SetOfSets if set_of_sets_class is None else set_of_sets_class
    @property
    def top(self) -> e: return self.set_of_sets_class() if not self else self[0]
    def is_empty(self) -> bool: return len(self) == 0
    def extract_min(self):
        """
        Returns :class:`~minorg.minweight_sc.Set` object with lowest weight.
        
        Returns:
            :class:`~minorg.minweight_sc.Set`
        """
        return min(self, key = lambda s:s.w)

## set class that tracks set name and set weight
class Set(set):
    def __init__(self, name, weight, elements):
        """
        Effectively a set object that tracks name and weight.
        
        Arguments:
            name (str): set name
            weight (float): set weight
            elements (iterable): elements in set
        """
        self.name = name
        self.weight = weight
        super().__init__(elements)
    def __repr__(self):
        return f"Set(name='{self.name}', {set(self)})"
        # return f"Set(name = {self.name}, weight = {self.weight}, length = {self.length})"
    ## hashable by name
    def __hash__(self):
        return hash(self.name)
    @property
    def w(self) -> float: return self.weight
    def copy(self):
        """
        Returns a shallow copy of self.
        """
        return self.__class__(self.name,
                              self.weight,
                              list(self))

class SetOfSets(set):
    def __init__(self, *Sets):
        """
        Arguments:
            Sets (iter): iter of Set objects
        """
        # self._sets = {s.name: s for s in Sets}
        super().__init__(Sets)
        # super().__init__(set().union(*Sets))
    def __sub__(self, other):
        output = super().__sub__(other)
        return self.__class__(*output)
    def __repr__(self):
        return (f"SetOfSets(weight={self.w}, num_sets={len(self)},"
                f" redundancy={self.redundancy}, num_elements={len(self.elements)})")
    def __eq__(self, other):
        return all(s in other for s in self) and all(s in self for s in other)
    def __hash__(self):
        return hash(self.sets)
    @property
    def sets(self): return tuple(self)
    @property
    def elements(self): return set().union(*self)
    @property
    def set_names(self): return [s.name for s in self.sets]
    @property
    def w(self): return float("Inf") if not self else sum(s.w for s in self.sets)
    @property
    def weight(self): return self.w
    @property
    def redundancy(self) -> int:
        """
        Returns sum of the number of times each element is present in more than one set.
        """
        return sum(len(s) for s in self.sets) - len(self.elements)
    def get(self, *names) -> list:
        """
        Get :class:`~minorg.minweight_sc.Set`s by set name.
        
        Returns:
            list of :class:`~minorg.minweight_sc.Set`
        """
        names = set(names)
        return [s for s in self if s.name in names]
    def union(self, *args, **kwargs):
        """
        Implements set.union method for this class.
        
        Returns:
            :class:`~minorg.minweight_sc.SetOfSets`
        """
        output = super().union(*args, **kwargs)
        return self.__class__(*output)
    def intersection(self, *args, **kwargs):
        """
        Implements set.intersection method for this class.
        
        Returns:
            :class:`~minorg.minweight_sc.SetOfSets`
        """
        output = super().intersection(*args, **kwargs)
        return self.__class__(*output)
    def symmetric_difference(self, *args, **kwargs):
        """
        Impelements set.symmetric_difference method for this class.
        
        Returns:
            :class:`~minorg.minweight_sc.Set`
        """
        output = super().symmetric_difference(*args, **kwargs)
        return self.__class__(*output)
    def remove(self, *Sets) -> None:
        """
        Remove one or more :class:`~minorg.minweight_sc.Set`.
        """
        self -= set(Sets)
        return
    def add(self, *Sets) -> None:
        """
        Add one or more :class:`~minorg.minweight_sc.Set`.
        """
        self.update(Sets)
        return
    def copy(self):
        """
        Create a shallow copy of self.
        
        Returns:
            :class:`~minorg.minweight_sc.SetOfSets`
        """
        return self.__class__(*self.sets)

def wc_ratio(s1, s2, low_coverage_penalty = 0) -> float:
    """
    Calculates ratio of weight of s1 to number of elements in s2 that can be covered by s1.
    (wc is shorthand for weight-cover, a name I came up with because the authors did not)
    
    Low coverage penalty was added to allow disincentivisation of ultra large sets full of small,
    non-overlapping, ultra low coverage sets.
    
    Arguments:
        s1 (set): first set
        s2 (set): second set
        low_coverage_penalty (float): if 0, no penalty for uncovered elements.
            Modifies output value into
            '<wc ratio> + (<low_coverage_penalty> * <wc ratio> * <uncovered>/<covered>)'
    """
    cover = s2.intersection(s1)
    if not cover: return float("Inf")
    wc = s1.w/len(cover)
    cov_penalty = len(s2-s1)*wc*low_coverage_penalty
    return wc + cov_penalty

def approx_min_SC(U, S, seed = None, low_coverage_penalty = 0) -> SetOfSets:
    """
    Approximate minimum set cover algorithm.
    
    Arguments:
        U (set): set of elements (targets) to cover
        S (:class:`~minorg.minweight_sc.SetOfSets`): SetOfSets (or child class) object 
            containing sets (gRNA coverage) for set cover
        seed (:class:`~minorg.minweight_sc.SetOfSets`): Set (or child class) object to add first
    
    Returns:
        :class:`~minorg.minweight_sc.SetOfSets`: Set cover solution. May also be child class of SetOfSets.
    """
    set_of_sets = type(S)
    U = set(U)
    if set(U) - S.elements:
        return False
    if seed:
        S = S.copy()
        C = set_of_sets(seed)
        U -= set(seed)
        S.remove(seed)
    else:
        C = set_of_sets()
    while U != set():
        s = min(S.sets,
                key = lambda s:(wc_ratio(s, U, low_coverage_penalty = low_coverage_penalty), -len(s)))
        C.add(s)
        U -= s
    return C

def argmin(set_of_sets) -> Set:
    """
    Returns Set with smallest weight.
    
    Arguments:
        set_of_sets (:class:`~minorg.minweight_sc.SetOfSets`): SetOfSets (or child class) object
    
    Returns:
        :class:`~minorg.minweight_sc.Set`
    """
    return min(set_of_sets, key = lambda s:s.w)

# output = []
def enum_approx_order_SC(U, S, num_iter = 100, seed = None, low_coverage_penalty = 0) -> list:
    """
    Minimum weight set cover algorithm.
    
    Algorithm described in: Ajami, Z. and Cohen, S. (2019). Enumerating Minimal Weight Set Covers. In Proceedings - International Conference on Data Engineering, 518-529
    
    Arguments:
        U (set): set of elements (targets) to cover
        S (:class:`~minorg.minweight_sc.SetOfSets`): SetOfSets (or child class) object 
            containing sets (gRNA coverage) for set cover
        num_iter (int): maximum number of iterations (default=100)
    
    Returns:
        list of :class:`~minorg.minweight_sc.SetOfSets`: list of SetOfSets object where each 
        is a set cover solution
    """
    ## get class of objects that are sets of sets (minimally must be SetOfSets, can be child classes)
    set_of_sets = type(S)
    output = []
    Q1, Q2 = Q(set_of_sets_class = set_of_sets), Q(set_of_sets_class = set_of_sets)
    C = approx_min_SC(U, S, seed = seed, low_coverage_penalty = low_coverage_penalty)
    if C:
        Q1.append(e1(C, set_of_sets(), S))
    ## limit number of iterations
    for i in range(num_iter):
        ## exit if queues are empty
        if (Q1.is_empty() and Q2.is_empty()): break
        ## if first entry of Q1 has lower (or equal) C weight to Q2's
        if Q1.top.w <= Q2.top.w:
            ## get entry with lowest C weight
            top_e = Q1.extract_min()
            C, Sf, So = top_e.C, top_e.Sf, top_e.So
            output.append(C)
            ## some algorithm things that I don't really understand
            if (So - C) != set():
                s = argmin(So - C)
                new_C = C.union({s})
                new_Sf = set_of_sets(s)
                new_So = set_of_sets(*(So - C - new_Sf))
                Q2.append(e2(new_C, new_Sf, new_So))
            for s in C - Sf:
                So = So - set_of_sets(s)
                U_prime = U - Sf.elements
                C_prime = approx_min_SC(U_prime, So, low_coverage_penalty = low_coverage_penalty)
                if C_prime is not False:
                    Q1.append(e1(C_prime.union(Sf), Sf, So))
                Sf = Sf.union({s})
            null = Q1.pop(0)
        ## else if first entry of Q2 has lower weight C
        else:
            ## get entry with lowest C weight
            top_e = Q2.extract_min()
            C, Sc, So = top_e.C, top_e.Sc, top_e.So
            output.append(C)
            if So != set():
                s_prime = argmin(So)
                new_C = C.union({s_prime})
                new_Sc = set_of_sets(s_prime)
                Q2.append(e2(new_C - Sc, new_Sc.copy(), So - new_Sc))
                Q2.append(e2(new_C, new_Sc.copy(), So - new_Sc))
            null = Q2.pop(0)
    return output

# make_examples_2()
# x = enum_approx_order_SC(U, S)
# [(setofsets.w, list(s.name for s in setofsets)) for setofsets in x[:10]]
# ## now we just gotta find a way to tie break this

# def n
