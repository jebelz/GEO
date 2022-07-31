"""This module does the heavy lifting of Einstein Summation Notation.
(It DOES NOT address any deep geometry in General Relativity: it is
just a typo-graphic analyzer).

It converts many forms of request into the correct method calls (or more).
Specific dotted requests:

T.iij  contracts the 1st 2 inidices

T.jik  transposes the 1st 2 indices

T.xzk fixes the 1st 2 indices at i=0, j=2, repsectively.

The general dotted request looks like, for example:

T.jmnxyikk

which tranposes j and i, contracts the k's, and then fixes x and y.

Finally, attribute requests (via getattr) that look like:

T.{ij}k

T.[ij]k

Form the tensor that is (anti)symmetric in the 1st 2 indices.

The work is done in the Albert class, with various helper functions.
"""
## \namespace geo.metric.euclid.einstein.albert Apply Einstein summation to
# Tensors.
import collections
import itertools
import operator

from ..euclid import ZOO
from . import AXES

## This is the only public interface.
__all__ = ("Albert",)


## These are fixed indices
FIX = AXES


## These are FREE indices, and the abstraction is leaky: you cannot use this
# modules on tensors with a rank larger than the length of this string (ikr:
# it's pretty long).
FREE = 'ijklmnopqrstuv'

## Symmetric Symbols
SYMBOLS = '{}'


## Asymetric Symbols
ASYMBOLS = '[]'


## Sign convention of symmetrizing symbols
SIGN = {ASYMBOLS[0]: operator.neg, SYMBOLS[0]: operator.pos}


## Open: Close pairs of symbols-> [] and {}.
PAIR = {key: chr(ord(key)+2) for key in SIGN}


# internal enumeration/flag
IS_FIXED = 1


# internal enumeration/flag
IS_CONTRACTING = 2


# internal enumeration/flag
IS_FREE = 0


## Classify Request by Index (very procedural, as there is a magic code).
# \param attr The attribute request
# \returns d A list of ::IS_FIXED, ::IS_CONTRACTING, ::IS_FREE
def step1a(attr):
    """Classify requests as IS_FIXES, IS_CONTRACTING, IS_FREE
    for each character in attr."""
    d = []
    for c in str(attr):
        if c in FIX:
            magic = IS_FIXED
        elif attr.count(c) == 2:
            magic = IS_CONTRACTING
        elif attr.count(c) == 1:
            magic = IS_FREE
        else:
            raise AttributeError(
                "can't deal with {} {}'s in {}".format(
                    attr.count(c), c, attr))
        d.append(magic)
    return d


## Convert non-free indices to free, and set free to None:
#  Since we MUST transpose before contraction/fixing, we have to
#  convert those to free indices, while holding free indices
# \param attr str of attribute requests
# \param dum The output of step1a()
# \returns new A list of new indices
def step1b(attr, dum):
    """new = step1b(attr, dum) """
    # 1st create a list to fill
    new = [None] * len(dum)
    # Address list by index
    for count, (dummy, d) in enumerate(zip(attr, dum)):
        if d != IS_FREE:
            # overwrite non-free index with ordered FREE index name
            new[count] = FREE[count]
    return new


## Rename free indices: since we used the free index name (maybe) in the
# renaming of non-free indices, we MUST chose the next available free
# index character
def step1c(attr, dum, new):
    """deque = step1c(attr, dum, new)"""
    # only so many letters between i and w.
    assert len(attr) <= len(FREE), "The Abstraction has Leaked!"
    # get available free indices
    free = collections.deque(i for i in FREE if i not in new)
    # get original free requests
    place = list(itertools.compress(attr,
                                    (item == IS_FREE for item in dum)))

    # keep a copy of them
    place_holder = place[:]
    # create a new target list for replacements
    replace = [None] * len(place)

    while place:
        # get minimum free reuest
        i = place.pop(place.index(min(place)))
        # replace it with min available
        replace[place_holder.index(i)] = free.popleft()
    return collections.deque(replace)


## Insert replacements into dummy requests-- now we have a complete
# set of indices to work with
def step1d(attr, dummy, new, rep):
    """str = step1c(attr, dum, new, deque)"""
    while rep:
        new[new.index(None)] = rep.popleft()
    return new


## Convert original attr request into a temporary request for transposition
# \param attr
# \returns tranposition requests.
def step1(attr):
    """Step 1: convert general dotted request into all free indices,
    with transposeition."""
    dum = step1a(attr)
    new = step1b(attr, dum)
    rep = step1c(attr, dum, new)
    return "".join(step1d(attr, dum, new, rep))


## Mixin for [] and {} tensor subscripts.
class Split(object):
    """This mixin handles non-dotted request, that is, anthing with
    [] or {}.

    It creates new sub requests and then adds of subtracts them (and
    normalizes them).
    """

    ## int --> rank
    def __int__(self):
        return self.tensor.rank

    ## Albert.attr
    def __str__(self):
        return self.attr

    ## Pass to str(self).
    def __contains__(self, x):
        return x in str(self)

    ## Am I splitting?
    # \returns bool If there is a split request ([] or {}).
    def issplit(self):
        """returns bool iff attribute splits tensor."""
        for left, right in PAIR.iteritems():
            if left in self:  # parse
                if right not in self:  # raise
                    msg = left + " is not balanced by a " + right
                    class GymnasticsError(UserWarning, ValueError):
                        pass
                    raise GymnasticsError(msg)
                return True
        return False
        
    ## Split on a key
    # \param key "[" or "}"
    # \returns tuple or (indent, transposed) indices
    # \throws AssertionError on key request.
    def split_pair(self, key):
        """ijk, ikj = i[jk].split('['), for example."""
        assert key in PAIR, "should not call split unless assertion."
        attr = str(self)
        start = attr.index(key)
        stop = attr.index(PAIR[key])
        swap = attr[start+1: stop]
        beg = attr[:start]
        end = attr[stop+1:]
        left = beg + swap + end
        right = beg + "".join(reversed(swap)) + end
        return left, right

    ## synonym, but why?
    split = split_pair
    
    ## Evaluate Bracket Request
    # \returns tuple (left, sign function, right)
    # \throws RuntimeError If bracket is not found (should never be reached)
    def bracket(self):
        """(left, +/-, right) indices (w/ sign of bracket)."""
        for key, op in SIGN.iteritems():
            if key in str(self):
                left, right = self.split(key)
                return left, op, right
        raise RuntimeError("failed to find sign. (How did you get here?)")

    ## Process any bracketed request
    # \return tensor The final version, after bifurcation...
    def process_brackets(self):
        """returns tensor w/ bracket request processed."""
        left, op, right = self.bracket()
        # this does T.left +/- T.right, e.g: [ij]--> T.ij - T.ji
        return (getattr(self.tensor, left) +
                op(getattr(self.tensor, right))) / 2.
    
    def process_brackets(self):
        if ASYMBOLS[0] in str(self):
            return self.dev(ASYMBOLS[0])
        return self.dev(SYMBOLS[0])

    ## This will do any symmetrization....
    def dev(self, key):
        from geo.metric.schur_weyl.monte import Sym
        assert key in PAIR, "code should not call split unless assertion."
        attr = str(self)
        start = attr.index(key)
        stop = attr.index(PAIR[key])
        # Get the size of the group
        indices = attr[start+1: stop]
        n = len(indices)

        beg = attr[:start]
        end = attr[stop+1:]

        attr_holder = beg + SYMBOLS + end

        # get the group of permuations
        Sn = Sym(n)

        # Compute target rank as input rank minus fixed indices
        target_rank = self.tensor.rank - sum(map(str(self).count, AXES))
        # intial null tensor of correct rank
        result = ZOO[target_rank](*list([0]*pow(3, target_rank)))

        # default to positive summation
        sign = 1

        # loop over all permutations
        for sigma in Sn:
            attr_sig = attr_holder.format("".join(sigma(indices)))
            if key == '[':
                # apply antsymmetric parity If needed
                sign = sigma.parity()
            # accumlate portion of tensor with new attr request
            result += sign * getattr(self.tensor, attr_sig)
        # normalize with total number of permutations
        result /= float(len(Sn))
        return result
    ## This will do any symmetrization....
    def dev(self, key):
        from geo.metric.schur_weyl.monte import Sym
        assert key in PAIR, "code should not call split unless assertion."
        attr = str(self)
        start = attr.index(key)
        stop = attr.index(PAIR[key])
        # Get the size of the group
        indices = attr[start+1: stop]
        n = len(indices)

        beg = attr[:start]
        end = attr[stop+1:]

        attr_holder = beg + SYMBOLS + end

        # get the group of permuations
        Sn = Sym(n)

        # Compute target rank as input rank minus fixed indices
        target_rank = self.tensor.rank - sum(map(str(self).count, AXES))

        ## (1), 2 for (anti)aymmetric
        parity_exp = 2 - (key == '[')

        ## Sum all permutations of indicies (with parity or not) starting
        # from a null tensor of target_rank, nomralized by n!
        return reduce(operator.add,
                      (pow(sigma.parity(), parity_exp) *
                      getattr(self.tensor,
                              attr_holder.format("".join(sigma(indices))))
                      for sigma in Sn),
                      ZOO[target_rank](
                          *list([0]*pow(3, target_rank)))) / len(Sn)
        

## Einstein Summation Notation Processor
class Albert(Split):
    """T = Albert(tensor, attr)()

    The actual instantiation is done by a package metaclass. It verifies
    that attr makes sense for tensor.

    Instaniations gets the tensor and the attribute request into 1 object.
    Calling it then does the work, returning the desired tensor result.

    A mess of a class. Ideally, we're going to self mutate as we
    move along from Transpose, Contraction, Fixation.

    Of course, the difficulty arises in processing multiple "features"
    See __iter__ for the solution.
    """

    ## Construct from a Tensor and an attribute request
    def __init__(self, tensor, attr):
        ## The Tensor (any rank)
        self.tensor = tensor
        ## The request.
        self.attr = attr

    ## Iterate through the steps and return the final answer.
    # \returns T The fully processed tensor
    def __call__(self):
        # each iteration progresses down the interpretation chain
        for item in self:
            pass
        # At then end, is the answer.
        return item

    ## Iterate each step
    # \yields process_brackets() Do any (anti)symmetric request and the leave.
    # \yields step1() That is, transpose all indices 1st.
    # \yields contract() Then do all contractions
    # \yields fixies() Then evaluate any fixed indices.
    def __iter__(self):
        if self.issplit():  # fork in processing tree
            # should spawn 2 new Alberts (which are completed and yielded)
            yield self.process_brackets()
        else:
            # step 1 is transpose
            T = transpose(self.tensor, self.attr)
            yield T, step1(self.attr)

            # step 2 is contract
            T, rat = contract(T, self.attr)
            yield T, rat

            # step 3 is evaluate fixies
            T = fixies(T, rat)
            yield T


## Transpose a Tensor according to alpahbetical indices
# \param T a rank N tensor, e.g: \f$ T_{ijklm} \f$
# \param attr length n-string of tensor indices, e.g.: \f$ ilkjm \f$
# \returns \f$ T_{ilkjm} \f$, for example.
# \throws AssertionError Iff rank(T) != len(attr)
def transpose(T, attr):
    """Process transposition requests."""
    assert T.rank == len(attr)
    # compute like-transposing attrs, with fixies and contractions free.
    attr = step1(attr)
    # If there is no transpose request, don't do it.
    if attr == FREE[:len(attr)]:
        return T
    # pass work off to tensor's transpose method, with numeric args pointing
    # to tensor indices.

    #TODO: RESOLVE THIS PROBLEM! WHICH WAY TO TRANSPOSES GO?????--Which one
    # Makes szyggy's work??

    indices = map(FREE.index, attr)

    ####return T.ntranspose(indices)  # this seems to make sense but...
    # this is the one that works and satisfies "test()".
    return T.transpose(*indices)

    

## Process fixed indices in a tensor
# \param T a rank N tensor, e.g: \f$ T_{ijklm} \f$
# \param attr length n-string of tensor indices with fixed
# indices, e.g.: \f$ ixyjk \f$
# \returns \f$ T_{i01jk} \f$, for example.
# \throws AssertionError Iff rank(T) != len(attr)
def fixies(T, attr):
    """Process fixed index requests."""
    assert T.rank == len(attr)
    # only apply If needed (x y or z is in attr).
    if any(c in FIX for c in attr):
        # Use Tensor_.e method to process fixies.
        return T.e(*[FIX.index(item) if item in FIX else
                     Ellipsis for item in attr])
    # or do nothing
    return T

## Process Index Contractions in a Tensor
# \param T a rank N tensor, e.g: \f$ T_{ijklm} \f$
# \param attr length n-string of tensor indices with
# contractions, e.g.: \f$ ijijk \f$
# \returns \f$ T_{ijik} \f$, for example.
# \throws AssertionError Iff rank(T) != len(attr)
def contract(T, attr):
    """Process contractions."""
    from collections import deque
    assert T.rank == len(attr)
    # get attr as a list
    lat = list(iter(attr))
    # now get a deque of paired indices
    pairs = deque(set([c for c in lat if c in FREE and lat.count(c) > 1]))

    # keep processing pairs
    while pairs:
        # pop left summation character
        cc = pairs.popleft()
        # get the indices's index
        i1 = lat.index(cc)
        # remove it from the list of attributes
        lat.pop(i1)
        # now find the second occurance of the summation index
        i2 = lat.index(cc)
        # and remove it from the list of attributes
        lat.pop(i2)
        # contract on the 2 indices (+1 correction is b/c 1st occurence was
        # removed from the list
        T = T.contract(i1, i2+1)

    # return index request, with pair-wise contractions removed.
    return T, "".join(lat)






## THIS MUST WORK, and einstein.py must make the last coulmn agree
def test():
    from geo.metric.euclid.syzygy import DELTA, DD
    from itertools import product

    f = lambda x: int(round(x))
    
    A = 'xyz'
    
    for (i, j, k, m), dikjm, dimjk in zip(product(A,A,A,A),DD.ikjm.iter(),DD.imjk.iter()):
        ik, jm, im, jk = i+k, j+m, i+m, j+k
        p1 = getattr(DELTA, ik)
        p2 = getattr(DELTA, jm)
        q1 = getattr(DELTA, im)
        q2 = getattr(DELTA, jk)
        p = p1*p2
        q = q1*q2
        print  i,j,k,m,':', \
               f(p), f(getattr(DD, ik+jm)), f(dikjm), '::', \
               f(q), f(getattr(DD, im+jk)), f(dimjk)
