# -*- coding: utf-8 -*-
"""Tensor Product Projecting Module for spherical tensors and
Quantum Mehcanical Angular Momentum.

Here is why we need it: Consider the rank-N tensor. It is the Cartesian
(outer) product of N basic vectors, e.g.:

xyzzz is rank 5.

There are 3**N combinations. Now we can replace (x, y, z) with a spherical
basis: (0, +, -). Big whoop. The point is, that we need to consider the
combinations of the basis vectors that are rotation eigenstates--in the
language of quantum mechanics that would be combo that are eigenstates
of TOTAL J and M:

|J, M>

So the prior:

xyzzz --> [|m1=1>+|m1=-1>][|m2=1>-|m2=-1>]|m3=0>|m4=0>|m5=0>/2

with j1 = j2 = j3 = j4 = j5 = 1 can be converted to a sum over J, M of:

a|J, M> + a'|J', M'> + ....

where each ket is a eigenstate of total spin and spin-projection (and also
a rank 5 tensor that transforms accordingly)

Those are special combinations of the N-adics, which are products of
states that are eigenstates of their individual rotation operators.

As we know, the transition from 2 states (dyads) to total rotation
eigenates is Clebsch Gordan coefficients. Well you just have to repeat
that to solve the N-adic, and this module helps.

So basis state dyadics look like this:

>>>k = Ket(j, m)   #  |j, m>

while a further combination looks like:

>>>for irrep in Ket(2, 0)*Ket(1, 1):
   ...:     print unicode(irrep)
√（1∕10) |1, 1⟩
-√（½) |2, 1⟩
√（⅖) |3, 1⟩


Like wise, the no brainer product is:

In [2]: Ket(1,1)**3
Out[2]: [Ket(3.0, 3)]

which says 3 aligned vectors is spin-3 with M=3.


It is by no means a full Racah Algebra-- but we're going with the Zen:
"practicality beats purity" It is a start.

The module does one step in process: it expands tensor product
representations into tensor sums of irreducible representations.

All possible combinations of a N-states of spin-j are:
>>>expand(N, j)


Finally:
natural_form_eigenstate(n, m, j) tells you how to construct the highest weight
(that is, symmetric and trace-free, aka natural form) state out of the
tensor prodoct of n spin-j basis states. Whhhattt? Here is an exmaple:

In [11]: for (m1,m2), ket in  racah.natural_form_eigenstates(2, 0).items():
...print "|{},{}>   -->".format(m1,m2), unicode(ket)
|-1,1>   --> √（⅔) |2, 0⟩
|0,0>   --> √（⅔) |2, 0⟩


combine that with:


In [44]: print (p*m) + (m*p)
[sqrt(4/3) | 0.0, 0> +
sqrt(2/3) | 2.0, 0>]

In [45]: print z*z
[sqrt(1/3) | 0.0, 0> +
sqrt(2/3) | 2.0, 0>]

and you see that |2, 0> = sqrt(3/2) S.zz (the trace free requirement
removes the J=0 state).


Also, the combinator adds any 2 quantum angular momenta (in this example,
2 spin-1's):

In [3]: print racah.Combinator(1, 1)
|0.0, -0.0> =  - sqrt(1/3)|0:0> + sqrt(1/3)|1:-1> + sqrt(1/3)|-1:1>
|1.0, 0.0> =  + sqrt(1/2)|1:-1> - sqrt(1/2)|-1:1>
|1.0, 1.0> =  - sqrt(1/2)|0:1> + sqrt(1/2)|1:0>
|1.0, -1.0> =  + sqrt(1/2)|0:-1> - sqrt(1/2)|-1:0>
|2.0, 0.0> =  + sqrt(2/3)|0:0> + sqrt(1/6)|1:-1> + sqrt(1/6)|-1:1>
|2.0, 1.0> =  + sqrt(1/2)|0:1> + sqrt(1/2)|1:0>
|2.0, 2.0> =  + sqrt(1)|1:1>
|2.0, -2.0> =  + sqrt(1)|-1:-1>
|2.0, -1.0> =  + sqrt(1/2)|0:-1> + sqrt(1/2)|-1:0>

where the left is |J, M> and the states on the right are |m1: m2> meaning
|j1, m1>|j2, m2>.


"""
## \namespace geo.metric.wigner.racah Representation Theory and Rotations
# for spherical tensors and Quantum Angular Momentum
import collections
import fractions
import functools
import itertools
import operator

from . import clebsch_gordan
from . import casimir

#__all__ = ('Ket', 'ITS')


HALF = fractions.Fraction(1, 2)

## Numerical Error which is considered Zero.
SMALL = 1.e-8


## List Unraveler (as a generator)
# \param x A list of lists of lists or of n
def iravel(x):
    """list = iravel(x)

    where x is unknown legnth and irregular depth list of list, or not.

    >>>iravel(4) --> 4

    >>iravel([4, [2, range(3)])) --> [4, 2, 0, 1, 2]
    """
    try:
        # is x iterable?
        len(x)
    except TypeError:
        # no it is not, yield it.
        yield x
    else:
        # yes it is, not iterate it
        for y in x:
            # and unravel what it yields.
            for z in ravel(y):
                # and yield it.
                yield z

def gravel(x):
    print 'gravel:', x
    from collections import Iterable            
    from itertools import imap         
    assert not isinstance(x, basestr)
    if isinstance(x, Iterable):
        for y in x:
            for z in gravel(y):
                yield z
    else:
        yield x

def ravel(x):
    """ravel"""
    return list(iravel(x))


## Whoops, it's possible that this is incoherent addition, as
# same l, m may apply to dIfferent slots (particles).
def coherent_add(*args):
    """Coherent Add."""
    kets = list(args)
    result = []
    while kets:
        # pop off a ket
        ket = kets.pop()
        # check to see If other jm-equivalent kets are in list
        while any([ket == k for k in kets]):
            # search for it
            for n in range(len(kets)):
                # find it
                if kets[n] == ket:
                    # pop it off and add it (coherently) to ket.
                    ket += kets.pop(n)
                    # start over looking for ket
                    break
                 # check next index...
        result.append(ket)
        continue
    # filter out alpha=0 items (it happens)
    return ITS(filter(lambda x: abs(x) > SMALL, result))


def coherent_add(*args):
    return reduce(iadd, args)


def iadd(x, y0):
    from copy import deepcopy
    result = deepcopy(list(x))
    y = deepcopy(list(y0))
    while y:
        d = y.pop()
        if d in result:
            idx = result.index(d)
            result[idx].alpha += d.alpha
        else:
            result.append(d)
    return ITS(filter(bool, result))


## Finish sort out the kets.
def _finish(kets):
    return ITS(sorted(iravel(kets)))


## Irrep Ket.
#@functools.total_ordering
class Ket(object):
    """|j, m> = Ket(j, m)"""

    ## Construct from degree, order, and phase-factor/weight
    # \param j degree
    # \param order
    # \param alpha =1 optional phase-factor/weight
    def __init__(self, j, m, alpha=1):
        ## Degree: total angular momentum quantum number
        self.j = j
        ## Order: Zed projection quantum number.
        self.m = m  # TODO: consider setter protection
        ## \f$ \alpha \f$ is an optional phase factor and/or scalar weight.
        self.alpha = alpha

    def __nonzero__(self):
        return abs(self.alpha) > SMALL
        
    def __pos__(self):
        return type(self)(self.j, self.m, alpha=+self.alpha)

    def __neg__(self):
        return type(self)(self.j, self.m, alpha=-self.alpha)

    ## \returns \f$ \delta_{jj'}\delta_{mm'}\alpha\alpha'^* \f$
    def __or__(self, other):
        return (self == other) * (self.alpha * other.alpha.conjugate())

    ## Conjugate???
    # \returns \f$ (-1)^j\alpha^*|j, -m\rangle \f$
    def __invert__(self):
        return (pow(1J, 2*self.j) * # TODO: CHECK J sign (and delete this comment)
                type(self)(self.j, -self.m, alpha=self.alpha.conjugate()))

    ## Positive integer power.
    # \param n for n in (1,2,...)
    # \returns \f$|j,m\rangle|j,m\rangle\cdots n{\mathrm -times}|j,m\rangle\f$
    def __pow__(self, n):
        return reduce(operator.mul, itertools.repeat(self, n))

    ## Convert Tensor Product to Tensor Sum
    # \param other \f$|j',m'\rangle\f$
    # \returns  \f$|j,m\rangle \otimes |j',m'\rangle =
    # \oplus_{J=|j-j'|, M=m+m'}^{j+j''}c_{jj'mm'JM}|J,M\rangle \f$
    def __mul__(self, other):
        if isinstance(other, type(self)):  # external poly
            return self.expand(other)
        try:
            return _finish(map(self.__mul__, other))
        except TypeError:
            return self.__rmul__(other)

    ## Generate valid J values
    # \param j \f$ j'\f$ for integer of 1/2 integer
    # \return list \f$ |j-j'| \le j le j+j' \f$
    def jrange(self, j):
        """[|j-j'|, |j-j'|+1, ..., j+j'-1, j+j']"""
        return list(clebsch_gordan.jrange(self.j, j))

    ## Generate Kets in Tensor Sum
    # \param other \f$|j'm'\rangle\f$
    # \kwd limit \f$=10^{-10}\f$ states below this abs are considered 0
    # \yields \f$ c_{jj'mm'JM}|J,M>\rangle
    def iexpand(self, other, limit=SMALL):
        """jm.iexpand(j'm', limit=SMALL)

        yields a*a'*c_(jj'm'JK)|JK>
        """
        #pylint: disable=E1102
        m = self.m + other.m
        for j in self.jrange(other.j):
            c_jm = self.clebsch_gordan(other.j, other.m, j, m)
            if c_jm:
                ket = c_jm*type(self)(j, m, alpha=(self.alpha * other.alpha))
                if abs(ket) > limit:
                    yield ket

    ## Project this Ket and another Ket onto JM.
    def clebsch_gordan(self, j, m, J, M):
        """self is j1m2...|j2,m2,J,M>"""
        return clebsch_gordan.clebsch_gordan(self.j, self.m, j, m, J, M)

    ## List of  Kets in Tensor Sum
    # \param other \f$|j'm'\rangle\f$
    # \kwd limit \f$=10^{-10}\f$ states below this abs are considered 0
    # \returns list \f$ c_{jj'mm'JM}|J,M>\rangle
    def expand(self, other, limit=1.e-10):
        """returns list(jm.iexpand(j'm'))"""
        return ITS(self.iexpand(other, limit=limit))

    ## Dilation, or multipliction across a list of kets
    def __rmul__(self, other):
        if isinstance(other, list):  # external poly
            return _finish([item*self for item in other])
        return self.dilation(other)

    ## Dilation
    def __div__(self, other):
        return self.dilation(1./other)

    ## Dilation
    # \param alpha \f$ \alpha'\f$
    # \returns \f$\alpha|jm\rangle \rightarrow \alpha'\alpha|jm\rangle\f$
    def dilation(self, alpha):
        """Dilate"""
        return type(self)(self.j, self.m, alpha=alpha*self.alpha)

    ## A nice ket string
    def __str__(self, root=" |{}, {}>"):
        # conver float int and  half int to non-float
        result = root.format(qn(self.j), qn(self.m))
        # nicefy sqrt of fractions
        if abs(self.alpha-1) > SMALL:  # format
            result = str(Nice(self.alpha)) + result
        # remove perfect sqaures
        for n in range(2, 100):
            new = result.replace("sqrt(1/{}".format(n**2), "1/{}".format(n))
            if new != result:
                new = new.replace(")", "")
            result = new
        return result

    def __unicode__(self):
        return str2unicode(self)

    def __repr__(self):
        name = type(self).__name__
        result = name + "({}, {}".format(self.j, self.m)
        if abs(self.alpha-1) > SMALL:  # format
            result += ", a={}".format(str(self.alpha))
        result += ")"
        return result

    ## Add (un)like kets (in)coherently
    def __add__(self, other):
        if self == other:
            return type(self)(self.j, self.m, alpha=(self.alpha+other.alpha))
        return ITS([self, other])

    def __sub__(self, other):
        return self + (-other)    

    ## 1st compare j-values, or try casimir.mcmp() on the m-values
    # \bugs overringd cmp breaks list.count (for real)
    def __cmp__(self, other, mcmp=casimir.mcmp):
        """Compare j or mcmp m."""
        return cmp(self.j, other.j) or mcmp(self.m, other.m)

    ## \f$ \alpha \f$
    def __complex__(self):
        """Normalization and phase factor."""
        return complex(self.alpha)

    ## \f$ \alpha \alpha^* \f$
    def __float__(self):
        return float((complex(self) * complex(self).conjugate()).real)

    ## \f$ |k| = \sqrt{\alpha\alpha^*} \f$
    def __abs__(self):
        return float(self) ** 0.5

    ## The ket's basis state.
    def hat(self):
        """Normalized (basis) Ket."""
        return type(self)(self.j, self.m)


## __I__ rreducible __T_ ensor __S__ um.
class ITS(list):
    """This list is an Irreducible Tensor Sum, so it's a list of
    Kets.
    """
    def __str__(self):
        return "[{}]".format(" +\n".join(map(str, self)))

    def __unicode__(self):
        return u"[{}]".format(u" +\n".join(map(unicode, self)))
    
    ## \f$ \sum_i{||\alpha_i||} \f$
    def __abs__(self):
        return sum(map(float, self)) ** 0.5

    ## Normalized Sum
    def hat(self):
        """Normalize by 1/abs"""
        return self / abs(self)

    ## Negate All Kets
    def __neg__(self):
        return type(self)(map(operator.neg, self))

    ## Posify All Kets
    def __pos__(self):
        return type(self)(map(operator.pos, self))

    ## Coherent addition:
    #  If |j'm'>
    def __add__(self, other):
        return coherent_add(self, other)

    def __sub__(self, other):
        return self + (-other)


    ## So here we need to deal with a another ITS, a Ket, or a
    # scale factor and whether we should auto simplify
    def __mul__(self, beta):
        # This is appropriate intimacy, as Ket multiplication is
        # handled by Ket.__rmul__
        print "its mul:", type(beta)
        if isinstance(beta, Ket):  # external poly
            print "pass to ket"
            return NotImplemented
        elif isinstance(beta, type(self)):
            print "nested its"
            result = []
            for ket in beta:
                result.extend(self*ket)
            return type(self)(result)
        print "dilation"
        return self.dilation(beta)

    ## This prevents list.__mul__ from concatenating, and invokes dilation.
    __rmul__ = __mul__

    def __div__(self, beta):
        return self.dilation(1./beta)

    def dilation(self, beta):
        """Dilate it."""
        return type(self)([beta*ket for ket in self])

    ## Project onto a ket
    # \param ket A ket: \f$ |j', m'\rangle \f$
    # \returns \f$ [\sum_{j, m}||{\langle j,m|j',m'\rangle||^2}]^
    # {\frac{1}{2}} \f$
    def __or__(self, ket):
        return sum([abs(item | ket)**2 for item in self]) ** 0.5


    ## Combine repeated Kets (linearly)
    # \kwd tol 1e-14 definition of zero
    # \returns None (This operators IN PLACE).
    def simplify(self, tol=1.e-14):
        """None = its.simlify(tol=1e-14) will simplify IN PLACE, with
        anything less than tol set to zero."""
        parts = []
        for ket in self:
            if self.count(ket) > 1:
                kets = []
                while ket in self:
                    kets.append(self.pop(self.index(ket)))
                total = sum(map(complex, kets))
                # TODO: quadrature or linear????
                if abs(total) > tol:
                    if total.imag == 0:  # type conversion
                        total = total.real
                    parts.append(Ket(ket.j, ket.m, alpha=total))
        if parts:
            self.extend(parts)
            
    ## Combine repeated Kets (linearly)
    # \kwd tol 1e-14 definition of zero
    # \returns ITS
    def simplified(self, tol=1.e-14, loops=1):
        """its.simlified(tol=1e-14) runs simplfy(tol=tol) and
        returns the simplified its."""
        for dum in xrange(0, min(10, max(0, loops), 1)):
            self.simplify(tol=tol)
            if len(self) == 1:
                return self[0]
        return self

                
## \f$ J_z \f$ Operator
# \param ket \f$ \alpha|j, m\rangle \f$
# \returns   \f$ m\alpha|j, m\rangle \f$
def Jz(ket):
    """Apply Jz to a ket"""
    return Ket(ket.j, ket.m, alpha=ket.m*ket.alpha)


## \f$ J_+ \f$ Raising Operator
# \param ket \f$ \alpha|j, m\rangle \f$
# \returns   \f$  \sqrt{j(j+1) - m(m+1)}\alpha|j, m+1\rangle \f$
def Jplus(ket):
    """Apply J+ to a ket"""
    return Ket(ket.j, ket.m+1,
               alpha=(clebsch_gordan.Jplus(ket.j, ket.m)*ket.alpha))


## \f$ J_- \f$ Lowering Operator
# \param ket \f$ \alpha|j, m\rangle \f$
# \returns   \f$  \sqrt{j(j+1) - m(m-1)}\alpha|j, m+1\rangle \f$
def Jminus(ket):
    """Apply J- to a ket"""
    return Ket(ket.j, ket.m-1,
               alpha=(clebsch_gordan.Jminus(ket.j, ket.m)*ket.alpha))


#u = Ket(0.5, 0.5)
#d = Ket(0.5, -0.5)

#uu = u*u
#ud = u*d
#du = d*u
#dd = d*d

#z = (ud+du)/2**0.5
#s = (ud-du)/2**0.5

#uuu = u*u*u
#uud = u*u*d
#udu = u*d*u
#udd = u*d*d
#duu = d*u*u
#dud = d*u*d
#ddu = d*d*u
#ddd = d*d*d


## Get a complete set of basis states for j.
# \param j
# \returns \f$ |j, m\rangle \f$ for all valid m.
def complete_set(j):
    """list = complete_set(j) --> [|j, m> for all m in (0,..,j, -j,.., -1)]"""
    return map(functools.partial(Ket, j), casimir.miter(j))


## Spinor Constructor
Spinor = functools.partial(Ket, HALF)


## Vector (boson) Constructor
Vector = functools.partial(Ket, 1)
Vector.__doc__ = "Vector(m) --> |1, m>"

## Spinor Basis
SPINOR = complete_set(HALF)


## Vector Basis
VECTOR = complete_set(1)


TENSOR = complete_set(2)


## Expand n-basis vectors.
def expand(n, j=1):
    """expand(n [,j=1]):

    yields (|m_0>, ..., |m_n>), [|J, M>, ...., |J', M'>]

    That is, it yields all combinations of n BASIS vectors,
    and their product expantion."""
    return ITS([reduce(operator.mul, item) for item in
                itertools.product(complete_set(j), repeat=n)])


## Get kets and their projections on to Irreps
def expantions(n, j=1):
    """see expantion"""
    return zip(itertools.product(complete_set(j), repeat=n), expand(n, j=1))


## Returns a list or sorted eigenstates
def eigenstates(n, j=1):
    """list = eigenstate(n [,j=1])

    sort all possible eigenstates for n compinations of j=1."""
    final_states = []
    for totals in expand(n, j=j):
        for item in (tot.hat() for tot in totals):
            if item not in final_states:
                final_states.append(item)
    return sorted(final_states)


## Solve m-projections that form a natural-form rank-n tensor
# \param n Rank of tensor
# \param m azimuth (order) parameter
# \param j=1 Defaults to Vector.
# \returns dictionary of unsymmetrized basis ket combinations.
def natural_form_basis_kets(n, m, j=1):
    """dict = natural_form_basis_kets(n, m)

    n is the rank of the  natural form tensor (made spin-j basis states)
    m is the TOTAL projection of z-angular momentum (order)
    j =1 (Vector default) defines the basis representation.

    dict is a dictionary of sorted m-combinations (tuple key) with lists of
    the ket combinations that have those combos.
    """
    result = collections.defaultdict(list)
    # loop over all spherical-dyad combinations
    #    for kets in itertools.imap(list, itertools.product(VECTOR, repeat=n)):
    for kets in itertools.product(complete_set(j), repeat=n):
        # get individual m values
        m_values = [ket.m for ket in kets]
        # filter m_1+...m_n == M (sum of m's is total M, which is 'm' param)
        if sum(m_values) == m:
            # now sort by m-states (since this is a symmetric combo)
            key = tuple(sorted(m_values))
            # put result in dictionary
            result[key].append(kets) 
    return result
            

## Solve total M states that form a natural-form rank-n tensor
# \param n Rank of tensor
# \param m azimuth (order) parameter
# \param j=1 Defaults to Vector
# \returns ITS of Kets that are |J=n, m=m>
def natural_form_eigenstates(n, m, j=1):
    """In development"""
    # prepare result
    result = {}
    # Get a nested list filtering on ket's total J
    f = functools.partial(filter, lambda dum: dum.j == n)
    # make a list multplier
    mul = functools.partial(reduce, operator.mul)
    # make a safe summer
    add = functools.partial(reduce, operator.add)
    # 1st get the basis ket combinations
    kets = natural_form_basis_kets(n, m, j=j)
    for key, value in kets.iteritems():
        # find total irreps built from individual kets.
        total_states = map(mul, value) 
        # now filter those states in heighest weight
        total_J_states = map(f, total_states)
        # sum up symmetric state combos (on a nested list)
        result[key] = add(map(add, total_J_states))
    return result


## Nicefy an irrational string
class Nice(object):
    """nice = Nice(0.70710678)

    >>>print nice
    sqrt(1/2)
    """
    
    ## Maximum denominator
    max_ = 2432
    
    def __init__(self, x):
        self.x = x

    def __float__(self):
        return float(self.x)

    def __abs__(self):
        return abs(float(self))
    
    ## Solve for < 1 case
    def sub_unit(self):
        for d in range(1, 1+self.max_):
            for n in range(1, d):
                if abs(operator.truediv(n, d) - float(self)**2) < 100*SMALL:
                    return n, d

    ## Solve for < 1 case
    def super_unit(self):
        data = (type(self)(1./float(self))).sub_unit()
        try:
            return list(reversed((data)))
        except TypeError:
            pass
        
    ## Find case and solve it.
    def __call__(self):
        if abs(float(self) - 1) < SMALL:
            return 1, 1
        elif float(self) > 1:
            return self.super_unit()
        else:
            return self.sub_unit()

    ## Solve and make a string for Clebsch-Gordan simplicity.
    def __str__(self):
        from math import sqrt
        try:
            n, d = self()
        except TypeError: # failed to find n and d
            return str(float(self))
        if float(self) < 0:
            s = '-'
        else:
            s = ''
        if d == 1:
            if int(sqrt(n)) == sqrt(n):
                return "{}{}".format(s, int(sqrt(n)))
            else:
                return "{}sqrt({})".format(s, n)
        else:
            if int(sqrt(n)) == sqrt(n) and int(sqrt(d)) == sqrt(d):
                return "{}{}/{}".format(s, int(sqrt(n)), int(sqrt(d)))
            else:
                return "{}sqrt({}/{})".format(s, n, d)

    def __unicode__(self):
        return str2unicode(str(self))
        


## Combine 2 angular momenta into total J, M.
class Combinator(object):
    """combo = Combinator(j1, j2)

    >>>print combo
    ...
    |J, M> = a|m1, m2> + b|m1', m2'> +...
    ...
    ...

    for all states.

    where |J, M> is combo[J, M].

    For exmaple:

    In [8]: c = racah.Combinator(0.5, 0.5)

    In [9]: print c
    |0.0, -0.0> =  + sqrt(1/2)|0.5:-0.5> - sqrt(1/2)|-0.5:0.5>
    |1.0, 0.0> =  + sqrt(1/2)|0.5:-0.5> + sqrt(1/2)|-0.5:0.5>
    |1.0, 1.0> =  + sqrt(1)|0.5:0.5>
    |1.0, -1.0> =  + sqrt(1)|-0.5:-0.5>


   In [12]: c = racah.Combinator(1,1)

    In [13]: print c
    |0.0, -0.0> =  - sqrt(1/3)|0:0> + sqrt(1/3)|1:-1> + sqrt(1/3)|-1:1>
    |1.0, 0.0> =  + sqrt(1/2)|1:-1> - sqrt(1/2)|-1:1>
    |1.0, 1.0> =  - sqrt(1/2)|0:1> + sqrt(1/2)|1:0>
    |1.0, -1.0> =  + sqrt(1/2)|0:-1> - sqrt(1/2)|-1:0>
    |2.0, 0.0> =  + sqrt(2/3)|0:0> + sqrt(1/6)|1:-1> + sqrt(1/6)|-1:1>
    |2.0, 1.0> =  + sqrt(1/2)|0:1> + sqrt(1/2)|1:0>
    |2.0, 2.0> =  + sqrt(1)|1:1>
    |2.0, -2.0> =  + sqrt(1)|-1:-1>
    |2.0, -1.0> =  + sqrt(1/2)|0:-1> + sqrt(1/2)|-1:0>
    """

    ## Construct from 2 spin quantum numbers
    def __init__(self, j1, j2):
        self.j1 = j1 if j1 <= j2 else j2
        self.j2 = j2 if j2 > j1 else j1

    ## Yield are pairs of $(m_1, m_2)$.
    def _iterms(self):
        return itertools.product(casimir.miter(self.j1),
                                 casimir.miter(self.j2))
    
    ## Yield are ket pairs $(|j_1, m_1\rangle ,|j_2, m_2\rangle)$.
    def _iterkets(self):
        for m1, m2 in self._iterms():
            k1 = Ket(self.j1, m1)
            k2 = Ket(self.j2, m2)
            yield k1, k2


    ## Yield all azimuthal numbers and ket _products_ in original basis:\n
    # $(m_1, m_2, |j_1, m_1\rangle |j_2, m_2\rangle)$.
    def _iterprods(self):
        for (m1, m2), (k1, k2) in itertools.izip(self._iterms(),
                                                 self._iterkets()):
            kk = k1*k2
            yield m1, m2, kk

    ## Yield all azimuthal numbers and ket _products_ in final basis:\n
    # $(m_1, m_2, [|J, m_1+m_2\rangle, \cdots  |J', m_1+m_2\rangle])$. 
    def expandprods(self):
        result = collections.defaultdict(list)
        for m1, m2, K in self._iterprods():
            #print "|{},{}>|{},{}> = {}".format(self.j1, m1, self.j2, m2, K)
            #print
            for k in K:
                J = float(k.j)
                M = float(k.m)
                alpha = complex(k)
                
                result[(J,M)].append([alpha, m1, m2])
        return result

    ## combinator[J, M] get all term of the form
    # $c|j_1, m_1\rangle|j_2, m_2\rangle$ in  $|J, M\rangle$
    def __getitem__(self, index):
        try:
            j, m = index
        except TypeError:
            raise IndexError("Use [J, M]")
        if not hasattr(self, 'prods'):
            self.prods = dict(self.expandprods())
        return self.prods.get(index)

    ## Yield all non-zero results from Combinator.__getitem__()
    def __iter__(self):
        for j in clebsch_gordan.jrange(self.j1, self.j2):
            for m in casimir.miter(j):
                yield j, m, self[j, m]

    ## Print all non-zero combinations from Combinator.__iter__()
    def __str__(self):
        result = []
        for count, (j, m, kets) in enumerate(self):
            front = "|{}, {}> = ".format(*map(qn, (j, m)))
            end = ""
            for alpha, m1, m2 in kets:
                if alpha.real < 0:
                    end += " - "
                else:
                    if count:
                        end += " + "
                    else:
                        end += " + "
                end += str(Nice(abs(alpha))) + "|{}:{}>".format(
                    *map(qn, (m1, m2)))

            result.append(front + end)
        return "\n".join(result)

    def __unicode__(self):
        return str2unicode(str(self))

        
    
## quantum number string function
# \param m str/int  int or 1/2-int quantum number (maybe as a float)
# \returns str nicer format 
def qn(m):
    """make a bar/ket string for l or m or j or jz"""
    try:
        m_int = int(m)
    except (TypeError, ValueError):
        print "qn fail on : ", str(m)
        return str(m)

    if m == m_int:
        return str(m_int)
    elif m < 0:
        return "-" + qn(-m)
    elif m < 1:
        return "1/2"
    else:
        n = int(round(m-HALF))
        return "{} 1/2".format(n)

## Convert standard Clebsch-Gordan format form ASCII to unicode
# \param x An ASCII string like "sqrt(1/3)|1, 1/2>"
# \param _maxdem=8 Maximum unicode vulgar fraction demoninator
# \returns unicode looking like $\f$\sqrt{\frac{1}{3}}|1, \frac 1 2 \rangle$
def str2unicode(x, _maxdem=8):
    """convert kets and CG factors to unicode, as much as possible"""
    from . import uc
    # start with a unicode version of x.__str__()
    result = unicode(str(x))
    # replace known unicode fractions
    for i in range(1, _maxdem+1):
        for j in range(2, _maxdem+1):
            try:
                x = uc.VULGAR[i]
            except KeyError:
                continue
            try:
                y = x[j]
            except KeyError: 
                continue
            for c in ",>):":
                org = "{}/{}{}".format(i, j, c)
                result = result.replace(org, y + c)
            org = "{}/{} ".format(i, j)
            result = result.replace(org, y)    
    # unicdoe bra-kets
    result = result.replace("<", uc.bra)
    result = result.replace(">", uc.ket)
    # remove space between quantum number int and 1/2 int
    for i in range(1, 100):
        result = result.replace("|{} ".format(i), "|{}".format(i))
    for key, value in uc.SYMBOLS.items():
        result = result.replace(key, value)
    return result
    

## Spin up
U = Ket(HALF, HALF)

## SPin down
D = Ket(HALF, -HALF)

## + Circular Polarization
P = Ket(1, 1)

## Transverse Polarization
Z = Ket(1, 0)

## - Circular Polarization
M = Ket(1, -1)


def timer():
    import time
    import collections
    RES = collections.namedtuple("result", "j time string")
    j = 0.5
    result = []
    try:
        while 1:
            a = Combinator(j, j)
            start = time.time()
            b = unicode(a)
            stop = time.time() - start
            res = RES(j, stop, b)
            result.append(res)
            print j
            j += 0.5
    except KeyboardInterrupt:
        return result
