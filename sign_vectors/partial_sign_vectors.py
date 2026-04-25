r"""
Auxiliary class for functions
=============================

There are several ways to define partial sign vectors::

    sage: from sign_vectors.partial_sign_vectors import *
    sage: from sign_vectors import *
    sage: partial_sign_vector("+/-+p0n")
    (+/-+p0n)
    sage: partial_sign_vector([[1], [0, 1], [-1, 0, 1], [1]])
    (+p*+)
    sage: v = vector([5, 2/5, -1, 0])
    sage: v = sign_vector(v)
    sage: partial_sign_vector(v)
    (++-0)

We construct the zero partial sign vector of a given length::

    sage: PartialSignVector.zero(6)
    (000000)

We construct the star partial sign vector of a given length::

    sage: PartialSignVector.star(6)
    (******)

There are different notions of support::

    sage: X = partial_sign_vector("+/**pn0--")
    sage: X.zero_support()
    [2, 3, 4, 5, 6]
    sage: X.positive_support()
    [0, 1, 2, 3, 4]
    sage: X.negative_support()
    [1, 2, 3, 5, 7, 8]

There are different notions of determined support::

    sage: X = partial_sign_vector("+/**pn0--")
    sage: X.determined_support()
    [0, 6, 7, 8]
    sage: X.determined_zero_support()
    [6]
    sage: X.determined_positive_support()
    [0]
    sage: X.determined_negative_support()
    [7, 8]

There is the also the undetermined support::

    sage: X.undetermined_support()
    [2, 3]

Next, we define two sign vectors and compose them::

    sage: X = partial_sign_vector("-0pn0")
    sage: Y = partial_sign_vector("0n/+0")
    sage: X.compose(Y)
    (-n//0)
    sage: Y.compose(X)
    (-n/+0)

The issubset function is a partial order on a set of partial sign vectors::

    sage: X = partial_sign_vector("-+00+")
    sage: Y = partial_sign_vector("-/p0+")
    sage: Z = partial_sign_vector("-+**0")
    sage: X.issubset(Y)
    True
    sage: Y.issubset(Z)
    False
    sage: X < Y
    True
    sage: X <= Z
    False
"""

#############################################################################
#  Copyright (C) 2026                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#          Arne Jenß (arne.jenss@uni-kassel.de)                             #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from __future__ import annotations
from types import NoneType

from sage.data_structures.bitset import FrozenBitset
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from . import SignVector
from itertools import permutations


def partial_sign_vector(iterable: list[list[int]] | str | SignVector | PartialSignVector) -> PartialSignVector:
    r"""
    Create a partial sign vector from a list, vector or string.

    INPUT:

    - ``iterable`` -- different inputs are accepted:

        - an iterable (e.g. a list or vector) of lists containing ``-1``, ``0``, ``1``.

        - a string consisting of ``-``, ``0``, ``+``, ``n``, ``/``, ``p``, ``*``.

    OUTPUT:

    Returns a partial sign vector.

    EXAMPLES::

        sage: from sign_vectors.partial_sign_vectors import *
        sage: from sign_vectors import *
        sage: partial_sign_vector([[1], [0], [-1, 1] , [-1, 0]])
        (+0/n)
        sage: v = vector([5, 0, -1, -2])
        sage: v = sign_vector(v)
        sage: partial_sign_vector(v)
        (+0--)

    We can also use a string to construct a partial sign vector::

        sage: partial_sign_vector("00--")
        (00--)
        sage: partial_sign_vector("++n*-0/p")
        (++n*-0/p)

    TESTS::

        sage: partial_sign_vector(partial_sign_vector("+-+"))
        (+-+)
    """
    if isinstance(iterable, PartialSignVector):
        return iterable
    if isinstance(iterable, SignVector):
        return PartialSignVector.from_sign_vector(iterable)
    if isinstance(iterable, str):
        return PartialSignVector.from_string(iterable)
    return PartialSignVector.from_iterable(iterable)


class PartialSignVector(SageObject):
    r"""A partial sign vector."""

    __slots__ = ("_negative_support","_zero_support", "_positive_support" )

    def __init__(self, negative_support: FrozenBitset, zero_support: FrozenBitset, positive_support: FrozenBitset) -> None:
        r"""
        Create a partial sign vector object.

        INPUT:

        - ``negative_support`` -- a ``FrozenBitset``
        - ``zero_support`` -- a ``FrozenBitset``
        - ``positive_support`` -- a ``FrozenBitset``

        .. NOTE::

            The union of``negative_support``, ``zero_support`` and ``positive_support``should be a full set.
            For efficiency, this is not checked.

        .. SEEALSO::

            - :func: `~partial_sign_vector`

        TESTS::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: PartialSignVector(FrozenBitset([1, 3], capacity=4), FrozenBitset([0], capacity=4), FrozenBitset([2,3], capacity=4))
            (0-+/)
        """
        self._negative_support = negative_support
        self._zero_support = zero_support
        self._positive_support = positive_support

    def _repr_(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return f"({self.to_string()})"

    def __getitem__(self, e: int | slice) -> int | SignVector:
        r"""
        Return the element at position ``e`` of the sign vector.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("0+*n")
            sage: X
            (0+*n)
            sage: X[0]
            [0]
            sage: X[1]
            [1]
            sage: X[2]
            [-1, 0, 1]
            sage: X[-1]
            [-1, 0]
            sage: X[1:3]
            (+*)

        TESTS::

            sage: X[-2]
            [-1, 0, 1]
            sage: X[100]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if isinstance(e, slice):
            # TODO improve
            return PartialSignVector.from_string(self.to_string()[e])
        if e >= self.length() or e < -self.length():
            raise IndexError("index out of range")
        if e < 0:
            e %= self.length()
        res = []
        if e in self._negative_support:
            res.append(-1)
        if e in self._zero_support:
            res.append(0)
        if e in self._positive_support:
            res.append(1)

        return res

    def to_string(self) -> str:
        r"""
        Return a string representation of this partial sign vector (without parentheses).

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("0+-")
            sage: X
            (0+-)
            sage: X.to_string()
            '0+-'

        Note that that `str` and `to_string` are different::

            sage: str(X)
            '(0+-)'
        """
        return "".join(self._index_to_string(e) for e in range(self.length()))

    def _index_to_string(self, index) -> str:
        match self[index]:
            case [-1]:
                return "-"
            case [0]:
                return "0"
            case [1]:
                return "+"
            case [-1, 0]:
                return "n"
            case [-1, 1]:
                return "/"
            case [0, 1]:
                return "p"
            case [-1, 0, 1]:
                return "*"

    def length(self) -> int:
        r"""
        Return the length of the partial sign vector.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("0+-")
            sage: X
            (0+-)
            sage: X.length()
            3
        """
        return self._negative_support.capacity()

    def __len__(self) -> int:
        r"""
        Return the length of this partial sign vector.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("0+-")
            sage: X
            (0+-)
            sage: len(X)
            3
        """
        return self.length()

    def negative_support(self) -> list[int]:
        r"""
        Return a list of indices where the partial sign vector can be negative.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-n+*/")
            sage: X
            (-n+*/)
            sage: X.negative_support()
            [0, 1, 3, 4]
        """
        return list(self._negative_support)

    def zero_support(self) -> list[int]:
        r"""
        Return a list of indices where the partial sign vector can be zero.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-n+*/")
            sage: X
            (-n+*/)
            sage: X.zero_support()
            [1, 3]
        """
        return list(self._zero_support)

    def positive_support(self) -> list[int]:
        r"""
        Return a list of indices where the partial sign vector can be positive.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-n+*/")
            sage: X
            (-n+*/)
            sage: X.positive_support()
            [2, 3, 4]
        """
        return list(self._positive_support)

    def cardinality(self) -> int:
        r"""
        Return the number of sign vectors packed in this partial sign vector .

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-n+*/")
            sage: X
            (-n+*/)
            sage: X.cardinality()
            12
            sage: len(X.unpack())
            12
        """
        return 2 ** (len(self) - len(self.determined_support()) - len(self.undetermined_support())) * 3 ** (len(self.undetermined_support()))

    def _determined_negative_support(self) -> FrozenBitset:
        return self._negative_support - (self._positive_support | self._zero_support)

    def _determined_zero_support(self) -> FrozenBitset:
        return self._zero_support - (self._negative_support | self._positive_support)

    def _determined_positive_support(self) -> FrozenBitset:
        return self._positive_support - (self._negative_support | self._zero_support)

    def _determined_support(self) -> FrozenBitset:
        return self._determined_negative_support() | self._determined_zero_support() | self._determined_positive_support()

    def _undetermined_support(self) -> FrozenBitset:
        return self._positive_support & self._negative_support & self._zero_support

    def determined_negative_support(self) -> list[int]:
        r"""
        Return a list of indices where the partial sign vector is negative.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-n+0/")
            sage: X
            (-n+0/)
            sage: X.determined_negative_support()
            [0]
        """
        return list(self._determined_negative_support())

    def determined_zero_support(self) -> list[int]:
        r"""
        Return a list of indices where the partial sign vector is zero.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-n+0/")
            sage: X
            (-n+0/)
            sage: X.determined_zero_support()
            [3]
        """
        return list(self._determined_zero_support())

    def determined_positive_support(self) -> list[int]:
        r"""
        Return a list of indices where the partial sign vector is positive.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-n+0/")
            sage: X
            (-n+0/)
            sage: X.determined_positive_support()
            [2]
        """
        return list(self._determined_positive_support())

    def determined_support(self) -> list[int]:
        r"""
        Return a list of indices where the partial sign vector is determined.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-n+0/")
            sage: X
            (-n+0/)
            sage: X.determined_support()
            [0, 2, 3]
        """
        return list(self._determined_support())

    def undetermined_support(self) -> list[int]:
        r"""
        Return a list of indices where the partial sign vector is undetermined.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-n+*/")
            sage: X
            (-n+*/)
            sage: X.undetermined_support()
            [3]
        """
        return list(self._undetermined_support())

    def list_from_positions(self, positions: list[int]) -> list[int]:
        r"""
        Return a list of components that are in the list of indices ``positions``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-+0*p")
            sage: X
            (-+0*p)
            sage: X.list_from_positions([0, 1, 4])
            [[-1], [1], [0, 1]]
        """
        return [self[e] for e in positions]

    def flip_signs(self, indices: list[int]) -> PartialSignVector:
        r"""
        Flips entries of given indices.

        INPUT:

        - ``indices`` -- list of indices

        OUTPUT:
        Returns a new sign vector. Components of ``indices`` are multiplied by ``-1``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-+p0+")
            sage: X
            (-+p0+)
            sage: X.flip_signs([0, 2, 3])
            (++n0+)
        """
        indices = FrozenBitset(indices) & (self._positive_support | self._negative_support)
        return self.__class__(self._negative_support ^ indices , self._zero_support, self._positive_support ^ indices)

    def set_sign(self, indices: list[int], sign:list[int] | str) -> PartialSignVector:
        r"""
        Set given entries to selected sign.

        INPUT:

        - ``indices`` -- list of indices
        - ``sign`` -- sign as integer list or string

        OUTPUT:
        Returns a new sign vector of same length. Components with indices in
        ``indices`` are set to ``-``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-++0+")
            sage: X
            (-++0+)
            sage: X.set_sign([0, 2, 3], "*")
            (*+**+)
            sage: X.set_sign([0, 2, 3], "n")
            (n+nn+)
            sage: X.set_sign([0, 2, 3], [0])
            (0+00+)

        """
        indices = FrozenBitset(indices)
        negative_support = self._negative_support
        zero_support = self._zero_support
        positive_support = self._positive_support

        if isinstance(sign, str):
            if sign in ["-", "n","/","*"]:
                negative_support = negative_support | indices
            else:
                negative_support = negative_support - indices

            if sign in ["0", "n","p","*"]:
                zero_support = zero_support | indices
            else:
                zero_support = zero_support - indices

            if sign in ["+", "p","/","*"]:
                positive_support = positive_support | indices
            else:
                positive_support = positive_support - indices

            return self.__class__(negative_support, zero_support, positive_support)

        if isinstance(sign, list):
            if -1 in sign:
                negative_support = negative_support | indices
            else:
                negative_support = negative_support - indices

            if 0 in sign:
                zero_support = zero_support | indices
            else:
                zero_support = zero_support - indices

            if 1 in sign:
                positive_support = positive_support | indices
            else:
                positive_support = positive_support - indices

            return self.__class__(negative_support, zero_support, positive_support)

    def delete_components(self, indices: list[int]) -> PartialSignVector:
        r"""
        Delete the given components from the partial sign vector.

        INPUT:

        - ``indices`` -- list of indices to delete

        OUTPUT:

        Returns a new partial sign vector with the specified components removed.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("+*0p-")
            sage: X
            (+*0p-)
            sage: X.delete_components([1, 3])
            (+0-)

        """
        return self.__class__.from_iterable(
            [self[i] for i in ~FrozenBitset(indices, capacity=self.length())]
        )

    def _connecting_elements(self, other: SignVector | PartialSignVector) -> FrozenBitset:
        if isinstance(other, SignVector):
            other = ExtendedSignVector.from_sign_vector(other)
        return (self._determined_positive_support() & other._determined_positive_support()) | (self._determined_negative_support() & other._determined_negative_support())

    def _separating_elements(self, other: SignVector | PartialSignVector) -> FrozenBitset:
        if isinstance(other, SignVector):
            other = ExtendedSignVector.from_sign_vector(other)
        return (self._determined_positive_support() & other._determined_negative_support()) | (self._determined_negative_support()  & other._determined_positive_support())

    def connecting_elements(self, other: SignVector | PartialSignVector) -> list[int]:
        r"""
        Compute the list of connecting elements of two (partial) sign vectors.

        INPUT:

        - ``other`` -- a sign vector or partial sign vector

        OUTPUT:
        List of elements ``e`` such that ``self[e] == -other[e] != 0``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("+*-p-")
            sage: X
            (+*-p-)
            sage: Y = partial_sign_vector("+--++")
            sage: Y
            (+--++)
            sage: X.connecting_elements(Y)
            [0, 2]
        """
        return list(self._connecting_elements(other))

    def separating_elements(self, other: SignVector | PartialSignVector) -> list[int]:
        r"""
        Compute the list of separating elements of two (partial) sign vectors.

        INPUT:

        - ``other`` -- a sign vector or partial sign vector

        OUTPUT:
        List of elements ``e`` such that ``self[e] == -other[e] != 0``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("+*-p-")
            sage: X
            (+*-p-)
            sage: Y = partial_sign_vector("+-+++")
            sage: Y
            (+-+++)
            sage: X.separating_elements(Y)
            [2, 4]
        """
        return list(self._separating_elements(other))

    def compose(self, other: SignVector | PartialSignVector) -> PartialSignVector:
        r"""
        Return the composition of two (partial) sign vectors.

        INPUT:

        - ``other`` -- a sign vector or partial sign vector

        OUTPUT:
        Composition of this partial sign vector with ``other``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("p*0")
            sage: X
            (p*0)
            sage: Y = partial_sign_vector("--0")
            sage: Y
            (--0)
            sage: X.compose(Y)
            (//0)
            sage: Y.compose(X)
            (--0)
            sage: X = partial_sign_vector("0p0+++---")
            sage: Y = partial_sign_vector("/+-0+-0+-")
            sage: X.compose(Y)
            (/+-+++---)
        """
        return self.__class__(
            self._negative_support | (other._negative_support & self._zero_support),
            self._zero_support & other._zero_support,
            self._positive_support | (other._positive_support & self._zero_support)
        )

    def is_orthogonal_to(self, other: SignVector | PartialSignVector) -> bool:
        r"""
        Return whether two (partial) sign vectors are orthogonal.

        INPUT:

        - ``other`` -- a sign vector or partial sign vector

        OUTPUT:
        - Return ``True`` if the (partial) sign vectors are orthogonal and ``False`` otherwise.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("n*0-+")
            sage: X
            (n*0-+)
            sage: X.is_orthogonal_to(X)
            False
            sage: X.is_orthogonal_to(partial_sign_vector("**0++"))
            True
            sage: X.is_orthogonal_to(partial_sign_vector("00+00"))
            True
            sage: X.is_orthogonal_to(partial_sign_vector("0*++0"))
            False
        """
        if isinstance(other, SignVector):
            other = ExtendedSignVector.from_sign_vector(other)

        if (~(self._determined_zero_support() | other._determined_zero_support())).isempty():
            return True

        return not(self._separating_elements(other).isempty() or self._connecting_elements(other).isempty())

    def issubset(self, other: PartialSignVector) -> bool:
        r"""
        subset relation of two partial sign vectors.

        .. NOTE::

            Alternatively, the operator ``<=`` can be used.
            Use ``>=``, ``<`` and ``>`` for the other relations.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("n*0-+")
            sage: X
            (n*0-+)
            sage: Y = partial_sign_vector("-+0-+")
            sage: Y
            (-+0-+)
            sage: Z = partial_sign_vector("-*0-+")
            sage: Z
            (-*0-+)
            sage: Y.issubset(X)
            True
            sage: X.issubset(X)
            True
            sage: Y.issubset(Z)
            True
            sage: Z.issubset(Y)
            False
            sage: Z.issubset(X)
            True
        """
        return self._positive_support.issubset(other._positive_support) and self._negative_support.issubset(other._negative_support) and self._zero_support.issubset(other._zero_support)

    def _intersection(self, other: PartialSignVector) -> PartialSignVector|NoneType:
        r"""
        Return the intersection of two partial sign vectors.

        INPUT:

        - ``other`` -- a partial sign vector

        OUTPUT:
        The intersection of this partial sign vector with ``other``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("n+0-+")
            sage: X
            (n+0-+)
            sage: Y = partial_sign_vector("-*0-+")
            sage: Y
            (-*0-+)
            sage: X.intersection(Y)
            (-+0-+)
            sage: X = partial_sign_vector("**-")
            sage: Y = partial_sign_vector("+++")
            sage: X.intersection(Y) is None
            True
        """
        n_support = self._negative_support & other._negative_support
        z_support = self._zero_support & other._zero_support
        p_support = self._positive_support & other._positive_support

        empty_support = ~(n_support | z_support | p_support)
        if empty_support.isempty():
            return self.__class__(n_support, z_support, p_support)
        return None

    def intersection(self, other: PartialSignVector|set[PartialSignVector]) -> PartialSignVector|NoneType:
        r"""
        Return the intersection of two or more partial sign vectors.

        INPUT:

        - ``other`` -- a partial sign vector

        OUTPUT:
        The intersection of this partial sign vector with ``other``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("n+0-+")
            sage: X
            (n+0-+)
            sage: Y = partial_sign_vector("-*0-+")
            sage: Y
            (-*0-+)
            sage: X.intersection(Y)
            (-+0-+)
        """
        if isinstance(other, PartialSignVector):
            return self._intersection(other)

        inter = self
        for o in other:
            inter = inter._intersection(o)
            if inter is None:
                return None
        return inter

    def _setminus(self, other: PartialSignVector) -> set[PartialSignVector]:
        r"""
        Return the set difference of two partial sign vectors.

        INPUT:

        - ``other`` -- a partial sign vector

        OUTPUT:

        The set difference of this partial sign vector with ``other``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("n+0-+")
            sage: X
            (n+0-+)
            sage: Y = partial_sign_vector("-*0-+")
            sage: Y
            (-*0-+)
            sage: Z = partial_sign_vector("*****")
            sage: X._setminus(Y)
            {(0+0-+)}
            sage: X._setminus(Z)
            set()
            sage: Z._setminus(X)
            {(+****), (*n***), (****n), (**/**), (***p*)}
            sage: Z.setminus(Y)
            {(**/**), (***p*), (p****), (****n)}
        """
        inter = self.intersection(other)
        if inter is None:
            return {self}
        if inter == self:
            return set()

        n_minus = self._negative_support - inter._negative_support
        z_minus = self._zero_support - inter._zero_support
        p_minus = self._positive_support - inter._positive_support

        minus_support = n_minus | z_minus | p_minus

        res = set()
        length = self.length()
        for i in minus_support:
            index = FrozenBitset([i], capacity=length)

            n_support = self._negative_support - index if i not in n_minus else self._negative_support
            z_support = self._zero_support - index if i not in z_minus else self._zero_support
            p_support = self._positive_support - index if i not in p_minus else self._positive_support
            res.add(self.__class__(n_support, z_support, p_support))

        return res

    def setminus(self, other: PartialSignVector | set[PartialSignVector]) -> set[PartialSignVector]:
        r"""
        Return the set difference of two or more partial sign vectors.

        INPUT:

        - ``other`` -- a partial sign vector

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("n+0-+")
            sage: X
            (n+0-+)
            sage: Y = partial_sign_vector("-*0-+")
            sage: Y
            (-*0-+)
            sage: Z = partial_sign_vector("*****")
            sage: X.setminus(Y)
            {(0+0-+)}
            sage: X.setminus(Z)
            set()
            sage: Z.setminus({X,Y})
            {(+****), (pn***), (p**p*), (p***n), (p*/**), (**/**), (***p*), (****n)}
        """
        if isinstance(other, PartialSignVector):
            return self._setminus(other)

        res = {self}
        for o in other:
            new_res = set()
            for r in res:
                new_res = new_res.union(r._setminus(o))
            res = new_res
        return  res

    def union(self, other: PartialSignVector) -> set[PartialSignVector]:
        r"""
        Return the disjoint union of two partial sign vectors.

        INPUT:

        - ``other`` -- a partial sign vector

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("n+0-+")
            sage: X
            (n+0-+)
            sage: Y = partial_sign_vector("-*0-+")
            sage: Y
            (-*0-+)
            sage: X.union(Y)
            {(-n0-+), (n+0-+)}
        """
        if self.issubset(other):
            return {other}
        if other.issubset(self):
            return {self}

        return {self}.union(other.setminus(self))

    def complement(self) -> set[PartialSignVector]:
        r"""
        Return the complement of this partial sign vector.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("n*0-+")
            sage: X
            (n*0-+)
            sage: C = X.complement()
            sage: C
            {(**/**), (+****), (***p*), (****n)}
            sage: CC = set().union(*(c.complement() for c in C))
            sage: CC
            {(****+), (***-*), (n****), (**0**)}
        """
        return PartialSignVector.star(self.length()).setminus(self)

    def __eq__(self, other: SignVector | PartialSignVector) -> bool:
        r"""
        Return whether this partial sign vector equals ``other``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("++0-")
            sage: X == X
            True
            sage: X == partial_sign_vector("0*++")
            False

        TESTS::

            sage: PartialSignVector.zero(3) == 0
            True
            sage: 0 == PartialSignVector.zero(3)
            True
        """
        if isinstance(other, (SignVector, ExtendedSignVector, PartialSignVector)):
            return (
                self._negative_support == other._negative_support
                and self._positive_support == other._positive_support
                and self._zero_support == other._zero_support
            )
        if other == 0:
            return not (self._positive_support | self._negative_support)
        return False

    def __le__(self, other: PartialSignVector) -> bool:
        r"""
        Return whether this partial sign vector is less or equal to ``other``.

        .. SEEALSO::

            - :meth: `issubset`

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-+00+")
            sage: X
            (-+00+)
            sage: Y = partial_sign_vector("-p**+")
            sage: Y
            (-p**+)
            sage: X <= Y
            True

        We can also use ``<=`` to compare a sign vector with ``0``::

            sage: partial_sign_vector("00--") <= 0
            False
            sage: PartialSignVector.zero(3) <= 0
            True
        """
        if isinstance(other, PartialSignVector):
            return self.issubset(other)
        if other == 0:
            return self == 0
        return False

    def __lt__(self, other: PartialSignVector) -> bool:
        r"""
        Return whether this partial sign vector is less than ``other``.

        .. SEEALSO::

            - :meth: `issubset`

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-+00+")
            sage: X
            (-+00+)
            sage: Y = partial_sign_vector("-p**+")
            sage: Y
            (-p**+)
            sage: X < Y
            True

        We can also use ``<`` to compare a sign vector with ``0``::

            sage: partial_sign_vector("00--") < 0
            False
            sage: PartialSignVector.zero(3) < 0
            False
        """
        return self != other and self <= other

    def __ge__(self, other: PartialSignVector) -> bool:
        r"""
        Return whether this partial sign vector is greater or equal to ``other``.

        .. SEEALSO::

            - :meth: `issubset`

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-+00+")
            sage: X
            (-+00+)
            sage: Y = partial_sign_vector("-p**+")
            sage: Y
            (-p**+)
            sage: Y >= X
            True

        We can also use ``>=`` to compare a sign vector with ``0``::

            sage: partial_sign_vector("00nn") >= 0
            True
            sage: PartialSignVector.zero(3) >= 0
            True
        """
        if isinstance(other, PartialSignVector):
            return other.issubset(self)
        if other == 0:
            return all( 0 in a for a in self)
        return False

    def __gt__(self, other: PartialSignVector) -> bool:
        r"""
        Return whether this partial sign vector is greater than ``other``.

        .. SEEALSO::

            - :meth: `conforms`

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("-+00+")
            sage: X
            (-+00+)
            sage: Y = partial_sign_vector("-p**+")
            sage: Y
            (-p**+)
            sage: Y > X
            True

        We can also use ``>`` to compare a sign vector with ``0``::

            sage: 0 >  partial_sign_vector("00--")
            False
            sage: PartialSignVector.zero(3) > 0
            False
        """
        return self != other and self >= other

    def __bool__(self) -> bool:
        return self != 0

    def __neg__(self) -> PartialSignVector:
        r"""
        Return this partial sign vectors multiplied by ``-1``.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("pn+-*/")
            sage: -X
            (np-+*/)
        """
        return self.__class__(self._positive_support, self._zero_support, self._negative_support)

    def __pos__(self) -> PartialSignVector:
        r"""
        Return this partial sign vectors multiplied by ``1`` which is a copy of this object.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("pn+-*/")
            sage: X
            (pn+-*/)
            sage: +X
            (pn+-*/)
        """
        return self.__class__(self._negative_support, self._zero_support, self._positive_support)

    def is_sign_vector(self) -> bool:
        r"""
        Return if this partial sign vector is a sign vector.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("pn+-")
            sage: X
            (pn+-)
            sage: X.is_sign_vector()
            False
            sage: X = partial_sign_vector("0++-")
            sage: X
            (0++-)
            sage: X.is_sign_vector()
            True
        """
        return len(self) == len(self.determined_support())

    def to_sign_vector(self, extended=False) -> bool:
        r"""Convert this partial sign vector to a sign vector."""
        if self.is_sign_vector():
            if extended:
                return ExtendedSignVector(self._positive_support, self._negative_support)
            return SignVector(self._positive_support, self._negative_support)
        raise TypeError("not a sign vector")

    def unpack(self, indices: list[int] = None) -> set[SignVector] | set[PartialSignVector]:
        r"""
        Return this partial sign vector as a set of (partial) sign vectors, where the given indices are determined.

        INPUT:

        - ``indices`` -- a list of indices, an integer or None

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("p*0")
            sage: X
            (p*0)
            sage: X.unpack([0])
            {(+*0), (0*0)}
            sage: X.unpack([1, 2])
            {(p00), (p+0), (p-0)}
            sage: X.unpack()
            {(000), (+-0), (0+0), (++0), (0-0), (+00)}
        """
        #TODO: parallelizing
        cast_sv = False
        res = set()
        res.add(self)
        temp = set()

        if indices is None:
            indices = range(len(self))
            cast_sv = True

        for e in indices:
            for psv in res:
                for s in psv[e]:
                    temp.add(psv.set_sign([e],[s]))
            res = temp
            temp = set()

        if cast_sv:
            return set(sv.to_sign_vector() for sv in res)

        return res

    def lower_closure(self) -> PartialSignVector:
        r"""
        Return the lower closure of this partial sign vector.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("pn+-*/")
            sage: X.lower_closure()
            (pnpn**)
        """
        return self.__class__(self._negative_support,
                               ~FrozenBitset([], capacity=self.length()),
                               self._positive_support)

    def upper_closure(self) -> PartialSignVector:
        r"""
        Return the upper closure of this partial sign vector.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("pn+-*/")
            sage: X.upper_closure()
            (**+-*/)
        """
        return self.__class__(self._negative_support | self._zero_support,
                               self._zero_support,
                               self._positive_support | self._zero_support)

    def closure(self) -> PartialSignVector:
       r"""
        Return the closure of this partial sign vector.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = partial_sign_vector("pn+-*/")
            sage: X.closure()
            (**pn**)
        """
       return self.__class__(self._negative_support | self._zero_support,
                               ~FrozenBitset([], capacity=self.length()),
                               self._positive_support | self._zero_support)

    def __hash__(self) -> int:
        r"""Return the hash value of this partial sign vector."""

        return hash((self._negative_support, self._zero_support, self._positive_support))

    @classmethod
    def zero(cls, length: int) -> PartialSignVector:
        r"""
        Create a zero partial sign vector of a given length.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: PartialSignVector.zero(5)
            (00000)
        """
        return cls(FrozenBitset([], capacity=length), ~FrozenBitset([], capacity=length), FrozenBitset([], capacity=length))

    @classmethod
    def star(cls, length: int) -> PartialSignVector:
        r"""
        Create a star partial sign vector of a given length.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: PartialSignVector.star(5)
            (*****)
        """
        return cls(~FrozenBitset([], capacity=length), ~FrozenBitset([], capacity=length), ~FrozenBitset([], capacity=length))

    @classmethod
    def from_string(cls, s: str) -> PartialSignVector:
        r"""
        Create a partial sign vector from a string.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: PartialSignVector.from_string("+*0+/")
            (+*0+/)
        """
        return cls.from_support(
            (pos for pos, t in enumerate(s) if t in ["-", "n","/","*"]),
            (pos for pos, t in enumerate(s) if t in ["0", "n","p","*"]),
            (pos for pos, t in enumerate(s) if t in ["+", "p","/","*"]),
            len(s)
        )

    @classmethod
    def from_iterable(cls, iterable) -> PartialSignVector:
        r"""
        Create a partial sign vector from an iterable.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: PartialSignVector.from_iterable([[-1, 0, 1], [0], [-1], [0, 1]])
            (*0-p)
        """
        negative_support = []
        zero_support = []
        positive_support = []
        length = 0
        for s in iterable:
            if -1 in s:
                negative_support.append(length)
            if 0 in s:
                zero_support.append(length)
            if 1 in s:
                positive_support.append(length)
            length += 1
        return cls.from_support(negative_support, zero_support, positive_support, length)

    @classmethod
    def from_support(cls, negative_support: list[int], zero_support: list[int], positive_support: list[int], length: int) -> PartialSignVector:
        r"""
        Return a partial sign vector that is given by lists representing negative  support, zero  support and positive support.

        INPUT:

        - ``negative_support`` -- a list of integers.
        - ``zero_support`` -- a list of integers.
        - ``positive_support`` -- a list of integers.
        - ``length`` -- a nonnegative integer

        OUTPUT:
        a partial sign vector

        .. NOTE::

            Every index from 0 to length-1 should occur at least once in one of the lists ``negative_support``, ``zero_support`` and  ``positive_support``.
            For efficiency, this is not checked.

        EXAMPLES::

            sage: from sign_vectors.partial_sign_vectors import *
            sage: PartialSignVector.from_support([1, 2, 4], [2, 3, 5], [0, 2, 5], 6)
            (+-*0-p)
        """
        return cls(FrozenBitset(negative_support, capacity=length), FrozenBitset(zero_support, capacity=length), FrozenBitset(positive_support, capacity=length))

    @classmethod
    def from_sign_vector(cls, sv: SignVector) -> PartialSignVector:
        r"""
        Return a partial sign vector that is give by a sign vector.

        INPUT:

        - ``sv`` -- a sign vector.

        .. NOTE::

            The list ``positive_support`` and ``negative_support`` should be disjoint.
            For efficiency, this is not checked.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = SignVector.from_support([1, 4], [2], 6)
            sage: PartialSignVector.from_sign_vector(X)
            (0+-0+0)
        """
        return cls(sv._negative_support, sv._zero_support, sv._positive_support)

class  ExtendedSignVector(SignVector):
    """Utility class for sign vectors with extra functionality."""
    def _determined_positive_support(self) -> FrozenBitset:
        return self._positive_support

    def _determined_negative_support(self) -> FrozenBitset:
        return self._negative_support

    def lower_closure(self) -> PartialSignVector:
        r"""
        Return the lower closure.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = ExtendedSignVector.from_sign_vector(sign_vector("+-00+-"))
            sage: X.lower_closure()
            (pn00pn)
        """
        return PartialSignVector(self._negative_support,
                               ~FrozenBitset([], capacity=self.length()),
                               self._positive_support)

    def upper_closure(self) -> PartialSignVector:
        r"""
        Return the upper closure.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = ExtendedSignVector.from_sign_vector(sign_vector("+-00+-"))
            sage: X.upper_closure()
            (+-**+-)
        """
        return PartialSignVector(self._negative_support | self._zero_support,
                               self._zero_support,
                               self._positive_support | self._zero_support)

    def closure(self) -> PartialSignVector:
       r"""
        Return the closure.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = ExtendedSignVector.from_sign_vector(sign_vector("+-00+-"))
            sage: X.closure()
            (pn**pn)
        """
       return PartialSignVector(self._negative_support | self._zero_support,
                               ~FrozenBitset([], capacity=self.length()),
                               self._positive_support | self._zero_support)

    def orthogonal_complement(self, other: PartialSignVector = None) -> list[PartialSignVector]:
        r"""
        Compute the orthogonal complement of the sign vector in a partial sign vector.

        If no argument is given, the complete orthogonal complement of the sign vector is returned.

        INPUT:

        - ``other`` -- partial sign vector (optional)

        OUTPUT:
        List of partial sign vectors.

        .. NOTE::

            The partial sign vector ``other`` should only contain the signs ``-``, ``0``, ``+`` and ``*``.
            For efficiency, this is not checked.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: from sign_vectors.partial_sign_vectors import *
            sage: X = ExtendedSignVector.from_sign_vector(sign_vector("++00-"))
            sage: X
            (++00-)
            sage: X.orthogonal_complement()
            [(00**0), (+-***), (+***+), (-+***), (*+**+), (-***-), (*-**-)]
            sage: Y = partial_sign_vector("**-+*")
            sage: Y
            (**-+*)
            sage: X.orthogonal_complement(Y)
            [(00-+0), (+--+*), (+*-++), (-+-+*), (*+-++), (-*-+-), (*--+-)]
            sage: Y = partial_sign_vector("****-")
            sage: Y
            (****-)
            sage: X.orthogonal_complement(Y)
            [(-***-), (*-**-)]
            sage: Y = partial_sign_vector("*+-+-")
            sage: Y
            (*+-+-)
            sage: X.orthogonal_complement(Y)
            [(-+-+-)]
        """
        if other is None:
            other = PartialSignVector.star(self.length())
        res = []
        if other._connecting_elements(self).isempty() and other._separating_elements(self).isempty():
            res.append(PartialSignVector(other._negative_support - self._support(),
                                         other._zero_support,
                                         other._positive_support - self._support()))

            for c, d in permutations(list(self._support() & other._undetermined_support()), 2):
                Z = other.set_sign([c],[self[c]])
                Z = Z.set_sign([d],[-self[d]])
                res.append(Z)

        elif not (other._connecting_elements(self).isempty() or other._separating_elements(self).isempty()):
            return [other]

        elif not other._connecting_elements(self).isempty():
            for d in list(self._support() & other._undetermined_support()):
                res.append(other.set_sign([d],[-self[d]]))

        elif not other._separating_elements(self).isempty():
            for c in list(self._support() & other._undetermined_support()):
                res.append(other.set_sign([c],[self[c]]))

        return res

    @classmethod
    def from_sign_vector(cls, sign_vector: SignVector) -> ExtendedSignVector:
        """Create an extended sign vector from a sign vector."""
        return cls(sign_vector._positive_support, sign_vector._negative_support)


def prune(iterable: set[PartialSignVector]) -> set[PartialSignVector]:
    r"""
    Remove all partial sign vectors that are subsets of other partial sign vectors in the iterable.

    INPUT:

    - ``iterable`` -- an iterable of partial sign vectors

    OUTPUT:
    Return a pruned set of partial sign vectors.

    EXAMPLES::

        sage: from sign_vectors.partial_sign_vectors import *
        sage: W = [partial_sign_vector("+**"), partial_sign_vector("-*+"), partial_sign_vector("-++"), partial_sign_vector("+-*")]
        sage: W
        [(+**), (-*+), (-++), (+-*)]
        sage: prune(W)
        {(-*+), (+**)}
        sage: prune([partial_sign_vector("000")])
        set()
    """
    if not iterable:
        return set()

    max_cardinality_length = 1
    same_cardinality_list = [set()]

    for x in iterable:
        while x.cardinality() >= max_cardinality_length:
            same_cardinality_list.append(set())
            max_cardinality_length += 1

        same_cardinality_list[x.cardinality() - 1].add(x)
    res = set()

    for i in range(max_cardinality_length - 1, 0, -1):
        for x in same_cardinality_list[i]:
            subset_flag = False
            for y in res:
                if x.issubset(y):
                    subset_flag = True
                    break
            if not subset_flag:
                res.add(x)

    return res
