r"""
Sign vectors
============

First, we load the functions from the package::

    sage: from sign_vectors import *

There are several ways to define sign vectors::

    sage: sign_vector([1, 0, -1, 1])
    (+0-+)
    sage: x = vector([5, 2/5, -1, 0])
    sage: sign_vector(x)
    (++-0)
    sage: sign_vector('++-+-00-')
    (++-+-00-)

We construct the zero sign vector of a given length::

    sage: zero_sign_vector(6)
    (000000)

It is also possible to generate random sign vectors::

    sage: random_sign_vector(7) # random
    (+-+00+-)

There are different notions of support::

    sage: X = sign_vector('+++000--')
    sage: X.support()
    [0, 1, 2, 6, 7]
    sage: X.zero_support()
    [3, 4, 5]
    sage: X.positive_support()
    [0, 1, 2]
    sage: X.negative_support()
    [6, 7]

Next, we define two sign vectors::

    sage: X = sign_vector([-1, 0, 1, -1, 0])
    sage: X
    (-0+-0)
    sage: Y = sign_vector([0, 1, 0, 1, 0])
    sage: Y
    (0+0+0)

We compose them::

    sage: X.compose(Y)
    (-++-0)
    sage: Y.compose(X)
    (-+++0)

Use the operator ``&`` as a shorthand for composition::

    sage: X & Y
    (-++-0)
    sage: Y & X
    (-+++0)

The conformal relation is a partial order on a set of sign vectors::

    sage: X = sign_vector([-1, 1, 0, 0, 1])
    sage: Y = sign_vector([-1, 1, 1, 0, 1])
    sage: Z = sign_vector([-1, 1, 1, -1, 1])
    sage: X
    (-+00+)
    sage: Y
    (-++0+)
    sage: Z
    (-++-+)
    sage: X.conforms(Y)
    True
    sage: X.conforms(X)
    True
    sage: Y.conforms(Z)
    True
    sage: Z.conforms(Y)
    False
    sage: X.conforms(Z)
    True
    sage: X < X
    False
    sage: X <= Y
    True
    sage: Y > Z
    False
    sage: Z >= X
    True

Similar as for real vectors, we define orthogonality for sign vectors.
First, we define some real vectors::

    sage: x = vector([0, 0, 1])
    sage: y = vector([1, -2, 0])
    sage: z = vector([2, 1, 2])

We compute some scalar products to investigate orthogonality of those vectors::

    sage: x.dot_product(y)
    0
    sage: y.dot_product(z)
    0
    sage: x.dot_product(z)
    2

Next, we define the corresponding sign vectors::

    sage: X = sign_vector(x)
    sage: X
    (00+)
    sage: Y = sign_vector(y)
    sage: Y
    (+-0)
    sage: Z = sign_vector(z)
    sage: Z
    (+++)

By definition, if two real vectors are orthogonal, then their corresponding
sign vectors are also orthogonal::

    sage: X.is_orthogonal_to(Y)
    True
    sage: Y.is_orthogonal_to(Z)
    True
    sage: X.is_orthogonal_to(Z)
    False
"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from __future__ import annotations
import warnings
from random import choices

from sage.data_structures.bitset import FrozenBitset
from sage.structure.sage_object import SageObject
from sage.symbolic.ring import SR


def sign_vector(iterable: list[int] | str | SignVector) -> SignVector:
    r"""
    Create a sign vector from a list, vector or string.

    INPUT:

    - ``iterable`` -- different inputs are accepted:

        - an iterable (e.g. a list or vector) of real values.
          Variables can also occur.

        - a string consisting of ``"-"``, ``"+"``, ``"0"``. Other characters are treated as ``"0"``.

    OUTPUT:

    Returns a sign vector.

    EXAMPLES::

        sage: from sign_vectors import *
        sage: sign_vector([5, 0, -1, -2])
        (+0--)
        sage: v = vector([5, 0, -1, -2])
        sage: sign_vector(v)
        (+0--)

    We can also use a string to construct a sign vector::

        sage: sign_vector("00--")
        (00--)
        sage: sign_vector('++-+-00-')
        (++-+-00-)

    Variables are supported to some extent::

        sage: var('a')
        a
        sage: v = vector([1, a, -1])
        sage: sign_vector(v)
        ...
        UserWarning: Cannot determine sign of symbolic expression, using ``0`` instead.
        (+0-)
        sage: assume(a > 0)
        sage: sign_vector(v)
        (++-)
        sage: forget()

    TESTS::

        sage: sign_vector(sign_vector("+-+"))
        (+-+)
    """
    if isinstance(iterable, SignVector):
        return iterable
    if isinstance(iterable, str):
        return SignVector.from_string(iterable)
    return SignVector.from_iterable(iterable)


def zero_sign_vector(length: int) -> SignVector:
    r"""
    Return the zero sign vector of a given length.

    INPUT:

    - ``length`` -- length

    EXAMPLES::

        sage: from sign_vectors import *
        sage: zero_sign_vector(4)
        (0000)
    """
    return SignVector.zero(length)


def random_sign_vector(length: int) -> SignVector:
    r"""
    Return a random sign vector of a given length.

    INPUT:

    - ``length`` -- length

    EXAMPLES::

        sage: from sign_vectors import *
        sage: random_sign_vector(5) # random
        (++-0-)

    TEST::

        sage: len(random_sign_vector(5))
        5
    """
    return SignVector.from_string("".join(choices("00+-", k=length)))


def sign_symbolic(value) -> int:
    r"""
    Return the sign of an expression. Supports symbolic expressions.

    OUTPUT:
    If the sign cannot be determined, a warning is shown and ``0`` is returned.

    EXAMPLES:

    For real numbers, the sign is determined::

        sage: from sign_vectors.sign_vectors import sign_symbolic
        sage: sign_symbolic(1)
        1
        sage: sign_symbolic(-1/2)
        -1
        sage: sign_symbolic(0)
        0

    For symbolic expressions, the sign is determined using assumptions::

        sage: var('a')
        a
        sage: sign_symbolic(a) # not tested
        ...
        UserWarning: Cannot determine sign of symbolic expression, using ``0`` instead.
        0
        sage: assume(a > 0)
        sage: sign_symbolic(a)
        1
        sage: forget()

    ::

        sage: var('a, b, c, d')
        (a, b, c, d)
        sage: assume(a > 0, b > 0, c > 0, d > 0)
        sage: sign_symbolic((a + b)*(c + d) - b*c)
        1
    """
    if value == 0:
        return 0
    if value > 0:
        return 1
    if value < 0:
        return -1

    expr = SR(value).simplify_full()
    if expr == 0:
        return 0
    if expr > 0:
        return 1
    if expr < 0:
        return -1

    warnings.warn("Cannot determine sign of symbolic expression, using ``0`` instead.")
    return 0


class SignVector(SageObject):
    r"""A sign vector."""

    __slots__ = ("_positive_support", "_negative_support")

    def __init__(self, positive_support: FrozenBitset, negative_support: FrozenBitset) -> None:
        r"""
        Create a sign vector object.

        INPUT:

        - ``positive_support`` -- a ``FrozenBitset`` representing the positive support.
        - ``negative_support`` -- a ``FrozenBitset`` representing the negative support.

        .. NOTE::

            The ``positive_support`` and ``negative_support`` should be disjoined.
            For efficiency, this is not checked when creating a sign vector object.

        .. SEEALSO::

            - :func: `~sign_vector`

        TESTS::

            sage: from sign_vectors import *
            sage: SignVector(FrozenBitset([1, 3], capacity=4), FrozenBitset([0], capacity=4))
            (-+0+)
        """
        self._positive_support = positive_support
        self._negative_support = negative_support

    def _repr_(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return f"({self.to_string()})"

    def to_string(self) -> str:
        r"""
        Return a string representation of this sign vector (without parentheses).

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector("0+-")
            sage: X
            (0+-)
            sage: X.to_string()
            '0+-'

        Note that that `str` and `to_string` are different::

            sage: str(X)
            '(0+-)'
        """
        return "".join(
            "+" if e in self._positive_support else ("-" if e in self._negative_support else "0")
            for e in range(self.length())
        )

    def length(self) -> int:
        r"""
        Return the length of the sign vector.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('0+-'); X
            (0+-)
            sage: X.length()
            3
        """
        return self._positive_support.capacity()

    def __len__(self) -> int:
        r"""
        Return the length of this sign vector.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('0+-'); X
            (0+-)
            sage: len(X)
            3
        """
        return self.length()

    def support(self) -> list[int]:
        r"""
        Return a list of indices where the sign vector is nonzero.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-0+-0'); X
            (-0+-0)
            sage: X.support()
            [0, 2, 3]
        """
        return list(self._support())

    def _support(self) -> FrozenBitset:
        return self._positive_support | self._negative_support

    def zero_support(self) -> list[int]:
        r"""
        Return a list of indices where the sign vector is zero.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-0+-0'); X
            (-0+-0)
            sage: X.zero_support()
            [1, 4]
        """
        return list(~self._support())

    def positive_support(self) -> list[int]:
        r"""
        Return a list of indices where the sign vector is positive.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-0+-0'); X
            (-0+-0)
            sage: X.positive_support()
            [2]
        """
        return list(self._positive_support)

    def negative_support(self) -> list[int]:
        r"""
        Return a list of indices where the sign vector is negative.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-0+-0'); X
            (-0+-0)
            sage: X.negative_support()
            [0, 3]
        """
        return list(self._negative_support)

    def __getitem__(self, e: int | slice) -> int | SignVector:
        r"""
        Return the element at position ``e`` of the sign vector.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('0++-'); X
            (0++-)
            sage: X[0]
            0
            sage: X[1]
            1
            sage: X[3]
            -1
            sage: X[-1]
            -1
            sage: X[1:3]
            (++)

        TESTS::

            sage: X[-2]
            1
            sage: X[100]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if isinstance(e, slice):
            # TODO improve
            return SignVector.from_string(self.to_string()[e])
        if e >= self.length() or e < -self.length():
            raise IndexError("index out of range")
        if e < 0:
            e %= self.length()
        if e in self._positive_support:
            return 1
        if e in self._negative_support:
            return -1
        return 0

    def list_from_positions(self, positions: list[int]) -> list[int]:
        r"""
        Return a list of components that are in the list of indices ``positions``.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-+00+'); X
            (-+00+)
            sage: X.list_from_positions([0, 1, 4])
            [-1, 1, 1]
        """
        return [self[e] for e in positions]

    def compose(self, other: SignVector) -> SignVector:
        r"""
        Return the composition of two sign vectors.

        INPUT:

        - ``other`` -- a sign vector

        OUTPUT:

        Composition of this sign vector with ``other``.

        .. NOTE::

            Alternatively, the operator ``&`` can be used.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('+00'); X
            (+00)
            sage: Y = sign_vector('--0'); Y
            (--0)
            sage: X.compose(Y)
            (+-0)
            sage: Y.compose(X)
            (--0)
            sage: Y & X
            (--0)
            sage: X = sign_vector('000+++---')
            sage: Y = sign_vector('0+-0+-0+-')
            sage: X.compose(Y)
            (0+-+++---)
        """
        return self.__class__(
            self._positive_support | (other._positive_support & ~self._negative_support),
            self._negative_support | (other._negative_support & ~self._positive_support)
        )

    def compose_harmonious(self, other: SignVector) -> SignVector:
        r"""
        Return the composition of two harmonious sign vectors.

        INPUT:

        - ``other`` -- a sign vector

        OUTPUT:

        Composition of this sign vector with ``other``.

        .. NOTE::

            This method is more efficient than - :meth:`compose`.
            However, it does not check whether the sign vectors are harmonious.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('+00')
            sage: X
            (+00)
            sage: Y = sign_vector('0-0')
            sage: Y
            (0-0)
            sage: X.compose_harmonious(Y)
            (+-0)
            sage: Y.compose_harmonious(X)
            (+-0)
        """
        return self.__class__(
            self._positive_support | other._positive_support,
            self._negative_support | other._negative_support
        )

    def __and__(self, other: SignVector) -> SignVector:
        r"""
        Return the composition of two sign vectors.

        .. SEEALSO::

            - :meth: `compose`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('+00')
            sage: X
            (+00)
            sage: Y = sign_vector('--0')
            sage: Y
            (--0)
            sage: X & Y
            (+-0)
            sage: Y & X
            (--0)
        """
        return self.compose(other)

    def __mul__(self, value: int) -> SignVector:
        r"""
        Multiplication with a scalar.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-+00+')
            sage: X
            (-+00+)
            sage: -1*X
            (+-00-)
            sage: 1*X
            (-+00+)
        """
        if value > 0:
            return +self
        if value < 0:
            return -self
        return zero_sign_vector(self.length())

    def __rmul__(self, value) -> SignVector:
        r"""
        Right multiplication with a scalar.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-+00+')
            sage: X
            (-+00+)
            sage: X*(-1)
            (+-00-)
            sage: X*1
            (-+00+)
        """
        return self * value

    def flip_signs(self, indices: list[int]) -> SignVector:
        r"""
        Flips entries of given indices.

        INPUT:

        - ``indices`` -- list of indices

        OUTPUT:
        Returns a new sign vector. Components of ``indices`` are multiplied by ``-1``.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-++0+')
            sage: X
            (-++0+)
            sage: X.flip_signs([0, 2, 3])
            (++-0+)

        ::

            sage: X = sign_vector('+-0+0')
            sage: X
            (+-0+0)
            sage: X.flip_signs([0, 3, 4])
            (--0-0)
        """
        indices = FrozenBitset(indices) & self._support()
        return self.__class__(self._positive_support ^ indices, self._negative_support ^ indices)

    def set_to_zero(self, indices: list[int]) -> SignVector:
        r"""
        Set given entries to zero.

        INPUT:

        - ``indices`` -- list of indices

        OUTPUT:
        Returns a new sign vector of same length. Components with indices in
        ``indices`` are set to ``0``.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-++0+')
            sage: X
            (-++0+)
            sage: X.set_to_zero([0, 2, 3])
            (0+00+)

        ::

            sage: X = sign_vector('+-+0+0')
            sage: X
            (+-+0+0)
            sage: X.set_to_zero([0, 1, 4])
            (00+000)
        """
        indices = FrozenBitset(indices)
        return self.__class__(self._positive_support - indices, self._negative_support - indices)

    def set_to_plus(self, indices: list[int]) -> SignVector:
        r"""
        Set given entries to plus.

        INPUT:

        - ``indices`` -- list of indices

        OUTPUT:
        Returns a new sign vector of same length. Components with indices in
        ``indices`` are set to ``+``.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-++0+')
            sage: X
            (-++0+)
            sage: X.set_to_plus([0, 2, 3])
            (+++++)

        ::

            sage: X = sign_vector('+-+0+0')
            sage: X
            (+-+0+0)
            sage: X.set_to_plus([0, 1, 4])
            (+++0+0)
        """
        indices = FrozenBitset(indices)
        return self.__class__(self._positive_support | indices, self._negative_support - indices)

    def set_to_minus(self, indices: list[int]) -> SignVector:
        r"""
        Set given entries to minus.

        INPUT:

        - ``indices`` -- list of indices

        OUTPUT:
        Returns a new sign vector of same length. Components with indices in
        ``indices`` are set to ``-``.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-++0+')
            sage: X
            (-++0+)
            sage: X.set_to_minus([0, 2, 3])
            (-+--+)

        ::

            sage: X = sign_vector('+-+0+0')
            sage: X
            (+-+0+0)
            sage: X.set_to_minus([0, 1, 4])
            (--+0-0)
        """
        indices = FrozenBitset(indices)
        return self.__class__(self._positive_support - indices, self._negative_support | indices)

    def delete_components(self, indices: list[int]) -> SignVector:
        r"""
        Delete the given components from the sign vector.

        INPUT:

        - ``indices`` -- list of indices to delete

        OUTPUT:

        Returns a new sign vector with the specified components removed.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('+0-')
            sage: X
            (+0-)
            sage: X.delete_components([1])
            (+-)

        ::

            sage: X = sign_vector('-+0+0')
            sage: X
            (-+0+0)
            sage: X.delete_components([0, 2])
            (++0)
            sage: X.delete_components([1])
            (-0+0)
        """
        return self.__class__.from_iterable(
            [self[i] for i in ~FrozenBitset(indices, capacity=self.length())]
        )

    def _connecting_elements(self, other: SignVector) -> FrozenBitset:
        return (self._positive_support & other._positive_support) | (self._negative_support & other._negative_support)

    def _separating_elements(self, other: SignVector) -> FrozenBitset:
        return (self._positive_support & other._negative_support) | (self._negative_support & other._positive_support)

    def separating_elements(self, other: SignVector) -> list[int]:
        r"""
        Compute the list of separating elements of two sign vectors.

        INPUT:

        - ``other`` -- sign vector

        OUTPUT:
        List of elements ``e`` such that ``self[e] == -other[e] != 0``.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('++00-')
            sage: X
            (++00-)
            sage: Y = sign_vector('+-+++')
            sage: Y
            (+-+++)
            sage: X.separating_elements(Y)
            [1, 4]
        """
        return list(self._separating_elements(other))

    def is_orthogonal_to(self, other: SignVector) -> bool:
        r"""
        Return whether two sign vectors are orthogonal.

        INPUT:

        - ``other`` -- a sign vector.

        OUTPUT:

        - Returns ``True`` if the sign vectors are orthogonal and ``False`` otherwise.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-+00+'); X
            (-+00+)
            sage: X.is_orthogonal_to(X)
            False
            sage: X.is_orthogonal_to(sign_vector('++000'))
            True
            sage: X.is_orthogonal_to(sign_vector('++00+'))
            True
            sage: X.is_orthogonal_to(sign_vector('00++0'))
            True
        """
        return not (self._separating_elements(other).isempty() ^ self._connecting_elements(other).isempty())

    def is_harmonious_to(self, other: SignVector) -> bool:
        r"""
        Check whether these two sign vectors are harmonious.

        INPUT:

        - ``other`` -- a sign vector or real vector

        OUTPUT:
        Returns true if there are no separating elements.
        Otherwise, false is returned.

        .. NOTE::

            Two sign vectors are harmonious if there is no component
            where one sign vector has ``+`` and the other has ``-``.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('++00-')
            sage: X
            (++00-)
            sage: Y = sign_vector('+-+++')
            sage: Y
            (+-+++)
            sage: X.is_harmonious_to(Y)
            False
            sage: sign_vector('0+00').is_harmonious_to(sign_vector('-+0+'))
            True
            sage: v = vector([1, 2/3, 0, -1, -1])
            sage: X.is_harmonious_to(v)
            True
            sage: Y.is_harmonious_to(v)
            False
        """
        if not isinstance(other, SignVector):
            other = sign_vector(other)

        return self._separating_elements(other).isempty()

    def disjoint_support(self, other: SignVector) -> bool:
        r"""
        Return whether these two sign vectors have disjoint support.

        INPUT:

        - ``other`` -- a sign vector or real vector

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('++00')
            sage: Y = sign_vector('0+0-')
            sage: Z = sign_vector('00--')
            sage: X.disjoint_support(Y)
            False
            sage: Y.disjoint_support(X)
            False
            sage: X.disjoint_support(Z)
            True
            sage: Z.disjoint_support(X)
            True
            sage: Y.disjoint_support(Z)
            False
        """
        if isinstance(other, SignVector):
            return self._support().isdisjoint(other._support())
        return set(self.support()).isdisjoint(other.support())

    def conforms(self, other: SignVector) -> bool:
        r"""
        Conformal relation of two sign vectors.

        .. NOTE::

            Alternatively, the operator ``<=`` can be used.
            Use ``>=``, ``<`` and ``>`` for the other relations.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-+00+'); X
            (-+00+)
            sage: Y = sign_vector('-++0+'); Y
            (-++0+)
            sage: Z = sign_vector('-++-+'); Z
            (-++-+)
            sage: X.conforms(Y)
            True
            sage: X.conforms(X)
            True
            sage: Y.conforms(Z)
            True
            sage: Z.conforms(Y)
            False
            sage: X.conforms(Z)
            True
        """
        return self._positive_support.issubset(other._positive_support) & self._negative_support.issubset(other._negative_support)

    def __eq__(self, other: SignVector) -> bool:
        r"""
        Return whether this sign vector is equal to ``other``.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector("++0-")
            sage: X == X
            True
            sage: X == sign_vector("00++")
            False

        TESTS::

            sage: zero_sign_vector(3) == 0
            True
            sage: 0 == zero_sign_vector(3)
            True
        """
        if isinstance(other, SignVector):
            return (
                self._negative_support == other._negative_support
                and self._positive_support == other._positive_support
            )
        if other == 0:
            return not self._support()
        return False

    def __le__(self, other: SignVector) -> bool:
        r"""
        Return whether this sign vector is less or equal to ``other``.

        .. SEEALSO::

            - :meth: `conforms`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-+00+'); X
            (-+00+)
            sage: Y = sign_vector('-++0+'); Y
            (-++0+)
            sage: X <= Y
            True

        We can also use ``<=`` to compare a sign vector with ``0``::

            sage: sign_vector('00--') <= 0
            True
            sage: sign_vector([1, 1, -1, 0]) <= 0
            False
            sage: 0 <= sign_vector([1, 1, 0, 0])
            True
            sage: zero_sign_vector(2) <= 0
            True
        """
        if isinstance(other, SignVector):
            return self.conforms(other)
        if other == 0:
            return all(a <= 0 for a in self)
        return False

    def __lt__(self, other: SignVector) -> bool:
        r"""
        Return whether this sign vector is less than ``other``.

        .. SEEALSO::

            - :meth: `conforms`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-+00+'); X
            (-+00+)
            sage: Y = sign_vector('-++0+'); Y
            (-++0+)
            sage: X < Y
            True

        We can also use ``<`` to compare a sign vector with ``0``::

            sage: sign_vector('00--') < 0
            True
            sage: sign_vector([1, 1, -1, 0]) < 0
            False
            sage: 0 < sign_vector([1, 1, 0, 0])
            True
            sage: zero_sign_vector(2) < 0
            False
        """
        return self != other and self <= other

    def __ge__(self, other: SignVector) -> bool:
        r"""
        Return whether this sign vector is greater or equal to ``other``.

        .. SEEALSO::

            - :meth: `conforms`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-+00+'); X
            (-+00+)
            sage: Y = sign_vector('-++0+'); Y
            (-++0+)
            sage: Y >= X
            True

        We can also use ``>=`` to compare a sign vector with ``0``::

            sage: sign_vector('00--') >= 0
            False
            sage: sign_vector([1, 1, -1, 0]) >= 0
            False
            sage: sign_vector([1, 1, 0, 0]) >= 0
            True
            sage: zero_sign_vector(2) >= 0
            True
        """
        if isinstance(other, SignVector):
            return other.conforms(self)
        if other == 0:
            return all(a >= 0 for a in self)
        return False

    def __gt__(self, other: SignVector) -> bool:
        r"""
        Return whether this sign vector is greater than ``other``.

        .. SEEALSO::

            - :meth: `conforms`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('-+00+'); X
            (-+00+)
            sage: Y =  sign_vector('-++0+'); Y
            (-++0+)
            sage: Y > X
            True

        We can also use ``>`` to compare a sign vector with ``0``::

            sage: 0 >  sign_vector('00--')
            True
            sage: sign_vector([1, 1, -1, 0]) > 0
            False
            sage: sign_vector([1, 1, 0, 0]) > 0
            True
            sage: zero_sign_vector(2) > 0
            False
        """
        return self != other and self >= other

    def __bool__(self) -> bool:
        return self != 0

    def __neg__(self) -> SignVector:
        r"""
        Return the sign vectors multiplied by ``-1``.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('0+-'); X
            (0+-)
            sage: -X
            (0-+)
        """
        return self.__class__(self._negative_support, self._positive_support)

    def __pos__(self) -> SignVector:
        r"""
        Return the sign vectors multiplied by ``1`` which is a copy of this sign vector.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector('0+-')
            sage: X
            (0+-)
            sage: +X
            (0+-)
        """
        return self.__class__(self._positive_support, self._negative_support)

    def is_vector(self) -> bool:
        r"""Return ``False`` since sign vectors are not vectors."""
        return False

    def __hash__(self) -> int:
        r"""Return the hash value of this sign vector."""

        return hash((self._positive_support, self._negative_support))

    @classmethod
    def zero(cls, length: int) -> SignVector:
        r"""
        Create a zero sign vector of a given length.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: SignVector.zero(5)
            (00000)
        """
        return cls(FrozenBitset([], capacity=length), FrozenBitset([], capacity=length))

    @classmethod
    def from_string(cls, s: str) -> SignVector:
        r"""
        Create a sign vector from a string.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: SignVector.from_string("+-0+0")
            (+-0+0)
        """
        return cls.from_support(
            (pos for pos, t in enumerate(s) if t == "+"),
            (pos for pos, t in enumerate(s) if t == "-"),
            len(s)
        )

    @classmethod
    def from_iterable(cls, iterable) -> SignVector:
        r"""
        Create a sign vector from an iterable.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: SignVector.from_iterable([5, 0, -1, -2])
            (+0--)
            sage: v = vector([5, 0, -1, 0, 8])
            sage: SignVector.from_iterable(v)
            (+0-0+)
        """
        positive_support = []
        negative_support = []
        length = 0
        for entry in iterable:
            sign_entry = sign_symbolic(entry)
            if sign_entry > 0:
                positive_support.append(length)
            elif sign_entry < 0:
                negative_support.append(length)
            length += 1
        return cls.from_support(positive_support, negative_support, length)

    @classmethod
    def from_support(cls, positive_support: list[int], negative_support: list[int], length: int) -> SignVector:
        r"""
        Return a sign vector that is given by lists representing positive support and negative  support.

        INPUT:

        - ``positive_support`` -- a list of integers.
        - ``negative_support`` -- a list of integers.
        - ``length`` -- a nonnegative integer

        OUTPUT:
        a sign vector

        .. NOTE::

            The list ``positive_support`` and ``negative_support`` should be disjoint.
            For efficiency, this is not checked.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: SignVector.from_support([1, 4], [2], 6)
            (0+-0+0)
        """
        return cls(FrozenBitset(positive_support, capacity=length), FrozenBitset(negative_support, capacity=length))
