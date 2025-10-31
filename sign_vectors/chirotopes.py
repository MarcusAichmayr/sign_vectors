r"""Chirotopes: An auxiliary module for oriented matroids."""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from bisect import bisect_left
from enum import IntEnum
from typing import Iterator

from sage.combinat.combination import Combinations

from . import sign_symbolic, SignVector


class Sign(IntEnum):
    r"""
    Class for chirotope entries.

    EXAMPLES::

        sage: from sign_vectors.chirotopes import Sign
        sage: Sign(1)
        +
        sage: Sign(-1)
        -
        sage: Sign(0)
        0
        sage: Sign(5)
        +
        sage: Sign(5).value
        1
        sage: -Sign(5)
        -
        sage: Sign("+")
        +
        sage: Sign("-")
        -
        sage: Sign("0")
        0
    """
    NEG = -1
    ZERO = 0
    POS = 1

    def __str__(self):
        if self.value > 0:
            return "+"
        if self.value < 0:
            return "-"
        return "0"

    def __repr__(self):
        return str(self)

    def __neg__(self):
        return Sign(-self.value) if self.value != 0 else Sign.ZERO

    @classmethod
    def _missing_(cls, value):
        if isinstance(value, str):
            if value == "+":
                return cls.POS
            if value == "-":
                return cls.NEG
            return cls.ZERO
        v = sign_symbolic(value)
        if v > 0:
            return cls.POS
        if v < 0:
            return cls.NEG
        return cls.ZERO


class Chirotope:
    r"""
    A chirotope of given rank and ground set size.

    EXAMPLES::

        sage: from sign_vectors import *
        sage: from sign_vectors.chirotopes import *
        sage: c = Chirotope.from_list([0, 0, -1, 1, 0, 1, -1, 1, -1, -1], 2, 5)
        sage: c
        Chirotope of rank 2 on ground set of size 5
        sage: c.entry((0, 3))
        -
        sage: c.entries()
        [0, 0, -, +, 0, +, -, +, -, -]
        sage: c.as_string()
        '00-+0+-+--'
        sage: c.oriented_matroid()
        Oriented matroid of dimension 1 with elements of size 5.
        sage: c.dual()
        Chirotope of rank 3 on ground set of size 5
        sage: c.dual().entries()
        [-, +, +, -, -, 0, -, -, 0, 0]

    We construct chirotopes from matrices::

        sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
        sage: c = Chirotope.from_matrix(M)
        sage: c.entry((0, 1))
        +
        sage: c.entries()
        [+, +, +, +, +, 0]
        sage: c.dual()
        Chirotope of rank 2 on ground set of size 4
        sage: c.dual().entries()
        [0, -, +, +, -, +]

    ::

        sage: M = matrix([[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1]])
        sage: c = Chirotope.from_matrix(M)
        sage: c
        Chirotope of rank 3 on ground set of size 4
        sage: c.entries()
        [+, +, -, +]
        sage: c.dual()
        Chirotope of rank 1 on ground set of size 4
        sage: c.dual().entries()
        [+, +, +, -]

    We construct chirotopes from circuits::

        sage: c = Chirotope.from_circuits([sign_vector("00+0"), sign_vector("000+")], 2, 4)
        sage: c.entry([0, 1])
        +
        sage: c.entry([0, 2])
        0
        sage: c.entry([1, 2])
        0
        sage: c = Chirotope.from_circuits([sign_vector("00+0"), sign_vector("000+")], 2, 4)
        sage: c.entries()
        [+, 0, 0, 0, 0, 0]
        sage: c.dual()
        Chirotope of rank 2 on ground set of size 4
        sage: c.dual().entries()
        [0, 0, 0, 0, 0, +]
        sage: c.as_string()
        '+00000'

    ::

        sage: c = Chirotope.from_circuits([sign_vector("++0--"), sign_vector("+0--0"), sign_vector("-0++0"), sign_vector("0--0+"), sign_vector("--0++"), sign_vector("0++0-")], 3, 5)
        sage: c.entries()
        [+, -, +, 0, -, +, +, 0, -, +]
        sage: c.dual().entries()
        [+, +, 0, -, +, +, 0, +, +, +]

    We construct chirotopes from cocircuits::

        sage: c = Chirotope.from_cocircuits([sign_vector("00+0"), sign_vector("000+")], 2, 4)
        sage: c.entries()
        [0, 0, 0, 0, 0, +]

    ::

        sage: c = Chirotope.from_cocircuits([sign_vector("00++"), sign_vector("0++0"), sign_vector("0+0-")], 2, 4)
        sage: c.entries()
        [0, 0, 0, +, +, +]
    """
    def __init__(self, rank: int, ground_set_size: int) -> None:
        self.rank = rank
        self.ground_set_size = ground_set_size
        self._chirotope_dict: dict[tuple[int], Sign] = {}

    def __repr__(self) -> str:
        r"""Return a string representation of the chirotope."""
        return f"Chirotope of rank {self.rank} on ground set of size {self.ground_set_size}"

    def as_string(self) -> str:
        r"""Represent the chirotope as a string."""
        return "".join(str(self.entry(tuple(rset))) for rset in Combinations(self.ground_set_size, self.rank))

    def entry(self, rset: list[int]) -> Sign:
        r"""Return the chirotope entry given by ``indices``."""
        if not isinstance(rset, tuple):
            rset = tuple(rset)
        if not self._has_entry(rset):
            self._set_entry(rset, self._compute_entry(rset))
        return self._chirotope_dict[tuple(rset)]

    def _has_entry(self, rset: list[int]) -> bool:
        return tuple(rset) in self._chirotope_dict

    def _set_entry(self, rset: tuple[int], value: Sign) -> None:
        self._chirotope_dict[rset] = value

    def _compute_entry(self, rset: list[int]) -> Sign:
        r"""Compute the chirotope entry given by ``indices``."""
        raise NotImplementedError

    def entries(self) -> list[Sign]:
        r"""Return all chirotope entries in lexicographic order."""
        return [self.entry(tuple(rset)) for rset in Combinations(self.ground_set_size, self.rank)]

    def _set_entries(self) -> None:
        for rset in Combinations(self.ground_set_size, self.rank):
            self.entry(rset)

    def __eq__(self, other: object) -> bool:
        r"""
        Return whether two chirotopes are equal.

        EXAMPLES::

            sage: from sign_vectors.chirotopes import *
            sage: c = Chirotope.from_matrix(random_matrix(ZZ, 3, 5))
            sage: d = Chirotope.from_list(c.entries(), c.rank, c.ground_set_size)
            sage: c == d
            True
        """
        if not isinstance(other, Chirotope):
            return NotImplemented
        return self.as_string() == other.as_string()

    def dual(self) -> "Chirotope":
        r"""Return the dual chirotope."""
        self._set_entries()
        dual_chirotope = Chirotope(self.ground_set_size - self.rank, self.ground_set_size)
        for rset, value in self._chirotope_dict.items():
            complement = tuple(i for i in range(self.ground_set_size) if i not in rset)
            inversions = sum(i < j for i in rset for j in complement)
            dual_chirotope._set_entry(complement, Sign(-value if inversions & 1 else value)) # check last bit
        return dual_chirotope

    def oriented_matroid(self) -> "OrientedMatroid":
        r"""Return the oriented matroid corresponding to the chirotope."""
        from sign_vectors import OrientedMatroid  # avoid circular import
        return OrientedMatroid.from_chirotope_class(self)

    @staticmethod
    def from_list(entries: list[Sign], rank: int, ground_set_size: int) -> "Chirotope":
        r"""Construct a chirotope from its entries."""
        return _ChirotopeFromEntries(entries, rank, ground_set_size)

    @staticmethod
    def from_matrix(matrix) -> "Chirotope":
        r"""Construct a chirotope from a matrix."""
        return _ChirotopeFromMatrix(matrix)

    @staticmethod
    def from_circuits(circuits: set[SignVector], rank: int, ground_set_size: int) -> "Chirotope":
        r"""Construct a chirotope from its circuits."""
        return _ChirotopeFromCircuits(circuits, rank, ground_set_size)

    @staticmethod
    def from_cocircuits(cocircuits: set[SignVector], rank: int, ground_set_size: int) -> "Chirotope":
        r"""Construct a chirotope from its cocircuits."""
        return _ChirotopeFromCocircuits(cocircuits, rank, ground_set_size)


class _ChirotopeFromEntries(Chirotope):
    def __init__(self, entries: list[Sign], rank: int, ground_set_size: int) -> None:
        super().__init__(rank, ground_set_size)
        for rset, value in zip(Combinations(ground_set_size, rank), entries):
            self._set_entry(tuple(rset), Sign(value))

    def _compute_entry(self, rset: list[int]) -> Sign:
        raise NotImplementedError("All entries are already set.")


class _ChirotopeFromMatrix(Chirotope):
    r"""
    A chirotope constructed from a matrix.

    .. NOTE::

        Pass a matrix with full rank.
        Otherwise, all entries of the chirotope will be zero.

    TESTS::

        sage: from sign_vectors.chirotopes import Chirotope
        sage: M = matrix(ZZ, 0, 3)
        sage: c = Chirotope.from_matrix(M)
        sage: c.entries()
        [+]

    ::

        sage: M = zero_matrix(ZZ, 1, 2)
        sage: c = Chirotope.from_matrix(M)
        sage: c.entries()
        [0, 0]
        sage: c.dual().entries()
        [0, 0]
    """
    def __init__(self, matrix) -> None:
        super().__init__(matrix.nrows(), matrix.ncols())
        self.matrix = matrix

    def _compute_entry(self, rset: list[int]) -> Sign:
        return Sign(self.matrix.matrix_from_columns(rset).det())


class _ChirotopeFromMinimalSupportElements(Chirotope):
    def __init__(self, faces: set[SignVector], rank: int, ground_set_size: int) -> None:
        super().__init__(rank, ground_set_size)
        self._faces_dict: dict[tuple[int], SignVector] = {}
        self._set_faces_dict(faces)

    def _set_faces_dict(self, faces: set[SignVector]) -> None:
        r"""
        Set up a dictionary mapping index sets to faces.

        .. NOTE::

        TESTS::

            sage: from sign_vectors import *
            sage: from sign_vectors.chirotopes import *
            sage: c = Chirotope.from_circuits([sign_vector("00+0"), sign_vector("000+")], 2, 4)
            sage: c._faces_dict
            {(0, 1, 2): (00+0), (0, 1, 3): (000+), (0, 2, 3): (000+), (1, 2, 3): (000+)}

        Note that the corresponding chirotopes are ``+00000``.
        Hence, ``(0, 2, 3)`` and ``(1, 2, 3)`` produce the zero face contrary to what is stored.
        However, this does not matter for the algorithm since these index sets connect two zero minors.

        We do the same for the cocircuits::

            sage: c = Chirotope.from_cocircuits([sign_vector("00+0"), sign_vector("000+")], 2, 4)
            sage: c._faces_dict
            {(0,): (000+), (1,): (000+), (2,): (000+), (3,): (00+0)}
        """
        checked_supports = set()
        for face in faces:
            support = frozenset(face.support())
            if support in checked_supports:
                continue
            checked_supports.add(support)
            for indices in self._corresponding_indices_for_face(support):
                self._faces_dict[indices] = face

    def _compute_entry(self, rset: list[int]) -> Sign:
        self._set_entries()
        return super().entry(rset)

    def _set_entries(self) -> None:
        self._set_other_nonzero_entries_from(self._set_zero_entries_and_first_nonzero())

    def _set_zero_entries_and_first_nonzero(self) -> tuple[int]:
        r"""
        Set all zero entries of the chirotope using (co)circuits and set the first nonzero entry to ``+``.
        """
        nonzero_rset = None
        for rset in Combinations(self.ground_set_size, self.rank):
            if self._is_entry_zero(rset):
                self._set_entry(tuple(rset), Sign.ZERO)
            elif nonzero_rset is None:
                nonzero_rset = tuple(rset)
                self._set_entry(nonzero_rset, Sign.POS)
        return nonzero_rset

    def _set_other_nonzero_entries_from(self, rset: tuple[int]) -> None:
        r"""
        Set all missing nonzero entries of the chirotope using the Grassmann-Plücker relations.
        """
        pending_entries: dict[tuple[int], tuple[int]] = {}
        while True:
            for adjacent_rset in self._get_adjacent_rsets(rset):
                if self._has_entry(adjacent_rset):
                    continue
                if adjacent_rset not in pending_entries:
                    pending_entries[adjacent_rset] = rset
            if not pending_entries:
                return

            adjacent_rset, rset = pending_entries.popitem()
            face_indices = self._connecting_face_indices(rset, adjacent_rset)
            face = self._faces_dict[face_indices]
            i, j = set(rset).symmetric_difference(adjacent_rset)
            adjacent_value = Sign(self.entry(rset) * face[i] * face[j])
            # set sign depending on positions of i and j in face_indices
            if (bisect_left(face_indices, i) + bisect_left(face_indices, j)) & 1:
                adjacent_value = -adjacent_value
            self._set_entry(adjacent_rset, adjacent_value)
            rset = adjacent_rset

    def _get_adjacent_rsets(self, rset: tuple[int]) -> Iterator[tuple[int]]:
        if self.rank == 0:
            return
        complement = [i for i in range(self.ground_set_size) if i not in rset]
        for subset in Combinations(rset, self.rank - 1):
            for i in complement:
                yield tuple(sorted(subset + [i]))

    def _corresponding_indices_for_face(self, support: set[int]) -> Iterator[tuple[int]]:
        r"""
        Generate all possible index sets that correspond to the face given by its support.

        For instance, ``{0, 1, 2}`` corresponds to the circuit ``(+++0)``.
        Also ``{0, 1, 2}`` would also correspond to the circuit ``(++00)``.

        On the other hand, ``{3}`` would correspond to the cocircuit ``(+++0)``.
        """
        raise NotImplementedError

    def _connecting_face_indices(self, rset1: tuple[int], rset2: tuple[int]) -> tuple[int]:
        raise NotImplementedError

    def _is_entry_zero(self, rset: list[int]) -> bool:
        r"""Return whether the entry of the chirotope given by ``rset`` is zero using (co)circuits."""
        raise NotImplementedError


class _ChirotopeFromCircuits(_ChirotopeFromMinimalSupportElements):
    r"""
    A chirotope constructed from its circuits.

    TESTS:

    When iterating over ``{7, 8}``, we obtain ``8`` first.
    Therefore, we need to sort indices::

        sage: from sign_vectors import *
        sage: from sign_vectors.chirotopes import Chirotope
        sage: M = matrix.ones(1, 9)
        sage: om = OrientedMatroid(M)
        sage: Chirotope.from_circuits(om.circuits(), 1, 9).entries()
        [+, +, +, +, +, +, +, +, +]
    """
    def _corresponding_indices_for_face(self, support: set[int]) -> Iterator[tuple[int]]:
        complement = set(range(self.ground_set_size)) - support
        for indices in Combinations(complement, self.rank + 1 - len(support)):
            yield tuple(sorted(support.union(indices)))

    def _connecting_face_indices(self, rset1: tuple[int], rset2: tuple[int]) -> tuple[int]:
        return tuple(sorted(set(rset1).union(rset2)))

    def _is_entry_zero(self, rset: list[int]) -> bool:
        for i in range(self.ground_set_size):
            if i in rset:
                continue
            face = self._faces_dict[tuple(sorted(rset + [i]))]
            if set(face.support()).issubset(rset):
                return True
        return False


class _ChirotopeFromCocircuits(_ChirotopeFromMinimalSupportElements):
    r"""A chirotope constructed from its cocircuits."""
    def _corresponding_indices_for_face(self, support: set[int]) -> Iterator[tuple[int]]:
        complement = set(range(self.ground_set_size)) - support
        for indices in Combinations(complement, self.rank - 1):
            yield tuple(sorted(indices))

    def _connecting_face_indices(self, rset1: tuple[int], rset2: tuple[int]) -> tuple[int]:
        return tuple(sorted(set(rset1).intersection(rset2)))

    def _is_entry_zero(self, rset: list[int]) -> bool:
        for i in rset:
            face = self._faces_dict[tuple(sorted(set(rset) - {i}))]
            if set(rset).issubset(face.zero_support()):
                return True
        return False
