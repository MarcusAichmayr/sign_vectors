"""Utility functions and other useful functions for working with oriented matroids"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from typing import Iterator

from . import SignVector


def are_parallel(iterable, component1, component2, return_ratio: bool = False):
    r"""
    Determine whether two components of sign vectors or vectors are parallel.

    INPUT:

    - ``iterable`` -- a list of sign vectors or vectors of length ``n``
    - ``component1`` -- an integer with ``0 <= component1 < n``
    - ``component2`` -- an integer with ``0  <= component2 < n``
    - ``return_ratio`` -- a boolean

    OUTPUT:

    Returns a boolean.
    If ``return_ratio`` is true, a tuple of a boolean and the ratio will be returned instead.

    .. NOTE::

        The elements ``component1`` and ``component2`` are parallel if there exists a ratio ``d`` such that
        ``v[component1] = d v[component2]`` for each ``v`` in ``iterable``.

    EXAMPLES::

        sage: from sign_vectors.utility import are_parallel
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("+-0+"), sign_vector("-0+0")]
        sage: L
        [(++0-), (+-0+), (-0+0)]
        sage: are_parallel(L, 0, 1)
        False
        sage: are_parallel(L, 1, 2)
        False
        sage: are_parallel(L, 1, 3)
        True

    Now, we consider some real vectors::

        sage: L = [vector([1, 1, 2, 3, 0, 0]), vector([-2, 1, -4, 3, 3, -17]), vector([0, 1, 0, 1, 0, 0])]
        sage: L
        [(1, 1, 2, 3, 0, 0), (-2, 1, -4, 3, 3, -17), (0, 1, 0, 1, 0, 0)]
        sage: are_parallel(L, 0, 2)
        True
        sage: are_parallel(L, 0, 1)
        False
        sage: are_parallel(L, 1, 3)
        False
        sage: are_parallel(L, 4, 5)
        True

    We can also return the ratio of the two components::

        sage: are_parallel(L, 0, 2, return_ratio=True)
        (True, 1/2)
        sage: are_parallel(L, 2, 0, return_ratio=True)
        (True, 2)
        sage: are_parallel(L, 0, 1, return_ratio=True)
        (False, None)

    Also works for matrices::

        sage: M = matrix([[0, 0, 1, -1, 0], [1, 0, 0, 0, 1], [1, 1, 1, 1, 1]])
        sage: are_parallel(M, 0, 4)
        True
        sage: are_parallel(M, 0, 1)
        False

    TESTS::

        sage: are_parallel([], 0, 1)
        True
        sage: are_parallel([], 0, 1, return_ratio=True)
        (True, 0)
    """
    ratio = 0

    for element in iterable:
        value1 = element[component1]
        value2 = element[component2]
        if ratio == 0:
            if value1 == 0 and value2 == 0:
                continue
            if value1 == 0 and value2 != 0:
                ratio = None
                break
            if value1 != 0 and value2 == 0:
                ratio = None
                break
            ratio = value1 / value2
        elif value1 != ratio * value2:
            ratio = None
            break
    if ratio is None:
        return (False, ratio) if return_ratio else False
    return (True, ratio) if return_ratio else True


def parallel_classes(iterable, length: int) -> list[set[int]]:
    r"""
    Compute the parallel classes of given sign vectors or vectors.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors or vectors with same length
    - ``length`` -- an integer ``n``

    OUTPUT:

    Returns a partition of ``[0, ..., n - 1]`` into parallel classes.

    .. NOTE::

        The elements ``component1`` and ``component2`` are parallel if there exists a ratio ``d`` such that
        ``X[component1] = d X[component2]`` for each ``X`` in ``iterable``.

    EXAMPLES::

        sage: from sign_vectors.utility import parallel_classes
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("+-0+"), sign_vector("-0+0")]
        sage: L
        [(++0-), (+-0+), (-0+0)]
        sage: parallel_classes(L, 4)
        [{0}, {1, 3}, {2}]

    Now, we compute the parallel classes of a list of real vectors::

        sage: L = [vector([1, 1, 2, 3, 0, 0]), vector([-2, 1, -4, 3, 3, -17]), vector([0, 1, 0, 1, 0, 0])]
        sage: L
        [(1, 1, 2, 3, 0, 0), (-2, 1, -4, 3, 3, -17), (0, 1, 0, 1, 0, 0)]
        sage: parallel_classes(L, 6)
        [{0, 2}, {1}, {3}, {4, 5}]

    Let us compute the parallel classes of the rows of a matrix::

        sage: M = matrix([[0, 0, 1, -2, 0], [1, 0, 0, 0, 1], [1, 1, -3, 6, 1]])
        sage: M
        [ 0  0  1 -2  0]
        [ 1  0  0  0  1]
        [ 1  1 -3  6  1]
        sage: parallel_classes(M, 5)
        [{0, 4}, {1}, {2, 3}]

    TESTS::

        sage: parallel_classes([], 5)
        [{0, 1, 2, 3, 4}]
    """
    result = []
    indices_to_check = set(range(length))

    while indices_to_check:
        component1 = indices_to_check.pop()
        parallel_class = {component1}
        for component2 in indices_to_check.copy():
            if are_parallel(iterable, component1, component2):
                parallel_class.add(component2)
                indices_to_check.remove(component2)
        result.append(parallel_class)
    return result


def classes_same_support(iterable) -> Iterator[set[SignVector]]:
    r"""
    Compute the classes with same support of given sign vectors.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors

    OUTPUT:
    A generator yielding sets of sign vectors with the same support.

    EXAMPLES::

        sage: from sign_vectors.utility import classes_same_support
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("+-0+"), sign_vector("-0+0")]
        sage: list(classes_same_support(L))
        [{(++0-), (+-0+)}, {(-0+0)}]
    """
    support_dict = {}
    for element in iterable:
        support = tuple(element.support())
        if support not in support_dict:
            support_dict[support] = {element}
        else:
            support_dict[support].add(element)
    yield from support_dict.values()
