from sage.all import *


def automorphism_group(g):
    auto_group = gap(g).AutomorphismGroup()
    return PermutationGroup(gap_group=auto_group.AsPermGroup())

def coprime_integers(n):
    return [m for m in range(1, n) if gcd(m, n) == 1]

def automorphisms_of_cyclic_group(n):
    """
    Returns all automorphisms of all cyclic group of order n
    """
    fs = [lambda x: x]

    ms = coprime_integers(n)
    for m in ms[1:]:
        f = {0: 0}
        for i in range(1, n):
            f[i] = (m*i) % n
        fs.append(lambda x, f=f: f[x])

    return fs

def automorphism_orbits(n):
    elems = set(list(range(n)))

    morphs = automorphisms_of_cyclic_group(n)

    groups = []
    mappings = {}

    while elems:
        elem = elems.pop()
        mappings[elem] = {elem: lambda x: x}
        group = [elem]
        for morph in morphs:
            image = morph(elem)
            if image in elems:
                group.append(image)
                elems.remove(image)
                mappings[elem][image] = morph

        groups.append(group)

    return groups, mappings

def inverse_automorphism(f, n):
    xs = range(n)
    ys = [f(x) for x in xs]
    mapping = {y: x for y, x in zip(ys, xs)}
    def finv(y):
        return mapping[y]
    return finv

def apply_automorphism_to_face(pi, phi):
    """
    Args:
        pi (list): The inequality/equality defining the face as a list of
            the coefficients, in order [b, a1, a2, ..., an]
        phi (func): A function mapping {1, ..., n} to itself
    """
    pi2 = pi.copy()
    n = len(pi) - 1
    for i in range(1, n+1):
        pi2[i] = pi[phi(i)]
    return pi2

def apply_automorphism_to_vertex(v, phi):
    """
    Args:
        v (list): A list of the coordinates of the vertex
        phi (func): A function mapping {1, ..., n} to itself
    """
    v2 = v.copy()
    n = len(v)
    for i in range(1, n+1):
        v2[i-1] = v[phi(i)-1]
    return v2
