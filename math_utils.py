from galois import GF2, Poly, egcd
import numpy as np

def random_perm_matrix(n: int) -> GF2:
    """Generate a random permutation matrix."""
    return GF2([[1 if i == loc_in_row else 0 for i in range(n)] for loc_in_row in np.random.permutation(n)])


def random_binary_inv_matrix(n: int) -> GF2:
    """Try to generate a random binary non-singular matrix, by randomly generating and testing. Prob ~28% for n -> inf"""
    while(True):
        candidate = np.random.randint(2, size=(n, n)).tolist()
        det2 = det_mod_2(candidate)
        if det2 == 1:
            return GF2(candidate)


def det_mod_2(matrix: list[list[int]]) -> int:
    """
    Calculate the determinant mod 2.

    Notes
    -----
    Solution found on [this](https://stackoverflow.com/a/43707133) stack overflow answer.
    We convert the matrix into an upper trianguar matrix, so that the product of the diagonal is the determinant.

    We use the follwing rules:
    1. If the first column is all 0 then the determinant is 0.
    2. If the first column starts with a 1 and is only 0 after, then the determinant is the same as the determinant of the matrix obtained by removing the first row and first column.
    3. Swapping two rows has no effect on the determinant mod 2.
    4. Replacing a row by a sum of that row and another row has no effect on the determinant.
    """
    n = len(matrix)

    # walk through the smaller (n-i)x(n-i) submatricies by ignoring the first rows and columns
    for i in range(n):

        # find first 1 in the current submatrix
        j = i
        while j < n and matrix[j][i] == 0:
            j += 1

        # column is all 0, so det is 0
        if j == n:
            return 0

        # swap the row with 1 with the first row, if it isn't already
        if i < j:
            matrix[i], matrix[j] = matrix[j], matrix[i]

        # clear rest of column
        for j in range(i+1, n):
            
            # if the row contains a 1, add the first row to it, to have a 0 in the first place
            if matrix[j][i] == 1:
                matrix[j] = [(a + b) % 2 for a, b in zip(matrix[i], matrix[j])]

    return 1


def inverse_mod_poly(poly: Poly, mod: Poly) -> Poly:
    """Use the extende deuclidian algorith to find the inverse of a given polynomial mod another."""
    return egcd(poly, mod)[1]


def split_poly(poly: Poly) -> (Poly, Poly):
    """Split a polynomial into an even and odd part according to [this](http://icit.zuj.edu.jo/icit11/PaperList/Papers/Information%20Security/526_Risse.pdf) paper."""
    coeffs = poly.coeffs[::-1]
    p_even = Poly([np.sqrt(even_term) for even_term in coeffs[0::2]][::-1], field=poly.field)
    p_odd = Poly([np.sqrt(odd_term) for odd_term in coeffs[1::2]][::-1], field=poly.field)
    return p_even, p_odd


def sqrt_mod_poly(poly: Poly, mod: Poly) -> Poly:
    """Calulate the square root of a polynomial mod another polynomial according to [this](http://icit.zuj.edu.jo/icit11/PaperList/Papers/Information%20Security/526_Risse.pdf) paper."""
    g0, g1 = split_poly(mod)
    t0, t1 = split_poly(poly)
    return (t0 + g0 * inverse_mod_poly(g1, mod) * t1) % mod


def norm(a: Poly, b: Poly) -> int:
    """Calculate the norm accoding to Bernstein in [this](https://cr.yp.to/codes/goppalist-20081107.pdf) paper."""
    return 2**((a**2+Poly.Identity(b.field)*b**2).degree)


def lattice_basis_reduction(poly: Poly, mod: Poly) -> (Poly, Poly):
    """Compute the lattice reduction similar to [this](https://digitalcommons.csbsju.edu/cgi/viewcontent.cgi?article=1019&context=honors_theses) paper, which in turn uses the explanation in [this one](https://cr.yp.to/codes/goppalist-20110303.pdf)."""
    t = mod.degree
    a, b = [], []
    q, r = divmod(mod, poly)
    a.append(mod-q*poly)
    b.append(-q)

    if norm(a[0], b[0]) > 2**t:
        q, r = divmod(poly, a[0])
        a.append(r)
        b.append(Poly.One(poly.field) - q * b[0])
    else:
        return a[0], b[0]
    
    i = 1
    while norm(a[i], b[i]) > 2**t:
        q, r = divmod(a[i-1], a[i])
        a.append(r)
        b.append(b[i-1] - q * b[i])
        i += 1
    
    return a[i], b[i]