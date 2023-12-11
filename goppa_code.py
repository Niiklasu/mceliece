import numpy as np
from galois import GF, GF2, Poly, irreducible_poly
from math_utils import sqrt_mod_poly, inverse_mod_poly, lattice_basis_reduction
"""
    Recommended [n, k, t] Goppa parameters according to https://eprint.iacr.org/2008/318.pdf
    80-Bit:  (1632, 1269, 33)  with m = 11
    128-Bit: (2960, 2288, 56)  with m = 12
    256-Bit: (6624, 5129, 118) with m = 13

    !!! Doesn't work with current Version that runs w/o list decoding

"""
class GoppaCode:
    def __init__(self, n: int, m: int, t: int):
        # Generate extention field 
        GF2m = GF(2**m)

        # support L is just list of elements, as this can be found in polynomial time anyway (i think, dont quote me on this)
        support = list(GF2m.elements)

        # Generate g as an irreducible polynomial, so no element is a root
        g_poly = irreducible_poly(n, t, method="random")

        # Generate Parity-Check Matrix
        X = GF2m([[g_poly.coeffs[row - col] if 0 <= row - col < t else 0 for col in range(t)] for row in range(t)])
        Y = GF2m([[support[col] ** row for col in range(n)] for row in range(t)])
        Z = GF2m([[g_poly(support[row])**-1 if col == row else 0 for col in range(n)] for row in range(n)])
        H = X @ Y @ Z

        # Convert each value to a column of the bin representator
        H_bin = GF2(np.concatenate([np.column_stack([elem.vector() for elem in row]) for row in H]))

        # Calculate the nullspace to find G
        G = H_bin.null_space()

        self.G = G
        self.g_poly = g_poly
        self.GF2m = GF2m
        self.supp = support


    def patterson_alg(self, received_message: np.ndarray) -> np.ndarray:
        """Patterson Algorithm to find and fix the errors. Mainly following the desription in [this](https://cr.yp.to/codes/goppalist-20081107.pdf) paper."""

        # generate the syndrome vector
        inverse_polys = [
            inverse_mod_poly(Poly.Str(f"x - {self.supp[i]}", field=self.GF2m), self.g_poly)
            for i in range(len(received_message))
            if received_message[i] == 1
        ]
        syndrome = Poly.Zero(self.GF2m)
        for poly in inverse_polys:
            syndrome += poly

        # Inverting the syndrome and adding x, then taking the square root.
        # Fails when s(x) == 0, but shouldn't (if i'm not mistaken) occur anyway, as we will always have an error
        syndrome_inv = inverse_mod_poly(syndrome, self.g_poly)    
        s = sqrt_mod_poly(syndrome_inv + Poly.Identity(self.GF2m), self.g_poly)

        # Use lattice basis reduction to find a = b * s mod g of least degrees.
        # This can also be done using a slightly modified extended euclidian algorithm.
        alpha, beta = lattice_basis_reduction(s, self.g_poly)
        
        # Combine the found vectors to a new vector sigma and fing the roots
        # im power-representation the exponent is the location of the error e.g. alpha^25 as a root means a error is in the 25th location
        sigma = (alpha**2 + Poly.Identity(self.GF2m) * beta**2)
        error_locations = np.array(sigma.roots())

        # fix errors
        for loc in error_locations:
            received_message[loc] ^= 1

        return received_message 
    

    def decode(self, received_message: np.ndarray) -> np.ndarray:
        """Decoding algorithm to decode a received cipher."""
        k = self.G.shape[0]

        # First remove the errors using the patterson algorithm
        fixed_message = self.patterson_alg(received_message)

        # Set up a system as (G^T | c^T)
        system_to_solve = np.append(self.G.T, fixed_message.reshape(len(fixed_message), 1), axis=1)  

        # Using Gaussian Elimination to create an identity matrix in all columns but the last
        I_m = GF2(system_to_solve).row_reduce(ncols=k)

        # Only looking at the first k rows, the system is now in form (I | m^T) so the message can easily be received
        return np.array(I_m[:k, -1].T)

