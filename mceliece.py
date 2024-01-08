from math_utils import random_binary_inv_matrix, random_perm_matrix
from goppa_code import GoppaCode, GF2
from dataclasses import dataclass
import random, numpy as np


@dataclass
class SecretKey:
    S_inv: GF2
    P_inv: GF2
    goppa: GoppaCode


class McEliece:
    def __init__(self, m: int, t: int):
        self.m = m
        self.t = t

        # Fix this for now, as other values require a different decoder for the goppa code
        self.n = 2**m
        self.k = self.n - t * m

    def generate_key_pair(self) -> (SecretKey, GF2):
        """Generate and return the secret and public key"""

        goppa_code = GoppaCode(self.n, self.m, self.t)
        S = random_binary_inv_matrix(self.k)
        P = random_perm_matrix(self.n)

        Gp = S @ goppa_code.G @ P

        P_inv = np.linalg.inv(P)
        S_inv = np.linalg.inv(S)

        return SecretKey(S_inv, P_inv, goppa_code), Gp

    def encrypt(self, message: np.ndarray, Gp: GF2) -> np.ndarray:
        """Encrypt a message with someones public key."""

        if len(message) != self.k:
            raise Exception(f"Wrong message length. It has to be {self.k} bits.")

        # generate codeword by multiplying message with generator matrix
        codeword = np.array(GF2(message) @ Gp)

        # generate t random errors by generating list of all possible error locations and shuffle
        error_locations = list(range(self.n))
        random.shuffle(error_locations)

        # apply errors
        for loc in error_locations[: self.t]:
            codeword[loc] ^= 1

        return codeword

    def decrypt(self, cipher: np.ndarray, sk: SecretKey) -> np.ndarray:
        """Decrypt a received cipher with a secret key."""

        if len(cipher) != self.n:
            raise Exception(f"Wrong cipher length. It has to be {self.n} bits.")

        unshuffled_cipher = GF2(cipher) @ sk.P_inv
        decoded_cipher = sk.goppa.decode(unshuffled_cipher)
        message = GF2(decoded_cipher) @ sk.S_inv
        return message


# EXAMPLE USAGE
# mc = McEliece(8, 22)
# sk, pk = mc.generate_key_pair()
# cipher = mc.encrypt(np.array([1] * mc.k), pk)
# message = mc.decrypt(cipher, sk)
# print(message)
