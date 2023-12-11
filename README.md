# McEliece Cryptosystem
The goal of this project is to offer a simple and understandable python implementation of the McEliece cryptosystem as well as the used Goppa Code. It is mainly aimed at students and others who have read papers and other theoretical explanations and need a visualization or just want to see how one might implement this cryptosystem.  

## Usage
You setup the cryptosystem by initializing a new instance of the McEliece class with the parameters of your choosing. In this instance, which acts as a single deployment of the system, one can use methods to:
- Generate new sets of private and public keys
- Encrypt a message with a public key
- Decrypt a cipher with a private key


## Explanation
I used the [galois](https://github.com/mhostetter/galois) library to do the heavy lifting of working in finite fields. If you are interested in how most of the methods work unter the hood, I can only refer the the incredible documentation of that library.

The start of the project is the `mceliece.py` file. You initialize a new cryptosystem with parameters $[m, t]$, instead of the usual $[n, k, t]$. This is because the cryptosystem only implements the default patterson algorithm and no other decoder, like [list decoding](https://cr.yp.to/codes/goppalist-20081107.pdf). Because of this, the $n$ is set as $2^m$ and $k$ is the minimal possible value and can be caluculated with the given values. 

The used Goppa Code and the patterson algorithm are implemented in the `goppa_code.py` file. I decided to go for the slower but more visual matrix implementation, instead of the polynomial approach that Bernstein explains [here](https://cr.yp.to/codes/goppalist-20081107.pdf). You can find sage code for that in [this](https://digitalcommons.csbsju.edu/cgi/viewcontent.cgi?article=1019&context=honors_theses) paper.

Lastly some math functions are implemented in the `mceliece.py` file. The source for most methods are given in the docstrings.

## Sources
My initial inspiration was [this](https://github.com/jkrauze/mceliece) repository, which already implements the McEliece cryptosystem in python. 

I mainly used [this](https://surface.syr.edu/cgi/viewcontent.cgi?article=1846&context=honors_capstone) paper to study the cryptosystem. The example given is very helpful, even though it contains some minor mistakes in the example calculated in section 2.6. If you want to verify the results, the X matrix of the parity-check matrix is computed with the goppa polynomial $g(x) = x^2 + \alpha^7x+1$, while the polynomial used for the Z matrix and any following computation is $g(x) = x^2+x+\alpha^3$.

A great resource for the Goppa Code and Patterson Algorithm is [this](https://cr.yp.to/codes/goppalist-20081107.pdf) paper. You can verify the results given in the example with this code. 

Lastly I took some help from [this](https://digitalcommons.csbsju.edu/cgi/viewcontent.cgi?article=1019&context=honors_theses) paper regarding some code in the Patterson Algorithm. 

## Disclaimer
I developed this project as part of my bachelor thesis. I am open to feedback and happy to try to fix any bugs, but currently do not plan on updating the code to release any new features.