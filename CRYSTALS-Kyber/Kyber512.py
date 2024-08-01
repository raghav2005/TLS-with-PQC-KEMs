import numpy as np
from numpy.fft import fft, ifft

# Kyber512 parameters
q = 3329  # prime modulus
n = 256   # polynomial degree (Kyber uses n=256)
k = 2     # dimension
eta1 = 3  # error distribution parameter for key generation
eta2 = 2  # error distribution parameter for encryption

# utility functions
def bit_reverse(a, n):
    log_n = n.bit_length() - 1
    return int(f'{a:0{log_n}b}'[::-1], 2)

def ntt(a, q): # number theoretic transform
    n = len(a)
    omega = np.array([pow(17, bit_reverse(i, n), q) for i in range(n)])
    return np.round(np.real(fft(a) * omega) % q).astype(int)

def intt(a, q): # inverse ntt
    n = len(a)
    inv_n = pow(n, q - 2, q)
    omega = np.array([pow(17, -bit_reverse(i, n), q) for i in range(n)])
    return np.round(np.real(ifft(a * omega) * inv_n) % q).astype(int)

def poly_add(a, b, q):
    return [(x + y) % q for x, y in zip(a, b)]

def poly_sub(a, b, q):
    return [(x - y) % q for x, y in zip(a, b)]

def poly_mul(a, b, q, n): # polynomial multiplication using ntt
    a_ntt = ntt(a, q)
    b_ntt = ntt(b, q)
    c_ntt = [(ai * bi) % q for ai, bi in zip(a_ntt, b_ntt)]
    res = intt(c_ntt, q)
    return res

def sample_gaussian(size, sigma):
    return np.random.normal(0, sigma, size).astype(int) % q

def sample_poly(n, eta):
    return sample_gaussian(n, eta)

def keygen(q, n, k, eta1):
    s = [sample_poly(n, eta1) for _ in range(k)]
    e = [sample_poly(n, eta1) for _ in range(k)]
    
    a = np.random.randint(0, q, (k, n))  # Random matrix
    t = [poly_add(ntt(ai, q), ntt(si, q), q) for ai, si in zip(a, s)]
    
    public_key = (a, t)
    private_key = s
    return public_key, private_key

def encapsulate(public_key, q, n, k, eta2):
    a, t = public_key
    r = [sample_poly(n, eta2) for _ in range(k)]
    e1 = [sample_poly(n, eta2) for _ in range(k)]
    e2 = sample_poly(n, eta2)
    
    u = [poly_add(poly_mul(ntt(ai, q), ntt(ri, q), q, n), e1i, q) for ai, ri, e1i in zip(a, r, e1)]
    v = poly_add(poly_mul(ntt(t[0], q), ntt(r[0], q), q, n), e2, q)
    
    m = np.random.randint(0, 2, n)
    v = poly_add(v, ntt(m, q), q)
    
    ciphertext = (u, v)
    shared_secret = m
    return ciphertext, shared_secret

def decapsulate(ciphertext, private_key, q, n, k):
    u, v = ciphertext
    s = private_key
    
    m_prime = poly_sub(v, poly_mul(ntt(u[0], q), ntt(s[0], q), q, n), q)
    shared_secret = intt(m_prime, q)
    
    return shared_secret

# example of how to use
public_key, private_key = keygen(q, n, k, eta1)
print("Public Key:", public_key)
print("Private Key:", private_key)

ciphertext, shared_secret_enc = encapsulate(public_key, q, n, k, eta2)
print("Ciphertext:", ciphertext)
print("Shared Secret (Encapsulation):", shared_secret_enc)

shared_secret_dec = decapsulate(ciphertext, private_key, q, n, k)
print("Shared Secret (Decapsulation):", shared_secret_dec)

assert np.array_equal(shared_secret_enc, shared_secret_dec), "Shared secrets do not match!"
print("Shared secrets match!")
