import numpy as np

# Kyber512 parameters
q = 3329  # prime modulus
n = 256   # polynomial degree
k = 2     # dimension
eta1 = 3  # error distribution parameter for key generation
eta2 = 2  # error distribution parameter for encryption

# NTT-related constants (specific to Kyber512)
def montgomery_reduce(a):
    u = a * 2285
    u &= (1 << 16) - 1
    u *= q
    a = a + u
    return a >> 16

def barrett_reduce(a):
    u = (a * 5) >> 16
    u *= q
    return a - u

zetas = [ # pre-computed NTT-related constants from Kyber512 spec
    2285,  2184,  467,   1202,  287,   848,   329,   171,   1861,  442, 
    424,   1069,  293,   1851,  1561,  1872,  1285,  1323,  1757,  768, 
    1637,  1188,  804,   1705,  1708,  1235,  1883,  1256,  216,   1402, 
    1817,  1452,  1281,  1544,  1842,  1517,  1925,  1701,  398,   939, 
    1846,  1534,  1391,  1907,  301,   1597,  2039,  1239,  2233,  1822, 
    952,   1295,  175,   1850,  3010,  303,   3032,  336,   1286,  1645, 
    2230,  374,   245,   2294,  1318,  356,   2244,  320,   1919,  1588, 
    328,   2204,  204,   144,   3216,  1871,  1550,  308,   3231,  2353, 
    1870,  3292,  1236,  2737,  2424,  2531,  1408,  3042,  1235,  303, 
    2126,  3029,  127,   2980,  289,   3058,  2224,  1051,  3237,  2810, 
    844,   1087,  986,   2963,  1675,  3114,  1633,  335,   1046,  3307, 
    1901,  188,   1210,  3006,  1695,  220,   785,   3276,  3273,  3257, 
    1861,  2251,  304,   1115,  2096,  720,   240,   3038,  2085,  170, 
    1190,  2141,  2592,  1956,  2883,  1915,  312,   2347,  108,   2922, 
    194,   1234,  1812,  3155,  2963,  1093,  1774,  470,   2412,  1850, 
    2618,  373,   1124,  2294,  104,   2756,  2837,  960,   2348,  1293, 
    1422,  2874,  2792,  2877,  2658,  1251,  1768,  325,   1704,  850, 
    2987,  1903,  2684,  1871,  1336,  1454,  345,   1762,  201,   286, 
    1508,  156,   2157,  927,   2926,  2816,  246,   1569,  2974,  1422, 
    3053,  1485,  1433,  3065,  3126,  2418,  1676,  2180,  2078,  1685, 
    1204,  2361,  727,   1816,  2846,  1610,  1814,  3087,  1750,  2143, 
    2686,  1760,  2292,  395,   2419,  482,   2035,  2665,  2925,  3085, 
    1427,  2693,  2173,  2933,  2037,  2486,  3222,  2910,  1458,  1179, 
    3186,  3138,  967,   2103,  1327,  2935,  3130,  3024,  1587,  2862, 
    286,   2491,  2815,  2462,  1115,  1195,  360,   278,   1600,  1504, 
    3194,  139,   2773,  3262,  1578,  2775,  3002,  1762,  2431,  2698, 
    485,   2802,  3068,  2689,  2278,  258,   2425,  2546,  2113,  2117
]

def ntt(a):
    n = len(a)
    for i in range(1, n, 2):
        zeta = zetas[i // 2]
        for j in range(n // 2):
            t = montgomery_reduce(zeta * a[j + n // 2])
            a[j + n // 2] = barrett_reduce(a[j] + 2 * q - t)
            a[j] = barrett_reduce(a[j] + t)
    return a

def intt(a):
    n = len(a)
    for i in range(n - 2, -1, -2):
        zeta = zetas[i // 2]
        for j in range(n // 2):
            t = a[j]
            a[j] = barrett_reduce(t + a[j + n // 2])
            a[j + n // 2] = t + 2 * q - a[j + n // 2]
            a[j + n // 2] = montgomery_reduce(zeta * a[j + n // 2])
    for i in range(n):
        a[i] = montgomery_reduce(a[i] * n)
    return a

def poly_add(a, b, q):
    return [(x + y) % q for x, y in zip(a, b)]

def poly_sub(a, b, q):
    return [(x - y) % q for x, y in zip(a, b)]

def poly_mul(a, b, q):
    a_ntt = ntt(a)
    b_ntt = ntt(b)
    c_ntt = [(ai * bi) % q for ai, bi in zip(a_ntt, b_ntt)]
    return intt(c_ntt)

def sample_gaussian(size, sigma):
    return np.random.normal(0, sigma, size).astype(int) % q

def sample_poly(n, eta):
    return sample_gaussian(n, eta)

def keygen(q, n, k, eta1):
    s = [sample_poly(n, eta1) for _ in range(k)]
    e = [sample_poly(n, eta1) for _ in range(k)]

    a = np.random.randint(0, q, (k, n)) # random matrix
    t = [poly_add(poly_mul(ai, si, q), ei, q) for ai, si, ei in zip(a, s, e)]
    
    public_key = (a, t)
    private_key = s
    return public_key, private_key

def encapsulate(public_key, q, n, k, eta2):
    a, t = public_key
    r = [sample_poly(n, eta2) for _ in range(k)]
    e1 = [sample_poly(n, eta2) for _ in range(k)]
    e2 = sample_poly(n, eta2)

    u = [poly_add(poly_mul(ai, ri, q), e1i, q) for ai, ri, e1i in zip(a, r, e1)]
    v = poly_add(poly_mul(t[0], r[0], q), e2, q)
    
    m = np.random.randint(0, 2, n)
    v = poly_add(v, m, q)
    
    ciphertext = (u, v)
    shared_secret = m
    return ciphertext, shared_secret

def decapsulate(ciphertext, private_key, q, n, k):
    u, v = ciphertext
    s = private_key
    
    m_prime = poly_sub(v, poly_mul(u[0], s[0], q), q)
    
    return m_prime

# examples of how to use
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
