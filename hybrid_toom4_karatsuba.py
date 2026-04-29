import random

def exGCD(a: int, n: int) -> int :
    i = 0
    r = [n, a]
    u = [0, 1]
    q = [0]
    while r[i+1] != 1:
        i = i + 1
        q.append(r[i-1]//r[i])
        r.append(r[i-1] - r[i]*q[i])
        u.append(u[i-1] - u[i]*q[i])
    return u[i+1]%n

def cut(f:list):
    for i in range(len(f)-1,-1,-1):
        if f[i] != 0:
            return [f[j] for j in range(0,i+1)]
    return []

def poly_sum_fixed_length(f:list, g:list, p:int):
    assert len(f) == len(g)
    if p == 0:
        return [f[i]+g[i] for i in range(0,len(f))]
    else:
        return [(f[i]+g[i])%p for i in range(0,len(f))]
    
def poly_diff_fixed_length(f:list, g:list, p:int):
    assert len(f) == len(g)
    if p == 0:
        return [f[i]-g[i] for i in range(0,len(f))]
    else:
        return [(f[i]-g[i])%p for i in range(0,len(f))]

def poly_add_mod(f:list[int], g:list[int], p:int) -> list[int]:
    assert len(f) == len(g)
    if p == 0:
        return [f[i]+g[i] for i in range(0,len(f))]
    else:
        return [(f[i]+g[i])%p for i in range(0,len(f))]
    
def poly_sub_mod(f:list, g:list, p:int) -> list[int]:
    assert len(f) == len(g)
    if p == 0:
        return [f[i]-g[i] for i in range(0,len(f))]
    else:
        return [(f[i]-g[i])%p for i in range(0,len(f))]

def poly_mlt(f:list, g:list, p:int):
    if (len(f) == 0) or (len(g) == 0):
        return[]
    n, m = len(f) - 1, len(g) - 1
    h = [0 for _ in range(0,n+m+1)]
    for k in range(0,n+m+1):
        for i in range(max([0,k-m]),min([k,n])+1):
            if p == 0:
                h[k] = h[k] + f[i]*g[k-i]
            else:
                h[k] = (h[k] + ((f[i]*g[k-i])%p))%p
    return cut(h)

def poly_shift(f:list, k:int) -> list:
    if not f:
        return []
    return [0] * k + f

def poly_mod_xn1(f: list, N: int, q: int) -> list:
    result = [0] * N
    for k in range(len(f)):
        c = f[k] % q
        if c == 0:
            continue
        if k < N:
            result[k] = (result[k] + c) % q
        else:
            idx = k - N
            result[idx] = (result[idx] - c) % q  
    return cut(result)

def compose_toom_blocks_mod(h, d, d1, q):
    """
    Compose shifted blocks.

    Given blocks h[i], this function computes
        sum_i h[i] x^{i*d1}
    modulo q.

    Each block is assumed to have length at most 2*d1 - 1.
    The returned list has length exactly 2*d - 1.
    """
    out = [0] * (2*d-1)
    block_len = 2*d1 - 1

    h0 = h[0] + [0]*(block_len - len(h[0]))
    for j in range(block_len):
        out[j] = h0[j]

    for i in range(1, len(h)):
        base = i * d1
        hi = h[i] + [0]*(block_len - len(h[i]))

        for j in range(d1):
            out[base + j] = (out[base + j] + hi[j]) % q

        for j in range(d1, 2*d1-1):
            out[base + j] = hi[j]

    return out

def Karatsuba(f: list[int], g: list[int], q: int) -> list[int]:
    assert len(f) == len(g)
    d = len(f)

    if d == 1:
        return [(f[0] * g[0]) % q]

    while d % 2 != 0:
        d += 1

    f = f + [0] * (d - len(f))
    g = g + [0] * (d - len(g))

    d1 = d // 2
    f0, f1 = f[:d1], f[d1:]
    g0, g1 = g[:d1], g[d1:]

    F1 = [(f0[i] + f1[i]) % q for i in range(d1)]
    G1 = [(g0[i] + g1[i]) % q for i in range(d1)]

    h0 = Karatsuba(f0, g0, q)
    h2 = Karatsuba(f1, g1, q)
    m1 = Karatsuba(F1, G1, q)
    h1 = poly_diff_fixed_length(poly_diff_fixed_length(m1, h0, q), h2, q)

    return compose_toom_blocks_mod([h0, h1, h2], d=d, d1=d1, q=q)

def toom4_evaluation(blocks, q):
    """
    Toom-4 evaluation for blocks [f0, f1, f2, f3].
    Returns:
        {
            0:   F(0),
            1:   F(1),
            -1:  F(-1),
            2:   F(2),
            -2:  F(-2),
            3:   F(3),
            "inf": F(infinity)
        }
    """
    f0, f1, f2, f3 = blocks

    a1 = poly_add_mod(f0, f2, q)
    a2 = poly_add_mod(f1, f3, q)
    a3_p = poly_add_mod(a1, a2, q)      # F(1)
    a3_m = poly_sub_mod(a1, a2, q)      # F(-1)

    a4 = poly_add_mod(f3, f3, q)        # 2 f3
    a5 = poly_add_mod(a4, f3, q)        # 3 f3

    a6_p = poly_add_mod(a3_p, a5, q)
    a6_m = poly_sub_mod(a3_m, a5, q)

    a7_p = poly_add_mod(a6_p, f2, q)
    a7_m = poly_add_mod(a6_m, f2, q)

    a8_p = poly_add_mod(a7_p, a7_p, q)
    a8_m = poly_add_mod(a7_m, a7_m, q)

    a9_p = poly_sub_mod(a8_p, f0, q)    # F(2)
    a9_m = poly_sub_mod(a8_m, f0, q)    # F(-2)

    a10 = poly_add_mod(a5, a5, q)
    a11 = poly_add_mod(poly_sub_mod(a9_p, f0, q), a10, q)
    a12 = poly_add_mod(a11, a11, q)
    a13 = poly_add_mod(a12, a3_m, q)    # F(3)

    return {
        0: f0,
        1: a3_p,
        -1: a3_m,
        2: a9_p,
        -2: a9_m,
        3: a13,
        "inf": f3,
    }

def addition_chain_multiples(poly, program, q):
    multiples = {1: poly.copy()}
    for c, a, b in program:
        multiples[c] = poly_sum_fixed_length(multiples[a], multiples[b], q)
    return multiples

def signed_sum_fixed_length(terms, q):
    """
    Compute a signed sum of fixed-length polynomials.

    terms = [(+1, poly), (-1, poly), ...]
    The first term must have sign +1.
    """
    assert len(terms) > 0
    assert terms[0][0] == 1

    out = terms[0][1]

    for sign, poly in terms[1:]:
        if sign == 1:
            out = poly_sum_fixed_length(out, poly, q)
        elif sign == -1:
            out = poly_diff_fixed_length(out, poly, q)
        else:
            raise ValueError("sign must be +1 or -1")

    return out

TOOM4_POINTS = [0, 1, -1, 2, -2, 3, "inf"]

TOOM4_CHAIN_PROGRAMS = {
    # H(0)
    0: [
        (2,1,1), (4,2,2), (8,4,4), (10,2,8), (20,10,10),
        (30,10,20), (40,20,20), (50,20,30), (100,50,50),
        (120,20,100), (150,50,100)
    ],
    # H(1)
    1: [
        (2,1,1), (4,2,2), (8,4,4), (10,2,8), (20,10,10),
        (30,10,20), (40,20,20), (70,30,40), (80,40,40),
        (120,40,80)
    ],
    # H(-1)
    -1: [
        (2,1,1), (4,2,2), (5,1,4), (10,5,5), (20,10,10),
        (40,20,20), (60,20,40), (80,40,40)
    ],
    # H(2)
    2: [
        (2,1,1), (4,2,2), (5,1,4), (10,5,5), (20,10,10),
        (30,10,20), (35,5,30)
    ],
    # H(-2)
    -2: [
        (2,1,1), (4,2,2), (5,1,4), (6,1,5)
    ],
    # H(3)
    3: [
        (2,1,1), (4,2,2), (5,1,4)
    ],
    # H(inf)
    "inf": [
        (2,1,1), (4,2,2), (8,4,4), (10,2,8), (20,10,10),
        (40,20,20), (80,40,40), (120,40,80), (240,120,120),
        (360,120,240), (480,120,360), (600,120,480),
        (1080,480,600), (1440,360,1080), (1800,360,1440)
    ],
}

def Toom4_Karatsuba_recursive(f:list[int], g:list[int], q:int, L:int, depth:int =0) -> list[int]:
    """
    Apply Toom-4 exactly until depth reaches L, then use Karatsuba.
    The output is scaled by 120^L.
    """
    assert len(f) == len(g)

    if depth == L:
        return Karatsuba(f,g,q)

    d = len(f)

    # Completely base case
    if d == 1:
        return [(f[0]*g[0])%q]
    
    # 0) Padding
    while (d%4) != 0:
            d = d + 1
    f = f + [0]*(d-len(f))
    g = g + [0]*(d-len(g))

    # 1) Splitting
    d1 = d//4
    f_blocks = [f[i*d1:(i+1)*d1] for i in range(4)]
    g_blocks = [g[i*d1:(i+1)*d1] for i in range(4)]

    # 2) Evaluation
    F_eval = toom4_evaluation(f_blocks, q)
    G_eval = toom4_evaluation(g_blocks, q)

    # 3) Recursive multiplication
    M = {
        alpha: Toom4_Karatsuba_recursive(F_eval[alpha], G_eval[alpha], q, L, depth+1)
        for alpha in TOOM4_POINTS
    } # [H(0),H(1),H(-1),H(2),H(-2),H(3),H(inf)]
    
    # 4) Interpolation: scalar multiples via fixed addition programs
    M_multiple = {
        alpha: addition_chain_multiples(M[alpha], TOOM4_CHAIN_PROGRAMS[alpha], q)
        for alpha in TOOM4_POINTS
    }
    
    h_primes = []
    
    # h0' = 120 M_0
    h_primes.append(M_multiple[0][120])
    
    # h1' = -40M_0 + 120M_1 - 60M_{-1} - 30M_2 + 6M_{-2} + 4M_3 - 1440M_inf
    h_primes.append(signed_sum_fixed_length([
        (+1, M_multiple[1][120]),
        (+1, M_multiple[-2][6]),
        (+1, M_multiple[3][4]),
        (-1, M_multiple[0][40]),
        (-1, M_multiple[-1][60]),
        (-1, M_multiple[2][30]),
        (-1, M_multiple["inf"][1440]),
    ], q))
    
    # h2' = -150M_0 + 80M_1 + 80M_{-1} - 5M_2 - 5M_{-2} + 480M_inf
    h_primes.append(signed_sum_fixed_length([
        (+1, M_multiple[1][80]),
        (+1, M_multiple[-1][80]),
        (+1, M_multiple["inf"][480]),
        (-1, M_multiple[0][150]),
        (-1, M_multiple[2][5]),
        (-1, M_multiple[-2][5]),
    ], q))
    
    # h3' = 50M_0 - 70M_1 - 5M_{-1} + 35M_2 - 5M_{-2} - 5M_3 + 1800M_inf
    h_primes.append(signed_sum_fixed_length([
        (+1, M_multiple[0][50]),
        (+1, M_multiple[2][35]),
        (+1, M_multiple["inf"][1800]),
        (-1, M_multiple[1][70]),
        (-1, M_multiple[-1][5]),
        (-1, M_multiple[-2][5]),
        (-1, M_multiple[3][5]),
    ], q))
    
    # h4' = 30M_0 - 20M_1 - 20M_{-1} + 5M_2 + 5M_{-2} - 600M_inf
    h_primes.append(signed_sum_fixed_length([
        (+1, M_multiple[0][30]),
        (+1, M_multiple[2][5]),
        (+1, M_multiple[-2][5]),
        (-1, M_multiple[1][20]),
        (-1, M_multiple[-1][20]),
        (-1, M_multiple["inf"][600]),
    ], q))
    
    # h5' = -10M_0 + 10M_1 + 5M_{-1} - 5M_2 - M_{-2} + M_3 - 360M_inf
    h_primes.append(signed_sum_fixed_length([
        (+1, M_multiple[1][10]),
        (+1, M_multiple[-1][5]),
        (+1, M_multiple[3][1]),
        (-1, M_multiple[0][10]),
        (-1, M_multiple[2][5]),
        (-1, M_multiple[-2][1]),
        (-1, M_multiple["inf"][360]),
    ], q))
    
    # h6' = 120 M_inf
    h_primes.append(M_multiple["inf"][120])

    out = compose_toom_blocks_mod(h_primes, d=d, d1=d1, q=q)
    return out

def Hybrid_Toom4_Karatsuba(L:int, f:list[int], g:list[int], u:int ,q:int) -> list[int]:
    """
    Toom-4/Karatsuba hybrid with exactly L Toom-4 levels.
    u should be the inverse of 120^L modulo q.
    """
    h = cut(Toom4_Karatsuba_recursive(f,g,q,L,depth=0))
    return [(u*h[i])%q for i in range(len(h))]

def Number_of_recurse_for_Toom4(d: int):
    if d == 1:
        return 0
    else:
        while (d%4) != 0:
            d = d + 1
        return 1 + Number_of_recurse_for_Toom4(d//4)

def main() -> None:
    """
    Small self-test and timing demo.  This block is executed only when this
    file is run directly, not when the file is imported as a module.
    """
    import time

    d, q = 512, 3329
    L = 2
    u = exGCD(pow(120, L, q), q)

    f = [random.randint(0, q - 1) for _ in range(d)]
    g = [random.randint(0, q - 1) for _ in range(d)]

    start = time.time()
    result_normal = poly_mlt(cut(f), cut(g), q)
    end1 = time.time()

    result_karatsuba = cut(Karatsuba(f, g, q))
    end2 = time.time()

    result_hybrid = Hybrid_Toom4_Karatsuba(L, f, g, u, q)
    end3 = time.time()

    print("time for    Normal =", end1 - start)
    print("time for Karatsuba =", end2 - end1)
    print("time for    Toom-4 =", end3 - end2)
    print("correctness", result_normal == result_karatsuba)
    print("correctness", result_normal == result_hybrid)


if __name__ == "__main__":
    main()
