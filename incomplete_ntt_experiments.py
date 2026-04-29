import random
import math
import time
from statistics import mean, median, stdev
import hybrid_toom4_karatsuba

def v2(m:int) -> int:
    """
    Return the 2-adic valuation of a given integer m.
    """
    c = 0
    while m % 2 == 0:
        c += 1
        m //= 2
    return c

def primitive_Nth_root_of_unity_mod_p(N:int, p:int, P:list[int]) -> int:
    """
    Return a primitive N-th root of unity modulo p.
    Assumes p is prime, N | (p-1), and P contains all prime divisors of p-1.
    """
    assert (p - 1) % N == 0
    Q = [(p-1)//P[i] for i in range(0,len(P))]
    for a in range(2,p):
        flag = 1
        for i in range(0,len(Q)):
            if pow(a,Q[i],p) == 1:
                flag = 0
                break
        if flag == 1:
            zeta = pow(a,(p-1)//N,p)
            assert pow(zeta,N,p) == 1
            if N > 1:
                assert pow(zeta,N//2,p) != 1
            return zeta
    raise ValueError("No primitive root found.")

def conv_prod(n:int, q:int, a:list[int], b:list[int]) -> list[int]:
    """
    Direct product modulo x^n+1.
    Used only for correctness checking.
    """
    assert len(a) == n and len(b) == n
    c = [0]*n
    for k in range(n):
        c[k] = (a[0]*b[k])%q
        for i in range(1,k+1):
            c[k] = (c[k] + (a[i]*b[k-i])%q)%q
        for i in range(k+1,n):
            c[k] = (c[k] - (a[i]*b[n+k-i])%q)%q
    return c

def exGCD(a:int, n:int) -> int:
    """
    Given coprime integers a and n,
    return an integer x with 1 =< x < n such that a*x = 1 (mod n).
    """
    i = 0
    r = [n,a]
    u = [0,1]
    while r[i+1] != 1:
        i = i+1
        q = r[i-1]//r[i]
        r.append(r[i-1]-q*r[i])
        u.append(u[i-1]-q*u[i])
    return (u[i+1])%n

def precompute_IncompleteNTT_data(N:int, q:int, t:int, zeta:int):
    assert (q - 1) % N == 0
    e=[[N//2]]
    for l in range(t):
        e_l=[0]*(2**(l+1))
        for i in range(2**l):
            e_l[2*i] = e[l][i]//2
            e_l[2*i+1] = e_l[2*i]+(N//2)
        e.append(e_l)
    zeta_powers_forward = [ [pow(zeta,e[l][i]//2,q) for i in range(2**l)] for l in range(t+1)]
    zeta_powers_inverse = [ [pow(zeta,(N//2)-(e[l][i]//2),q) for i in range(2**l)] for l in range(t)]
    zeta_powers_components = [ pow(zeta,e[t][i],q) for i in range(2**t)]
    inv_2t = exGCD(2**t,q)
    return {
        "forward": zeta_powers_forward,
        "inverse": zeta_powers_inverse,
        "components": zeta_powers_components,
        "inv_2t": inv_2t,
    }
    
def IncompleteNTT(n:int, t:int, q:int, zeta_powers, f:list[int]) -> list[int]:
    """
    In-place-style incomplete NTT.
    """
    length = n//2
    for l in range(t):
        i = 0
        for s in range(0,n,2*length):
            zeta = zeta_powers[l][i]
            i += 1
            for j in range(s,s+length):
                tmp = (zeta*f[j+length])%q
                f[j+length] = (f[j] - tmp)%q
                f[j] = (f[j] + tmp)%q
        length //= 2
    return f

def IncompleteINTT(n:int, t:int, q:int, zeta_powers, inv_2t:int, h:list[int]) -> list[int]:
    """
    In-place-style incomplete inverse NTT.
    """
    length = n//(2**t)
    for l in range(t-1,-1,-1):
        i = 0
        for s in range(0,n,2*length):
            zeta = zeta_powers[l][i]
            i += 1
            for j in range(s,s+length):
                tmp = h[j]
                h[j] = (tmp + h[j+length])%q
                h[j+length] = (zeta*((h[j+length]-tmp)%q))%q
        length *= 2
    return [(inv_2t*h[i])%q for i in range(n)]

# ============================================================
# Component multiplications modulo x^d - zeta
# ============================================================

def component_mlt_schoolbook(d:int, q:int, zeta:int, a:list[int], b:list[int]) -> list[int]:
    """
    Multiply two degree < d polynomials modulo x^d - zeta by schoolbook.
    """
    assert len(a) == d and len(b) == d
    c = [0]*d
    for k in range(d):
        c[k] = (a[0]*b[k])%q
        for i in range(1,k+1):
            c[k] = (c[k] + (a[i]*b[k-i])%q)%q
        if zeta == q - 1:
            for i in range(k+1, d):
                c[k] = (c[k] - (a[i]*b[d+k-i]) % q) % q
        elif zeta == 1:
            for i in range(k+1, d):
                c[k] = (c[k] + (a[i]*b[d+k-i]) % q) % q
        else:
            for i in range(k+1, d):
                c[k] = (c[k] + zeta*((a[i]*b[d+k-i]) % q) % q) % q
    return c

def component_mlt_hybrid(d:int, L:int, q:int, inv_120L:int, zeta:int, a:list[int], b:list[int]) -> list[int]:
    """
    Multiply modulo x^d - zeta.
    The ordinary product is computed by Hybrid_Toom4_Karatsuba.Hybrid_Toom4.
    """
    assert len(a) == d and len(b) == d
    prod = hybrid_toom4_karatsuba.Hybrid_Toom4_Karatsuba(L,a,b,inv_120L,q)
    # prod1 = Hybrid_Toom4_Karatsuba.poly_mlt(a,b,q)
    # assert np.all(prod == prod1)
    prod = prod + [0 for _ in range(len(prod),2*d-1)]
    c = prod[:d]
    if zeta == q-1:
        for k in range(d-1):
            c[k] = (c[k] - prod[k+d])%q
    elif zeta == 1:
        for k in range(d-1):
            c[k] = (c[k] + prod[k+d])%q
    else:
        for k in range(d-1):
            c[k] = (c[k] + (zeta*prod[k+d])%q)%q
    return c

def conv_prod_IncompleteNTT_hybrid(
    n: int,
    q: int,
    t: int,
    bottom_deg: int,
    L: int,
    inv_120L: int,
    ntt_data,
    f: list[int],
    g: list[int],
    method: str
) -> list[int]:
    """
    method:
      - "hybrid": component multiplication by Hybrid_Toom4 with depth L
      - "schoolbook": component multiplication by schoolbook
    """
    f_ntt = IncompleteNTT(n,t,q,ntt_data["forward"],f[:])
    g_ntt = IncompleteNTT(n,t,q,ntt_data["forward"],g[:])
    prod = []
    for j in range(2**t):
        a = f_ntt[bottom_deg*j: bottom_deg*(j+1)]
        b = g_ntt[bottom_deg*j: bottom_deg*(j+1)]
        if method == "hybrid":
            c = component_mlt_hybrid(bottom_deg,L,q,inv_120L,ntt_data["components"][j],a,b)
        elif method == "schoolbook":
            c = component_mlt_schoolbook(bottom_deg,q,ntt_data["components"][j],a,b)
        else:
            raise ValueError(f"Unknown method: {method}")
        prod.extend(c)
    prod = IncompleteINTT(n,t,q,ntt_data["inverse"],ntt_data["inv_2t"],prod)
    return prod

# ============================================================
# Theoretical L_opt
# ============================================================

def Hybrid_Toom4_complexity(d: int, L:int) -> tuple[int,int]:
    """
    Same theoretical model as in the paper:
    ordinary product in F_q[x], not including reduction modulo x^d-zeta.
    """
    assert d > 0 and (d & (d - 1) == 0)
    assert L >=0 and 4**L <= d
    base = d//(4**L)
    Lambda = int(math.log2(base))
    T_M = (7**L)*(3**Lambda) 
    if L >0:
        T_M += 2*d-1
    T_A = (7**L)*8*(3**Lambda- base) + 72*(7**L)*base - 72*d
    return T_M, T_A

def L_max_for_d(d: int) -> int:
    assert d > 0 and (d & (d - 1) == 0)
    return int(math.log2(d)) // 2

def theoretical_L_opt(d: int, w: float = 0.2) -> tuple[int, float]:
    L_max = L_max_for_d(d)
    vals = []
    for L in range(L_max + 1):
        T_M, T_A = Hybrid_Toom4_complexity(d, L)
        vals.append((T_M + w*T_A,L))
    cost, L = min(vals)
    return L, cost

def print_theoretical_L_table(d_values: list[int], w: float = 0.2):
    print("\nTheoretical costs for Toom-4/Karatsuba hybrid:")
    print("d, L, T_M, T_A, C_w")
    Lopt_dict = {}

    for d in d_values:
        L_max = L_max_for_d(d)
        vals = []
        print(f"\nd = {d}")
        for L in range(L_max + 1):
            T_M, T_A = Hybrid_Toom4_complexity(d, L)
            Cw = T_M + w*T_A
            vals.append((Cw, L))
            print(f"  L={L}: T_M={T_M}, T_A={T_A}, C_w={Cw:.1f}")

        Cw_opt, L_opt = min(vals)
        Lopt_dict[d] = L_opt
        print(f"  ==> L_opt({d}) = {L_opt}, C_w(d,L_opt) = {Cw_opt:.1f}")

    return Lopt_dict

def IncompleteNTT_hybrid_complexity(n:int, t:int, L:int, w:float = 0.2):
    T_M, T_A = Hybrid_Toom4_complexity(n//(2**t), L)
    C_M = 3*((t*n)//2) + (2**t)*T_M
    if t >0:
        C_M += 2*n - 2**t
    C_A = (3*t+1)*n - 2**t + (2**t)*T_A
    C_w = C_M + w*C_A
    return C_M, C_A, C_w

# ============================================================
# Benchmark runner
# ============================================================

def benchmark_one_parameter_set(
    n: int,
    q: int,
    q_minus_1_factors: list[int],
    Lopt_dict: dict[int, int],
    trials: int = 100,
    w: float = 0.2,
    seed: int = 20260426,
    measure_schoolbook: bool = True,
):
    assert n > 0 and (n & (n - 1) == 0)
    random.seed(seed + n + q)
    lambda_n = int(math.log2(n))
    t_max = min(v2(q-1)-1,lambda_n)
    t = t_max
    d = n // (2**t)
    L_max = L_max_for_d(d)
    if d in Lopt_dict:
        L_opt = Lopt_dict[d]
        T_M_opt, T_A_opt = Hybrid_Toom4_complexity(d, L_opt)
        L_opt_cost = T_M_opt + w*T_A_opt
    else:
        L_opt, L_opt_cost = theoretical_L_opt(d, w=w)
    CM0, CA0, Cw0 = IncompleteNTT_hybrid_complexity(n, t, 0, w=w)
    CMmax, CAmax, Cwmax = IncompleteNTT_hybrid_complexity(n, t, L_max, w=w)
    CMopt, CAopt, Cwopt = IncompleteNTT_hybrid_complexity(n, t, L_opt, w=w)
        
    N = 2**(t+1)
    zeta = primitive_Nth_root_of_unity_mod_p(N,q,q_minus_1_factors)
    ntt_data = precompute_IncompleteNTT_data(N, q, t, zeta)

    print("=" * 72)
    print(f"n={n}, q={q}, v2(q-1)={v2(q-1)}, t_max={t_max}, d={d}")
    print(f"L_max={L_max}, L_opt={L_opt}, Cw_opt={L_opt_cost:.1f}")
    print(f"Cw: Karatsuba={Cw0:.1f}, Toom-4={Cwmax:.1f}, Hybrid={Cwopt:.1f}")
    print(f"N = 2^(t+1) = {N}")
    print(f"primitive N-th root zeta = {zeta}")
    print("-" * 72)

    strategies = []
    if measure_schoolbook:
        strategies.append(("Schoolbook", "schoolbook", -1))

    strategies.append(("Karatsuba", "hybrid", 0))
    strategies.append(("Toom-4", "hybrid", L_max))

    if L_opt != 0:
        strategies.append(("Hybrid", "hybrid", L_opt))

    # Generate random inputs and references once, then reuse them for all methods.
    inputs = []
    references = []

    for _ in range(trials):
        f = [random.randrange(q) for _ in range(n)]
        g = [random.randrange(q) for _ in range(n)]
        ref = conv_prod(n, q, f, g)
        inputs.append((f, g))
        references.append(ref)

    # Precompute inv_120L for each strategy.
    inv_120L_by_label = {}
    for label, method, L in strategies:
        if L == -1:
            inv_120L_by_label[label] = 1
        else:
            inv_120L_by_label[label] = exGCD(pow(120, L, q), q)

    # Warm-up: run each strategy a few times without measuring.
    warmup_trials = min(5, trials)
    for (f, g), ref in zip(inputs[:warmup_trials], references[:warmup_trials]):
        for label, method, L in strategies:
            out = conv_prod_IncompleteNTT_hybrid(
                n, q, t, d, L,
                inv_120L_by_label[label],
                ntt_data, f, g, method
            )
            assert out == ref, f"Correctness failed during warm-up for {label}."

    # Main benchmark: shuffle the order of strategies for each input.
    times_by_label = {label: [] for label, _, _ in strategies}

    for (f, g), ref in zip(inputs, references):
        order = strategies[:]
        random.shuffle(order)

        for label, method, L in order:
            start = time.perf_counter()
            out = conv_prod_IncompleteNTT_hybrid(
                n, q, t, d, L,
                inv_120L_by_label[label],
                ntt_data, f, g, method
            )
            elapsed = time.perf_counter() - start

            assert out == ref, f"Correctness failed for {label}."
            times_by_label[label].append(elapsed)

    results = {}
    for label, method, L in strategies:
        times = times_by_label[label]
        results[label] = {
            "L": L,
            "mean_ms": 1000.0 * mean(times),
            "median_ms": 1000.0 * median(times),
            "std_ms": 1000.0 * stdev(times) if len(times) >= 2 else 0.0,
        }

    # If L_opt=0, Hybrid is identical to Karatsuba, so we do not measure it separately.
    if L_opt == 0:
        results["Hybrid"] = results["Karatsuba"].copy()
        results["Hybrid"]["L"] = 0

    # Header
    header = ["method", "L", "mean", "median", "std"]
    labels = list(results.keys())
    print(f"{'metric':<10}", end="")
    for label in labels:
        print(f"{label:>12}", end="")
    print()
    
    # L
    print(f"{'L':<10}", end="")
    for label in labels:
        L = results[label]["L"]
        L_str = "-" if L == -1 else str(L)
        print(f"{L_str:>12}", end="")
    print()
    
    # mean
    print(f"{'mean':<10}", end="")
    for label in labels:
        print(f"{results[label]['mean_ms']:12.3f}", end="")
    print()

    # median
    print(f"{'median':<10}", end="")
    for label in labels:
        print(f"{results[label]['median_ms']:12.3f}", end="")
    print()
    
    # std
    print(f"{'std':<10}", end="")
    for label in labels:
        print(f"{results[label]['std_ms']:12.3f}", end="")
    print()

    print()

    return {
        "n": n,
        "q": q,
        "q-1 prime factors": q_minus_1_factors,
        "t_max": t_max,
        "t": t,
        "d": d,
        "L_max": L_max,
        "L_opt": L_opt,
        "theory": {
            "Karatsuba": Cw0,
            "Toom-4": Cwmax,
            "Hybrid": Cwopt,
        },
        "results": results,
    }

def main():
    d_values = [64, 128, 256, 512, 1024]
    Lopt_dict = print_theoretical_L_table(d_values, w=0.2)
    print()
    params = [
    # (n, q, prime factors of q-1)
    # 8380402 = 2 * 1583 * 2647
    # 8380380 = 2^2 * 3 * 5 * 197 * 709
    # 8380248 = 2^3 * 3 * 349177
    # 8380368 = 2^4 * 3^3 * 19 * 1021

    (256,  8380403, [2, 1583, 2647]),
    (256,  8380381, [2, 3, 5, 197, 709]),

    (512,  8380403, [2, 1583, 2647]),
    (512,  8380381, [2, 3, 5, 197, 709]),
    (512,  8380249, [2, 3, 349177]),

    (1024, 8380403, [2, 1583, 2647]),
    (1024, 8380381, [2, 3, 5, 197, 709]),
    (1024, 8380249, [2, 3, 349177]),
    (1024, 8380369, [2, 3, 19, 1021]),
    ]

    all_results = []

    for n, q, q_minus_1_factors in params:
        res = benchmark_one_parameter_set(
            n=n,
            q=q,
            q_minus_1_factors=q_minus_1_factors,
            Lopt_dict=Lopt_dict,
            trials=100,
            w=0.2,
            seed=20260426,
            measure_schoolbook=True,
        )
        all_results.append(res)
    
    print("\nLaTeX theoretical C_w(n;t,L) summary:")
    print("n & t & d & Karatsuba & Toom-4 & Hybrid \\\\")
    for res in all_results:
        th = res["theory"]
        print(
            f"{res['n']} & {res['t_max']} & {res['d']} & "
            f"{th['Karatsuba']:.1f} & {th['Toom-4']:.1f} & "
            f"{th['Hybrid']:.1f} ($L_{{\\rm opt}}={res['L_opt']}$) \\\\"
        )

    print("\nLaTeX summary:")
    print("n & q & t_max & d & Karatsuba & Toom-4 & Hybrid \\\\")
    for res in all_results:
        r = res["results"]

        def fmt(label):
            return f"{r[label]['median_ms']:.3f}"
        
        print(
            f"{res['n']} & {res['q']} & {res['t_max']} & {res['d']} & "
            f"{fmt('Karatsuba')} & {fmt('Toom-4')} & {fmt('Hybrid')} \\\\"
            )

if __name__ == "__main__":
    main()