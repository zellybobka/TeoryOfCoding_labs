import math
PRIM = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011}
def gf_tables(m):
    n = 2**m - 1
    prim = PRIM[m]
    exp = [0] * (2 * n)
    log = [0] * 2**m
    exp[0] = 1
    log[0] = -1
    log[1] = 0
    for i in range(1, n):
        x = exp[i - 1] << 1
        if x & (1 << m):
            x ^= prim
        exp[i] = x & ((1 << m) - 1)
        log[exp[i]] = i
    for i in range(n, 2 * n):
        exp[i] = exp[i - n]
    return exp, log, n
def gf_mul(a, b, exp, log, n):
    if a == 0 or b == 0: return 0
    return exp[(log[a] + log[b]) % n]
def gf_inv(a, exp, log, n):
    if a == 0: raise ZeroDivisionError("Обратный элемент для 0 не существует")
    return exp[-log[a]]
def trim(p):
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def pmul2(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                if b: r[i + j] ^= 1
    return trim(r)
def pdiv2(dividend, divisor):
    dd = trim(dividend[:])
    dv = trim(divisor[:])
    q = [0] * max(0, len(dd) - len(dv) + 1)
    while len(dd) >= len(dv) and dd != [0]:
        sh = len(dd) - len(dv)
        q[sh] = 1
        for i, c in enumerate(dv):
            if c: dd[i + sh] ^= 1
        dd = trim(dd)
    return trim(q), trim(dd)
def pstr(p):
    t = [f"{'1' if i==0 else ('x' if i==1 else f'x^{i}')}" for i, c in enumerate(p) if c]
    return "0" if not t else " + ".join(reversed(t))
def pmul_gf(p, q, exp, log, n):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                if b: r[i + j] ^= gf_mul(a, b, exp, log, n)
    return r
def coset(i, n):
    c, x = [], i % n
    while x not in c:
        c.append(x)
        x = (2 * x) % n
    return c
def minimal_poly(cos, exp, log, n):
    poly = [1]
    for e in cos:
        a = exp[e]
        poly = pmul_gf(poly, [a, 1], exp, log, n)
    for c in poly:
        if c not in (0, 1): raise ValueError("Минимальный многочлен не в GF(2).")
    return trim(poly)
def build_bch(n, d):
    m = int(round(math.log2(n + 1)))
    exp, log, _ = gf_tables(m)
    used = set()
    g = [1]
    info = []
    for i in range(1, d):
        if i in used: continue
        cs = coset(i, n)
        used.update(cs)
        mp = minimal_poly(cs, exp, log, n)
        g = pmul2(g, mp)
        info.append((i, cs, mp))
    xn1 = [1] + [0] * (n - 1) + [1]
    h, rem = pdiv2(xn1, g)
    if rem != [0]: raise RuntimeError("Ошибка: (x^n-1)/g(x) делится с остатком.")
    r = len(g) - 1
    k = n - r
    G = [[0] * n for _ in range(k)]
    for i in range(k):
        for j, c in enumerate(g):
            if c: G[i][i + j] = 1
    hs = list(reversed(h))
    H = [[0] * n for _ in range(r)]
    for i in range(r):
        for j, c in enumerate(hs):
            if c: H[i][i + j] = 1
            
    return exp, log, g, h, G, H, k, r, info
def encode(msg, g, n, k):
    r = n - k
    shifted = [0] * r + msg[:]
    rem = pdiv2(shifted, g)[1]
    c = shifted[:]
    for i in range(len(rem)): c[i] ^= rem[i]
    return c[:n]
def eval_bin_poly(bits, a, exp, log, n):
    res, p = 0, 1
    for b in bits:
        if b: res ^= p
        p = gf_mul(p, a, exp, log, n)
    return res
def syndromes(rbits, t, exp, log, n):
    return [eval_bin_poly(rbits, exp[i], exp, log, n) for i in range(1, 2 * t + 1)]
def solve_gf(A, b, exp, log, n):
    v = len(A)
    M = [A[i][:] + [b[i]] for i in range(v)]
    for col in range(v):
        piv = next((row for row in range(col, v) if M[row][col] != 0), None)
        if piv is None: return None
        M[col], M[piv] = M[piv], M[col]
        invp = gf_inv(M[col][col], exp, log, n)
        for j in range(col, v + 1): M[col][j] = gf_mul(M[col][j], invp, exp, log, n)
        for row in range(v):
            if row == col: continue
            f = M[row][col]
            if f:
                for j in range(col, v + 1): M[row][j] ^= gf_mul(f, M[col][j], exp, log, n)
    return [M[i][v] for i in range(v)]

def pgz_decode(rbits, n, d, exp, log):
    t = (d - 1) // 2
    if t == 0: return rbits[:], [], []
    S = syndromes(rbits, t, exp, log, n)
    if all(s == 0 for s in S): return rbits[:], [], S
    sigma = None
    for v in range(t, 0, -1):
        A, b = [], []
        ok = True
        for krow in range(1, v + 1):
            row = [S[v + krow - j - 1] for j in range(1, v + 1)]
            if any(not (1 <= (v + krow - j) <= 2 * t) for j in range(1, v + 1)) or \
               not (1 <= (v + krow) <= 2 * t): ok = False; break
            A.append(row)
            b.append(S[v + krow - 1])
        if not ok: continue
        sol = solve_gf(A, b, exp, log, n)
        if sol is not None:
            sigma = [1] + sol[:]
            break
    if sigma is None: return rbits[:], [], S
    def sigma_eval(z):
        acc, p = 0, 1
        for c in sigma:
            if c: acc ^= gf_mul(c, p, exp, log, n)
            p = gf_mul(p, z, exp, log, n)
        return acc
    err_pos = [j for j in range(n) if sigma_eval(exp[(n - j) % n]) == 0]
    corrected = rbits[:]
    for j in err_pos: corrected[j] ^= 1
    return corrected, err_pos, S
def read_int(prompt, okset=None):
    while True:
        try:
            x = int(input(prompt).strip())
            if okset and x not in okset: print(f"Можно только: {sorted(okset)}"); continue
            return x
        except: print("Нужно целое число.")
def read_bits(prompt, L):
    while True:
        s = input(prompt).strip().replace(" ", "")
        if len(s) != L or any(ch not in "01" for ch in s): print(f"Нужно ровно {L} бит (0/1)."); continue
        return [int(ch) for ch in s]
def print_mat(M, name):
    print(f"\n{name} ({len(M)}x{len(M[0])}):")
    for row in M: print(" ".join(map(str, row)))
if __name__ == "__main__":
    print("=== BCH-код (n=3/7/15) + кодирование + PGZ-декодирование ===")
    n = read_int("Введите n (3, 7 или 15): ", {3, 7, 15})
    d = read_int(f"Введите d (2..{n}): ")
    if d < 2 or d > n or d > 15: raise SystemExit("Некорректное d.")
    exp, log, g, h, G, H, k, r, info = build_bch(n, d)
    print(f"\nПараметры кода: n={n}, k={k}, r={r}, d(designed)={d}, t={(d-1)//2}")
    print(f"g(x): {pstr(g)}\nh(x): {pstr(h)}")
    print("\nЦиклотомические классы:")
    for i0, cs, mp in info: print(f"  i={i0}: {cs}  =>  M_{i0}(x) = {pstr(mp)}")
    print_mat(G, "G")
    print_mat(H, "H")
    print("\n--- КОДИРОВАНИЕ ---")
    msg = read_bits(f"Введите сообщение m (k={k} бит): ", k)
    code = encode(msg, g, n, k)
    print("\nКодовое слово:")
    print("  c0..c_{n-1}:", "".join(map(str, code)))
    print("  c_{n-1}..c0:", "".join(map(str, reversed(code))))
    print(f"Проверочные: {code[:r]}, Сообщение: {code[r:]}.")
    print("\n--- ДЕКОДИРОВАНИЕ (PGZ) ---")
    recv = read_bits(f"Принятое слово (n={n} бит): ", n)
    corr, pos, S = pgz_decode(recv, n, d, exp, log)
    print("\nСиндромы (S1..S_{2t}):")
    if S: print("  " + "  ".join(f"S{i+1}={S[i]}" for i in range(len(S))))
    else: print("  t=0.")
    if pos: print("\nПозиции ошибок:", pos)
    else: print("\nОшибок не найдено ИЛИ PGZ не справился.")
    print("\nИсправленное слово:")
    print("  c0..c_{n-1}:", "".join(map(str, corr)))
    print("  c_{n-1}..c0:", "".join(map(str, reversed(corr))))
    msg_hat = corr[r:r+k]
    print("\nДекодированное сообщение (k бит):", "".join(map(str, msg_hat)))
