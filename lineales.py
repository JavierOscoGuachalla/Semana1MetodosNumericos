"""SIMTELE-RISK · Semana 1: errores transversales y sistemas lineales (Jacobi, GS, LU, CG, referencia)."""
import numpy as np

def errores(x, xref):
    e = x - xref
    nref = np.linalg.norm(xref)
    ea = np.linalg.norm(e)
    er = ea / nref if nref > 0 else float('nan')   # sin dividir por cero
    return ea, er, 100 * er

def residual(A, x, b): return np.linalg.norm(b - A @ x)

def jacobi(A, b, tol=1e-10, maxit=10000):
    x = np.zeros_like(b); D = np.diag(A); R = A - np.diagflat(D); hist = []
    for k in range(1, maxit + 1):
        x = (b - R @ x) / D; r = residual(A, x, b); hist.append((k, x.copy(), r))
        if r < tol * np.linalg.norm(b): return x, hist, True
    return x, hist, False

def gauss_seidel(A, b, tol=1e-10, maxit=10000):
    n = len(b); x = np.zeros(n); hist = []
    for k in range(1, maxit + 1):
        for i in range(n):
            x[i] = (b[i] - A[i, :i] @ x[:i] - A[i, i+1:] @ x[i+1:]) / A[i, i]
        r = residual(A, x, b); hist.append((k, x.copy(), r))
        if r < tol * np.linalg.norm(b): return x, hist, True
    return x, hist, False

def lu_doolittle(A):
    n = len(A); L = np.eye(n); U = A.astype(float).copy()   # sin pivoteo (A SPD)
    for k in range(n - 1):
        for i in range(k + 1, n):
            L[i, k] = U[i, k] / U[k, k]; U[i, k:] -= L[i, k] * U[k, k:]
    return L, U

def lu_solve(A, b):
    L, U = lu_doolittle(A); n = len(b); y = np.zeros(n); x = np.zeros(n)
    for i in range(n): y[i] = b[i] - L[i, :i] @ y[:i]
    for i in range(n - 1, -1, -1): x[i] = (y[i] - U[i, i+1:] @ x[i+1:]) / U[i, i]
    return x, L, U

def gradiente_conjugado(A, b, tol=1e-10, maxit=10000):
    x = np.zeros_like(b); r = b - A @ x; p = r.copy(); rs = r @ r; hist = []
    for k in range(1, maxit + 1):
        Ap = A @ p; a = rs / (p @ Ap); x = x + a * p; r = r - a * Ap
        rn = r @ r; hist.append((k, x.copy(), np.sqrt(rn)))
        if np.sqrt(rn) < tol * np.linalg.norm(b): return x, hist, True
        p = r + (rn / rs) * p; rs = rn
    return x, hist, False

def referencia(A, b): return np.linalg.solve(A, b)

# --- Matrices de los tres escenarios (todas SPD) ---
def caso_normal():
    A = np.array([[4,1,0,0],[1,4,1,0],[0,1,4,1],[0,0,1,4]], float)
    xs = np.array([1,2,3,4], float); return A, A @ xs, xs
def caso_estres(n=10):
    A = 2.2*np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
    xs = np.arange(1, n+1, dtype=float); return A, A @ xs, xs
def caso_mal(n=5):
    A = np.array([[1/(i+j+1) for j in range(n)] for i in range(n)]); xs = np.ones(n)
    return A, A @ xs, xs   # Hilbert: SPD, kappa ~ 4.8e5
