"""Genera docs/avance_semana1.pdf (fórmulas, ejemplo manual, matrices, tablas de error)."""
import sys, os; sys.path.insert(0, os.path.dirname(__file__))
from lineales import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4; from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle; from reportlab.lib.units import cm
from reportlab.platypus import *; from reportlab.pdfbase import pdfmetrics; from reportlab.pdfbase.ttfonts import TTFont
from reportlab.graphics.shapes import Drawing, Rect, String, Line
OUT = os.path.join(os.path.dirname(__file__), "..", "docs"); os.makedirs(OUT, exist_ok=True)
FD = "/usr/share/fonts/truetype/dejavu/"
pdfmetrics.registerFont(TTFont("DJ", FD+"DejaVuSans.ttf")); pdfmetrics.registerFont(TTFont("DJB", FD+"DejaVuSans-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DJM", FD+"DejaVuSansMono.ttf"))
pdfmetrics.registerFontFamily("DJ", normal="DJ", bold="DJB", italic="DJ", boldItalic="DJB")
B = ParagraphStyle("b", fontName="DJ", fontSize=9, leading=12.5, spaceAfter=5)
H1 = ParagraphStyle("h1", parent=B, fontName="DJB", fontSize=13, leading=16, spaceBefore=10, spaceAfter=6, textColor=colors.HexColor("#1F3A5F"))
H2 = ParagraphStyle("h2", parent=B, fontName="DJB", fontSize=10.5, spaceBefore=6, textColor=colors.HexColor("#1F3A5F"))
F = ParagraphStyle("f", parent=B, fontName="DJM", fontSize=8.5, leading=12, leftIndent=12, backColor=colors.HexColor("#F1F4F8"), borderPadding=4, spaceBefore=3, spaceAfter=8)
SM = ParagraphStyle("sm", parent=B, fontSize=7.5, leading=9.5, spaceAfter=0)
def P(t, s=B): return Paragraph(t, s)
def tabla(data, cw=None, fs=7.5, hdr=True):
    data = [[Paragraph(str(c), ParagraphStyle("c", parent=SM, fontSize=fs, leading=fs+2, fontName="DJB" if (hdr and i == 0) else "DJ")) for c in r] for i, r in enumerate(data)]
    t = Table(data, colWidths=cw, repeatRows=1)
    t.setStyle(TableStyle([("GRID", (0,0), (-1,-1), .4, colors.grey), ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#DCE6F2")), ("VALIGN", (0,0), (-1,-1), "MIDDLE"), ("TOPPADDING", (0,0), (-1,-1), 2), ("BOTTOMPADDING", (0,0), (-1,-1), 2)]))
    return t
def mat(M, d=4): return "[" + "; ".join("  ".join(f"{v:.{d}g}" for v in r) for r in M) + "]"
def vec(v, d=6): return "[" + ", ".join(f"{u:.{d}g}" for u in v) + "]"
def e(x): return f"{x:.3e}"

casos = [("NORMAL", caso_normal()), ("ESTRÉS", caso_estres()), ("MAL CONDICIONADO", caso_mal())]
res = {}
for nm, (A, b, xs) in casos:
    r = {"kappa": np.linalg.cond(A), "lmin": np.linalg.eigvalsh(A).min()}
    D = np.diag(np.diag(A)); BJ = np.eye(len(A)) - np.linalg.solve(D, A); r["rhoJ"] = max(abs(np.linalg.eigvals(BJ)))
    mi = 5000 if nm == "MAL CONDICIONADO" else 20000
    for mn, fn in [("Jacobi", jacobi), ("Gauss-Seidel", gauss_seidel), ("Grad. conjugado", gradiente_conjugado)]:
        with np.errstate(all="ignore"): x, h, c = fn(A, b, maxit=mi if mn != "Gauss-Seidel" else 200000)
        r[mn] = (x, h, c)
    r["LU"] = lu_solve(A, b)
    xr = referencia(A, b); r["ref"] = xr
    rng = np.random.default_rng(0); d = rng.standard_normal(len(b)); d *= 1e-6*np.linalg.norm(b)/np.linalg.norm(d)
    xp = referencia(A, b + d); r["sens"] = (np.linalg.norm(d)/np.linalg.norm(b), np.linalg.norm(xp-xr)/np.linalg.norm(xr))
    res[nm] = r

# ---- gráfico de convergencia ----
fig, ax = plt.subplots(1, 3, figsize=(10, 2.9))
for a, (nm, _) in zip(ax, casos):
    for mn, col in [("Jacobi", "tab:red"), ("Gauss-Seidel", "tab:blue"), ("Grad. conjugado", "tab:green")]:
        h = res[nm][mn][1][:600]; rr = [t[2] for t in h]
        rr = [v if np.isfinite(v) else np.nan for v in rr]
        a.semilogy(range(1, len(rr)+1), rr, label=mn, color=col)
    a.set_title(nm, fontsize=9); a.set_xlabel("iteración", fontsize=8); a.grid(alpha=.3); a.tick_params(labelsize=7)
ax[0].set_ylabel("‖b − A x_k‖₂", fontsize=8); ax[0].legend(fontsize=7); plt.tight_layout(); plt.savefig(os.path.join(OUT, "convergencia.png"), dpi=170); plt.close()

S = []
S += [P("SIMTELE-RISK · Avance de la Semana 1", ParagraphStyle("t", parent=H1, fontSize=18, leading=22)),
      P("Desafío Integrador Final — Métodos Numéricos · Modelos, teoría de errores y sistemas lineales"),
      P("Integrantes: ____________________ &nbsp;&nbsp; Repositorio: github.com/____/simtele-risk &nbsp;&nbsp; Fecha: ____/____/2026", SM),
      P("<i>Contexto académico: los escenarios son hipotéticos y no describen fallas verificadas de ninguna empresa.</i>", SM), Spacer(1, 6)]
S += [P("1. Planteamiento de los cuatro casos y relación modelo–método", H1),
 tabla([["Módulo", "Modelo matemático (variables, unidades)", "Métodos asignados", "Referencia analítica"],
  ["1 Detención", "dx/dt = v, dv/dt = −a(v,t). x [m], v [m/s], a [m/s²], t [s]. Hip.: freno constante (caso base) o variable (ampliación).", "RK2, RK4; Trapecio, Simpson, Romberg sobre v(t); Bisección, Newton, Secante para v(t)=0", "v=v₀−at; x=v₀t−at²/2; t_det=v₀/a; d=v₀²/(2a)"],
  ["2 Deterioro", "D(t)=D₀+αt+βt², D adimensional [0–1], t [meses]. Umbral D*. Ajuste de (D₀,α,β) por mínimos cuadrados: (VᵀV)c=Vᵀy (SPD).", "Bisección, Newton, Secante sobre D(t)−D*=0; Jacobi, GS, LU, CG para (VᵀV)c=Vᵀy; integración opcional de D′(t)", "t* = (−α+√(α²−4β(D₀−D*)))/(2β), raíz positiva"],
  ["3 Oscilaciones", "x′=v; v′=(F(t)−cv−kx)/m. m [kg], c [N·s/m], k [N/m], F [N]. Régimen subamortiguado (ζ&lt;1).", "RK2, RK4; Trapecio, Simpson, Romberg sobre ∫x² dt; sist. lineales para calibrar k, c", "x(t)=e^(−ζω₀t)[C₁cos ω_d t + C₂ sin ω_d t]"],
  ["4 Sensores", "Señal y(t) muestreada en nodos t_i con ruido ε~N(0,σ²); reconstrucción en puntos ocultos.", "Lagrange, Newton (dif. divididas), spline cúbico; sist. lineales de calibración; raíz de cruce de umbral", "y=sin t, y=e^(−t) (error real exacto)"]],
  cw=[2*cm, 5.6*cm, 5.2*cm, 4.2*cm]),
 Spacer(1, 6), P("<b>Escenarios definidos por módulo</b> (parámetros nominales sintéticos; la fuente de cada parámetro es «supuesto del equipo» salvo indicación):", B),
 tabla([["Módulo", "NORMAL", "ESTRÉS", "MAL CONDICIONADO (subexperimento numérico)"],
  ["1", "v₀=5 m/s, a=1 m/s², h=0.05 s → t_det=5 s, d=12.5 m", "v₀=6.5 m/s, a=0.6 m/s², reacción 1.5 s; h=0.5 s vs 0.01 s", "v(t)=v₀(1−t/T)²: raíz doble, Newton converge linealmente y no hay cambio de signo para bisección"],
  ["2", "D₀=0.10, α=0.02, β=0.001, D*=1 → t*≈21.6 meses", "β×3, α±10 %, mantenimiento retrasado 6 meses", "Ajuste con matriz de Hilbert (κ≈4.8×10⁵) y discriminante≈0 (raíz casi doble)"],
  ["3", "m=1, k=4, c=0.4 (ζ=0.1), F=0, x₀=0.1; h=0.01", "F=F₀sin(ωt) con ω≈ω₀, h=0.1 (RK2 vs RK4), T mayor", "c²≈4mk (ζ→1, raíces características casi iguales); calibración (k,c) casi singular"],
  ["4", "y=sin t, 9 nodos equiespaciados en [0,2π], 3 puntos ocultos", "Ruido σ=0.05; 5 vs 15 nodos", "20 nodos equiespaciados de e^(−t) con ruido 10⁻³ (Vandermonde muy mal condicionada / Runge)"]],
  cw=[1.2*cm, 5*cm, 5*cm, 5.8*cm])]

S += [P("2. Teoría de errores transversal", H1),
 P("Sea x* la referencia (exacta o numérica de alta precisión justificada) y x̃ la aproximación. Se calculan siempre:"),
 P("Error absoluto: E<sub>a</sub> = ‖x̃ − x*‖₂ &nbsp;&nbsp;|&nbsp;&nbsp; Error relativo: E<sub>r</sub> = E<sub>a</sub>/‖x*‖₂ (solo si ‖x*‖≠0; si no, se informa solo E<sub>a</sub>) &nbsp;&nbsp;|&nbsp;&nbsp; Porcentual: E<sub>%</sub> = 100·E<sub>r</sub>", F),
 P("Residual: r = b − A x̃ ; se usa también para raíces: |f(x̃)|. <b>Residuo pequeño no implica error pequeño:</b> ‖x̃−x*‖/‖x*‖ ≤ κ(A)·‖r‖/‖b‖ con κ(A)=λ<sub>máx</sub>/λ<sub>mín</sub> para A SPD."),
 P("Sensibilidad: perturbando b→b+δb, ‖δx‖/‖x‖ ≤ κ(A)·‖δb‖/‖b‖. Se compara el cociente observado contra el factor κ(A) (cota superior). Criterio de parada de iterativos: ‖b−Ax<sub>k</sub>‖₂ &lt; 10⁻¹⁰‖b‖₂ o k=k<sub>máx</sub>; se distingue <b>no convergencia</b> (ρ≥1 o k<sub>máx</sub>) de <b>mal condicionamiento</b> (convergencia con error grande respecto del residuo)."),
 P("3. Sistemas lineales: fórmulas y condiciones de convergencia", H1)]
S += [tabla([["Método", "Fórmula", "Condición / observación"],
 ["Jacobi", "x<sub>i</sub><super>(k+1)</super> = (b<sub>i</sub> − Σ<sub>j≠i</sub> a<sub>ij</sub>x<sub>j</sub><super>(k)</super>)/a<sub>ii</sub>", "Converge ⇔ ρ(B<sub>J</sub>)&lt;1, B<sub>J</sub>=I−D⁻¹A. Suficiente: dominancia diagonal estricta."],
 ["Gauss-Seidel", "x<sub>i</sub><super>(k+1)</super> = (b<sub>i</sub> − Σ<sub>j&lt;i</sub> a<sub>ij</sub>x<sub>j</sub><super>(k+1)</super> − Σ<sub>j&gt;i</sub> a<sub>ij</sub>x<sub>j</sub><super>(k)</super>)/a<sub>ii</sub>", "Converge siempre si A es SPD (o diag. dominante). Típicamente ≈2× más rápido que Jacobi."],
 ["LU (Doolittle)", "A=LU, l<sub>ik</sub>=u<sub>ik</sub>/u<sub>kk</sub>; Ly=b (sustitución adelante), Ux=y (atrás)", "Directo, O(n³/3). Sin pivoteo es válido para A SPD (pivotes&gt;0)."],
 ["Grad. conjugado", "α<sub>k</sub>=r<sub>k</sub>ᵀr<sub>k</sub>/p<sub>k</sub>ᵀAp<sub>k</sub>; x<sub>k+1</sub>=x<sub>k</sub>+α<sub>k</sub>p<sub>k</sub>; r<sub>k+1</sub>=r<sub>k</sub>−α<sub>k</sub>Ap<sub>k</sub>; β<sub>k</sub>=r<sub>k+1</sub>ᵀr<sub>k+1</sub>/r<sub>k</sub>ᵀr<sub>k</sub>; p<sub>k+1</sub>=r<sub>k+1</sub>+β<sub>k</sub>p<sub>k</sub>", "Exige A SPD (todas las matrices usadas lo son: λ<sub>mín</sub>&gt;0). En aritmética exacta termina en ≤n pasos; cota: ‖e<sub>k</sub>‖<sub>A</sub>≤2((√κ−1)/(√κ+1))<super>k</super>‖e<sub>0</sub>‖<sub>A</sub>."],
 ["Referencia", "x* = A⁻¹b (LAPACK, numpy.linalg.solve) y x* conocido por construcción (b=Ax*)", "Se verifica que ‖x<sub>ref</sub>−x*‖ ≈ ε<sub>máq</sub>·κ."]], cw=[2.6*cm, 8.4*cm, 6*cm])]

# ---- ejemplo manual ----
A, b, xs = caso_normal(); L, U = lu_doolittle(A); xcg, hcg, _ = gradiente_conjugado(A, b)
S += [P("4. Ejemplo manual (escenario NORMAL, n=4, x₀=0)", H1),
 P(f"A = {mat(A,2)}, &nbsp; x* = {vec(xs)}, &nbsp; b = Ax* = {vec(b)}. Autovalor mínimo λ<sub>mín</sub>={res['NORMAL']['lmin']:.4f}&gt;0 ⇒ SPD."),
 P("<b>Jacobi, k=1:</b> x<sup>(1)</sup> = D⁻¹b = (6/4, 12/4, 18/4, 19/4) = (1.5, 3, 4.5, 4.75). <b>k=2:</b> x₁=(6−3)/4=0.75; x₂=(12−1.5−4.5)/4=1.5; x₃=(18−3−4.75)/4=2.5625; x₄=(19−4.5)/4=3.625.".replace("<sup>","<super>").replace("</sup>","</super>")),
 P("<b>Gauss-Seidel, k=1:</b> x₁=6/4=1.5; x₂=(12−1.5)/4=2.625; x₃=(18−2.625)/4=3.84375; x₄=(19−3.84375)/4=3.7890625. Ya usa valores nuevos: error tras 1 paso menor que Jacobi."),
 P("<b>LU:</b> l₂₁=1/4=0.25, u₂₂=4−0.25=3.75; l₃₂=1/3.75=0.2667, u₃₃=4−0.2667=3.7333; l₄₃=1/3.7333=0.2679, u₄₄=4−0.2679=3.7321. Resultado del código:"),
 P(f"L = {mat(L,4)}<br/>U = {mat(U,4)}<br/>Ly=b → y={vec(np.linalg.solve(L,b),5)} ; Ux=y → x={vec(res['NORMAL']['LU'][0],6)}", F),
 P(f"<b>Gradiente conjugado, k=1:</b> r₀=b=(6,12,18,19), p₀=r₀; r₀ᵀr₀={b@b:.0f}; Ap₀={vec(A@b)}; p₀ᵀAp₀={b@(A@b):.0f}; α₀={(b@b)/(b@(A@b)):.6f}; x₁=α₀p₀={vec((b@b)/(b@(A@b))*b,6)}. Tras k=4 (=n) el residuo es {hcg[3][2]:.2e}: terminación finita."),
 P("5. Matrices y soluciones de referencia (SPD)", H1)]
rows = [["Escenario", "n", "A (resumen)", "x* (conocido)", "κ₂(A)", "λ<sub>mín</sub>", "ρ(B<sub>J</sub>)"]]
desc = {"NORMAL": "tridiag(1,4,1)", "ESTRÉS": "tridiag(−1,2.2,−1)", "MAL CONDICIONADO": "Hilbert: a<sub>ij</sub>=1/(i+j−1)"}
for nm, (A_, b_, xs_) in casos:
    r = res[nm]; rows.append([nm, len(A_), desc[nm], vec(xs_, 3) if len(xs_) < 6 else "(1,2,…,10)", e(r["kappa"]), e(r["lmin"]), f"{r['rhoJ']:.4f}"])
S += [tabla(rows, cw=[2.6*cm, .8*cm, 3.6*cm, 3.2*cm, 2.3*cm, 2.2*cm, 2*cm]),
 P("Para NORMAL: b = (6, 12, 18, 19). La referencia directa reproduce x* con error ~10⁻¹⁶. En MAL CONDICIONADO ρ(B<sub>J</sub>)&gt;1: Jacobi <b>diverge</b>, fenómeno distinto del mal condicionamiento (que afecta a todos los métodos).", SM),
 P("6. Tablas de error por método y escenario", H1)]
rows = [["Escenario", "Método", "Iter.", "Convergió", "E<sub>a</sub>", "E<sub>r</sub>", "E<sub>%</sub>", "‖r‖₂"]]
for nm, (A_, b_, xs_) in casos:
    r = res[nm]
    for mn in ["Jacobi", "Gauss-Seidel", "Grad. conjugado"]:
        x, h, c = r[mn]
        if not np.all(np.isfinite(x)): rows.append([nm, mn, f">{len(h)}", "No (diverge)", "∞/NaN", "∞/NaN", "∞/NaN", "∞/NaN"]); continue
        ea, er, ep = errores(x, xs_); rows.append([nm, mn, len(h), "Sí" if c else "No", e(ea), e(er), e(ep), e(residual(A_, x, b_))])
    x = r["LU"][0]; ea, er, ep = errores(x, xs_); rows.append([nm, "LU", "directo", "Sí", e(ea), e(er), e(ep), e(residual(A_, x, b_))])
    ea, er, ep = errores(r["ref"], xs_); rows.append([nm, "Referencia (solve)", "directo", "Sí", e(ea), e(er), e(ep), e(residual(A_, r["ref"], b_))])
S += [tabla(rows, cw=[2.6*cm, 3*cm, 1.4*cm, 2*cm, 2*cm, 2*cm, 2*cm, 2*cm]), Spacer(1, 6),
 P("Sensibilidad ante perturbación relativa de b (≈10⁻⁶, semilla fija 0):"),
 tabla([["Escenario", "‖δb‖/‖b‖", "‖δx‖/‖x‖ observado", "Cota κ·‖δb‖/‖b‖", "Amplificación observada"]] +
   [[nm, e(res[nm]["sens"][0]), e(res[nm]["sens"][1]), e(res[nm]["kappa"]*res[nm]["sens"][0]), f"{res[nm]['sens'][1]/res[nm]['sens'][0]:.2f}"] for nm, _ in casos], cw=[3.5*cm, 3*cm, 3.8*cm, 3.5*cm, 3.5*cm]),
 Spacer(1, 6), Image(os.path.join(OUT, "convergencia.png"), width=17*cm, height=4.9*cm),
 P("Figura 1. Residual vs iteración (primeras 600 iteraciones; escala log).", SM)]

S += [P("7. Conclusiones de la semana", H1),
 P("<b>Conclusión numérica.</b> En NORMAL y ESTRÉS los cuatro métodos coinciden con la referencia; CG converge en ≤n iteraciones, GS necesita ≈½ de las iteraciones de Jacobi y LU entrega precisión de máquina sin iterar. En el escenario mal condicionado (Hilbert, κ≈4.8×10⁵) Jacobi diverge (ρ&gt;1), Gauss-Seidel necesita más de 10⁵ iteraciones y, aun con residual ~10⁻¹⁰, queda con error relativo ~3×10⁻⁵ (residuo ≠ error); CG y LU conservan entre 10⁻¹³ y 10⁻¹¹ de error relativo, coherente con la pérdida de ≈log₁₀κ≈5.7 cifras. La amplificación observada de perturbaciones está acotada por κ."),
 P("<b>Relación con el caso.</b> Estos solvers calibrarán los parámetros del Módulo 2 (ajuste D(t)) y, de forma auxiliar, los de los Módulos 3 y 4. Un error alto en los coeficientes se propaga al tiempo estimado de inspección t*, por lo que una calibración mal condicionada no debe interpretarse como una predicción confiable del deterioro. Nada de esto describe fallas reales: son escenarios sintéticos."),
 P("8. Modelo de datos inicial (DER)", H1)]

def box(d, x, y, w, name, fields):
    h = 20 + 10*len(fields); d.add(Rect(x, y-h, w, h, fillColor=colors.HexColor("#F1F4F8"), strokeColor=colors.HexColor("#1F3A5F")))
    d.add(Rect(x, y-14, w, 14, fillColor=colors.HexColor("#1F3A5F"), strokeColor=colors.HexColor("#1F3A5F")))
    d.add(String(x+4, y-10.5, name, fontName="DJB", fontSize=7.5, fillColor=colors.white))
    for i, f in enumerate(fields): d.add(String(x+4, y-24-10*i, f, fontName="DJ", fontSize=6.3))
    return (x, y, w, h)
d = Drawing(480, 330); W = 135
c1 = box(d, 5, 325, W, "caso_estudio", ["PK id", "nombre", "modelo_version"])
c2 = box(d, 172, 325, W, "escenario", ["PK id", "FK caso_id", "tipo (NORMAL/ESTRES/MAL_COND)", "parametros_json", "fuente_datos"])
c3 = box(d, 340, 325, W, "simulacion", ["PK id", "FK escenario_id", "fecha", "equipo (opcional)", "estado"])
c4 = box(d, 340, 205, W, "ejecucion_metodo", ["PK id", "FK simulacion_id", "metodo", "tolerancia, paso, max_iter", "convergio, tiempo_ms"])
c5 = box(d, 172, 205, W, "resultado_resumen", ["PK id", "FK ejecucion_id", "referencia", "valor_json, metrica_error_json", "interpretacion"])
c6 = box(d, 340, 95, W, "iteracion_resultado", ["PK id", "FK ejecucion_id", "numero", "aproximacion_json", "residual, error_abs, error_rel"])
def ln(x1, y1, x2, y2, a, b_):
    d.add(Line(x1, y1, x2, y2, strokeColor=colors.HexColor("#1F3A5F"))); d.add(String(x1+(4 if x2 >= x1 else -10), y1+3, a, fontSize=7)); d.add(String(x2-(10 if x2 >= x1 else -4), y2+3, b_, fontSize=7))
ln(5+W, 295, 172, 295, "1", "N"); ln(172+W, 295, 340, 295, "1", "N")
d.add(Line(407, 325-14-10*5-4, 407, 205, strokeColor=colors.HexColor("#1F3A5F"))); d.add(String(410, 200+0, "1", fontSize=7)); d.add(String(410, 325-14-10*5-14, "N", fontSize=7)) if False else None
d.add(String(410, 250, "1 — N", fontSize=7))
ln(340, 170, 172+W, 170, "1", "N") if False else d.add(Line(340, 170, 172+W, 170, strokeColor=colors.HexColor("#1F3A5F")))
d.add(String(318, 173, "1", fontSize=7)); d.add(String(176+W-8, 173, "N", fontSize=7))
d.add(Line(407, 205-14-10*5-4, 407, 95, strokeColor=colors.HexColor("#1F3A5F"))); d.add(String(410, 125, "1 — N", fontSize=7))
S += [d, P("Figura 2. DER inicial (notación pata de gallo simplificada: 1 — N). Modelo relacional y script SQL en <font name='DJM'>db/schema.sql</font>.", SM),
 P("Prototipo web: <font name='DJM'>index.html</font> con cuatro módulos visibles (Detención, Deterioro, Oscilaciones, Sensores), selector NORMAL/ESTRÉS/MAL CONDICIONADO, ejecución de los cuatro solvers lineales con tabla de errores, gráfico de convergencia e historial (en este prototipo se guarda en el navegador; en la Semana 2 se migra a la base SQL del DER)."),
 P("Cronograma siguiente: S2 → raíces, interpolación e integración + persistencia SQL; S3 → EDO (RK2/RK4), integración de módulos, informe y prueba ejecutar→guardar→recuperar→repetir.", SM)]
def pie(c, doc): c.setFont("DJ", 7); c.drawString(2*cm, 1.2*cm, "SIMTELE-RISK | Proyecto final de Métodos Numéricos | Avance Semana 1"); c.drawRightString(A4[0]-2*cm, 1.2*cm, f"{doc.page}")
SimpleDocTemplate(os.path.join(OUT, "avance_semana1.pdf"), pagesize=A4, leftMargin=2*cm, rightMargin=2*cm, topMargin=1.8*cm, bottomMargin=2*cm, title="SIMTELE-RISK Avance Semana 1").build(S, onFirstPage=pie, onLaterPages=pie)
print("ok")
