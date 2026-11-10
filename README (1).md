# SIMTELE-RISK

Desafío Integrador Final de **Métodos Numéricos**: aplicación web con cuatro casos de estudio inspirados en un sistema de transporte por cable (La Paz). Los escenarios son **hipotéticos** y no describen fallas verificadas de ninguna empresa.

## Módulos
1. **Detención de emergencia** — RK2/RK4, integración (Trapecio, Simpson, Romberg), raíces.
2. **Deterioro y mantenimiento** — raíces, sistemas lineales (Jacobi, Gauss-Seidel, LU, CG).
3. **Oscilaciones por viento** — RK2/RK4, integración.
4. **Reconstrucción de sensores** — Lagrange, Newton, splines cúbicos.

Cada módulo corre los escenarios **NORMAL, ESTRÉS y MAL CONDICIONADO** y aplica teoría de errores (absoluto, relativo, porcentual, residual, sensibilidad).

## Estado (Semana 1)
- [x] Prototipo web con 4 módulos visibles (`index.html`)
- [x] Teoría de errores + sistemas lineales (JS en el prototipo, Python en `python/lineales.py`)
- [x] DER y SQL (`db/DER.md`, `db/schema.sql`)
- [x] Avance PDF (`docs/avance_semana1.pdf`, generado con `python/generar_avance.py`)
- [ ] Semana 2: raíces, interpolación, integración, persistencia SQL
- [ ] Semana 3: EDO, integración final, prueba ejecutar → guardar → recuperar → repetir

## Cómo ejecutar
```bash
# Prototipo web: abrir index.html en el navegador (o)
python3 -m http.server 8000      # y visitar http://localhost:8000

# Regenerar el PDF de avance (requiere numpy, matplotlib, reportlab)
pip install numpy matplotlib reportlab
python3 python/generar_avance.py
```

## Estructura
```
index.html            prototipo web (4 módulos, solvers lineales, historial)
python/lineales.py    Jacobi, Gauss-Seidel, LU, CG, referencia, errores
python/generar_avance.py
db/schema.sql         modelo relacional
db/DER.md             diagrama entidad–relación (Mermaid)
docs/                 avance PDF
```

## Nota de persistencia
En el prototipo el historial se guarda en `localStorage` (sobrevive a recargas). En la Semana 2 se migra a la base SQL definida en `db/schema.sql`.

## Equipo
Integrantes: _completar_
