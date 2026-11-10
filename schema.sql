-- SIMTELE-RISK · Modelo relacional inicial (SQLite / PostgreSQL compatible)
PRAGMA foreign_keys = ON;

CREATE TABLE caso_estudio (
  id             INTEGER PRIMARY KEY,
  nombre         TEXT NOT NULL UNIQUE,          -- Detención, Deterioro, Oscilaciones, Sensores
  modelo_version TEXT NOT NULL DEFAULT '1.0'
);

CREATE TABLE escenario (
  id             INTEGER PRIMARY KEY,
  caso_id        INTEGER NOT NULL REFERENCES caso_estudio(id) ON DELETE CASCADE,
  tipo           TEXT NOT NULL CHECK (tipo IN ('NORMAL','ESTRES','MAL_CONDICIONADO')),
  parametros_json TEXT NOT NULL,                -- parámetros con unidades, reproducibles
  fuente_datos   TEXT NOT NULL DEFAULT 'supuesto del equipo',
  UNIQUE (caso_id, tipo)
);

CREATE TABLE simulacion (
  id           INTEGER PRIMARY KEY,
  escenario_id INTEGER NOT NULL REFERENCES escenario(id) ON DELETE CASCADE,
  fecha        TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
  equipo       TEXT,
  estado       TEXT NOT NULL DEFAULT 'COMPLETA' CHECK (estado IN ('PENDIENTE','COMPLETA','ERROR'))
);

CREATE TABLE ejecucion_metodo (
  id            INTEGER PRIMARY KEY,
  simulacion_id INTEGER NOT NULL REFERENCES simulacion(id) ON DELETE CASCADE,
  metodo        TEXT NOT NULL,                  -- JACOBI, GAUSS_SEIDEL, LU, CG, REFERENCIA, ...
  tolerancia    REAL,
  paso          REAL,
  max_iter      INTEGER,
  convergio     INTEGER NOT NULL CHECK (convergio IN (0,1)),
  tiempo_ms     REAL
);

CREATE TABLE iteracion_resultado (
  id             INTEGER PRIMARY KEY,
  ejecucion_id   INTEGER NOT NULL REFERENCES ejecucion_metodo(id) ON DELETE CASCADE,
  numero         INTEGER NOT NULL,
  aproximacion_json TEXT NOT NULL,
  residual       REAL,
  error_abs      REAL,
  error_rel      REAL,                          -- NULL si la referencia es 0 (sin dividir por cero)
  UNIQUE (ejecucion_id, numero)
);

CREATE TABLE resultado_resumen (
  id             INTEGER PRIMARY KEY,
  ejecucion_id   INTEGER NOT NULL REFERENCES ejecucion_metodo(id) ON DELETE CASCADE,
  referencia     TEXT,
  valor_json     TEXT,
  metrica_error_json TEXT,
  interpretacion TEXT
);

CREATE INDEX idx_esc_caso   ON escenario(caso_id);
CREATE INDEX idx_sim_esc    ON simulacion(escenario_id);
CREATE INDEX idx_ej_sim     ON ejecucion_metodo(simulacion_id);
CREATE INDEX idx_it_ej      ON iteracion_resultado(ejecucion_id);

INSERT INTO caso_estudio(nombre) VALUES ('Detención'),('Deterioro'),('Oscilaciones'),('Sensores');
