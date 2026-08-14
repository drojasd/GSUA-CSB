# MIGRATION — tarea pendiente para el agente que abra este proyecto

**A ti, agente (Claude/Codex) que estás leyendo esto:** este proyecto tiene
divergencia real entre home y office sin reconciliar. No es solo un aviso —
es una tarea: tu trabajo es **investigar y responder cómo fusionar
correctamente**, y luego ejecutarlo.

## Contexto

Home (`DESKTOP-25V4FD3`) tuvo una falla de RAM confirmada y salió a
reparación. La copia de office venía sincronizándose solo por archivos
(MATLAB Drive, symlink `~/MATLAB-Drive/Toolbox` →
`octopus-large/matlab-drive-source-of-truth-20260718T153305Z/Toolbox`), sin
`.git` real. Se clonó el historial real de GitHub en office y se comparó
contra el contenido en disco: **124 diferencias**. Se protegió el estado de
office con un commit propio en vez de forzar un merge.

- **Home/GitHub** (rama `master`): commit `8d8ea20`.
- **Office**: commit local (snapshot de disco tras conectar el historial).

## Qué está en conflicto

- **Falta todo un paquete `python/` completo** en el disco de office
  (`python/src/gsua_csb/`, `python/examples/`, `python/paper_experiments/`,
  `python/README.md`, `python/USERGUIDE.md`, `python/pyproject.toml`) que
  sí está en el historial de git de home. Esto es un puerto/adaptación
  Python del toolbox MATLAB — verifica si nunca se sincronizó a office o si
  se eliminó deliberadamente ahí.
- También aparecen como faltantes en office: `LICENSE`, `README.md`,
  `docs/gsua_userguide.html`, `docs/index.html`, y varios archivos bajo
  `tests/` (gsuaCostfNanTolerantTest.m, gsuaNoisefloor*Test.m, etc.) —
  revisa si esto es real o un artefacto de cómo MATLAB Drive sincroniza
  (algunos patrones de archivo pueden estar excluidos de su sync).
- **Funciones MATLAB modificadas de forma distinta en cada lado:**
  `Functions/gsua_costf.m`, `gsua_costfMulti.m`, `gsua_depth.m`,
  `gsua_likecost.m`, `gsua_likelihood.m`, `gsua_oatr2.m`, `gsua_pe.m`,
  `gsua_rcostf.m`, `gsua_timer.m`, `Examples/user_dependent.m`,
  `progressbar/progressbar.m`. Estas SÍ son código activo — no fusiones a
  ciegas.
- `.gitignore` también difiere (aparece como eliminado en office) — revisa
  si el `.gitignore` de home necesita actualizarse para reflejar qué debe
  excluirse de MATLAB Drive vs qué debe excluirse de git.

## Cómo investigar (pistas de dónde se envió qué)

1. Compara `git log --oneline -- Functions/gsua_pe.m` (y las demás
   funciones en conflicto) en ambos hosts — los mensajes de commit indican
   qué bug se corrigió o qué feature se agregó en cada versión.
2. Busca si hay documentación de sincronización específica de este
   proyecto (revisa `SSH_PROJECT_COORDINATION.md` y
   `OCTOPUS_TRANSFER.md`, que existen en la raíz de `~/MATLAB-Drive` en
   office — pueden explicar el protocolo usado para este toolbox
   específicamente).
3. Para el paquete `python/` faltante: confirma con Daniel si es trabajo
   activo que debe existir en office, o si fue un experimento solo de home
   que no necesita replicarse.
4. Pregunta a Daniel si tienes dudas sobre cuál versión de una función es
   la correcta antes de fusionar código que otros(as) puedan estar usando.

## Qué hacer al terminar

1. Resuelve la divergencia con criterio (funciones: comparar lógica real,
   no solo tomar la más reciente por fecha; `python/`: decidir si
   transferir o descartar).
2. Completa el merge y pushea a `origin/master`.
3. **Actualiza este archivo** con qué se decidió y por qué, o bórralo si ya
   no aplica.

_Generado 2026-08-14 como parte de la migración pre-reparación de hardware._
