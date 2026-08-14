# Migración a office PC — LEE ESTO

**Contexto (2026-08-14):** el PC de casa (`DESKTOP-25V4FD3`) tiene una falla de
hardware en RAM confirmada (Windows Memory Diagnostic, 2026-07-27) que ha
causado 5 pantallazos azules en 9 días. Va a salir a reparación pronto. Todo
el trabajo activo debe quedar accesible desde el PC de la oficina
(`octopus-office`, Linux, Tailscale IP 100.111.171.1:2222) antes de eso.

## Estado de este proyecto (GSUA-CSB toolbox)

- Repositorio git, remote `https://github.com/drojasd/GSUA-CSB.git`, rama
  `master` limpia (sin cambios sin confirmar) en este PC.
- **RESUELTO 2026-08-14:** office ya tiene `.git` real (clonado desde
  GitHub). Se comparó contra el contenido sincronizado por MATLAB Drive:
  124 diferencias reales — falta todo un paquete `python/` (puerto
  `gsua_csb`, ejemplos, experimentos de paper) que sí existe en el
  historial de git, y varias funciones `Functions/*.m` difieren del
  código pusheado. Se protegió el estado actual del disco de office con un
  commit nuevo (sin borrar ni sobrescribir nada) — reconciliar cuál versión
  de cada función es la vigente queda pendiente de revisión manual.

## Instrucciones para cualquier agente (Claude/Codex) trabajando aquí

1. Este proyecto vive bajo una carpeta con sufijo
   `.source-of-truth-20260718T095207` — es intencional, resultado de un fix
   de reubicación de MATLAB Drive (ver memoria del proyecto). No renombrar
   ni "limpiar" ese sufijo sin entender por qué existe.
2. **Commitea y haz push** — no dependas de la sincronización de MATLAB
   Drive para preservar historial, esa sincronización no incluye `.git`.
3. **Si trabajas en office:** confirma si ya existe `.git` en la copia local
   antes de asumir que tienes historial de commits.

_Este archivo fue generado automáticamente como parte de la migración pre-reparación._
