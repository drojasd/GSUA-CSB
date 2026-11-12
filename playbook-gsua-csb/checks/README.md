# Checks GSUA-CSB

Verificación ejecutable del entorno y de las convenciones del paquete.

```matlab
cd playbook-gsua-csb/checks
ok = run_checks();      % true si todo pasa
```

`run_checks.m` comprueba:
1. La toolbox (`Functions/`) y las funciones clave están en el path.
2. El invariante de **factores fijos**: un factor con rango degenerado se reinserta en su posición y el
   modelo recibe el vector completo (ver `../references/lineamientos.md §1`).
3. El desempaque **CI-primero** de un modelo del curso vía `gsua_eval` (la condición inicial llega a la
   posición correcta; ver §2).

Para verificar el contrato de un modelo de usuario propio, usa la ayuda del curso
`course/00_Escribir_Modelos/matlab/verificar_modelo.m` (firma, forma/orientación de la salida, cobertura
del dominio, NaN/Inf). En Python, la verificación de paridad numérica está en
`course/00_Escribir_Modelos/pruebas/paridad.py`.
