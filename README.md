# Ajuste riguroso de un diodo

Este repositorio contiene `ajuste_diodo_riguroso.py`, un script para ajustar curvas I-V de diodos comparando modelos físicos alternativos.

## Por qué no basta con `V = nVt ln(I/Is + 1) + I Rs + I Rp`

La resistencia serie y la resistencia paralela no se incorporan sumando dos caídas `I*R` a la tensión del diodo:

- `Rs` está en serie y reduce la tensión interna del diodo: `Vd = V - I*Rs`.
- `Rp` está en paralelo con la juntura, por lo que aporta una corriente de fuga `Vd/Rp`.
- El modelo completo se plantea como una ecuación implícita:

```text
I = Is * (exp((V - I*Rs)/(n*Vt)) - 1) + (V - I*Rs)/Rp + Ioff
```

Por eso el script ajusta `I(V)` en lugar de forzar un `V(I)` incorrecto para todos los casos.

## Modelos comparados

El script compara automáticamente:

1. `M0`: offset instrumental.
2. `M1`: Shockley ideal.
3. `M2`: Shockley con resistencia serie `Rs`.
4. `M3`: Shockley con resistencia paralela `Rp`.
5. `M4`: Shockley con `Rs` y `Rp`.

Los modelos se ordenan por BIC, para penalizar parámetros extra que no estén justificados por los datos.

## Uso

El CSV debe tener, por defecto, columnas `Voltaje` y `Corriente`; la corriente se interpreta en `mA` para ser compatible con el script original.

```bash
python ajuste_diodo_riguroso.py --csv datos.csv
```

Opciones útiles:

```bash
python ajuste_diodo_riguroso.py \
  --csv datos.csv \
  --voltage-col Voltaje \
  --current-col Corriente \
  --current-unit mA \
  --plot ajuste_diodo_riguroso.png
```

Si solo quieres imprimir parámetros y métricas sin abrir/generar figuras:

```bash
python ajuste_diodo_riguroso.py --csv datos.csv --no-plot
```

## Dependencias

Instala las dependencias con:

```bash
python -m pip install -r requirements.txt
```
