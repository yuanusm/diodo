# Ajuste riguroso de un diodo

Este repositorio contiene `ajuste_diodo_riguroso.py`, un ejecutable directo para ajustar curvas I-V de diodos comparando modelos físicos alternativos.

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

El ejecutable compara automáticamente:

1. `M0`: offset instrumental.
2. `M1`: Shockley ideal.
3. `M2`: Shockley con resistencia serie `Rs`.
4. `M3`: Shockley con resistencia paralela `Rp`.
5. `M4`: Shockley con `Rs` y `Rp`.

Los modelos se ordenan por BIC, para penalizar parámetros extra que no estén justificados por los datos.

## Uso directo

Coloca un archivo `datos.csv` en esta carpeta. Debe tener, por defecto, columnas `Voltaje` y `Corriente`; la corriente se interpreta en `mA` para ser compatible con el script original.

Ejecuta:

```bash
./ajuste_diodo_riguroso.py
```

También puedes ejecutarlo con Python si tu sistema no respeta el shebang:

```bash
python ajuste_diodo_riguroso.py
```

El resultado se imprime en terminal, el CSV teórico se guarda como `valores_teoricos_modelos_20uA_100uA.csv` y el gráfico se guarda como `ajuste_diodo_riguroso.png`.

## Si necesitas cambiar archivo, columnas o unidades

Este archivo ya no expone argumentos de línea de comandos. Para mantenerlo como ejecutable simple, edita directamente la sección `CONFIGURACIÓN DEL USUARIO` al inicio de `ajuste_diodo_riguroso.py`:

```python
CSV_PATH = Path("datos.csv")
VOLTAGE_COL = "Voltaje"
CURRENT_COL = "Corriente"
CURRENT_UNIT = "mA"
PLOT_PATH = Path("ajuste_diodo_riguroso.png")
GENERAR_GRAFICO = True
EXPORT_THEORY_CSV = Path("valores_teoricos_modelos_20uA_100uA.csv")
CORRIENTE_MIN_EXPORT_A = 20e-6
CORRIENTE_MAX_EXPORT_A = 100e-6
PUNTOS_EXPORT = 200
```

## CSV teórico exportado

Al finalizar el ajuste, el ejecutable exporta un CSV con una malla de corrientes objetivo entre `20 µA` y `100 µA`. Para cada corriente, incluye el voltaje teórico invertido para cada modelo ajustado (`M1` a `M4`); el modelo `M0` de offset instrumental queda como `NaN` porque no define una curva I-V invertible.

Columnas principales:

- `Corriente_objetivo_A`
- `Corriente_objetivo_uA`
- `m0_Voltaje_teorico_V`
- `m1_Voltaje_teorico_V`
- `m2_Voltaje_teorico_V`
- `m3_Voltaje_teorico_V`
- `m4_Voltaje_teorico_V`

## Dependencias

Instala las dependencias con:

```bash
python -m pip install -r requirements.txt
```
