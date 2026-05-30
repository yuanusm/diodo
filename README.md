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

El resultado se imprime en terminal y el gráfico se guarda como `ajuste_diodo_riguroso.png`.

## Si necesitas cambiar archivo, columnas o unidades

Este archivo ya no expone argumentos de línea de comandos. Para mantenerlo como ejecutable simple, edita directamente la sección `CONFIGURACIÓN DEL USUARIO` al inicio de `ajuste_diodo_riguroso.py`:

```python
CSV_PATH = Path("datos.csv")
VOLTAGE_COL = "Voltaje"
CURRENT_COL = "Corriente"
CURRENT_UNIT = "mA"
PLOT_PATH = Path("ajuste_diodo_riguroso.png")
GENERAR_GRAFICO = True
```

## Dependencias

Instala las dependencias con:

```bash
python -m pip install -r requirements.txt
```
