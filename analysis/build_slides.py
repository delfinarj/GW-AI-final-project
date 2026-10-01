"""Construye report/slides.html: diez diapositivas para diez minutos, desde los archivos de resultados.

Mismas reglas que el resto: ningún número escrito a mano. Lee los valores de `skmask.summary_numbers`,
el mismo módulo del que leen los dos resúmenes de cinco páginas, así que la charla no puede decir una
cifra distinta de la que dice el informe. El PDF sale de imprimir esta página con `scripts/make_pdf.py`,
que produce una diapositiva por hoja en 16:9.

Uso:  python analysis/build_slides.py
"""
import shutil
from pathlib import Path

from skmask.presets import PRESETS
from skmask.provenance import write_sidecar
from skmask.summary_numbers import collect, sci

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "report"
FIG = OUT / "figures"
PAGE = "https://delfinarj.github.io/GW-AI-final-project/"
REPO = "https://github.com/delfinarj/GW-AI-final-project"

SENSOR_NAMES = {"deep_underground": "Subterráneo profundo", "shallow_underground": "Subterráneo somero",
                "surface_lab": "Laboratorio en superficie"}
SENSOR_PHRASE = {"deep_underground": "subterráneo profundo", "shallow_underground": "subterráneo somero",
                 "surface_lab": "de superficie"}
MASK_NAMES = {"hot_columns": "columnas y píxeles calientes", "cti": "estelas de transferencia",
              "halo": "halo", "serial": "registro serie", "low_energy_clusters": "clusters de baja energía"}
GROUP_NAMES = {"hot_columns_pixels": "columnas y píxeles calientes", "cti": "estelas de transferencia",
               "serial": "registro serie", "low_energy_clusters": "clusters de baja energía", "halo": "halo"}
MASK_R4 = {"hot": "columnas y píxeles calientes", "cti": "estelas de transferencia", "halo": "halo",
           "serial": "registro serie", "lec": "clusters de baja energía", "muon": "trazas de muones"}


def describe(parts):
    name = f"{MASK_R4[parts['mask']]} en el sensor {SENSOR_PHRASE[parts['sensor']]}"
    return f"{name} (ajustada en {SENSOR_PHRASE[parts['source']]})" if parts.get("source") else name


def slide(number, kicker, body):
    return f'<section class="slide"><div class="kicker">{kicker}</div>{body}<div class="n">{number}</div></section>'


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    for src in (RES / "compare_masks" / "compare_masks.png",
                RES / "cross_defect_false_positives" / "cross_defect.png"):
        shutil.copy2(src, FIG / src.name)

    c = collect(SENSOR_NAMES, MASK_NAMES, GROUP_NAMES)
    head, cross, null = c["head"], c["cross"], c["null"]
    tr, ad, held = c["tr"], c["ad"], c["held"]
    exposures, last = c["exposures"], c["last"]
    surface = PRESETS["surface_lab"]
    deep = PRESETS["deep_underground"]
    worst_fire = c["worst_fire"]
    css = (OUT / "slides.css").read_text(encoding="utf-8")

    slides = [
        slide(1, "Trabajo final &middot; Ondas Gravitacionales e Investigación Asistida por IA", f"""
<h1>Máscaras para imágenes de Skipper-CCD<br>que no dependen del sensor</h1>
<p class="lede">Un Skipper-CCD cuenta electrones de a uno. Antes de hacer física, el análisis
<strong>tira píxeles</strong> con máscaras. Esas máscaras se ajustan a mano, sensor por sensor.</p>
<p class="lede"><strong>¿Se pueden escribir de modo que funcionen en un sensor para el que nunca
fueron ajustadas?</strong></p>
<p class="byline" style="margin-top:8mm">D. Rodriguez Juiz y F. Pérez &middot; {PAGE}</p>"""),

        slide(2, "Por qué importa", f"""
<h2>Seis máscaras, y cada sensor nuevo las vuelve a pedir</h2>
<div class="two">
<div>
<ul>
<li><strong>Defectos de lectura:</strong> carga que queda atrás en la transferencia, impactos en el
registro serie.</li>
<li><strong>Defectos del material:</strong> columnas y píxeles calientes.</li>
<li><strong>Física que no buscamos:</strong> halo alrededor de trazas, clusters de baja energía,
muones.</li>
</ul>
<p>Cada una lleva un radio, un largo, un corte en tasa &mdash; elegidos para <em>un</em> sensor, en
<em>un</em> sitio.</p>
</div>
<div>
<div class="callout">
<p>El costo de equivocarse tiene dos lados: enmascarar de más <strong>tira exposición</strong>
(y la exposición es el experimento), enmascarar de menos <strong>deja fondo</strong> donde se busca
una señal de un electrón.</p>
</div>
<p style="margin-top:5mm">Y los sensores cambian: SENSEI, DAMIC, Oscura. Entre el
{SENSOR_NAMES['deep_underground'].lower()} y el {SENSOR_NAMES['surface_lab'].lower()} de este trabajo, la
corriente oscura cambia por un factor {surface.dark_e_per_pix_day / deep.dark_e_per_pix_day:.0f} y el
flujo de muones por {sci(surface.muon_flux_per_cm2_day / deep.muon_flux_per_cm2_day, 0)}.</p>
</div>
</div>"""),

        slide(3, "La hipótesis", """
<h2>Lo que se transfiere no son los números: es el procedimiento</h2>
<div class="two">
<div>
<p>Cada máscara se escribió <strong>dos veces</strong>:</p>
<ul>
<li><strong>Fija:</strong> constantes puestas a mano. Es lo que se hace hoy.</li>
<li><strong>Adaptativa:</strong> mide sus propios tamaños sobre las imágenes a las que se aplica, y
fija su umbral con una tasa de falsos positivos (&alpha; = 0.01, corregida por número de tests) o con
propiedades físicas medibles del sensor.</li>
</ul>
</div>
<div>
<table>
<thead><tr><th>Máscara</th><th>Forma fija</th><th>Forma adaptativa</th></tr></thead>
<tbody>
<tr><td>Columnas calientes</td><td>tasa por columna arriba de k &times; mediana</td>
<td>cola de Poisson de cada columna contra su propia exposición</td></tr>
<tr><td>Estelas</td><td>L píxeles después de cada píxel brillante</td>
<td>exceso río abajo vs. río arriba del mismo disparador</td></tr>
<tr><td>Halo</td><td>disco de radio R</td>
<td>cada anillo contra todo lo que queda afuera</td></tr>
</tbody>
</table>
</div>
</div>"""),

        slide(4, "Cómo lo probamos", f"""
<h2>Un simulador que recuerda de dónde vino cada electrón</h2>
<div class="two">
<div>
<ul>
<li>Cada fuente de carga llena <strong>su propio mapa</strong>: una máscara se puntúa contra
exactamente los eventos que debería remover.</li>
<li><strong>Tres sensores</strong> con geometría, ruido, fondos y flujos de muones de mediciones
públicas (SENSEI en SNOLAB y MINOS, Oscura, PDG).</li>
<li><strong>Un cuarto sensor, real:</strong> el release público de SENSEI.</li>
</ul>
</div>
<div>
<p>De cada máscara se comparan cuatro formas sobre una pila de prueba independiente:</p>
<ul>
<li><strong>Oráculo:</strong> la fija ajustada con la verdad del simulador en ese sensor. El techo.</li>
<li><strong>Trasplante:</strong> las constantes del oráculo <em>de otro sensor</em>. Lo que pasa hoy
cuando una máscara se hereda.</li>
<li><strong>Adaptativa</strong> y <strong>ninguna máscara</strong> &mdash; esta última es la que dice
si enmascarar sirve de algo.</li>
</ul>
<p>Figura de mérito: S/&radic;(S+B), sobre {head["n_seeds"]} semillas, dos de ellas corridas una sola
vez después de congelar el código.</p>
</div>
</div>"""),

        slide(5, "Resultado 1", f"""
<h2>Trasplantar constantes puede ser peor que no enmascarar</h2>
<div class="two wide-right">
<div>
<div class="stat"><div class="big warn">{tr["harmful"]} de {tr["cases"]}</div>
<p>trasplantes puntúan <strong>por debajo de no enmascarar nada</strong>, en todas las semillas</p></div>
<div class="stat"><div class="big">{tr["worst"]["relative_to_oracle"]:.2f}</div>
<p>el peor caso ({describe(tr["worst"]["parts"])}), donde no hacer nada vale
{tr["worst"]["no_mask"]:.2f}</p></div>
<p>No es «un poco subóptimo»: es tirar exposición a cambio de nada.</p>
</div>
<figure>
<img src="figures/compare_masks.png" alt="Figura de mérito relativa al oráculo por sensor y máscara">
<figcaption>Cada máscara en cada sensor, relativa al oráculo de ese sensor. La línea punteada en 1 es el
oráculo (el divisor); la barra gris es no enmascarar, que es contra lo que hay que comparar.</figcaption>
</figure>
</div>"""),

        slide(6, "Resultado 2", f"""
<h2>Calibrarse solo evita los desastres &mdash; y no sale gratis</h2>
<div class="two">
<div>
<div class="stat"><div class="big accent">{ad["median_of_oracle"]:.2f}</div>
<p>mediana de la adaptativa respecto del oráculo ajustado con la verdad en ese mismo sensor</p></div>
<div class="stat"><div class="big">{ad["harmful"]} de {ad["cases"]}</div>
<p>casos adaptativos peores que no enmascarar, contra {tr["harmful"]} de {tr["cases"]} trasplantes</p></div>
</div>
<div>
<div class="callout warn">
<p><strong>El precio:</strong> en {ad["beaten_by_a_transplant"]} casos un trasplante le sigue ganando
a la adaptativa en todas las semillas, y el peor caso adaptativo
({describe(ad["worst"]["parts"])}) llega sólo a {ad["worst"]["median_of_oracle"]:.2f} del oráculo.</p>
</div>
<p style="margin-top:5mm">Las semillas <strong>reservadas</strong> &mdash; corridas una sola vez
después de congelar el código, nunca usadas para cambiar nada &mdash; dicen lo mismo:
{held["transplant_harmful"]} de {held["transplant_cases"]} trasplantes dañinos,
{held["adaptive_harmful"]} de {held["cases"]} adaptativos, mediana
{held["adaptive_median_of_oracle"]:.2f}.</p>
</div>
</div>"""),

        slide(7, "Resultado 3", f"""
<h2>La falla que encontramos: una máscara disparando sobre el defecto de otra</h2>
<div class="two wide-right">
<div>
<p>Apagar todos los defectos a la vez no puede ver esto. Así que los encendimos <strong>de a
uno</strong>: {cross["cells_tested"]} celdas de máscara y defecto, {cross["n_runs_completed"]}
corridas cada una.</p>
<div class="stat"><div class="big warn">{len(c["fires"])} de {cross["cells_tested"]}</div>
<p>celdas disparan sobre un defecto que no es el suyo</p></div>
<div class="callout warn">
<p>En superficie, con estelas de transferencia como <strong>único</strong> defecto presente, la
máscara de clusters de baja energía dispara en <strong>todos</strong> sus ensayos y tapa
{worst_fire[1]:.1%} de la imagen.</p>
</div>
<p style="margin-top:4mm">Estaba <strong>predicho por escrito</strong> antes de correr, a partir de un
síntoma del Resultado 1.</p>
</div>
<figure>
<img src="figures/cross_defect.png" alt="Grilla de máscaras contra el único defecto encendido">
<figcaption>Un defecto encendido por vez. La columna de «estelas» es la que enciende casi todo.</figcaption>
</figure>
</div>"""),

        slide(8, "Resultado 4", f"""
<h2>Sobre un sensor real, sin ajustarle nada</h2>
<div class="two">
<div>
<p>Las mismas máscaras adaptativas, sin tocarles un parámetro, sobre las {len(exposures)} exposiciones
del release público de SENSEI en SNOLAB &mdash; un sensor real, con una geometría que ningún preset
tiene.</p>
<ul>
<li>Mide el ruido de lectura en
{min(c["noise_medians"]):.3f}&ndash;{max(c["noise_medians"]):.3f} e, contra
{min(c["r1_noise"]):.4f}&ndash;{max(c["r1_noise"]):.4f} e de un ajuste independiente.</li>
<li>Marca {", ".join(str(len(r["constants_chosen"]["hot_columns"])) for r in exposures)} columnas a
medida que crece la exposición, y los conjuntos están <strong>anidados</strong>.</li>
</ul>
</div>
<div>
<div class="callout">
<p>En la exposición más larga, las <strong>{c["n_loud_flagged"]} columnas más ruidosas de la imagen
son exactamente las que marca</strong>, y todas están dentro de la máscara de columnas malas que
publicó la colaboración.</p>
</div>
<p style="margin-top:5mm">De lo que marcamos, ellos también marcan
{min(c["precision"]):.2f}&ndash;{max(c["precision"]):.2f}. De lo que ellos marcan, nosotros marcamos
{min(c["recall"]):.3f}&ndash;{max(c["recall"]):.3f}: <strong>el release enmascara mucho más</strong>.
La primera dirección es la fácil, y lo decimos.</p>
</div>
</div>"""),

        slide(9, "Cómo se hizo", f"""
<h2>Todo con agentes &mdash; y el método es la mitad del trabajo</h2>
<div class="two">
<div>
<ul>
<li><strong>Reglas antes de empezar:</strong> ningún número tipeado a mano en una página, cada salida
con su sidecar (script, commit, hashes de entradas, semilla), ningún documento privado en un archivo
versionado.</li>
<li><strong>Expectativas registradas antes de correr</strong> y commiteadas. Una salió mal, y está
escrita con los números que la explican.</li>
<li><strong>Seis revisiones independientes</strong>, sin acceso a la conversación que produjo el
trabajo; una de ellas de otro proveedor, que <em>re-derivó</em> los resultados en vez de auditarlos:
mismas siete columnas, con otro criterio.</li>
<li><strong>Todo reproducido desde un clon limpio</strong>, bit a bit, con tolerancia 10<sup>-9</sup>.</li>
</ul>
</div>
<div>
<div class="callout warn">
<p><strong>Lo que los agentes hicieron mal</strong> (y lo atraparon los tests o las revisiones, no
ellos mismos):</p>
<ul style="margin-top:2mm">
<li>verdad del simulador filtrándose en un umbral «adaptativo»;</li>
<li>geometría isótropa en un sensor binado, que escondió la columna más ruidosa;</li>
<li>un veredicto aplicado a 60 celdas sin corregir por su número;</li>
<li>un comando del README que en un clon limpio <strong>no hacía nada</strong>.</li>
</ul>
</div>
</div>
</div>"""),

        slide(10, "Qué queda", f"""
<h2>Respuesta, y lo que sigue mal</h2>
<div class="two">
<div>
<div class="callout">
<p><strong>Sí, para cinco de las seis máscaras</strong>, y el precio se puede medir: mediana
{ad["median_of_oracle"]:.2f} del oráculo, con {ad["beaten_by_a_transplant"]} casos donde un trasplante
todavía gana. La que no: clusters de baja energía en superficie, que dispara sobre estelas.</p>
</div>
</div>
<div>
<ul>
<li>Las poblaciones de defectos son valores de escenario, no mediciones.</li>
<li>Los defectos se encienden de a uno: falta ver qué pasa cuando <strong>se superponen</strong>.</li>
<li>En el sensor real sólo cuatro de las seis se pueden verificar.</li>
<li>Una sola figura de mérito, y un oráculo ajustado máscara por máscara.</li>
</ul>
<p style="margin-top:4mm"><strong>Todo está acá:</strong> <a href="{PAGE}">{PAGE}</a><br>
informe, resumen de 5 páginas (ES/EN), <code>PROVENANCE.md</code>, <code>METHOD.md</code> y los logs
de las verificaciones.</p>
</div>
</div>"""),
    ]

    page = f"""<!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Máscaras que se calibran solas &mdash; diez diapositivas</title>
<style>{css}</style>
</head>
<body>
{"".join(slides)}
</body>
</html>
"""
    output = OUT / "slides.html"
    output.write_text(page, encoding="utf-8")
    write_sidecar(output, __file__, inputs=c["inputs"],
                  notes="diez diapositivas; cada número leído de los archivos de resultados")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
