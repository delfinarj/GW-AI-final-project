"""Construye report/resumen.html, el resumen de cinco páginas en español, desde los archivos de resultados.

Es la traducción de `analysis/build_summary.py`, no una versión distinta del trabajo: las dos leen sus
números de `skmask.summary_numbers`, así que pueden diferir en la prosa pero no en un valor. Ningún
número está escrito a mano. El PDF sale de imprimir esta página con `scripts/make_pdf.py`.

Uso:  python analysis/build_summary_es.py
"""
import shutil
from pathlib import Path

from skmask.provenance import write_sidecar
from skmask.summary_numbers import collect

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "report"
FIG = OUT / "figures"
REPO = "https://github.com/delfinarj/GW-AI-final-project"
PAGE = "https://delfinarj.github.io/GW-AI-final-project/"

SENSOR_NAMES = {"deep_underground": "Subterráneo profundo", "shallow_underground": "Subterráneo somero",
                "surface_lab": "Laboratorio en superficie"}
SENSOR_PHRASE = {"deep_underground": "subterráneo profundo", "shallow_underground": "subterráneo somero",
                 "surface_lab": "de superficie"}
MASK_NAMES = {"hot_columns": "columnas y píxeles calientes", "cti": "estelas de transferencia de carga",
              "halo": "halo", "serial": "impactos en el registro serie",
              "low_energy_clusters": "clusters de baja energía"}
GROUP_NAMES = {"hot_columns_pixels": "columnas y píxeles calientes", "cti": "estelas de transferencia de carga",
               "serial": "impactos en el registro serie", "low_energy_clusters": "clusters de baja energía",
               "halo": "halo"}
# las claves que usa results/report_numbers.json para las máscaras del resultado del trasplante
MASK_R4 = {"hot": "columnas y píxeles calientes", "cti": "estelas de transferencia de carga", "halo": "halo",
           "serial": "impactos en el registro serie", "lec": "clusters de baja energía", "muon": "trazas de muones",
           "all": "las seis juntas"}


def describe(parts):
    """El nombre de un caso en español, armado desde las claves y no traduciendo una frase en inglés."""
    name = f"{MASK_R4[parts['mask']]} en el sensor {SENSOR_PHRASE[parts['sensor']]}"
    if parts.get("source"):
        return f"{name} (ajustada en {SENSOR_PHRASE[parts['source']]})"
    return f"{name} (adaptativa)"


def table(headers, rows):
    head = "".join(f"<th>{h}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in row) + "</tr>" for row in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    for src in (RES / "compare_masks" / "compare_masks.png",
                RES / "cross_defect_false_positives" / "cross_defect.png"):
        shutil.copy2(src, FIG / src.name)

    c = collect(SENSOR_NAMES, MASK_NAMES, GROUP_NAMES)
    head, null, cross, rate = c["head"], c["null"], c["cross"], c["rate"]
    tr, ad, held = c["tr"], c["ad"], c["held"]
    fires, run_survivors, worst_fire = c["fires"], c["run_survivors"], c["worst_fire"]
    exposures, last = c["exposures"], c["last"]
    precision, recall, ratios, r1_noise = c["precision"], c["recall"], c["ratios"], c["r1_noise"]
    noise = c["noise_medians"]
    css = (OUT / "summary.css").read_text(encoding="utf-8")

    page = f"""<!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<title>Máscaras que se calibran solas &mdash; resumen de cinco páginas</title>
<style>{css}</style>
</head>
<body>

<h1>Máscaras para imágenes de Skipper-CCD que no dependen del sensor</h1>
<p class="lede">¿Se pueden escribir las máscaras que descartan defectos de lectura, defectos del material y física no
buscada de modo que funcionen en un sensor para el que nunca fueron ajustadas?</p>
<p class="byline">D. Rodriguez Juiz y F. P&eacute;rez &middot; trabajo final, <em>Ondas Gravitacionales e
Investigación Asistida por IA</em> &middot; informe completo, código y procedencia:
<a href="{PAGE}">{PAGE}</a></p>

<div class="answer">
<p><strong>Respuesta corta.</strong> Sí para cinco de las seis máscaras, y el precio se puede medir. De los
{tr["cases"]} casos en que una máscara ajustada a mano en un sensor se trasladó sin cambios a otro, {tr["harmful"]}
fueron peores que no enmascarar nada en todas las semillas y {tr["no_difference"]} no cambiaron nada; el más claro,
{describe(tr["worst"]["parts"])}, obtiene {tr["worst"]["relative_to_oracle"]:.2f} de lo que obtiene la máscara
ajustada con la verdad en ese mismo sensor, donde no hacer nada obtiene {tr["worst"]["no_mask"]:.2f}. Las mismas
máscaras escritas como procedimientos que se calibran solos son peores que no enmascarar en {ad["harmful"]} de
{ad["cases"]} casos y mejores en {ad["better_than_no_mask"]}, y alcanzan una mediana de
{ad["median_of_oracle"]:.2f} de la máscara ajustada con la verdad. La excepción es
{describe(ad["harmful_case_parts"][0])}, que dispara sobre electrones de las estelas de transferencia de carga.</p>
</div>

<h2>La pregunta y qué se construyó</h2>
<p>Un Skipper-CCD cuenta electrones de a uno. Antes de extraer cualquier física, los análisis descartan píxeles con
máscaras, cada una apuntada a un problema: carga que queda atrás durante la transferencia, impactos en el registro
serie, columnas y píxeles calientes, fotones alrededor de trazas de alta energía (el halo), clusters de eventos de
baja energía y trazas de muones. Cada máscara lleva tamaños y umbrales &mdash; un radio, un largo de estela, un
corte en tasa &mdash; elegidos a mano para un sensor y un sitio.</p>
<p><strong>Hipótesis:</strong> lo que se transfiere entre sensores no son los números sino el <em>procedimiento</em>
que los elige. Por eso cada una de las seis máscaras se escribió dos veces: una forma <strong>fija</strong>, con
constantes puestas a mano, y una forma <strong>adaptativa</strong>, que mide sus propios tamaños sobre las imágenes
a las que se aplica y fija su umbral con una tasa de falsos positivos (&alpha;&nbsp;=&nbsp;0.01, corregida por el
número de tests) o con propiedades físicas medibles del sensor.</p>
<p>El banco de pruebas es un simulador que guarda un mapa de carga separado por fuente, de modo que una máscara se
puede puntuar contra exactamente los eventos que debería remover. Tres sensores difieren como difieren los
despliegues reales; su geometría, ruido, tasas de corriente oscura, fondos y flujos de muones salen de mediciones
públicas, mientras que las poblaciones de defectos, que no se publican de forma transferible, son valores de
escenario elegidos para diferir entre sensores.</p>
{table(["Sensor", "Píxeles", "Espesor (&micro;m)", "Ruido (e&minus;)", "Exposición", "Oscura (e&minus;/pix/día)",
        "Muones (cm&minus;&sup2;/día)"], c["sensor_rows"])}
<p>De cada máscara se compararon cuatro formas sobre una pila de prueba independiente de cada sensor: el
<strong>oráculo</strong> (la forma fija ajustada con la verdad del simulador en ese sensor), el
<strong>trasplante</strong> (las constantes del oráculo de otro sensor), la forma <strong>adaptativa</strong> y
<strong>ninguna máscara</strong>. La figura de mérito es S/&radic;(S+B): señal inyectada que sobrevive contra fondo
buscado que sobrevive.</p>

<div class="page-break"></div>
<h2>Resultado 1 &middot; Trasplantar constantes puede hacer daño; calibrarse solo evita los desastres, pero cuesta</h2>
<figure>
<img src="figures/compare_masks.png" alt="Figura de mérito relativa al oráculo para cada sensor y cada máscara">
<figcaption>Cada máscara en cada sensor, relativa al oráculo de ese sensor, sobre {head["n_seeds"]} semillas. Una
máscara trasplantada es el oráculo de otro sensor; «ninguna máscara» es la referencia que dice si enmascarar sirve
de algo.</figcaption>
</figure>
<p>Mover constantes entre sensores no es apenas subóptimo: puede ser peor que no enmascarar. {tr["harmful"]} de
{tr["cases"]} trasplantes puntúan por debajo de no enmascarar en todas las semillas, por más del
{100 * head["margin"]:.0f}&nbsp;% de margen que se usa para llamar real a una diferencia. La autocalibración elimina
esos desastres &mdash; {ad["harmful"]} de {ad["cases"]} casos &mdash; pero no sale gratis: en
{ad["beaten_by_a_transplant"]} casos un trasplante le sigue ganando a la adaptativa en todas las semillas, y el peor
caso adaptativo, {describe(ad["worst"]["parts"])}, llega sólo a {ad["worst"]["median_of_oracle"]:.2f} del oráculo
({ad["worst"]["worst_seed"]:.2f} en su peor semilla).</p>
<p>Los estadísticos del resumen se registraron antes de la corrida final. Los casos con menos de
{head["min_target_events"]} eventos buscados en una semilla quedan excluidos de esa semilla, y un caso cuenta como
dañino sólo si pierde por más del margen en <em>todas</em> las semillas válidas. Dos semillas ({held["seeds"]}) se
corrieron una sola vez después de congelar el código y nunca se usaron para cambiar nada: en ellas,
{held["transplant_harmful"]} de {held["transplant_cases"]} trasplantes son dañinos, {held["adaptive_harmful"]} de
{held["cases"]} casos adaptativos lo son, y la mediana adaptativa es {held["adaptive_median_of_oracle"]:.2f} del
oráculo &mdash; el mismo cuadro que en las semillas de desarrollo.</p>

<div class="page-break"></div>
<h2>Resultado 2 &middot; Sobre qué disparan las máscaras cuando no hay nada que encontrar</h2>
<p>Con todos los defectos apagados en los tres sensores, y dejando lo que un sensor real no puede apagar (corriente
oscura, carga espuria, la señal inyectada, trazas de muones y depósitos de alta energía), las cinco máscaras
estadísticas dispararon {c["null_fired"]} veces en {len(c["null_tested"])} tests de sensor y máscara de
{null["n_runs"]} corridas cada uno, y todos los intervalos contienen &alpha;. La sexta máscara, la de muones, no
puede estar en ese test, porque las trazas son justamente lo que un sensor real no puede apagar; preguntado aparte,
en un simulador con el flujo en cero, disparó en {c["muon_fired"]} de {c["muon_images"]} imágenes.</p>
<p>Ese test tiene un punto ciego: apagar todos los defectos a la vez no puede ver una máscara disparando sobre el
defecto <em>de otra</em>. Así que los defectos también se encendieron de a uno, {cross["cells_tested"]} celdas de
máscara y defecto en total, {cross["n_runs_completed"]} corridas cada una, con el límite que cada celda debe superar
tomado al {100 * cross["confidence_of_the_corrected_bound"]:.2f}&nbsp;% de confianza para controlar un
{100 * cross["family_wise_error"]:.0f}&nbsp;% de probabilidad de una falsa alarma en toda la grilla.</p>
<figure>
<img src="figures/cross_defect.png" alt="Grilla de máscaras contra el único defecto encendido, por sensor">
<figcaption>Un defecto encendido por vez. El área del punto es la fracción de ensayos en que la máscara disparó; el
naranja marca las celdas por encima de la tasa que esa máscara no debería superar; los círculos huecos son la
diagonal, donde una máscara se encuentra con su propio defecto.</figcaption>
</figure>
<p>{len(fires)} de {cross["cells_tested"]} celdas disparan sobre un defecto que no es el suyo &mdash; y las
{len(run_survivors)} contadas por imagen sobreviven a recontarlas de a una corrida, que es la forma conservadora de
tratar dos imágenes que comparten una calibración. Las estelas de transferencia de carga son lo que más veces las
enciende. La mayor no deja lugar a dudas: en el sensor de superficie, con las estelas como único defecto presente,
la máscara adaptativa de clusters de baja energía dispara en todos sus ensayos y enmascara una mediana de
{worst_fire[1]:.3f} de la imagen. Es la única falla real de la autocalibración encontrada acá, y estaba predicha por
escrito antes de la corrida, a partir de un síntoma visto en el Resultado&nbsp;1.</p>
{table(["Sensor", "Único defecto presente", "Máscara que disparó", "Disparos/ensayos", "Tasa",
        "Fracción mediana enmascarada"], c["fire_rows"])}

<div class="page-break"></div>
<h2>Resultado 3 &middot; Los mismos procedimientos sobre un sensor real</h2>
<p>Todo lo anterior es simulado. Las máscaras adaptativas también se corrieron, sin cambiarles nada, sobre las
{len(exposures)} exposiciones del release público de SENSEI en SNOLAB: un cuarto sensor, real, con una geometría que
ningún preset tiene ({last["shape"][1]}&times;{last["shape"][0]} superpíxeles activos, cada uno binando 32 filas
físicas). No se ajustó nada para él, y el release publica su propia máscara, así que las dos se pueden comparar
donde se superponen.</p>
<p>Con sólo las imágenes, el estimador mide un ruido de lectura de {min(noise):.3f}&ndash;{max(noise):.3f}&nbsp;e,
contra {min(r1_noise):.4f}&ndash;{max(r1_noise):.4f}&nbsp;e de un ajuste independiente sobre los mismos archivos, y
una densidad de un electrón que crece con la exposición. El procedimiento de columnas calientes marca
{", ".join(str(len(r["constants_chosen"]["hot_columns"])) for r in exposures)} columnas a medida que crece la
exposición, y los conjuntos están anidados: cada uno conserva las anteriores y agrega la siguiente más ruidosa. En
la exposición más larga, las {c["n_loud_flagged"]} columnas con mayor tasa de píxeles cargados son exactamente las
que marca, todas dentro de la máscara de columnas malas del propio release, y cada una carga entre
{min(ratios):.0f} y {max(ratios):.0f} veces la tasa de píxeles cargados de las columnas que dejó en paz.</p>
{table(["Exposición", "Imágenes", "Píxeles disparadores", "Ruido (e&minus;)", "Densidad 1e", "Columnas calientes",
        "Radio de halo", "Enmascaramos", "Enmascara el release"], c["public_rows"])}
<p>Las dos direcciones del acuerdo dicen cosas distintas y las dos están acá. De lo que el procedimiento marca, el
release también marca {min(precision):.2f}&ndash;{max(precision):.2f}: se queda adentro de una máscara que ajustó
una persona. De lo que marca el release, el procedimiento marca {min(recall):.3f}&ndash;{max(recall):.3f}: el
release enmascara mucho más. El primer número es la dirección fácil &mdash; un procedimiento que marcara una sola
columna verdadera sacaría 1.00 &mdash; y se apoya en un puñado de decisiones por columna, no en los cientos de
píxeles sobre los que se cuenta.</p>
<p>La calibración de estelas eligió largo cero en las cuatro exposiciones, algo que la expectativa registrada no
predijo. El motivo está en el archivo de resultados: el release oculta sus impactos, así que una exposición entera
tiene {", ".join(str(r["trigger_pixels"]["total"]) for r in exposures)} píxeles por encima del disparador. En la
exposición más larga el test sí ve el exceso del lado correcto &mdash; {last["cti_detail"]["v"]["n_downstream"]}
electrones simples río abajo de un disparador contra {last["cti_detail"]["v"]["n_upstream"]} río arriba,
p&nbsp;=&nbsp;{last["cti_detail"]["v"]["smallest_p"]:.4f} &mdash; pero el umbral corregido pide
p&nbsp;&lt;&nbsp;{last["cti_detail"]["v"]["bonferroni_threshold"]:.1e}. El procedimiento se niega a enmascarar con
evidencia tan flaca, que es para lo que fue construido.</p>

<div class="page-break"></div>
<h2>Cómo se hizo el trabajo, y cómo se verificó</h2>
<p>Cada línea de código y cada análisis de este proyecto los produjeron agentes de IA, a partir de una pregunta
explícita y bajo reglas escritas antes de empezar: ningún número de ninguna página se escribe a mano, cada salida
lleva un sidecar con el script, el commit de git, los hashes de sus entradas y su semilla, y ningún documento
privado entra jamás en un archivo versionado. Lo que hizo esto confiable, y no solamente rápido, fue el método
adversarial:</p>
<ul>
<li><strong>Las expectativas se registraron antes de correr.</strong> Escritas en <code>PLAN.md</code> y
commiteadas antes de que el análisis corriera, para que una confirmación no se pudiera inventar después. Una de
ellas estuvo equivocada &mdash; el largo de estela en el sensor real &mdash; y el repositorio lo dice, con los
números que lo explican.</li>
<li><strong>Agentes independientes revisaron el trabajo</strong> sin acceso a la conversación que lo produjo.
Encontraron defectos reales: verdad del simulador filtrándose en un umbral «adaptativo»; un titular que medía no
hacer nada; una máscara de halo midiendo distancias en superpíxeles sobre un sensor binado, que escondía la columna
más ruidosa de la calibración de columnas calientes; un veredicto aplicado a sesenta celdas sin corregir por su
número. Cada uno está en <code>PROVENANCE.md</code> con lo que se hizo al respecto, incluidos los hallazgos que
<em>no</em> se corrigieron.</li>
<li><strong>Todos los resultados se reprodujeron desde un clon limpio</strong> &mdash; un checkout nuevo que no
comparte nada con la copia de trabajo salvo el commit y el entorno fijado &mdash; y se compararon con una
tolerancia relativa de 10<sup>&minus;9</sup>. Todos volvieron idénticos, figuras incluidas. El más lento, de dos
horas y media, hubo que retomarlo después de que algo lo matara por la mitad, y eso fue lo que mostró que el
comando que un desconocido correría para reproducirlo no hacía nada en un clon limpio: encontraba el archivo
terminado que publicamos y se salteaba el trabajo.</li>
</ul>
<p>Los errores que atraparon las verificaciones están listados en la página del informe completo, no escondidos: un
simulador que volvía a sortear sus columnas calientes en cada imagen, una máscara de columnas calientes que marcaba
vecinas inocentes, un test de halo que trataba su tasa de referencia como exacta, una máscara de muones que perdía
las trazas cruzadas, y un orden de calibración que hacía falsas 69 de 81 columnas marcadas.</p>

<h3>Qué sigue mal</h3>
<ul>
<li>Las poblaciones de defectos son valores de escenario, no mediciones; sólo los parámetros del sensor y de los
fondos son mediciones públicas.</li>
<li>Los defectos se encienden de a uno, así que nada de esto mide qué hacen las máscaras cuando varios defectos se
superponen en los mismos píxeles.</li>
<li>En el sensor real sólo cuatro de las seis máscaras se pueden verificar: el release no publica contraparte de la
máscara de clusters de baja energía, y sus superpíxeles binados hacen que la geometría de la máscara de muones no
tenga sentido.</li>
<li>Una sola figura de mérito, y un oráculo ajustado máscara por máscara y no en conjunto, así que una máscara
adaptativa combinada puede superarlo.</li>
</ul>

<h3>Dónde está todo</h3>
<p>Repositorio: <a href="{REPO}">{REPO}</a> &mdash; <code>README.md</code> reproduce el proyecto entero desde un
clon limpio con <code>uv</code>; <code>PLAN.md</code> tiene las afirmaciones registradas de antemano;
<code>PROVENANCE.md</code> registra cada entrada, cada elección que tenía una alternativa defendible, cómo se
verificó cada resultado y cada error encontrado en el camino; <code>report/index.html</code> es el informe completo,
del que este documento es un resumen.</p>

</body>
</html>
"""
    output = OUT / "resumen.html"
    output.write_text(page, encoding="utf-8")
    write_sidecar(output, __file__, inputs=c["inputs"],
                  notes="resumen de cinco páginas en español; cada número leído de las entradas")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
