# Reproduce the initial OTA's finite sensing divider

This example keeps the initial six-unit OTA and adds two 1-MΩ sensing
resistors. Its input is 0.3 V through 10 kΩ; its output drives 1 MΩ and 5 pF.
The supply is 1 V at 26.85 °C, with the original finite bias resistor and
saved zero-raw MOS units. It is a gain-two connection with finite error and
loading, distinct from the unity buffer's required function.

Read [the operating connection](https://ds54e.github.io/analog-design-notes/learn/04-differential-pair-and-ota.html#sense-low-while-driving-a-higher-output)
and [the precision argument](https://ds54e.github.io/analog-design-notes/learn/05-feedback-and-response.html#a-finite-sensing-divider-changes-both-the-target-and-its-errors).
Download [learning-feedback-divider-v1.0.zip](https://ds54e.github.io/analog-design-notes/downloads/learning-feedback-divider-v1.0.zip).
The [evidence record](evidence.json) binds the selected inputs and actual checks.

## Inputs, definitions and exposure

The [exact SPICE body](circuits/D1-half-divider.cir) lists all MOS D/G/S/B
terminals and zero-volt current sensors. Five units are 4/0.24 µm and the
reference unit is 1/0.40 µm: 5.2 µm² total drawn MOS area. Rbias is
273692.63541213644 Ω. Rupper connects output to the sensing node; Rlower
connects that node to ground. The negative input reaches the sensing node
through the retained zero-volt injection source. No selected-D1 mirror
correction, geometry change, redraw or rebias is applied.

The [frozen corrected plan](plan.json) declares exactly one OP and one
361-point AC sweep from 1 kHz to 1 GHz. The [observation index](data/index.json)
retains its actual start, circuit/deck, request/realization and native-output
identities. Physical/model readback and output completeness were checked
before applying the unchanged [measurement script](analyze.py).
AC magnitude 1 normalizes a small-signal derivative; it does not apply a
physical 1-V large-signal excursion to the 1-V circuit.

The prior uses three exact members of the already published
[learning-day-v1.0.zip](https://ds54e.github.io/analog-design-notes/downloads/learning-day-v1.0.zip), SHA256
`ddfbad7a06c9c45a94fa4f0379a421816f04d3fa1c56bd21ec8a09ac99f13b15`.
[source.json](data/source.json) identifies its configuration, original D1-5p
circuit and complete observations. Only the known central OP/AC are used.
The original source schema name is retained as provenance; the new divider
plan and analysis have their own identities.

### A rejected input representation is retained

The [original input plan](original-input-plan.json) redundantly listed the
two fixed divider resistors as per-analysis overrides. The finite acquisition
adapter accepts that field only for Rload/Cload. It rejected the request
**before creating a native deck, attempt directory or simulator process**.
No scientific target was exposed by that invocation.

The [correction record](input-correction.json) preserves the original plan
identity and rejection. Removing only those redundant overrides kept both
resistors in the exact, bound circuit. The circuit bytes, six physical units,
source conditions, prior forecast, measurement script and all tolerances
remained unchanged. A new input/freeze destination preceded the **one completed
native condition**. There was no native retry or outcome-driven extension.
No new transient, noise, Local cohort or regulator condition was acquired.

## Recompute the selected observations

Extract the new ZIP in a clean directory. It contains `example/`, with the
exact scripts, compact arrays, source hashes and retained notices. Use Python
3.9 or later and the established pinned environment:

~~~bash
python3 -m venv .venv
.venv/bin/python -m pip install -r example/requirements.txt
curl --fail --location https://ds54e.github.io/analog-design-notes/downloads/learning-day-v1.0.zip --output learning-day-v1.0.zip
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/prior.py --source-zip learning-day-v1.0.zip --output "$PWD/prior-001"
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/analyze.py --source-zip learning-day-v1.0.zip --output "$PWD/analysis-001"
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/network.py --source-zip learning-day-v1.0.zip --output "$PWD/network-001"
~~~

The five expected files are `prior.json`, `summary.json`,
`operating-point.csv`, `network.json` and `error-accounting.csv`. Both CSV
writers use LF separators. Existing output directories are refused.
The [data manifest](data/manifest.json) binds the selected files. This one
condition uses a worked table; it adds no graph or schematic generator.
Previously published graph packages remain available.

Actual resistor parameters and independent DC V/I checks cover all five
resistors. The saved MOS and node KCL checks include the sensing node.
The capacitor is checked over 100 kHz–10 MHz from output current balance:

```text
Iupper = (Vout - Vfeedback) / Rupper
Icap = -i(VtM2d) - i(VtM4d) - Vout/Rload - Iupper
C = Icap / (j*2*pi*f*Vout)
```

The bandwidth definition is the first downward **−3.000-dB** crossing relative
to the observed 1-kHz gain, interpolated in dB/log frequency. It remains the
bridge definition; it is not the earlier studies' half-power interpolation
or a measured return ratio.

## Keep the three reasoning steps distinct

Before exposure, the old-OP affine current equations predicted 0.580819 V,
local gain 1.921755 and bandwidth 2.169364 MHz. They retained the known
terminal-current intercept and held the old gm/gmb/gds fixed. Independent
six-node and reduced three-node constructions agree on the local transfer.

The new observation is **0.579784 V**, gain **1.906249** and bandwidth
**2.259680 MHz**. Retain the −1.035-mV output, −0.807% gain and +4.163%
bandwidth discrepancies. The ideal divider values 0.6 V and gain 2 are
reference values, not retrospectively adjusted pass thresholds.

The separate new-OP calculation gives gain 1.902202 and a one-external-
capacitor bandwidth of 2.205055 MHz. Its remaining discrepancies are −0.212%
and −2.417%. It omits intrinsic charging and gate/body-current derivatives.
The additional 1-kHz gate-path explanation uses **already measured** input
transfers and the external capacitor; it is labeled retrospective, not an
independent forecast or a whole-frequency-range validation.

Exact stationary KCL accounts for output error as +0.046895 mV from the
positive-input source term, −17.519228 mV from finite input difference and
−2.743496 mV from sensing-gate loading. Their sum is −20.215829 mV.
They are correlated quantities at one operating point, not separate random
errors. The divider consumes 0.168079 µW; total VDD current is 16.016208 µA.
All signed MOS/resistor/port powers close within the declared residual check.

## Replay the known native condition

Install **ngspice 47** separately. The site build does not install or run
SPICE. Obtain the public APM source at the exact v5.0.0 pin:

~~~bash
git clone --branch v5.0.0 --depth 1 https://github.com/ds54e/analog-process-models.git apm-v5
git -C apm-v5 rev-parse HEAD
# Expected: 381517fda5107fabf98af7801d5a5103f38e230c
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B example/replay.py --source-zip learning-day-v1.0.zip --apm "$PWD/apm-v5" --binary /usr/local/bin/ngspice --output "$PWD/native-replay-001"
~~~

The [post-exposure replay plan](replay-plan.json) fixes the known condition and
comparison tolerances. Relocation creates new request/realization identities;
the six saved physical values remain identical. Native APM checks actual
model/W/L/m/nf/DELVTO/MULU0 readback. The replay compares all selected terminal
OP/AC observations and independently checks resistor V/I and capacitor current.
The derivative-based saved calculation uses the checked original OP parameters.

The output directory retains its request, realization, native logs and
`result.json`. This is a reproducibility check using public inputs on the same
host, not fresh confirmation or independent scientific replication. It does
not establish a voltage range, step/overload behavior, noise performance,
silicon yield or an integrated LDO. A regulator's pass device adds another
inversion and a different receiving port; its feedback connection must be
derived separately.

The checked clean run used Python 3.9.25, the listed pins and an anonymously
retrieved APM v5.0.0 checkout. All five calculated JSON/CSV files were
byte-identical. One OP and 361 AC rows reproduced every selected terminal
voltage/current exactly. Nineteen execution/input payloads and sixteen
license notices matched the selected package. These are reproducibility
checks of this known nominal condition, with the limits above.
