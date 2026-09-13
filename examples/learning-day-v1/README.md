# Reproduce the learning-day supplement

This **learning-day-v1.0** package accompanies learning edition **v1.1**.
It adds a worked sizing/mirror reconstruction and a separately planned
initial-D1 receiving-port example. The five study ZIPs and
[learning-path-v1.0](https://ds54e.github.io/analog-design-notes/studies/learning-path-v1/reproduce.html) remain unchanged.

Download [learning-day-v1.0.zip](https://ds54e.github.io/analog-design-notes/downloads/learning-day-v1.0.zip), or
inspect the [public example files](https://github.com/ds54e/analog-design-notes/tree/main/examples/learning-day-v1).
The ZIP has one `example/` directory. Start with the
[English route](https://ds54e.github.io/analog-design-notes/learn/index.html) or [Japanese overview](https://ds54e.github.io/analog-design-notes/ja/learning-path.html)
for the circuit operation and design judgments.

## What can be recalculated from saved observations?

| Selection | Inputs and role | Recalculation |
| --- | --- | --- |
| Twelve feasibility points | Four nominal grounded-source 2-µm devices, three gate voltages; known first-study development observations. | Current efficiency, projected total width and optimistic loaded gain. |
| R37/R75 first and final designs | Four saved OP/AC cases; original finite source/load and unit inventory. | Current/voltage constraints, actual first/final gain, local-derivative explanation and resource cost. |
| Preparation mirror | One saved nonzero-raw SOURCE_TRANSFER_HYPOTHESIS realization; all 51 output-clamp points. | Copying error versus voltage and a retrospective illustrative ±1% grid comparison. |
| A/B/C headroom | Exact arrays copied from the original learning package with existing identities. | Actual N-input and P-input voltage budgets, not a new common-mode sweep. |
| Initial-D1 receiving port | New same-core 5/25-pF functional attempts, one 25-pF timestep refinement, six C-fixture OP clamps. | DC error, signal bandwidth, current/headroom, independent-endpoint settling and capacitor charge. |

The [design configuration](data/design-inputs.json),
[headroom configuration](data/headroom-inputs.json) and
[new interface configuration](data/interface-inputs.json) provide exact circuit
hashes, original identities, column names, recipes, physical-unit records and
verified parameter readbacks. Named NPZ files contain complete columns and
original float64 points of the selected analyses. There are 1,671 rows in
the nine selected preparation/design cases and 18,337 rows in the nine new
interface cases; the reused A/B/C arrays retain their original selection.
These are finite explanatory records, not independent samples to pool.

The feasibility circuits fix drain at 0.5 V and source/body at zero. Their
current/derivative width scaling is an approximation, not a native observation
at every estimated width. Individual supported dimensions remain W=1–4 µm,
L=0.12–0.40 µm. The R37/R75 design bodies retain their original 10-mV
development sine source, but only OP/AC are selected here. The larger-signal
distortion failure and the later 5-mV study identity remain distinct.

The mirror's two saved raw vectors are nonzero, under its original
SOURCE_TRANSFER_HYPOTHESIS profile. They are not relabeled nominal or a
qualified Local population. Its source circuit is an unmodified APM example
with the original Apache-2.0 header. The new code and numerical analysis do
not relicense that circuit or the separately retrieved model inputs.

## Regenerate the summaries and five new graphs

Use Python 3.9 or later. The checked scientific dependencies are NumPy 1.26.4,
SciPy 1.13.1 and Matplotlib 3.9.4. In a new directory after extracting the ZIP:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r example/requirements.txt
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/design.py --output "$PWD/design-001"
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/interface.py --output "$PWD/interface-001"
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/plot_design.py --output "$PWD/design-figures-001"
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/plot_interface.py --output "$PWD/interface-figures-001"
```

No SPICE is called. Each command refuses an existing output directory.
Compare the generated [design summary](design-summary.json),
[interface summary](interface-summary.json), [feasibility CSV](feasibility.csv),
[mirror CSV](mirror.csv) and [clamp CSV](clamps.csv) with the supplied versions.
The two design and three interface SVGs appear in chapters 2–5; the scripts
also produce PNG inspection copies. The top-level data
[manifest](data/manifest.json) and narrower calculation manifests expose file
identities. No private archive, database or runtime record is needed.

Expected selected values include R37's grounded-source width estimate
8.378069 µm versus its solved total width 8.841660 µm, and the mirror's
−0.434%/+2.505% copying error at 0.30/0.80 V. The new 25-pF D1 has the
same DC outputs as 5 pF at the three selected input values, bandwidth
0.847733 MHz, and 2-ns-step rising/falling entry brackets
889.3–891.3 / 878.9–880.9 ns. More digits in the JSON support calculation
comparison; they do not represent silicon precision.

## New receiving-port definitions and exposure

The [finite plan](receiving-port-plan.json) was prepared and committed before
the new target acquisition, with the same six saved zero-raw ARTIFICIAL units.
The [pre-exposure note](receiving-port-estimate-note.json) corrects insignificant
excess digits transcribed in the known training bandwidth; the original plan
bytes and scientific conditions are retained. No sample was drawn or trimmed.
The later archived selected-D1 25-pF transfer is excluded.

Functional D1 retains 1 V, 26.85 °C, Rbias=273692.63541213644 Ω, a 10-kΩ
source, and 1 MΩ at the output. Only the external capacitor changes from
5 to 25 pF. The original PWL is 0.25→0.35→0.25 V, with 12.5-ns edges.
Separate OP solves at 0.25/0.30/0.35 V establish the stationary endpoints.
The ±0.983220652-mV settling band is 1% of the actual output step, with time
measured from each input's 50% crossing and observation through 4.5/8.5 µs.
A single predeclared 25-pF repeat halves maximum timestep from 2 to 1 ns.

AC uses the original learning bridge's −3.000-dB crossing relative to 1 kHz,
interpolated in dB/log frequency. The older ADS half-power/magnitude definition
is not overwritten. The new OP endpoint at 0.25 V differs from the older
DC-sweep endpoint by about 3.55 nV; the old array and definition remain intact.

Signed MOS currents point into each terminal. The functional output-capacitor
current is `-i(VtM2d)-i(VtM4d)-v(out)/Rload+i(Vfeedback)`. The feedback sensor
points from `feedback` toward `out`, so its sign is explicit. Integrating this
current and comparing with `Cload*(Vout(t)-Vout(0))` checks the charge identity;
intrinsic MOS charge is already included in the terminal currents. The C-fixture
clamps keep the input common mode and differential voltage fixed while a
separate output voltage source takes the imbalance. They are not the functional
buffer driven over the entire 0.05–0.70-V diagnostic span.

All nine planned research attempts completed with usable finite outputs.
Their observations are now known. Source/current/headroom loss is retained
as a result; no bias tuning, retries, redrawing, extended cohort or relaxed
criterion was used to repair the example. Numerical consistency and a chosen
physical performance requirement remain separate judgments.

## Replay the new 25-pF case from public model inputs

Install **ngspice 47** separately. Retrieve the public pinned model repository:

```bash
git clone --branch v5.0.0 --depth 1 https://github.com/ds54e/analog-process-models.git apm-v5
git -C apm-v5 rev-parse HEAD
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B example/replay_interface.py --apm "$PWD/apm-v5" --binary /usr/local/bin/ngspice --output "$PWD/replay-25p-001"
```

The required APM commit is `381517fda5107fabf98af7801d5a5103f38e230c`.
Change `--binary` to your installed ngspice 47 location. APM's own Apache-2.0
code and FreePDK45 model notices remain in that separately retrieved checkout.
No device-model file is redistributed inside this package.

This one finite replay uses the exact new 25-pF circuit and saved physical
units through APM's public API. It runs three OPs, 361 central AC points and
the 8.5-µs transient (4,332 rows on the checked host). It compares native
outputs with the supplied observations, requiring differences below 1 µV
for OP/transient and 0.001 V/V for complex unit-stimulus AC. Adaptive transient
results are linearly compared at the saved time coordinates; endpoints,
maximum timestep and complete horizon are checked independently.

APM verifies actual model identity, W/L/m/nf and DELVTO/MULU0 before/applied/after
each analysis. This public replay reads the three resistor values from measured
DC voltage/current ratios. It reads capacitance through complex AC output
current/voltage and transient charge. The private acquisition additionally
recorded direct native passive parameter vectors, which the public API does
not expose. It does not substitute a copied expected parameter for a readback.

See the [selected evidence receipt](evidence.json) for the clean input/graph/native
checks actually included in the release. Same-host replay with an existing
ngspice 47 installation is not a clean-OS or all-platform qualification.

## Scope and review limits

These examples explain pinned-model nominal behavior and a selected historical
source-transfer realization. They do not establish manufacturing yield,
calibrated passive/spatial statistics, full pass-gate charge, arbitrary overload
recovery, a regulator return ratio, layout/PEX or silicon qualification. The
old self-bias equality-check failure remains unresolved and outside this work.
The same Codex agent wrote and checked the material; there is no claim of prior
human scientific approval or independent replication. Native acquisition stays
outside the static Pages build, and every earlier release remains retrievable.
