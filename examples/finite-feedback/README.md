# Reproduce finite feedback and its input cost

Download the [finite feedback 1.0 example](https://ds54e.github.io/analog-design-notes/downloads/feedback-study-v1.0.zip),
or use `examples/finite-feedback` in the public Notes repository. The package
contains eight exact SPICE bodies, 12,504 saved MOS unit records, compact
nominal/development/confirmation data and graph/statistic/replay code. The
public APM dependency is retrieved separately under its original terms.

The commands below were actually executed in `/tmp/adn-feedback-clean-001`,
using only selected public inputs, a fresh public model clone and a fresh Python
environment. The host was Linux x86_64, Python 3.9.25, with ngspice 47 already
installed at `/usr/local/bin/ngspice`. This checks clean inputs/Python on the
same host; it does not qualify a clean OS or every platform. Use new output
destinations for another replay. SPICE acquisition is outside the site build.

## Obtain the pinned public dependencies

```bash
git -c credential.helper= -c http.extraHeader= clone --depth 1 --branch v5.0.0 https://github.com/ds54e/analog-process-models.git /tmp/adn-feedback-clean-001/apm-v5
python3 -m venv /tmp/adn-feedback-clean-001/.venv
/tmp/adn-feedback-clean-001/.venv/bin/python -m pip install --disable-pip-version-check --no-cache-dir --index-url https://pypi.org/simple -r /tmp/adn-feedback-clean-001/example/requirements.txt
```

The APM tag resolves to `381517fda5107fabf98af7801d5a5103f38e230c`.
The replay checks that commit, the actual model hashes and public data hashes.
Obtain ngspice 47 separately if necessary from the
[official project](https://ngspice.sourceforge.io/). No private repository or
GitHub authentication is required for these dependencies or example inputs.

## Regenerate the compact statistics and figures

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/adn-feedback-clean-001/.venv/bin/python -B /tmp/adn-feedback-clean-001/example/statistics.py --data /tmp/adn-feedback-clean-001/example/data --output /tmp/adn-feedback-clean-001/recomputed-statistics.json
MPLCONFIGDIR=/tmp/adn-feedback-clean-001/mplconfig PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/adn-feedback-clean-001/.venv/bin/python -B /tmp/adn-feedback-clean-001/example/plot.py --data /tmp/adn-feedback-clean-001/example/data --output /tmp/adn-feedback-clean-001/regenerated-plots
```

All eleven SVGs matched the source bytes. The complete statistics JSON also
matched: eight candidate/condition summaries, eight frozen paired intervals,
and ten reconstructed source-role variance terms. The five common functional
checks were recomputed from scalars. These operations do not rerun the full
4,168-condition cohort or its reference calibration. Compact dynamic plotting
traces use interpolation grids; the saved scientific metrics came from the
complete native observations.

## Replay the four examples selected before target exposure

All four use coordinate 2000, its saved UID/geometry/DELVTO/ln(MULU0), and
the frozen circuit. No sampling, reference remapping, trimming, rebiasing or
recentering is performed.

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/adn-feedback-clean-001/.venv/bin/python -B /tmp/adn-feedback-clean-001/example/replay.py --apm /tmp/adn-feedback-clean-001/apm-v5 --output /tmp/adn-feedback-clean-001/replay-F75L-nominal-001 --candidate F75L --index 2000 --condition nominal
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/adn-feedback-clean-001/.venv/bin/python -B /tmp/adn-feedback-clean-001/example/replay.py --apm /tmp/adn-feedback-clean-001/apm-v5 --output /tmp/adn-feedback-clean-001/replay-F75L-supply-001 --candidate F75L --index 2000 --condition supply_097
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/adn-feedback-clean-001/.venv/bin/python -B /tmp/adn-feedback-clean-001/example/replay.py --apm /tmp/adn-feedback-clean-001/apm-v5 --output /tmp/adn-feedback-clean-001/replay-R75-signal20m-001 --candidate R75 --index 2000 --condition signal_20m
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/adn-feedback-clean-001/.venv/bin/python -B /tmp/adn-feedback-clean-001/example/replay.py --apm /tmp/adn-feedback-clean-001/apm-v5 --output /tmp/adn-feedback-clean-001/replay-F75L-signal20m-001 --candidate F75L --index 2000 --condition signal_20m
```

| Example | Output (V) | AC gain (V/V) | DC delivered current (µA) | 20-mV THD (%) |
| --- | ---: | ---: | ---: | ---: |
| F75L nominal | 0.5226815247957144 | −6.014289590680717 | 76.33383202878405 | — |
| F75L at 0.97 V | 0.4845597330136635 | −5.9742987788816855 | 72.07466172662554 | — |
| R75, 20-mV sine | 0.4970156264947820 | −6.021233234361601 | 75.44767722212904 | 4.129288749650584 |
| F75L, 20-mV sine | 0.5226815247939318 | −6.014289590682177 | 76.33383202878794 | 0.7683049674521053 |

The two OP/AC-only examples matched all three original scalars exactly.
The instrumented sine bodies differ from the main-cohort OP/AC readings by
at most 1.8 pV at the output; all four waveform metrics matched their original
instrumented observations exactly. R75 still fails the 2% distortion limit;
F75L passes the joint checks in this selected case. Reproduction success
preserves a physical failure when that is the observed result.

Pinned APM verifies actual model/W/L/m/nf/raw readbacks around each analysis.
The example also checks terminal domains, MOS KCL, frequency/time endpoints
and actual resistor voltage/current ratios. The 20-mV examples check output
transient KCL and infer capacitance from measured charge and voltage change:
both agree with 1 pF within 4.90 ppm. That charge check is distinct from
the private acquisition's direct capacitance-parameter readback.

The public recipe API emits `tran step stop`; the private recipe additionally
spells out the same maximum step. For these inputs the documented default is
the smaller of step and stop/50, giving the identical 100-ns bound. The replay
checks the actual grid and eight-cycle harmonic result. See the
[ngspice transient-analysis documentation](https://nmg.gitlab.io/ngspice-manual/analysesandoutputcontrol_batchmode/analyses/tran_transientanalysis.html).
The signal bodies retain the exact zero-volt capacitor sensor; public replays
filter private native-parameter diagnostic vectors unsupported by the public API.

## Evidence and limits

- [Selected scientific and clean-replay evidence](evidence.json)
- [Frozen confirmation policy](confirmation-plan.json), [function/resources](contract-v1.json) and [measurements](measurements-v1.json)
- [Main confirmation observations](data/confirmation.csv) and [saved physical units](data/realizations.csv)
- [Preselected 20-mV observations](data/selected-signal.csv)
- [Nominal design attempts](data/design-attempts.csv) and [development observations](data/development.csv)
- [Nominal dynamics and refinements](data/dynamics.csv)
- [Source transformations](data/source-transforms.csv) and [file manifest](data/manifest.json)

The nominal and .97-V settings were known; this confirmation uses fresh physical
coordinates under the documented retained-exposure check. Its intervals are
conditional on the APM source-transfer hypothesis. Six preselected waveform
coordinates do not estimate a dynamic population rate. Global process,
calibrated passive/spatial/temperature statistics, layout parasitics, real input
buffers and upstream supply implementations remain outside this package.
Neither manufacturing yield nor silicon qualification follows from these data.
The earlier [gain study](https://ds54e.github.io/analog-design-notes/studies/gain-operating-point-accuracy/),
[passive study](https://ds54e.github.io/analog-design-notes/studies/passive-sensitivity/) and [allocation study](https://ds54e.github.io/analog-design-notes/studies/active-allocation/)
retain their original definitions, URLs and immutable downloads.
