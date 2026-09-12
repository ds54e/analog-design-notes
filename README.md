# Analog Design Notes

Read the [research site](https://ds54e.github.io/analog-design-notes/) or the
[Japanese overview](https://ds54e.github.io/analog-design-notes/ja/gain-operating-point-accuracy.html).

The first study compares gain, output operating point and untrimmed accuracy
in four finite amplifiers. It includes a frozen 464-coordinate paired study,
3,712 condition observations, exact circuit bodies, compact data and saved
Local physical parameters. Model limitations and valid failures are retained.

Start with the [example README](examples/gain-operating-point-accuracy/README.md)
to regenerate graphs or replay selected circuits. Native replay uses separately
retrieved APM v5.0.0 and ngspice 47. The static site build never runs SPICE.

Original prose/plots/data: CC BY 4.0. Original code/circuit examples: MIT.
See [LICENSE.md](LICENSE.md) and [third-party notices](THIRD_PARTY_NOTICES.md).
This first release has internal AI checks, without prior human or peer review.

The editable canonical manuscripts are maintained in the private research
repository. This repository receives an explicit selected publication set.
Reader corrections belong in Issues; they are applied to that canonical source
before export. The immutable research tag is
`gain-operating-point-accuracy-v1.0`.

Maintainers render with Quarto **1.10.18**:

```bash
python3 tools/build_site.py
quarto render site
python3 tools/check_site.py --site site/_site
```

The `Publish Pages` workflow is dispatched on `main` after a selected release
is ready. It builds only this public repository, records its exact commit in
every page and `revision.json`, and deploys the Pages artifact. No private
repository, model acquisition, credentials file or research runtime is bundled.
