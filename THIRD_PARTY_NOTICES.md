# Third-party notices

The original study prose, numerical observations, finite circuit bodies and
example code were prepared for Analog Design Notes. No historical private
circuit archive, third-party figure or model package is redistributed here.

The native example retrieves [APM v5.0.0](https://github.com/ds54e/analog-process-models/tree/381517fda5107fabf98af7801d5a5103f38e230c),
commit `381517fda5107fabf98af7801d5a5103f38e230c`. APM's original code is Apache
2.0; its FreePDK45 model inputs retain their upstream notices and terms in that
checkout. This repository calls APM's public APIs and model wrappers; it does
not include the APM implementation or model files. ngspice 47 is separately
installed and retains its own BSD/other component notices.

Plots are generated with Matplotlib 3.9.4 and use DejaVu glyphs. The applicable
[DejaVu/Bitstream Vera font notice](LICENSES/DejaVu.txt) is retained. External
Python dependencies listed in the example are obtained separately, with their
own licenses.

The site uses Quarto 1.10.18. Its [MIT license](LICENSES/Quarto-MIT.txt) and
[copyright notice](LICENSES/Quarto-COPYRIGHT.txt) are retained. The generated
runtime assets retain their upstream headers; additional component licenses
and retrieval sources are listed in [the asset notices](LICENSES/ASSET_SOURCES.json).
MathJax 4.1.3 is loaded from the public jsDelivr npm endpoint, under its Apache
2.0 terms. It retrieves its matching New Computer Modern web fonts and speech
resources from that public CDN. Quarto also requests Cloudflare's public
cdnjs ES6 polyfill endpoint. These separately retrieved components retain
their own upstream terms. They render the page and are not used for scientific
acquisition. Ordinary body text uses system fonts.

The guides link to MIT OpenCourseWare, Analog Devices' resistance tutorial and the pinned APM source decision as
technical references. No source illustrations or extended quoted passages
are copied. Those linked works retain their own terms.
