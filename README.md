# pySimFrac

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.cageo.2024.105665-blue.svg)](https://doi.org/10.1016/j.cageo.2024.105665)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue)](https://opensource.org/licenses/BSD-3-Clause)

**pySimFrac** is a Python library for generating and analyzing 3D synthetic rough fracture surfaces. It supports stochastic generation using multiple methods and provides built-in tools for geostatistical analysis and integration with pore-scale and discrete fracture network (DFN) simulators.

🔗 **GitHub**: [https://github.com/lanl/pySimFrac](https://github.com/lanl/pySimFrac)  
📖 **Paper**: Guiltinan et al., *Computers & Geosciences* (2024) – [DOI: 10.1016/j.cageo.2024.105665](https://doi.org/10.1016/j.cageo.2024.105665)  
📚 **Docs**: [https://lanl.github.io/pySimFrac](https://lanl.github.io/pySimFrac)

---

# Documentation 
Pysimfrac documentation can be found [here](https://lanl.github.io/pySimFrac/)

    https://lanl.github.io/pySimFrac/

Or in the pysimfrac.pdf document in the main directory 



## 🔍 Features

- Three generation methods: **Spectral**, **Gaussian**, **Box**
- Generate realistic 3D fracture geometries with periodic boundary conditions
- Geostatistical analysis tools: moments, autocorrelation, PDFs, variograms
- Import real surface profilometry data
- Seamless integration with:
  - [`dfnWorks`](https://github.com/lanl/dfnWorks) – DFN modeling
  - [`MF-LBM`](https://github.com/lanl/MF-LBM) – multiphase lattice Boltzmann
- Visualize surface topographies, aperture fields, and more

---

## 📦 Installation

### ✅ Requirements

- Python ≥ 3.8
- `pip`, `setuptools ≥ 64`, and `wheel`

Install required libraries:

```bash
pip install numpy scipy matplotlib seaborn scikit-gstat vedo
```

### 🔧 Build and Install from Source

```bash
git clone https://github.com/lanl/pySimFrac.git
cd pySimFrac/src

# Install build tool
pip install build

# Build and install
python -m build
pip install dist/pysimfrac-1.1-py3-none-any.whl
```

---

## 🧪 Quick Start Example

```python
from pysimfrac.general.simFrac import SimFrac

# Create a fracture using the spectral method
fracture = SimFrac(h=0.01, lx=3, ly=1, method="spectral", units="mm")

# Set generation parameters
fracture.params["H"]["value"] = 0.5
fracture.params["mean-aperture"]["value"] = 0.5
fracture.params["roughness"]["value"] = 0.5

# Generate and visualize
fracture.create_fracture()
fracture.plot_surfaces()
fracture.plot_aperture_field()
```

---

## 🧠 Scientific Background

pySimFrac enables generation of synthetic fractures using:
- **Spectral methods**: self-affine fractal surfaces with tunable Hurst exponent, anisotropy, and mismatch
- **Convolution methods**: Gaussian and box kernel smoothing for smooth, anisotropic surfaces

It includes statistical tools to analyze surfaces and apertures:
- Autocorrelation functions
- Probability density functions
- Variograms using SciKit-GStat

pySimFrac supports integration with:
- **MP-LBM** for single-phase flow
- **dfnWorks** for embedding variable aperture surfaces into DFNs

See full details in the [Computers & Geosciences paper](https://doi.org/10.1016/j.cageo.2024.105665).

---

## 📘 Documentation & Examples

- 📚 Full documentation: [https://lanl.github.io/pySimFrac](https://lanl.github.io/pySimFrac)
- 📓 Example notebooks: [GitHub/examples](https://github.com/lanl/pySimFrac/tree/main/examples)

---

## 📜 Citation

If you use `pySimFrac` in your research, please cite:

> **Guiltinan, E.**, Santos, J.E., Purswani, P., Hyman, J.D. (2024).  
> *pySimFrac: A Python library for synthetic fracture generation and analysis*.  
> Computers & Geosciences, 191, 105665.  
> DOI: [10.1016/j.cageo.2024.105665](https://doi.org/10.1016/j.cageo.2024.105665)

**BibTeX**:
```bibtex
@article{guiltinan2024pysimfrac,
  title={pySimFrac: A Python library for synthetic fracture generation and analysis},
  author={Guiltinan, Eric and Santos, Javier E. and Purswani, Prakash and Hyman, Jeffrey D.},
  journal={Computers \& Geosciences},
  volume={191},
  pages={105665},
  year={2024},
  publisher={Elsevier},
  doi={10.1016/j.cageo.2024.105665}
}
```

---

## 👩‍💻 Authors

- Jeffrey Hyman – jhyman@lanl.gov  
- Prakash Purswani – ppurswani@lanl.gov  
- Eric Guiltinan – eric.guiltinan@lanl.gov  
- Javier Santos – jesantos@lanl.gov

---

## 📄 License

This project is licensed under the **BSD 3-Clause License**. See the [LICENSE](https://github.com/lanl/pySimFrac/blob/main/LICENSE) file for details.


# Open-Source License

This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


---

## 📞 Contact

For questions or technical issues, contact:
- Eric Guiltinan – eric.guiltinan@lanl.gov
- Jeffrey Hyman – jhyman@lanl.gov




# Open-Source License

This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

