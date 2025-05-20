# Vortex-Legendre Method for Digital Holographic Microscopy (DHM)

This repository contains the implementation of the **Vortex-Legendre Method** for phase aberration correction in Digital Holographic Microscopy (DHM) of intricate bio and non-biological samples. The code supports the efficient correction of low- and high-order aberrations in off-axis DHM, enabling high-speed, accurate quantitative phase imaging (QPI).

---

## Overview

Digital Holographic Microscopy (DHM) is a powerful technique for label-free imaging, particularly useful in biomedical applications. However, the accuracy of DHM measurements often depends on compensating for phase aberrations caused by optical setups and sample properties. However, when imaging intricate biological samples, such as those involving high spatial frecuency objects or rapidly changing contrast details, conventional phase compensation method usually fail. Therefore, this method is specially suited to properly compensate and image these kinds of samples.

This repository implements the **Vortex-Legendre Method**, a two-step computational approach that:
1. Utilizes **numerical optical vortices** for precise sub-pixel localization of diffraction orders, enabling robust tilt aberration correction.
2. Applies **Legendre Polynomial Fitting (LPF)** to address residual higher-order aberrations, ensuring fully compensated phase maps.

---

## Features

- **Suited for intricate sample imaging**: Specially suited for intricate samples where conventional fail to properly compensate.
- **Fast Processing**: Achieves phase aberration correction in ~2 second for conventional format holograms on a standard laptop.
- **High Accuracy**: Validated with USAF and star phase targets, yielding phase thickness measurements within a 5.3% error margin.
- **Generalizable**: Works for both telecentric and non-telecentric configurations and with complex biological samples.
- **Streamlined Workflow**: Eliminates the need for multi-shot acquisitions, iterative procedures, and additional optical components.

---

## Why It Is Useful

The Vortex-Legendre Method addresses common challenges in DHM by:
- Reducing computational overhead, enabling real-time imaging.
- Improving measurement accuracy with robust compensation of both low- and high-order aberrations for complex bio and non-biological samples.
- Enhancing QPI capabilities for a broad range of applications, including biomedical imaging, material analysis, and metrology.

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/sophiaresearchlaboratory/VortexDHM.git
2. Install dependencies (e.g., MATLAB 2024b or compatible versions).

Open the provided .m files in MATLAB to explore and run the implementation.

---

## How to Use

The following is a step-by-step guide to processing holograms using the Vortex-Legendre method and visualizing the corrected phase maps. 
Once you download the files, it is important to save them in the same folder.
### 1. Preparation of recorded holograms
The Vortex-Legendre method requires recorded holograms, which can come from various sources, such as USAF targets, target stars, or biological samples. 
These holograms can be of telecentric or non-telecentric configuration.
Example holograms can be found in the 'holograms' folder. This folder includes four examples:
-	USAF (Telecentric configuration, "usaf.bmp")
-	StarTarget (Telecentric configuration, "starTarget.bmp")
-	Red Blood Cells (No telecentric configuration, 30 mm offset, "redBlood30mm.tiff")
-	Star target (No telecentric configuration, 40 mm offset, "-4cm_20x_star.tiff")
### 2. Using the Vortex-Legendre Method
The core of the Vortex-Legendre method is implemented in a special function in the repository. 
The repository contains a main script (VortexLegendreCompensation_main), whose execution performs the process of load hologram, processing and visualizing the phase map compensation. 
Before running this script, the hologram to be compensated must be specified. 
Consider the physical parameters with which the hologram was recorded (these parameters are specified for the example holograms in the code such as: wavelength - 632.8 nm and pixel size - 3.75 um). 
Once these parameters have been verified, you can run the main script directly.
The following functions are executed:

- **functions_evaluation** This function allows the selection of the order (+1 DO).

- **Vortex compensation** This algorithm allows to find the frequencies for the tilt compensation using the Hilbert 2D transform (function hilbertTransform2D) and the residual theorem using vortex interpolation. At the end of these functions, the tilt aberration compensated phase map can be displayed.

- **Square Legendre fitting**. This algorithm calculates the first 15th square Legendre polynomials. Then, the phase compensation is performed by calculating the coefficients corresponding to each polynomial. Finally, the final compensation is performed. 

### 3. Displaying the compensated phase map
After the hologram is processed, you can visualize four figures including the spectrum of the hologram, the filtered spectrum hologram, the corrected tilt phase map and the complete map compensated with a Legendre polynomial. 

---

## Authors
Karina Ortega-Sánchez - Universidad Politécnica de Tulancingo - karina.ortega2415006@upt.edu.mx 

Rene Restrepo - Universidad EAFIT - rrestre6@eafit.edu.co

Alfonso Padilla-Vivanco - Universidad Politécnica de Tulancingo

Raúl Castañeda - Universidad EAFIT

Ana Doblas - University of Massachusetts Dartmouth - adoblas@umassd.edu

Carlos Trujillo - Universidad EAFIT - catrujilla@eafit.edu.co

---

## Citation
If you use this code in your research, please cite the corresponding publication:

Karina Ortega-Sánchez, Rene Restrepo, Raúl Castañeda, Alfonso Padilla-Vivanco, Ana Doblas, Carlos Trujillo, Intricate Quantitative Phase Imaging via Vortex-Legendre High-Order Phase Compensation. Submitted to Optics and Lasers in Engineering (2025). DOI: [Not ready yet]

---

## License
This code is distributed under the MIT License. See the LICENSE file for details.

---

## Acknowledgments
We thank the Optics and Photonics lab at EAFIT for their valuable discussions and foundational contributions to the development of this work.

---

## Contact
For questions or further assistance, please reach out to: Carlos Trujillo, René Restrepo and Karina Ortega-Sanchez.
Email: catrujilla@eafit.edu.co, rrestre6@eafit.edu.co, karina.ortega2415006@upt.edu.mx
