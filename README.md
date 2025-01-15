# Vortex-Legendre Method for Digital Holographic Microscopy (DHM)

This repository contains the implementation of the **Vortex-Legendre Method** for phase aberration correction in Digital Holographic Microscopy (DHM). The code supports the efficient correction of low- and high-order aberrations in off-axis DHM, enabling high-speed, accurate quantitative phase imaging (QPI).

---

## Overview

Digital Holographic Microscopy (DHM) is a powerful technique for label-free imaging, particularly useful in biomedical applications. However, the accuracy of DHM measurements often depends on compensating for phase aberrations caused by optical setups and sample properties.

This repository implements the **Vortex-Legendre Method**, a two-step computational approach that:
1. Utilizes **numerical optical vortices** for precise sub-pixel localization of diffraction orders, enabling robust tilt aberration correction.
2. Applies **Legendre Polynomial Fitting (LPF)** to address residual higher-order aberrations, ensuring fully compensated phase maps.

---

## Features

- **Fast Processing**: Achieves phase aberration correction in ~1 second for 960×1280 px holograms on a standard laptop.
- **High Accuracy**: Validated with USAF and star phase targets, yielding phase thickness measurements within a 5.3% error margin.
- **Generalizable**: Works for both telecentric and non-telecentric configurations and with biological samples.
- **Streamlined Workflow**: Eliminates the need for multi-shot acquisitions, iterative procedures, and additional optical components.

---

## Why It Is Useful

The Vortex-Legendre Method addresses common challenges in DHM by:
- Reducing computational overhead, enabling real-time imaging.
- Improving measurement accuracy with robust compensation of both low- and high-order aberrations.
- Enhancing QPI capabilities for a broad range of applications, including biomedical imaging, material analysis, and metrology.

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/sophiaresearchlaboratory/VortexDHM.git
