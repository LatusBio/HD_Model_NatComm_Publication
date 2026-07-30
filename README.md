# Huntington's Disease Therapeutic Model

[![DOI](https://img.shields.io/badge/DOI-10.64898%2F2026.01.06.697909-blue)](https://www.biorxiv.org/content/10.64898/2026.01.06.697909v1)

An interactive web application for simulating somatic CAG repeat expansion in Huntington's disease and modeling the potential therapeutic impact of interventions targeting somatic instability.

**Live app:** [www.LatusHDmodel.com](www.latushdmodel.com)

## Overview

Huntington's disease (HD) is caused by an expanded CAG trinucleotide repeat in the *HTT* gene. While the inherited (germline) repeat length determines disease risk, **somatic expansion** of the CAG repeat in striatal medium spiny neurons (MSNs) over a patient's lifetime is increasingly recognized as a key driver of disease onset and progression. This application allows users to simulate this process and explore how therapeutic interventions that reduce somatic expansion rate could delay neuronal dysfunction and HD milestones.

The model implements the two-phase linear expansion rate framework reported by [Handsaker et al. (2025)](https://doi.org/10.1016/j.cell.2024.11.038), with extensions for therapeutic intervention modeling described in [Simpson and Ranum et al. (2026)](https://www.biorxiv.org/content/10.64898/2026.01.06.697909v1).

## Features

- **Two simulation engines:**
  - **Handsaker CTMC model (recommended):** Continuous-time Markov chain simulation using matrix exponentiation, parameterized by expansion probability (p_exp = 0.676) and two-phase linear mutation rates from Handsaker et al.
  - **Expansion-to-contraction ratio model:** Event-sampling approach with configurable expansion:contraction bias (default 6:1)

- **Therapeutic intervention modeling:**
  - Configurable age at intervention
  - Adjustable MSH3 knockdown percentage (reduction in somatic expansion rate)
  - Adjustable cell transduction rate (fraction of neurons targeted)
  - Side-by-side comparison of treated vs. untreated trajectories

- **Outputs:**
  - Animated visualization of CAG repeat expansion across individual neurons over time
  - Mean CAG repeat length trajectories (treated vs. untreated)
  - Healthy cell fraction over time with HD-ISS stage progression
  - All plots and underlying data downloadable (PNG, GIF, CSV)

- **Configurable parameters:**
  - Germline CAG repeat length (36–56)
  - Number of MSNs to simulate (100–3,000)
  - Healthy cell CAG threshold (default 150, based on transcriptional dysregulation onset)
  - Maximum simulated lifetime

## Safety & Disclaimers

This application includes safeguards appropriate for a tool that displays disease progression projections:

- **Gated entry:** A disclaimer overlay requires users to acknowledge that this is a research/educational tool before accessing the application
- **Persistent crisis support bar:** Always-visible footer with direct links to the 988 Suicide & Crisis Lifeline, HDSA Help Line (1-800-345-HDSA), and Crisis Text Line
- **Contextual framing:** Caveats are placed directly below each visualization explaining that projections are population-level estimates and individual trajectories will vary
- **Therapeutic benefit disclaimer:** Clearly states that modeled therapeutic benefit is hypothetical and no approved therapy has demonstrated these specific outcomes
- **Genetic counselor referral:** Links to HDSA Find Help, the NSGC Find a Genetic Counselor directory, and HDSA Centers of Excellence

## Installation

### Requirements

- R (≥ 4.0)
- Required R packages (installed automatically on first run):

```r
shiny, ggplot2, cowplot, viridis, R6, reshape2, gifski, av,
RCurl, httr, remotes, shinyWidgets, gganimate, expm
```

### Running locally

```bash
git clone https://github.com/LatusBio/HD_Model_NatComm_Publication.git
cd HD_Model_NatComm_Publication
```

Open R or RStudio and run:

```r
library(shiny)
runApp(".")
```

The app will open in your default web browser.

### Data files

The following data files are required and should be in the app directory:

| File | Description |
|------|-------------|
| `HD_ISS_Thresholds_ByAge.csv` | HD Integrated Staging System thresholds by age |
| `UHDRS_Total_Functional_Capacity_AIscraped_5yearRange_individualYears.csv` | UHDRS Total Functional Capacity reference data |
| `meanMotorScore_AIscraped_5yearRange_individualYears.csv` | Mean motor score reference data |

## Model Details

### Two-phase linear expansion rate

The somatic expansion rate follows a two-phase linear model (Handsaker et al., 2025):

```
expansion_rate(CAG) = r1 × max(CAG - T1, 0) + r2 × max(CAG - T2, 0)
```

Where:
- **Phase 1:** r1 = 3.51% per year for CAG > T1 (T1 = 33.5)
- **Phase 2:** r2 = 57.6% per year for CAG > T2 (T2 = 72.2 in CTMC model, 55 in event-sampling model)

This creates a slow initial expansion followed by rapid acceleration once repeats exceed ~55–72 CAG.

### Handsaker CTMC implementation

The recommended simulation engine constructs a continuous-time Markov chain (CTMC) rate matrix Q where each CAG state has:
- Expansion rate = total_mutation_rate × p_exp
- Contraction rate = total_mutation_rate × (1 − p_exp)

The one-year transition probability matrix is computed as P = exp(Q) using matrix exponentiation. Each neuron's CAG length is propagated year-by-year by sampling from the appropriate row of P.

Default p_exp = 0.676 ± 0.011 (~2.1:1 expansion:contraction ratio), from Handsaker et al. Supplementary Fig. SN4.4.

### Therapeutic intervention

At the specified intervention age, a fraction of neurons (transduction rate) receive the therapeutic modification. In treated neurons, the mutation rate is reduced by the MSH3 knockdown percentage, producing a modified rate matrix Q_treated with lower expansion rates. Treated and untreated neurons are then simulated independently for the remainder of the lifetime.

## Citation

If you use this tool in your research, please cite:

```
Simpson et al. (2026). [Title]. bioRxiv.
https://doi.org/10.64898/2026.01.06.697909
```

And the foundational expansion rate model:

```
Handsaker et al. (2025). Long-read sequencing of Huntington disease
brains reveals complex, repeat-length-dependent somatic CAG repeat
expansion and DNA damage. Cell.
https://doi.org/10.1016/j.cell.2024.11.038
```

## Contact

Developed and maintained by [Latus Bio](https://www.latusbio.com).

Correspondence: paul.ranum@latusbio.com

## License

[Specify your license here]
