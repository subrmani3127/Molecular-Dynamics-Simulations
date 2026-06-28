Biomolecular Dynamics & Conformational Analysis Workflow

This repository contains a collection of computational pipelines, analysis protocols, and automated scripts designed for long-timescale molecular dynamics (MD) simulations and biophysical characterization of protein-ligand systems.

The framework integrates physics-based MD trajectories generated via GROMACS with data-driven statistical mechanics, dimensionality reduction, and structural analysis tools to quantify macromolecular conformational landscapes and binding kinetics.
📂 Repository Structure & Key Components
1. MD Simulation & Core Workflows

    Simulation Setup in GROMACS: Complete configuration files (.mdp), topology generation scripts, and workflows for system solvation, neutralization, energy minimization, and equilibration (NVT/NPT) leading up to production MD runs.

    TCL Scripting: Custom VMD (Visual Molecular Dynamics) core scripts automated for structural alignments, trajectory processing, dynamic selection handling, and frame-by-frame data extraction.

2. Dimensionality Reduction & Conformational Landscapes

    Principal Component Analysis (PCA): Scripts to project high-dimensional trajectory coordinates onto low-dimensional collective variables (essential dynamics). Used to identify dominant global motions, calculate covariance matrices, and visualize free energy landscapes (FEL).

    Markov State Models (MSMs): Framework for clustering structural ensembles into microstates, constructing transition probability matrices, and calculating kinetic macrostates to uncover long-timescale biophysical processes from shorter parallel MD runs.

3. Advanced Biophysical & Thermodynamic Analysis

    Minimum Distance Distribution Function (MDDF) & Kirkwood-Buff Integrals (KBI): Mathematical and spatial distribution tools used to quantify local solution theory and structural coordination of solvents/cosolutes around biomolecular domains.

    Schlitter Entropy & Secondary Structure Propensity: Implementations to estimate configurational entropy from covariance matrices (Schlitter formula) alongside time-resolved secondary structure tracking to evaluate protein loop flexibility and stability.

    Water Analysis in Protein-Ligand Systems: Targeted analytical scripts to investigate active-site hydration dynamics, water bridging networks, hydrogen bonding lifetimes, and active-site occupancy.

4. Kinetics & Binding Lifetime Scripts

    residence_time.py / residence_time.ipynb: Automated Python script and interactive Jupyter notebook designed to calculate ligand residence times, unbinding kinetics, and continuous correlation functions from trajectory data.

🚀 Getting Started
Prerequisites

Ensure your local or high-performance computing (HPC) environment has the following core tools installed:

    MD Engines: GROMACS (2020+ recommended)

    Visualization & Scripting: VMD (with tclsh)

    Python Stack: Python 3.8+, NumPy, SciPy, MDAnalysis, MDTraj, PyEMMA (for MSMs), Scikit-learn (for PCA), Jupyter Lab/Notebook

Sample Usage: Computing Residence Time

To execute the standalone script for tracking ligand binding kinetics:
Bash

python residence_time.py --trajectory prod_run.xtc --topology sys.tpr --ligand LIG
