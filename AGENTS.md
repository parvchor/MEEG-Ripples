# Codex Usage Guide

This repository contains MATLAB scripts for analyzing MEEG ripple events. The raw datasets are very large and are **not** included in the repository. They reside on our VNC-based supercomputer and must be accessed there.

## Research goal
These scripts support analysis of high-frequency ripple events in intracranial recordings. They compute ripple statistics and generate diagnostic plots for a cohort of subjects.

## Dataset location
The BIDS-formatted recordings live on the VNC server at a path similar to:

```
/space/seh10/2/halgdev/projects/pchordiya/Cogigate/BIDS MEEG Data/COG_MEEG_EXP1_BIDS_SAMPLE
```

## Environment setup
1. Start MATLAB (R2021b or later is typical).
2. Add FieldTrip to the path and run `ft_defaults`:
   ```matlab
   addpath('/home/bqrosen/matlab/fieldtrip-20210825');
   ft_defaults;
   ```
3. Ensure the repository root is on the MATLAB path so that the `analysis`, `computation`, `plotting`, and `tools` folders are visible.

## Directory structure
- `analysis/` – High level scripts orchestrating ripple analysis and wrappers.
- `computation/` – Functions for event detection and statistic calculations.
- `plotting/` – Scripts and utilities for generating figures.
- `tools/` – Helper functions such as spike detection utilities located in `tools/LFP_Tools`.

## Entry points
Typical analysis begins with scripts in `analysis/`, for example `AnalyzeRipple.m` or `CortRipple_wrapper_MGH.m`. These scripts expect the dataset path above and may save figures and intermediate files within subject-specific directories.

## Usage notes
These scripts are normally executed on the VNC server due to dataset size. They rely on absolute paths (see above) and may require substantial RAM. Adapt path variables if running in a different environment.
