# sPHENIX-style p+p fast simulation

These scripts generate RHIC p+p collisions at sqrt(s) = 200 GeV and apply a fast, parametrized EMCal response. They do not run the full sPHENIX Geant4 detector simulation.

## Environment setup

Activate the existing Conda environment:

```bash
conda activate pythia
cd "/Users/sijanregmi/Documents/ChatGPT/PYTHIA"
```

The required packages are `pythia8mc` and `pyhepmc`. If needed:

```bash
python -m pip install pythia8mc pyhepmc
```

## 1. Generate p+p events

Run the PYTHIA8 generator:

```bash
python pp_rhic_200gev.py
```

Default output:

```text
pp_rhic_200gev.hepmc
```

Generate a larger sample with a reproducible seed:

```bash
python pp_rhic_200gev.py \
  --events 100000 \
  --seed 12345 \
  --output pp_rhic_200gev.hepmc
```

The generator uses sPHENIX-style minimum-bias settings: p+p beams, 200 GeV center-of-mass energy, inelastic soft QCD, CTEQ6L PDFs, `pTHatMin = 0.5 GeV`, and the sPHENIX K-factor setting.

## 2. Basic EMCal smearing

`fast_emcal_sim.py` reads the HepMC file and applies a simple response to final photons, electrons, and positrons:

```bash
python fast_emcal_sim.py \
  pp_rhic_200gev.hepmc \
  --output pp_rhic_200gev_emcal.csv
```

It applies EMCal acceptance, efficiency, energy smearing, and angular smearing. The default energy resolution is:

```text
sigma(E)/E = 3.0% + 15.4%/sqrt(E[GeV])
```

The output is a CSV file containing truth and smeared kinematics.

## 3. Detector-effects EMCal model

`fast_emcal_det_effects.py` is the more detailed, separate fast model. It adds approximate photon conversions, nearby-particle merging, tower thresholds, cluster thresholds, electronic noise, and cluster-level energy smearing:

```bash
python fast_emcal_det_effects.py \
  pp_rhic_200gev.hepmc \
  --output pp_rhic_200gev_emcal_clusters.csv
```

Useful options:

```bash
python fast_emcal_det_effects.py \
  pp_rhic_200gev.hepmc \
  --merge-dr 0.04 \
  --material-x0 0.30 \
  --tower-threshold 0.10 \
  --cluster-threshold 0.20 \
  --noise 0.02 \
  --output pp_rhic_200gev_emcal_clusters.csv
```

The output CSV contains cluster energy and direction, number of constituents, photon/electron/positron content, converted-photon count, truth cluster energy, and sampled noise.

## Script summary

| File | Purpose |
| --- | --- |
| `pp_rhic_200gev.py` | Generate PYTHIA8 p+p events and HepMC3 output |
| `fast_emcal_sim.py` | Basic final-particle EMCal smearing |
| `fast_emcal_det_effects.py` | Extended fast EMCal effects and clustering |

## Important limitations

These are analysis-level parameterizations, not a substitute for full sPHENIX Geant4. They do not reproduce detailed material geometry, shower shapes, conversion kinematics, electronics, calibrated clustering, pileup, or detector correlations. Validate any physics result against sPHENIX simulation or performance data before using it for publication.
