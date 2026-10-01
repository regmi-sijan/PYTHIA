"""Fast, parametrized sPHENIX-like EMCal response for HepMC3 events.

This is an analysis-level model for final photons and e+/e-. It does not
simulate Geant4 showers, conversions, clustering, or electronics.
"""

import argparse
import csv
import math
import random

import pyhepmc


def resolution_fraction(energy_gev: float, stochastic: float, constant: float) -> float:
    """Return sigma(E)/E using a stochastic-plus-constant model."""
    return math.sqrt((stochastic / math.sqrt(energy_gev)) ** 2 + constant**2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="?", default="pp_rhic_200gev.hepmc")
    parser.add_argument("-o", "--output", default="pp_rhic_200gev_emcal.csv")
    parser.add_argument("--seed", type=int, default=24680)
    parser.add_argument("--eta-max", type=float, default=1.1)
    parser.add_argument("--efficiency", type=float, default=1.0)
    parser.add_argument("--sigma-eta", type=float, default=0.005)
    parser.add_argument("--sigma-phi", type=float, default=0.005)
    args = parser.parse_args()

    if not 0.0 <= args.efficiency <= 1.0:
        raise ValueError("--efficiency must be between 0 and 1")

    rng = random.Random(args.seed)
    n_events = 0
    n_candidates = 0
    n_accepted = 0

    with pyhepmc.io.ReaderAscii(args.input) as reader, open(args.output, "w", newline="") as output:
        columns = [
            "event", "barcode", "pid", "particle",
            "truth_e_gev", "smeared_e_gev", "truth_pt_gev",
            "smeared_pt_gev", "truth_eta", "smeared_eta",
            "truth_phi", "smeared_phi", "energy_sigma_gev",
        ]
        writer = csv.DictWriter(output, fieldnames=columns)
        writer.writeheader()

        while not reader.failed():
            event = reader.read()
            if reader.failed():
                break
            n_events += 1

            for particle in event.particles:
                if particle.status != 1 or particle.pid not in (22, 11, -11):
                    continue
                n_candidates += 1

                momentum = particle.momentum
                truth_e = max(momentum.e, 1e-9)
                truth_eta = momentum.eta()
                if abs(truth_eta) > args.eta_max or rng.random() > args.efficiency:
                    continue

                # sPHENIX EMCal prototype: 3.0% constant and 15.4%/sqrt(E).
                frac_sigma = resolution_fraction(truth_e, 0.154, 0.030)
                sigma_e = frac_sigma * truth_e
                smeared_e = max(0.0, rng.gauss(truth_e, sigma_e))
                smeared_eta = rng.gauss(truth_eta, args.sigma_eta)
                smeared_phi = rng.gauss(momentum.phi(), args.sigma_phi)
                smeared_pt = smeared_e / math.cosh(smeared_eta)

                writer.writerow({
                    "event": n_events,
                    "barcode": particle.id,
                    "pid": particle.pid,
                    "particle": {22: "photon", 11: "electron", -11: "positron"}[particle.pid],
                    "truth_e_gev": truth_e,
                    "smeared_e_gev": smeared_e,
                    "truth_pt_gev": momentum.pt(),
                    "smeared_pt_gev": smeared_pt,
                    "truth_eta": truth_eta,
                    "smeared_eta": smeared_eta,
                    "truth_phi": momentum.phi(),
                    "smeared_phi": smeared_phi,
                    "energy_sigma_gev": sigma_e,
                })
                n_accepted += 1

    print(f"Read {n_events} events")
    print(f"EMCal candidates: {n_candidates}; accepted: {n_accepted}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
