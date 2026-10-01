"""Advanced fast EMCal response for sPHENIX-like studies.

Separate from fast_emcal_sim.py: adds approximate photon conversions, tower-like
clustering/overlap, noise, thresholds, and cluster-level energy smearing.
This is a tunable analysis model, not a Geant4 replacement.
"""

import argparse
import csv
import math
import random

import pyhepmc


def delta_phi(a: float, b: float) -> float:
    return (a - b + math.pi) % (2.0 * math.pi) - math.pi


def delta_r(a_eta: float, a_phi: float, b_eta: float, b_phi: float) -> float:
    return math.hypot(a_eta - b_eta, delta_phi(a_phi, b_phi))


def frac_resolution(energy: float) -> float:
    return math.sqrt((0.154 / math.sqrt(max(energy, 1e-9))) ** 2 + 0.030**2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="?", default="pp_rhic_200gev.hepmc")
    parser.add_argument("-o", "--output", default="pp_rhic_200gev_emcal_clusters.csv")
    parser.add_argument("--seed", type=int, default=97531)
    parser.add_argument("--eta-max", type=float, default=1.1)
    parser.add_argument("--merge-dr", type=float, default=0.04, help="cluster merge radius")
    parser.add_argument("--material-x0", type=float, default=0.30, help="effective pre-EMCal material in X0")
    parser.add_argument("--efficiency", type=float, default=1.0)
    parser.add_argument("--tower-threshold", type=float, default=0.10, help="GeV")
    parser.add_argument("--cluster-threshold", type=float, default=0.20, help="GeV")
    parser.add_argument("--noise", type=float, default=0.02, help="GeV RMS per cluster")
    args = parser.parse_args()

    if args.material_x0 < 0 or not 0 <= args.efficiency <= 1:
        raise ValueError("material-x0 must be nonnegative and efficiency must be in [0, 1]")

    rng = random.Random(args.seed)
    n_events = 0
    n_particles = 0
    n_clusters = 0

    columns = [
        "event", "cluster", "energy_gev", "pt_gev", "eta", "phi",
        "n_constituents", "n_photons", "n_electrons", "n_positrons",
        "n_converted_photons", "truth_energy_gev", "noise_gev",
    ]

    with pyhepmc.io.ReaderAscii(args.input) as reader, open(args.output, "w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=columns)
        writer.writeheader()

        while not reader.failed():
            event = reader.read()
            if reader.failed():
                break
            n_events += 1
            deposits = []

            for particle in event.particles:
                if particle.status != 1 or particle.pid not in (22, 11, -11):
                    continue
                n_particles += 1
                mom = particle.momentum
                energy = max(mom.e, 0.0)
                eta = mom.eta()
                phi = mom.phi()
                if abs(eta) > args.eta_max or rng.random() > args.efficiency:
                    continue
                if energy < args.tower_threshold:
                    continue

                converted = False
                if particle.pid == 22:
                    # High-energy photon conversion probability through t X0.
                    converted = rng.random() < 1.0 - math.exp(-7.0 * args.material_x0 / 9.0)

                deposits.append({
                    "energy": energy, "eta": eta, "phi": phi,
                    "pid": particle.pid, "converted": converted,
                })

            clusters = []
            for deposit in sorted(deposits, key=lambda d: d["energy"], reverse=True):
                target = None
                for cluster in clusters:
                    if delta_r(deposit["eta"], deposit["phi"], cluster["eta"], cluster["phi"]) < args.merge_dr:
                        target = cluster
                        break
                if target is None:
                    clusters.append({
                        "energy": deposit["energy"],
                        "eta": deposit["eta"],
                        "phi": deposit["phi"],
                        "deposits": [deposit],
                    })
                else:
                    old_e = target["energy"]
                    new_e = old_e + deposit["energy"]
                    target["eta"] = (old_e * target["eta"] + deposit["energy"] * deposit["eta"]) / new_e
                    target["phi"] = target["phi"] + deposit["energy"] * delta_phi(deposit["phi"], target["phi"]) / new_e
                    target["energy"] = new_e
                    target["deposits"].append(deposit)

            for cluster_number, cluster in enumerate(clusters):
                truth_energy = cluster["energy"]
                noise = rng.gauss(0.0, args.noise)
                energy = max(0.0, rng.gauss(truth_energy, truth_energy * frac_resolution(truth_energy)) + noise)
                if energy < args.cluster_threshold:
                    continue
                deps = cluster["deposits"]
                nphotons = sum(d["pid"] == 22 for d in deps)
                nelectrons = sum(d["pid"] == 11 for d in deps)
                npositrons = sum(d["pid"] == -11 for d in deps)
                nconverted = sum(d["converted"] for d in deps)
                writer.writerow({
                    "event": n_events, "cluster": cluster_number,
                    "energy_gev": energy, "pt_gev": energy / math.cosh(cluster["eta"]),
                    "eta": cluster["eta"], "phi": cluster["phi"],
                    "n_constituents": len(deps), "n_photons": nphotons,
                    "n_electrons": nelectrons, "n_positrons": npositrons,
                    "n_converted_photons": nconverted,
                    "truth_energy_gev": truth_energy, "noise_gev": noise,
                })
                n_clusters += 1

    print(f"Read {n_events} events; input EM particles = {n_particles}")
    print(f"Output EMCal clusters = {n_clusters}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
