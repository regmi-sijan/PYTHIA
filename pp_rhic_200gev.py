"""Standalone sPHENIX-style PYTHIA8 minimum-bias p+p generator."""

import argparse
import os

# Use the XML data bundled with pythia8mc, not an older system installation.
os.environ.pop("PYTHIA8DATA", None)

import pyhepmc
import pythia8mc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--events", type=int, default=10)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("-o", "--output", default="pp_rhic_200gev.hepmc")
    args = parser.parse_args()

    if args.events < 1:
        raise ValueError("--events must be positive")
    if not 1 <= args.seed <= 900000000:
        raise ValueError("--seed must be between 1 and 900000000")

    pythia = pythia8mc.Pythia()

    # Match sPHENIX generators/PHPythia8/phpythia8_mb.cfg.
    commands = [
        "Beams:idA = 2212",
        "Beams:idB = 2212",
        "Beams:eCM = 200.",
        "Init:showChangedSettings = on",
        "Next:numberShowInfo = 0",
        "PDF:pSet = 7",                    # CTEQ6L
        "SoftQCD:inelastic = on",         # minimum-bias inelastic sample
        "MultipartonInteractions:Kfactor = 0.5",
        "PhaseSpace:pTHatMin = 0.5",
        "Random:setSeed = on",
        f"Random:seed = {args.seed}",
    ]
    for command in commands:
        if not pythia.readString(command):
            raise RuntimeError(f"PYTHIA rejected setting: {command}")

    if not pythia.init():
        raise RuntimeError("PYTHIA initialization failed")

    # Preserve the full decay ancestry while exporting the event to HepMC3.
    converter = pythia8mc.Pythia8ToFilteredHepMC3()
    generated = 0
    with pyhepmc.io.WriterAscii(args.output) as writer:
        event_number = 0
        while generated < args.events:
            if not pythia.next():
                print("PYTHIA event generation failed; retrying")
                continue

            event_number += 1
            generated += 1
            hepmc_event = pyhepmc.GenEvent()
            converter.fill_next_event(hepmc_event, event_number, pythia=pythia)
            writer.write_event(hepmc_event)

            final_state = [p for p in pythia.event if p.isFinal()]
            charged_final_state = [p for p in final_state if p.isCharged()]
            print(
                f"Event {generated:6d}: final = {len(final_state):3d}, "
                f"charged final = {len(charged_final_state):3d}"
            )

    print(f"Wrote {generated} events to {args.output}")
    pythia.stat()


if __name__ == "__main__":
    main()
