#!/usr/bin/env python3
"""
mace_server.py — MACE-MP-0 energy + forces server

Called by atlas.exe in --once mode:
    python mace_server.py --once request.json response.json

Request JSON:  {"N":13, "species":["Au","Au",...], "pos":[x0,y0,z0,...]}
Response JSON: {"energy":-42.3, "forces":[fx0,fy0,fz0,...]}
             or {"error":"message"}

Units: eV and eV/Å (atlas.exe converts to K internally)
"""
import sys
import json
import argparse
import numpy as np
import os

def check_mace_installed():
    try:
        import mace
        return True
    except ImportError:
        return False

def get_calculator(model="medium"):
    """Load MACE-MP-0 calculator. Downloads model on first use (~50MB)."""
    try:
        from mace.calculators import mace_mp
        calc = mace_mp(model=model, dispersion=False, default_dtype="float64", device="cpu")
        return calc
    except Exception as e:
        raise RuntimeError(f"Failed to load MACE-MP-0: {e}\n"
                           f"Install with: pip install mace-torch")

def evaluate(calc, N, species, pos_flat):
    """Evaluate energy and forces for a cluster."""
    from ase import Atoms
    import numpy as np

    pos = np.array(pos_flat).reshape(N, 3)
    # Add 30Å vacuum box so periodic boundary conditions don't interact
    box_size = max(np.max(pos) - np.min(pos) + 30.0, 40.0)
    atoms = Atoms(
        symbols=species,
        positions=pos + box_size/2,  # centre in box
        cell=[box_size, box_size, box_size],
        pbc=False
    )
    atoms.calc = calc
    energy = atoms.get_potential_energy()    # eV
    forces = atoms.get_forces()              # eV/Å, shape (N,3)
    return float(energy), forces.flatten().tolist()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--once", action="store_true",
                        help="Single-shot mode: read req_file, write res_file")
    parser.add_argument("req_file",  nargs="?", default="_atlas_mace_req.json")
    parser.add_argument("res_file",  nargs="?", default="_atlas_mace_res.json")
    parser.add_argument("--model",   default="medium",
                        help="MACE model: small|medium|large (default medium)")
    args = parser.parse_args()

    # Check dependencies
    if not check_mace_installed():
        err = json.dumps({"error": "mace-torch not installed. Run: pip install mace-torch"})
        if args.once:
            with open(args.res_file, "w") as f: f.write(err)
        else:
            print(err, flush=True)
        sys.exit(1)

    # Load model
    try:
        calc = get_calculator(args.model)
    except Exception as e:
        err = json.dumps({"error": str(e)})
        if args.once:
            with open(args.res_file, "w") as f: f.write(err)
        else:
            print(err, flush=True)
        sys.exit(1)

    if args.once:
        # Single-shot mode: read request from file, write response to file
        try:
            with open(args.req_file) as f:
                req = json.load(f)
            if req.get("quit"):
                sys.exit(0)
            N       = req["N"]
            species = req["species"]
            pos     = req["pos"]
            energy, forces = evaluate(calc, N, species, pos)
            result = {"energy": energy, "forces": forces}
            with open(args.res_file, "w") as f:
                json.dump(result, f)
        except Exception as e:
            with open(args.res_file, "w") as f:
                json.dump({"error": str(e)}, f)
        sys.exit(0)

    else:
        # Stdin/stdout streaming mode (for future persistent server use)
        for line in sys.stdin:
            line = line.strip()
            if not line: continue
            try:
                req = json.loads(line)
                if req.get("quit"): break
                N       = req["N"]
                species = req["species"]
                pos     = req["pos"]
                energy, forces = evaluate(calc, N, species, pos)
                print(json.dumps({"energy": energy, "forces": forces}), flush=True)
            except Exception as e:
                print(json.dumps({"error": str(e)}), flush=True)

if __name__ == "__main__":
    main()
