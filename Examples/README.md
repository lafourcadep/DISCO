# Examples

| Path | Contents |
|---|---|
| `*.yaml` | Hand-written DISCO inputs (quadrupoles, rectangular loops, cubic and orthorhombic boxes; `fcc_loops.yaml`: FCC copper; `fcc_faulted_loops.yaml`: Frank and Shockley loops in copper; `fcc_oriented_loops.yaml`: oriented box; `polycrystal_cu.yaml`: 12-grain copper polycrystal; `polycrystal_cu_loops.yaml`: loops in grains, clipped at the GBs; `fcc_twin_lamella.yaml`, `polycrystal_cu_twins.yaml`: twin lamellae with twin-aware loops) |
| `generators/` | Input generators, each with a JSON ground-truth sidecar: single loops, quadrupoles (one or many per sample), prismatic dipoles (`make_prismatic_dipoles.py`, one or N), BCC prismatic loops, FCC prismatic and Frank loops (`make_fcc_loops.py`), dislocation-populated polycrystals (`make_polycrystal_loops.py`) |
| `in.lammps.example` | LAMMPS run template (needs `Ta1_Ravelo_2013.eam.alloy` from the NIST Interatomic Potentials Repository) |

See the documentation: [generators](https://lafourcadep.github.io/DISCO/generators/),
[quadrupoles](https://lafourcadep.github.io/DISCO/quadrupoles/),
[MD workflow](https://lafourcadep.github.io/DISCO/md_workflow/).
