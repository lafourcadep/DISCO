# Examples

| Path | Contents |
|---|---|
| `*.yaml` | Hand-written DISCO inputs (quadrupoles, rectangular loops, cubic and orthorhombic boxes; `fcc_loops.yaml`: FCC copper; `fcc_faulted_loops.yaml`: Frank and Shockley loops in copper; `fcc_oriented_loops.yaml`: oriented box; `polycrystal_cu.yaml`: 12-grain copper polycrystal) |
| `generators/` | Input generators: single loops, quadrupoles (one or many per sample), prismatic dipoles (one, or N at random positions), BCC prismatic loops, FCC prismatic and Frank loops (`make_fcc_loops.py`), each with a JSON ground-truth sidecar |
| `replicate_quadrupole/` | Scales and tiles the reference quadrupole (`quadrupole_circle.yaml`), with random shifts |
| `in.lammps.example` | LAMMPS run template (needs `Ta1_Ravelo_2013.eam.alloy` from the NIST Interatomic Potentials Repository) |

See the documentation: [generators](https://lafourcadep.github.io/DISCO/generators/),
[quadrupoles](https://lafourcadep.github.io/DISCO/quadrupoles/),
[MD workflow](https://lafourcadep.github.io/DISCO/md_workflow/).
