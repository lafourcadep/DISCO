# Examples

| Path | Contents |
|---|---|
| `*.yaml` | Hand-written DISCO inputs (quadrupoles, rectangular loops, cubic and orthorhombic boxes) |
| `generators/` | Input generators: single loops, quadrupoles, prismatic dipoles, each with a JSON ground-truth sidecar |
| `replicate_quadrupole/` | Scales and tiles the reference quadrupole (`quadrupole_circle.yaml`), with random shifts |
| `in.lammps.example` | LAMMPS run template (needs `Ta1_Ravelo_2013.eam.alloy` from the NIST Interatomic Potentials Repository) |

See the documentation: [generators](https://lafourcadep.github.io/DISCO/generators/),
[quadrupoles](https://lafourcadep.github.io/DISCO/quadrupoles/),
[MD workflow](https://lafourcadep.github.io/DISCO/md_workflow/).
