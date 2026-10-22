# MD workflow

DISCO writes an unrelaxed sample: the atoms carry the elastic displacement of the loops. A typical
next step is to relax it, and possibly load it, with molecular dynamics, then extract the
dislocation lines.

## LAMMPS

`Examples/in.lammps.example` is a template run for BCC tantalum. It:

1. reads `Geo_40x40x40_with_loop.lmp` (`Output_format: lammps`);
2. sets the Ta EAM potential `Ta1_Ravelo_2013.eam.alloy`;
3. gives velocities at 500 K;
4. compresses x by 5 % in one step at constant volume;
5. runs 50 000 NVT steps at 500 K, with an atom dump every 100 steps in `dir.dumps/`.

Edit the `read_data` path, the deformation and the thermostat to suit your study.

!!! info "Potential file"
    The EAM potential is not distributed with DISCO. Download `Ta1_Ravelo_2013.eam.alloy`
    (R. Ravelo et al., *Phys. Rev. B* **88**, 134101 (2013)) from the
    [NIST Interatomic Potentials Repository](https://www.ctcms.nist.gov/potentials/) and place
    it next to the LAMMPS input.

## Dislocation extraction with OVITO

`Code/analysis_dxa.py` runs OVITO's Dislocation Analysis (DXA, BCC) on every `dump.*.dump` file in
the current directory. It prints the dislocation density (total line length / cell volume), then
writes the lines in two formats for every frame:

- `dislocations.*.ca`: the OVITO crystal-analysis format;
- `dislocations.*.vtk`: the ParaView format.

```bash
cd dir.dumps
ovitos ../analysis_dxa.py      # or: python analysis_dxa.py with the `ovito` Python module
```

The `.ca` files can then be analysed further, for example with DANCE, to classify segments by
Burgers vector and slip system and to detect junctions.
