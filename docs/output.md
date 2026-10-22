# Output files

| `Output_format` | Files written |
|---|---|
| `xyz` (default) | `Geo_<Nx>x<Ny>x<Nz>_with_loop.xyz` |
| `lammps` | `Geo_<Nx>x<Ny>x<Nz>_with_loop.lmp` |
| `both` | both of the above |
| `none` | nothing (timing runs) |

Files are written in the current directory. The sample has 2 · Nx · Ny · Nz atoms.

## Extended XYZ

Extended XYZ with a `Lattice="..."` header holding the 9 components of the box matrix, listed
column by column (one box vector after the other). OVITO reads it directly.

## LAMMPS data file

`atomic` atom style, readable with `read_data`. If the box is triclinic the file has an
`xy xz yz` tilt line. Atom positions are expressed in the LAMMPS right-handed frame
obtained by a QR (Gram–Schmidt) decomposition of the box vectors.

Atoms are wrapped back into the periodic box after the displacement is added.
