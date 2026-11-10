# Theory

## Loops as sums of loopons

A dislocation loop of Burgers vector **b** bounding a planar surface S produces a displacement field
that is the sum of the fields of infinitesimal loops tiling S. DISCO tiles each loop surface with a
regular grid of step `Loopon_size` centred on the loop centre. The grid is square for circles,
rectangles and polygons, and oblique for parallelograms. Every grid node inside the loop contour is a
*loopon*: an infinitesimal loop that carries the area dS of one grid cell. For a circle, nodes
exactly on the contour are excluded; for a rectangle or a polygon, they are included.

![Loop tiled by loopons](assets/kroupa_loopons.png){ width="320" align="right" }

For each loop:

- m is the unit Burgers direction (b = |b| m, |b| = a0 |m_Miller| / 2);
- n is the unit loop-plane normal;
- dS = `Loopon_size`² (square grid) or `Loopon_size`² sin(Dir_1, Dir_2) (oblique grid);
- the swept area is A = N_loopons · dS. It tends to the exact contour area as `Loopon_size` → 0;
  for polygons the log prints both.

The contour only decides which grid nodes are loopons. Every shape then goes through the same
Kroupa sum, so a new shape needs no change to the kernels.

## Kroupa displacement field

The displacement at an atom due to one loopon, with ρ = r_atom − r_loopon and ρ = |ρ|, is the
isotropic-elasticity field of an infinitesimal loop (Kroupa, 1966):

$$
u_i = -\frac{|b|\, dS}{8\pi(1-\nu)}\,\frac{1}{\rho^2}
\left[\frac{1-2\nu}{\rho}\left(n_i\,(\mathbf m\cdot\boldsymbol\rho) + m_i\,(\mathbf n\cdot\boldsymbol\rho) - \rho_i\,(\mathbf m\cdot\mathbf n)\right)
+ \frac{3\,\rho_i\,(\mathbf m\cdot\boldsymbol\rho)(\mathbf n\cdot\boldsymbol\rho)}{\rho^3}\right]
$$

with Poisson's ratio ν = 0.35 (set in `Code/kroupa.hpp`). An atom closer than 1e-10 Å to a loopon
gets no contribution from it.

The displacement of an atom is the sum over all loopons of all loops, and over the 27 periodic
images (the box and its 26 neighbours) of each loopon. Atoms are then wrapped back into the box.
The cost is therefore proportional to atoms × total loopons × 27 (see [CPU or GPU](backends.md#speed)).

!!! note "Loopon size"
    The loopon step sets how finely the loop contour is resolved. The default used by the
    generators, 0.4 a0, is well below the interatomic distance.

## Prismatic loops: adding or removing matter

For a loop with b · n ≠ 0, the Volterra cut opens a gap or makes its lips overlap by |b · n̂|.
The displacement field cannot create or destroy matter. Without atoms added or removed, the
configuration is only a strained perfect crystal, and it relaxes back.
`Code/prismatic.cpp` follows Atomsk's procedure (`opt_disloc.f90`, loop case):
it takes the atoms in a slab of thickness |b · n̂| inside the contour, then deletes them
(vacancy, m · n > 0) or duplicates them (interstitial, m · n < 0). Each copy gets the displacement
u(x) + sign(s) **b**, so it lands on the other side of the cut. Since **b** is a lattice
translation, originals and copies fill the gap with perfect crystal.
See [Input file](input.md#prismatic-and-mixed-loops-b-n-0).

The matter step wraps the Kroupa sum and does not replace it:

1. **Slab selection** (host). The slab atoms are those between the two mid-planes that bound
   |b · n̂| and inside the contour, with the same contour test as the loopon grid. All slabs are
   selected on the perfect lattice.
2. **Vacancy loops**: the slab atoms are deleted.
3. **Kroupa sum** (GPU or CPU): one call for all loops, glide and prismatic alike. Each loopon
   carries the full **b**, including its component along n (the m · n terms of the formula above).
4. **Interstitial loops**: each copy is appended at the site of its source atom, with the
   displacement already computed for that atom plus sign(s) **b**. No kernel is run for the copies.

## Loop placement

The cut surface of each loop must lie between two atomic planes, not on one. The requested centre is
moved along n onto the nearest mid-plane between dense planes of the (hkl) = n family, with spacing
d = a0/|n| (h + k + l even) or a0/(2|n|) (odd). See [Input file](input.md#loop-centre-and-cut-plane).

## Box deformation

Inserting a shear loop of area A displaces the material above its plane by b relative to the material below.
In a periodic box this corresponds to a homogeneous plastic distortion. `Code/loop_deform.cpp`
computes the matching box deformation.

!!! warning "Currently not applied"
    The routine is implemented and unit-tested (`make test`), but the driver currently does not apply
    it: in `Code/main.cpp` the call to `computeU` is commented out and U is set to the identity.
    The box is therefore the undeformed orthorhombic box, and the log prints `U − I = 0`.

### Plastic distortion

For loops k with Burgers vector b⁽ᵏ⁾, unit normal n⁽ᵏ⁾ and swept area A⁽ᵏ⁾ in a box of volume V:

$$
\beta_{ij} = \frac{1}{V}\sum_k b_i^{(k)}\, n_j^{(k)}\, A^{(k)}
$$

β is not symmetric in general. Its symmetric part is the plastic strain, and its antisymmetric part
the plastic lattice rotation.

### Pure stretch

The deformation gradient is F = I + β. Its polar decomposition F = R U is computed from the
eigendecomposition of C = Fᵀ F = V diag(λ) Vᵀ, which gives U = V diag(√λ) Vᵀ. Only the symmetric,
positive-definite stretch U is applied, which removes the rigid rotation R:

- box matrix (columns = box vectors): H′ = U H;
- atom positions: r′ = U r.

For pure shear loops (b · n = 0) the volume is conserved, det U = 1. The code warns if the volume
changes by more than 0.01 %.

## References

- F. Kroupa, "Dislocation dipoles and dislocation loops", *J. Phys. Colloques* **27** (1966) C3-154.
- J. P. Hirth and J. Lothe, *Theory of Dislocations*, 2nd ed., Wiley (1982).
- E. B. Tadmor and R. E. Miller, *Modeling Materials*, Cambridge University Press (2011).
