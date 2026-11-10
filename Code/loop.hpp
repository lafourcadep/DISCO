#ifndef LOOP_HPP
#define LOOP_HPP

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>

// Shape of the dislocation loop surface
enum class LoopShape { CIRCULAR, RECTANGULAR, PARALLELOGRAM, POLYGON };

class loop {

private:

  // All lengths in Angstrom
  LoopShape _shape;
  double _radius;      // Circular loop radius
  double _half_m;      // Rectangular half-width along u (m by default); parallelogram Half_1
  double _half_p;      // Rectangular half-width along v (p by default); parallelogram Half_2
  double _r_corner;    // Rounded-corner radius for rectangular / polygon (0 = sharp)
  // Polygon: vertices in the orthonormal in-plane basis (_u, _v), counter-clockwise
  // about n, and the inner polygon offset inward by _r_corner (the rounded contour is
  // the set of points within _r_corner of it). _half_m / _half_p = bounding box.
  std::vector<double> _pu, _pv;
  std::vector<double> _qu, _qv;
  std::vector<double> _nu, _nv;   // inward unit normal of side i (P_i -> P_i+1)

  std::vector<double> _center;   // Loop centre (Angstroms, already on a mid-plane)
  std::vector<double> _m;        // Normalised glide direction (Burgers vector direction)
  std::vector<double> _n;        // Normalised normal to the glide plane
  std::vector<double> _p;        // p = m × n  (in-plane, perpendicular to m)
  std::vector<double> _u;        // Loopon grid axis of Half_m (default m; rectangular Side;
                                 // parallelogram Dir_1)
  std::vector<double> _v;        // Loopon grid axis of Half_p: v = u × n (default p);
                                 // parallelogram Dir_2 (oblique grid)
  std::vector<double> _burgers;  // Burgers vector magnitude per Cartesian direction
  std::vector<double> _compress; // Compression coefficients
  double _a0;                    // Lattice parameter (Angstroms)
  double _da;                    // Uniform grid spacing between loopons (Angstroms)
  std::vector< std::vector<double> > _loopons; // Loopon positions (Angstroms)
  double _dS;                    // Loopon surface element (Ang²)

public:

  // All lengths in Angstrom. `center` should come from snapToMidPlane().
  // da : loopon grid step (same for every loop)

  // Circular loop constructor
  loop(double da, double radius,
       std::vector<double> center, std::vector<double> m, std::vector<double> n,
       double latticeParameter, std::vector<double> compress);

  // Rectangular loop constructor
  // half_m   : half-length along m (screw sides)
  // half_p   : half-length along p (edge  sides)
  // corner_r : corner rounding radius (0 = sharp right angle),
  //            clamped internally to min(half_m, half_p)
  // side     : in-plane direction of the half_m axis (empty = m, i.e. screw
  //            sides along m); the half_p axis is side × n. b is always m.
  loop(double da, double half_m, double half_p, double corner_r,
       std::vector<double> center, std::vector<double> m, std::vector<double> n,
       double latticeParameter, std::vector<double> compress,
       std::vector<double> side = {});

  // Parallelogram constructor: sides along the in-plane directions dir1, dir2
  // (unit vectors taken, need not be orthogonal), half-lengths half1, half2.
  // The loopon grid is oblique along dir1, dir2, so each cell is a small
  // parallelogram of area da^2 sin(angle) and the tiles follow both sides.
  loop(double da, double half1, double half2,
       std::vector<double> dir1, std::vector<double> dir2,
       std::vector<double> center, std::vector<double> m, std::vector<double> n,
       double latticeParameter, std::vector<double> compress);

  // Convex polygon constructor
  // vertices : in-plane offsets from the centre (Angstrom, 3 components each, any
  //            order of traversal, strictly convex; the component along n is dropped)
  // corner_r : corner rounding radius (0 = sharp), clamped to the largest radius
  //            that keeps every side
  // axis     : in-plane direction of the loopon grid axis u (empty = default)
  loop(double da, const std::vector< std::vector<double> >& vertices, double corner_r,
       std::vector<double> center, std::vector<double> m, std::vector<double> n,
       double latticeParameter, std::vector<double> compress,
       std::vector<double> axis = {});

  void generateLoopons();
  void generateSurface();

  // In-plane point (u, v) along the grid axes _u, _v (oblique for a parallelogram)
  // inside the loop contour, with the boundary rules of generateLoopons().
  bool insideUV(double u, double v) const;
  // rel = point - centre (Angstrom, already minimum-image): its projection on the
  // loop plane lies inside the contour (the component along n is ignored).
  bool contains(const double rel[3]) const;

  std::vector<double> getCenter()           const;
  std::vector<double> getM()                const;
  std::vector<double> getN()                const;
  std::vector<double> getP()                const;
  std::vector<double> getOneLoopon(int nb)  const;
  std::vector<double> getBurgers()          const;
  int                 getNumberOfLoopons()  const;
  double              getSurface()          const;
  LoopShape           getShape()            const;
  double              getCornerRadius()     const;
  int                 getNumberOfVertices() const;   // polygon only
  double              getContourArea()      const;   // exact area of the (rounded) contour

  // Contour as a closed polyline for the line-integral field (barnett.hpp):
  // absolute positions (Angstrom), counter-clockwise about n, last vertex joined
  // to the first. Straight sides are exact; circles and rounded corners are cut
  // into chords of sagitta <= eps (a circle's chord radius keeps its area).
  std::vector< std::vector<double> > polyline(double eps) const;
};

// Vertices (in-plane offsets from the centre, Angstrom) of a regular polygon with
// `sides` sides and circumradius R in the plane normal to n. The first vertex lies
// along vertexDir (in-plane; empty = the default grid axis u of the other shapes).
std::vector< std::vector<double> > regularPolygon(int sides, double R,
                                                  std::vector<double> vertexDir,
                                                  const std::vector<double>& m,
                                                  const std::vector<double>& n);

// Move `center` (Angstrom, box-centred frame) along n onto the nearest plane
// halfway between two atomic planes of the (hkl) = n family of the BCC
// supercell built by generateBCC(). In-plane components are kept.
// n must be integer Miller indices; `spacing` returns the interplanar distance.
std::vector<double> snapToMidPlane(const std::vector<double>& center,
                                   const std::vector<double>& n,
                                   double a0, int nx, int ny, int nz,
                                   double& spacing);

#endif
