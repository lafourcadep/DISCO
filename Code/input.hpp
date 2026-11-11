#ifndef INPUT_HPP
#define INPUT_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <set>
#include <yaml-cpp/yaml.h>

#include "lattice.hpp"
#include "polycrystal.hpp"

class input {

private:
  std::string _name;
  int    _numberOfLoops;
  int    _nx, _ny, _nz;
  double _a0;
  Lattice _lattice;           // host lattice: BCC (default) | FCC
  Orientation _orient;        // crystal directions of the box axes (identity if absent)
  PolycrystalSpec _poly;      // polycrystal: block (inactive if absent)
  std::string _element;       // species name written to the output files (default Ta)
  double _mass;               // atomic mass in the LAMMPS data file (Dalton)
  double _looponSize;  // Loopon grid step, same for every loop (Angstrom); 0 if unused
  std::string _method;        // "line" (closed-form line integral, default) | "loopons" (Kroupa sum)
  double _segmentTolerance;   // line: max chord sagitta of curved contours (Angstrom)
  std::string _outputFormat;  // "xyz" | "lammps" | "both" | "none"  (default: "xyz")

  // Per-loop geometry
  std::vector< std::vector<double> > _listCenters;
  std::vector< std::vector<double> > _listM;
  std::vector< std::vector<double> > _listN;
  // All lengths in Angstrom; the centre is relative to the box centre.
  std::vector<std::string>           _listShapes;       // "circular" | "rectangular" | "parallelogram" | "polygon" | "hexagonal"
  std::vector<double>                _listRadius;       // circular: radius; polygon: circumradius
  std::vector<double>                _listHalfM;        // rectangular: half-length along m; parallelogram: Half_1
  std::vector<double>                _listHalfP;        // rectangular: half-length along p; parallelogram: Half_2
  std::vector<double>                _listCornerRadius; // rectangular, polygon: corner radius
  std::vector< std::vector<double> > _listSide;         // rectangular: Half_m axis (empty = m);
                                                        // parallelogram: Dir_1; polygon: Vertex_dir
  std::vector< std::vector<double> > _listDir2;         // parallelogram: Dir_2
  std::vector<int>                   _listSides;        // polygon: number of sides (0 = Vertices)
  std::vector< std::vector< std::vector<double> > > _listVertices; // polygon: explicit vertices
  std::vector<std::string>           _listMatter;       // "auto" | "none" (loops with b.n != 0)

public:

  input(std::string inputFile);
  void readInputFile();

  double              getA0();
  Lattice             getLattice();
  const Orientation&  getOrientation() const { return _orient; }
  const PolycrystalSpec& getPolycrystal() const { return _poly; }
  std::vector<double> getBoxLengths();      // Angstrom: N_i * a0, or N_i * period along X/Y/Z
  std::string         getElement();
  double              getMass();
  std::vector<int>    getNumberOfCells();
  int                 getNumberOfLoops();
  int                 getNumberOfAtoms();   // perfect lattice, before any loop
  std::string         getOutputFileName();
  std::string         getOutputFormat();
  double              getLooponSize();
  std::string         getMethod();
  double              getSegmentTolerance();

  // Per-loop accessors
  std::vector<double> getCenter(int dislo);
  std::vector<double> getM(int dislo);
  std::vector<double> getN(int dislo);
  std::string         getShape(int dislo);        // see _listShapes
  double              getRadius(int dislo);       // circular: radius; polygon: circumradius
  double              getHalfM(int dislo);        // rectangular only
  double              getHalfP(int dislo);        // rectangular only
  double              getCornerRadius(int dislo); // rectangular, polygon
  std::vector<double> getSide(int dislo);         // rectangular: Half_m axis (empty = m);
                                                  // polygon: Vertex_dir (empty = default)
  std::vector<double> getDir1(int dislo);         // parallelogram only
  std::vector<double> getDir2(int dislo);         // parallelogram only
  std::string         getMatter(int dislo);       // "auto" or "none"
  int                 getSides(int dislo);        // polygon: regular, 0 = explicit Vertices
  std::vector< std::vector<double> > getVertices(int dislo);   // polygon, explicit
};

#endif
