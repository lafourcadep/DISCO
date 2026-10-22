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

class input {

private:
  std::string _name;
  int    _numberOfLoops;
  int    _nx, _ny, _nz;
  double _a0;
  double _looponSize;  // Loopon grid step, same for every loop (Angstrom)
  std::string _outputFormat;  // "xyz" | "lammps" | "both" | "none"  (default: "xyz")

  // Per-loop geometry
  std::vector< std::vector<double> > _listCenters;
  std::vector< std::vector<double> > _listM;
  std::vector< std::vector<double> > _listN;
  // All lengths in Angstrom; the centre is relative to the box centre.
  std::vector<std::string>           _listShapes;       // "circular" | "rectangular" | "parallelogram"
  std::vector<double>                _listRadius;       // circular: radius
  std::vector<double>                _listHalfM;        // rectangular: half-length along m; parallelogram: Half_1
  std::vector<double>                _listHalfP;        // rectangular: half-length along p; parallelogram: Half_2
  std::vector<double>                _listCornerRadius; // rectangular: corner radius
  std::vector< std::vector<double> > _listSide;         // rectangular: Half_m axis (empty = m);
                                                        // parallelogram: Dir_1
  std::vector< std::vector<double> > _listDir2;         // parallelogram: Dir_2

public:

  input(std::string inputFile);
  void readInputFile();

  double              getA0();
  std::vector<int>    getNumberOfCells();
  int                 getNumberOfLoops();
  int                 getNumberOfAtoms();
  std::string         getOutputFileName();
  std::string         getOutputFormat();
  double              getLooponSize();

  // Per-loop accessors
  std::vector<double> getCenter(int dislo);
  std::vector<double> getM(int dislo);
  std::vector<double> getN(int dislo);
  std::string         getShape(int dislo);        // "circular" or "rectangular"
  double              getRadius(int dislo);       // circular only
  double              getHalfM(int dislo);        // rectangular only
  double              getHalfP(int dislo);        // rectangular only
  double              getCornerRadius(int dislo); // rectangular only
  std::vector<double> getSide(int dislo);         // rectangular only; empty = along m
  std::vector<double> getDir1(int dislo);         // parallelogram only
  std::vector<double> getDir2(int dislo);         // parallelogram only
};

#endif
