#ifndef INPUT_CPP
#define INPUT_CPP

#include "input.hpp"

using namespace std;

input::input(string name)
  : _name(name), _numberOfLoops(0),
    _nx(10), _ny(10), _nz(10),
    _a0(3.3024), _looponSize(0.),
    _outputFormat("xyz")
{}

static void fail(const string& where, const string& msg)
{
  cerr << "ERROR in input file (" << where << "): " << msg << "\n";
  abort();
}

// Unknown keys are errors, so a typo cannot silently fall back to a default.
static void checkKeys(const YAML::Node& block, const set<string>& allowed, const string& where)
{
  for (const auto& kv : block) {
    const string k = kv.first.as<string>();
    if (!allowed.count(k)) fail(where, "unknown key '" + k + "'");
  }
}

static YAML::Node require(const YAML::Node& block, const string& key, const string& where)
{
  const YAML::Node node = block[key];
  if (!node) fail(where, "missing key '" + key + "'");
  return node;
}

static double readPositive(const YAML::Node& block, const string& key, const string& where)
{
  const double v = require(block, key, where).as<double>();
  if (v <= 0.) fail(where, "'" + key + "' must be > 0");
  return v;
}

static vector<double> readVec3(const YAML::Node& block, const string& key, const string& where)
{
  const YAML::Node v = require(block, key, where);
  if (!v.IsSequence() || v.size() != 3) fail(where, "'" + key + "' must be a list of 3 numbers");
  return { v[0].as<double>(), v[1].as<double>(), v[2].as<double>() };
}

void input::readInputFile()
{
  ifstream stream(_name.c_str());
  if (!stream.good()) {
    cerr << "Unable to open file " << _name << "\n";
    abort();
  }
  string first;
  stream >> first;
  if (first == "a0=") {
    cerr << "ERROR: old 'key= value' input format.\n"
         << "  Inputs are now YAML (a 'general:' block and a 'dislocations:' list).\n";
    abort();
  }
  stream.close();
  cout << "Reading parameters from " << _name << "\n";

  try {
    const YAML::Node root = YAML::LoadFile(_name);
    checkKeys(root, {"general", "dislocations"}, "top level");

    const string head = "general";
    const YAML::Node g = require(root, head, "top level");
    checkKeys(g, {"a0", "Nx", "Ny", "Nz", "Loopon_size", "Output_format"}, head);
    _a0         = readPositive(g, "a0", head);
    _nx         = require(g, "Nx", head).as<int>();
    _ny         = require(g, "Ny", head).as<int>();
    _nz         = require(g, "Nz", head).as<int>();
    _looponSize = readPositive(g, "Loopon_size", head);
    if (g["Output_format"]) _outputFormat = g["Output_format"].as<string>();

    const YAML::Node dislos = require(root, "dislocations", "top level");
    if (!dislos.IsSequence() || dislos.size() == 0)
      fail("dislocations", "must be a non-empty list");
    _numberOfLoops = dislos.size();

    cout << "  a0             = " << _a0 << " Ang\n"
         << "  Cells          = " << _nx << " x " << _ny << " x " << _nz << "\n"
         << "  Loopon_size    = " << _looponSize << " Ang\n"
         << "  Loops          = " << _numberOfLoops << "\n";
    if (g["Output_format"])
      cout << "  Output_format  = " << _outputFormat << "\n";

    for (int d = 0; d < _numberOfLoops; d++) {
      const string where = "dislocation " + to_string(d + 1);
      const YAML::Node dl = dislos[d];
      cout << "\n  Dislocation " << d+1 << "\n";

      const string shape = require(dl, "Shape", where).as<string>();
      cout << "    Shape          = " << shape << "\n";

      double radius = 0., half_m = 0., half_p = 0., corner_r = 0.;
      vector<double> side;   // empty = rectangle sides along m and m x n; parallelogram Dir_1
      vector<double> dir2;   // parallelogram Dir_2
      if (shape == "circular") {
        checkKeys(dl, {"Shape", "Radius", "Center", "m", "n"}, where);
        radius = readPositive(dl, "Radius", where);
        cout << "    Radius         = " << radius << " Ang\n";
      } else if (shape == "rectangular") {
        checkKeys(dl, {"Shape", "Half_m", "Half_p", "Corner_radius", "Side", "Center", "m", "n"}, where);
        half_m   = readPositive(dl, "Half_m", where);
        half_p   = readPositive(dl, "Half_p", where);
        corner_r = require(dl, "Corner_radius", where).as<double>();  // 0 = sharp corners
        cout << "    Half_m         = " << half_m   << " Ang\n"
             << "    Half_p         = " << half_p   << " Ang\n"
             << "    Corner_radius  = " << corner_r << " Ang\n";
        if (dl["Side"]) {   // optional: direction of the Half_m axis (default m)
          side = readVec3(dl, "Side", where);
          cout << "    Side           = " << side[0] << " " << side[1] << " " << side[2] << "\n";
        }
      } else if (shape == "parallelogram") {
        checkKeys(dl, {"Shape", "Half_1", "Half_2", "Dir_1", "Dir_2", "Center", "m", "n"}, where);
        half_m = readPositive(dl, "Half_1", where);
        half_p = readPositive(dl, "Half_2", where);
        side   = readVec3(dl, "Dir_1", where);
        dir2   = readVec3(dl, "Dir_2", where);
        cout << "    Half_1         = " << half_m << " Ang  along Dir_1 = "
             << side[0] << " " << side[1] << " " << side[2] << "\n"
             << "    Half_2         = " << half_p << " Ang  along Dir_2 = "
             << dir2[0] << " " << dir2[1] << " " << dir2[2] << "\n";
      } else {
        fail(where, "unknown shape '" + shape + "' (circular | rectangular | parallelogram)");
      }

      const vector<double> center = readVec3(dl, "Center", where);
      const vector<double> m      = readVec3(dl, "m", where);
      const vector<double> n      = readVec3(dl, "n", where);
      cout << "    Center         = " << center[0] << " " << center[1] << " " << center[2] << " Ang\n"
           << "    m              = " << m[0] << " " << m[1] << " " << m[2] << "\n"
           << "    n              = " << n[0] << " " << n[1] << " " << n[2] << "\n";
      // In-plane directions: Side (rectangular), Dir_1 / Dir_2 (parallelogram)
      const double ln = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
      const string k1 = (shape == "parallelogram") ? "Dir_1" : "Side";
      for (const auto& kv : { make_pair(k1, side), make_pair(string("Dir_2"), dir2) }) {
        const vector<double>& v = kv.second;
        if (v.empty()) continue;
        const double lv = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        if (lv == 0.) fail(where, "'" + kv.first + "' is zero");
        if (fabs(v[0]*n[0] + v[1]*n[1] + v[2]*n[2]) > 1e-9 * lv * ln)
          fail(where, "'" + kv.first + "' must lie in the loop plane (" + kv.first + " . n = 0)");
      }
      if (!dir2.empty()) {
        const double c[3] = { side[1]*dir2[2] - side[2]*dir2[1],
                              side[2]*dir2[0] - side[0]*dir2[2],
                              side[0]*dir2[1] - side[1]*dir2[0] };
        const double l1 = sqrt(side[0]*side[0] + side[1]*side[1] + side[2]*side[2]);
        const double l2 = sqrt(dir2[0]*dir2[0] + dir2[1]*dir2[1] + dir2[2]*dir2[2]);
        if (sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) < 1e-6 * l1 * l2)
          fail(where, "'Dir_1' and 'Dir_2' are parallel");
      }

      _listCenters.push_back(center);
      _listM.push_back(m);
      _listN.push_back(n);
      _listShapes.push_back(shape);
      _listRadius.push_back(radius);
      _listHalfM.push_back(half_m);
      _listHalfP.push_back(half_p);
      _listCornerRadius.push_back(corner_r);
      _listSide.push_back(side);
      _listDir2.push_back(dir2);
    }
  } catch (const YAML::Exception& e) {
    cerr << "ERROR in input file " << _name << ": " << e.what() << "\n";
    abort();
  }
}

double input::getA0()         { return _a0; }
string input::getOutputFormat() { return _outputFormat; }
double input::getLooponSize() { return _looponSize; }

vector<int> input::getNumberOfCells() { return { _nx, _ny, _nz }; }
int         input::getNumberOfLoops() { return _numberOfLoops; }
int         input::getNumberOfAtoms() { return 2 * _nx * _ny * _nz; }

string input::getOutputFileName()
{
  return "Geo_" + to_string(_nx) + "x" + to_string(_ny) + "x" + to_string(_nz)
         + "_with_loop.xyz";
}

vector<double> input::getCenter(int d)       { return _listCenters[d]; }
vector<double> input::getM(int d)            { return _listM[d]; }
vector<double> input::getN(int d)            { return _listN[d]; }
string         input::getShape(int d)        { return _listShapes[d]; }
double         input::getRadius(int d)       { return _listRadius[d]; }
double         input::getHalfM(int d)        { return _listHalfM[d]; }
double         input::getHalfP(int d)        { return _listHalfP[d]; }
double         input::getCornerRadius(int d) { return _listCornerRadius[d]; }
vector<double> input::getSide(int d)         { return _listSide[d]; }
vector<double> input::getDir1(int d)         { return _listSide[d]; }
vector<double> input::getDir2(int d)         { return _listDir2[d]; }

#endif
