#ifndef INPUT_CPP
#define INPUT_CPP

#include "input.hpp"
#include <stdexcept>

using namespace std;

input::input(string name)
  : _name(name), _numberOfLoops(0),
    _nx(10), _ny(10), _nz(10),
    _a0(3.3024), _lattice(Lattice::BCC), _element("Ta"), _mass(180.9479), _looponSize(0.), _method("line"), _segmentTolerance(0.01),
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

// A number, or an exact fraction written "p/q" (e.g. m: [2/3, 2/3, 2/3] for a
// Frank loop): 0.666667 is not accurate enough for the plane-count checks.
static double readNumber(const YAML::Node& node, const string& key, const string& where)
{
  const string txt = node.as<string>();
  const size_t slash = txt.find('/');
  if (slash == string::npos) return node.as<double>();
  try {
    size_t endP = 0, endQ = 0;
    const string ps = txt.substr(0, slash), qs = txt.substr(slash + 1);
    const double p = stod(ps, &endP), q = stod(qs, &endQ);
    if (endP != ps.size() || endQ != qs.size() || q == 0.) throw invalid_argument(txt);
    return p / q;
  } catch (const exception&) {
    fail(where, "'" + key + "': cannot read '" + txt + "' as a number or a fraction p/q");
  }
  return 0.;
}

static vector<double> readVec3(const YAML::Node& block, const string& key, const string& where)
{
  const YAML::Node v = require(block, key, where);
  if (!v.IsSequence() || v.size() != 3) fail(where, "'" + key + "' must be a list of 3 numbers");
  return { readNumber(v[0], key, where), readNumber(v[1], key, where), readNumber(v[2], key, where) };
}

// Orientation: {X: [..], Y: [..], Z: [..]}: integer crystal directions of the
// lab axes, orthogonal and right-handed
static Orientation parseOrientation(const YAML::Node& o, const string& where)
{
  Orientation out;
  if (!o.IsMap()) fail(where, "'Orientation' must be a map {X: [..], Y: [..], Z: [..]}");
  checkKeys(o, {"X", "Y", "Z"}, where + ": Orientation");
  const char* ax[3] = { "X", "Y", "Z" };
  for (int a = 0; a < 3; a++) {
    const vector<double> v = readVec3(o, ax[a], where + ": Orientation");
    double nn = 0.;
    for (int k = 0; k < 3; k++) {
      const long h = lround(v[k]);
      if (fabs(v[k] - h) > 1e-9) fail(where, string("Orientation ") + ax[a] + " must be integer Miller indices");
      out.dir[a][k] = h;
      nn += double(h) * h;
    }
    if (nn == 0.) fail(where, string("Orientation ") + ax[a] + " is zero");
    for (int k = 0; k < 3; k++) out.R[a][k] = out.dir[a][k] / sqrt(nn);
  }
  const long (*d)[3] = out.dir;
  for (int a = 0; a < 3; a++) {
    const int b = (a + 1) % 3;
    if (d[a][0]*d[b][0] + d[a][1]*d[b][1] + d[a][2]*d[b][2] != 0)
      fail(where, string("Orientation ") + ax[a] + " and " + ax[b] + " are not orthogonal");
  }
  const long cz[3] = { d[0][1]*d[1][2] - d[0][2]*d[1][1], d[0][2]*d[1][0] - d[0][0]*d[1][2],
                       d[0][0]*d[1][1] - d[0][1]*d[1][0] };
  if (cz[0]*d[2][0] + cz[1]*d[2][1] + cz[2]*d[2][2] <= 0)
    fail(where, "Orientation is left-handed: Z must be along X x Y");
  out.active = true;
  return out;
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
    checkKeys(root, {"general", "dislocations", "polycrystal", "twins"}, "top level");

    const string head = "general";
    const YAML::Node g = require(root, head, "top level");
    checkKeys(g, {"a0", "Nx", "Ny", "Nz", "Loopon_size", "Output_format", "Method",
                  "Segment_tolerance", "Lattice", "Element", "Mass", "Orientation"}, head);
    _a0         = readPositive(g, "a0", head);
    if (g["Lattice"] && !parseLattice(g["Lattice"].as<string>(), _lattice))
      fail(head, "'Lattice' must be BCC or FCC");
    if (g["Element"]) {
      _element = g["Element"].as<string>();
      if (_element.empty() || _element.find_first_of(" \t") != string::npos)
        fail(head, "'Element' must be one word (e.g. Ta, Cu)");
    }
    if (g["Mass"]) _mass = readPositive(g, "Mass", head);
    else {
      _mass = standardAtomicMass(_element);
      if (_mass <= 0.) fail(head, "no standard mass known for Element '" + _element + "': give 'Mass' (Dalton)");
    }
    // Element and Lattice are independent keys: flag a combination that is not
    // the element's room-temperature structure (allowed, e.g. for model studies)
    if (!g["Element"] && _lattice != Lattice::BCC)
      cerr << "WARNING: Lattice: " << latticeName(_lattice) << " without 'Element': the species"
           << " defaults to Ta (BCC); set Element (and Mass) for your metal\n";
    else if (const ElementData* e = findElement(_element))
      if (e->lattice != _lattice)
        cerr << "WARNING: Element " << _element << " is " << latticeName(e->lattice)
             << " at room temperature, but Lattice: " << latticeName(_lattice) << " was requested\n";
    // Orientation: crystal directions (integer) of the box axes X, Y, Z; orthogonal,
    // right-handed. Nx, Ny, Nz then count lattice periods along X, Y, Z.
    if (g["Orientation"]) _orient = parseOrientation(g["Orientation"], head);
    if (root["polycrystal"]) {
      // polycrystal: the box and the grains replace Nx, Ny, Nz and Orientation
      if (g["Nx"] || g["Ny"] || g["Nz"] || g["Orientation"])
        fail(head, "with a 'polycrystal' block, give the box there (Box), not Nx/Ny/Nz/Orientation");
      const string ph = "polycrystal";
      const YAML::Node pc = root[ph];
      checkKeys(pc, {"Box", "Grains", "Seed", "Grain_list", "GB_min_distance"}, ph);
      _poly.active = true;
      const vector<double> box = readVec3(pc, "Box", ph);
      for (int k = 0; k < 3; k++) {
        if (box[k] <= 0.) fail(ph, "'Box' lengths must be > 0");
        _poly.box[k] = box[k];
      }
      if (pc["Seed"]) _poly.seed = pc["Seed"].as<uint64_t>();
      if (pc["GB_min_distance"]) _poly.gbMinDistance = readPositive(pc, "GB_min_distance", ph);
      if (pc["Grains"] && pc["Grain_list"]) fail(ph, "give either 'Grains' (random) or 'Grain_list'");
      if (pc["Grains"]) {
        _poly.nGrains = pc["Grains"].as<int>();
        if (_poly.nGrains < 1) fail(ph, "'Grains' must be >= 1");
      } else if (pc["Grain_list"]) {
        const YAML::Node gl = pc["Grain_list"];
        if (!gl.IsSequence() || gl.size() == 0) fail(ph, "'Grain_list' must be a non-empty list");
        for (size_t i = 0; i < gl.size(); i++) {
          const string gw = ph + ": grain " + to_string(i + 1);
          const YAML::Node gi = gl[i];
          checkKeys(gi, {"Center", "Orientation", "Euler"}, gw);
          Grain gr;
          const bool hasC = (bool)gi["Center"];
          if (hasC) {
            const vector<double> c = readVec3(gi, "Center", gw);
            for (int k = 0; k < 3; k++) gr.seed[k] = c[k];
          }
          if (gi["Orientation"] && gi["Euler"]) fail(gw, "give either 'Orientation' or 'Euler'");
          const bool hasO = gi["Orientation"] || gi["Euler"];
          if (gi["Orientation"]) gr.o = parseOrientation(gi["Orientation"], gw);
          if (gi["Euler"]) {
            const vector<double> e = readVec3(gi, "Euler", gw);
            gr.o = orientationFromEuler(e[0], e[1], e[2]);
          }
          _poly.list.push_back(gr);
          _poly.hasCenter.push_back(hasC);
          _poly.hasOrientation.push_back(hasO);
        }
      } else
        fail(ph, "missing 'Grains' (number of random grains) or 'Grain_list'");
      _nx = _ny = _nz = 0;
    } else {
      _nx         = require(g, "Nx", head).as<int>();
      _ny         = require(g, "Ny", head).as<int>();
      _nz         = require(g, "Nz", head).as<int>();
    }
    if (g["Method"]) _method = g["Method"].as<string>();
    if (_method != "loopons" && _method != "line")
      fail(head, "'Method' must be loopons or line");
    if (_method == "loopons") {
      _looponSize = readPositive(g, "Loopon_size", head);
      if (g["Segment_tolerance"]) fail(head, "'Segment_tolerance' needs Method: line");
    } else {
      if (g["Loopon_size"]) _looponSize = readPositive(g, "Loopon_size", head);   // kept, unused
      if (g["Segment_tolerance"]) _segmentTolerance = readPositive(g, "Segment_tolerance", head);
    }
    if (g["Output_format"]) _outputFormat = g["Output_format"].as<string>();

    // twins: coherent twin lamellae (twin.hpp), in grains or in an oriented single crystal
    if (root["twins"]) {
      const YAML::Node tl = root["twins"];
      if (!tl.IsSequence() || tl.size() == 0) fail("twins", "must be a non-empty list");
      const int nG = _poly.list.empty() ? _poly.nGrains : (int)_poly.list.size();
      for (size_t i = 0; i < tl.size(); i++) {
        const string tw = "twin " + to_string(i + 1);
        const YAML::Node ti = tl[i];
        checkKeys(ti, {"Grain", "n", "Center", "Thickness"}, tw);
        TwinSpec sp;
        if (_poly.active) {
          const int gk = require(ti, "Grain", tw).as<int>();
          if (gk < 1 || gk > nG) fail(tw, "'Grain' must be between 1 and " + to_string(nG));
          sp.grain = gk - 1;
        } else if (ti["Grain"])
          fail(tw, "'Grain' needs a polycrystal block");
        sp.n = readVec3(ti, "n", tw);
        if (!isTwinPlane(_lattice, sp.n))
          fail(tw, string("'n' must be a twin plane: ") + (_lattice == Lattice::FCC ? "{111} in FCC" : "{112} in BCC"));
        sp.center = readVec3(ti, "Center", tw);
        sp.thickness = readPositive(ti, "Thickness", tw);
        _twins.push_back(sp);
      }
    }

    // dislocations: required, except for a polycrystal or a twinned crystal
    YAML::Node dislos;
    if (_poly.active) {
      // loops in grains (Grain: k): optional for a polycrystal
      dislos = root["dislocations"] ? root["dislocations"] : YAML::Node(YAML::NodeType::Sequence);
      if (!dislos.IsSequence()) fail("dislocations", "must be a list");
    } else if (!_twins.empty() && !root["dislocations"]) {
      dislos = YAML::Node(YAML::NodeType::Sequence);     // twins only
    } else {
      dislos = require(root, "dislocations", "top level");
      if (!dislos.IsSequence() || dislos.size() == 0)
        fail("dislocations", "must be a non-empty list");
    }
    _numberOfLoops = dislos.size();

    cout << "  Lattice        = " << latticeName(_lattice) << "  (" << _element << ", "
         << _mass << " Da)\n"
         << "  a0             = " << _a0 << " Ang\n"
         << (_poly.active ? "  Polycrystal    = " : _orient.active ? "  Periods        = " : "  Cells          = ");
    if (_poly.active)
      cout << (_poly.list.empty() ? _poly.nGrains : (int)_poly.list.size()) << " grains, box "
           << _poly.box[0] << " x " << _poly.box[1] << " x " << _poly.box[2] << " Ang, seed " << _poly.seed
           << ", GB_min_distance " << _poly.gbMinDistance << " NN\n";
    else
      cout << _nx << " x " << _ny << " x " << _nz << "\n";
    cout
         << "  Method         = " << _method
         << (_method == "line" ? "  (line integral, closed form)" : "  (Kroupa loopon sum)") << "\n";
    if (_method == "loopons")
      cout << "  Loopon_size    = " << _looponSize << " Ang\n";
    else {
      cout << "  Segment_tol.   = " << _segmentTolerance << " Ang\n";
      if (_looponSize > 0.)
        cout << "  Loopon_size    = " << _looponSize << " Ang  (unused with Method: line)\n";
    }
    cout
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
      int sides = 0;         // regular polygon; 0 = explicit Vertices
      vector< vector<double> > vertices;
      if (shape == "circular") {
        checkKeys(dl, {"Shape", "Radius", "Center", "m", "n", "Matter", "Grain", "Clip"}, where);
        radius = readPositive(dl, "Radius", where);
        cout << "    Radius         = " << radius << " Ang\n";
      } else if (shape == "rectangular") {
        checkKeys(dl, {"Shape", "Half_m", "Half_p", "Corner_radius", "Side", "Center", "m", "n", "Matter", "Grain", "Clip"}, where);
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
        checkKeys(dl, {"Shape", "Half_1", "Half_2", "Dir_1", "Dir_2", "Center", "m", "n", "Matter", "Grain", "Clip"}, where);
        half_m = readPositive(dl, "Half_1", where);
        half_p = readPositive(dl, "Half_2", where);
        side   = readVec3(dl, "Dir_1", where);
        dir2   = readVec3(dl, "Dir_2", where);
        cout << "    Half_1         = " << half_m << " Ang  along Dir_1 = "
             << side[0] << " " << side[1] << " " << side[2] << "\n"
             << "    Half_2         = " << half_p << " Ang  along Dir_2 = "
             << dir2[0] << " " << dir2[1] << " " << dir2[2] << "\n";
      } else if (shape == "polygon" || shape == "hexagonal") {
        // Regular polygon (Sides, Circumradius, optional Vertex_dir) or, for
        // "polygon" only, an explicit convex vertex list (Vertices)
        if (shape == "hexagonal")
          checkKeys(dl, {"Shape", "Circumradius", "Vertex_dir", "Corner_radius",
                         "Center", "m", "n", "Matter", "Grain", "Clip"}, where);
        else
          checkKeys(dl, {"Shape", "Sides", "Circumradius", "Vertex_dir", "Vertices", "Corner_radius",
                         "Center", "m", "n", "Matter", "Grain", "Clip"}, where);
        if (dl["Vertices"]) {
          if (dl["Sides"] || dl["Circumradius"] || dl["Vertex_dir"])
            fail(where, "'Vertices' excludes 'Sides', 'Circumradius' and 'Vertex_dir'");
          const YAML::Node vl = dl["Vertices"];
          if (!vl.IsSequence() || vl.size() < 3)
            fail(where, "'Vertices' must be a list of at least 3 points [x, y, z]");
          for (size_t k = 0; k < vl.size(); k++) {
            if (!vl[k].IsSequence() || vl[k].size() != 3)
              fail(where, "'Vertices' must be a list of at least 3 points [x, y, z]");
            vertices.push_back({ vl[k][0].as<double>(), vl[k][1].as<double>(), vl[k][2].as<double>() });
          }
          cout << "    Vertices       = " << vertices.size() << " (Ang, from the centre)\n";
        } else {
          sides = (shape == "hexagonal") ? 6 : require(dl, "Sides", where).as<int>();
          if (sides < 3) fail(where, "'Sides' must be >= 3");
          radius = readPositive(dl, "Circumradius", where);
          cout << "    Sides          = " << sides << "\n"
               << "    Circumradius   = " << radius << " Ang\n";
          if (dl["Vertex_dir"]) {   // optional: centre -> first vertex (default grid axis)
            side = readVec3(dl, "Vertex_dir", where);
            cout << "    Vertex_dir     = " << side[0] << " " << side[1] << " " << side[2] << "\n";
          }
        }
        if (dl["Corner_radius"]) {
          corner_r = dl["Corner_radius"].as<double>();
          if (corner_r < 0.) fail(where, "'Corner_radius' must be >= 0");
        }
        cout << "    Corner_radius  = " << corner_r << " Ang\n";
      } else {
        fail(where, "unknown shape '" + shape
                    + "' (circular | rectangular | parallelogram | polygon | hexagonal)");
      }

      const vector<double> center = readVec3(dl, "Center", where);
      const vector<double> m      = readVec3(dl, "m", where);
      const vector<double> n      = readVec3(dl, "n", where);
      cout << "    Center         = " << center[0] << " " << center[1] << " " << center[2] << " Ang\n"
           << "    m              = " << m[0] << " " << m[1] << " " << m[2] << "\n"
           << "    n              = " << n[0] << " " << n[1] << " " << n[2] << "\n";
      // In-plane directions: Side (rectangular), Dir_1 / Dir_2 (parallelogram)
      const double ln = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
      const string k1 = (shape == "parallelogram") ? "Dir_1"
                      : (shape == "polygon" || shape == "hexagonal") ? "Vertex_dir" : "Side";
      for (const auto& kv : { make_pair(k1, side), make_pair(string("Dir_2"), dir2) }) {
        const vector<double>& v = kv.second;
        if (v.empty()) continue;
        const double lv = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        if (lv == 0.) fail(where, "'" + kv.first + "' is zero");
        if (fabs(v[0]*n[0] + v[1]*n[1] + v[2]*n[2]) > 1e-9 * lv * ln)
          fail(where, "'" + kv.first + "' must lie in the loop plane (" + kv.first + " . n = 0)");
      }
      // Vertices: in-plane offsets from the centre (Angstrom), tolerance 1e-6 Ang
      for (const auto& v : vertices)
        if (fabs(v[0]*n[0] + v[1]*n[1] + v[2]*n[2]) > 1e-6 * ln)
          fail(where, "every point of 'Vertices' must lie in the loop plane (vertex . n = 0)");
      if (!dir2.empty()) {
        const double c[3] = { side[1]*dir2[2] - side[2]*dir2[1],
                              side[2]*dir2[0] - side[0]*dir2[2],
                              side[0]*dir2[1] - side[1]*dir2[0] };
        const double l1 = sqrt(side[0]*side[0] + side[1]*side[1] + side[2]*side[2]);
        const double l2 = sqrt(dir2[0]*dir2[0] + dir2[1]*dir2[1] + dir2[2]*dir2[2]);
        if (sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) < 1e-6 * l1 * l2)
          fail(where, "'Dir_1' and 'Dir_2' are parallel");
      }

      // Matter: auto (default) adds or removes the atoms of a loop with b.n != 0;
      // none applies the displacement field only (no atom added or removed)
      string matter = "auto";
      if (dl["Matter"]) matter = dl["Matter"].as<string>();
      if (matter != "auto" && matter != "none")
        fail(where, "'Matter' must be auto or none");
      if (fabs(m[0]*n[0] + m[1]*n[1] + m[2]*n[2]) > 1e-9)
        cout << "    Matter         = " << matter << "  (b.n != 0)\n";

      _listCenters.push_back(center);
      _listMatter.push_back(matter);

      // Grain (polycrystal only, 1-based): the loop's crystal frame and the cell it is clipped to
      const int nG = _poly.list.empty() ? _poly.nGrains : (int)_poly.list.size();
      if (_poly.active) {
        const int gk = require(dl, "Grain", where).as<int>();
        if (gk < 1 || gk > nG) fail(where, "'Grain' must be between 1 and " + to_string(nG));
        _listGrain.push_back(gk - 1);
        string clip = "auto";
        if (dl["Clip"]) clip = dl["Clip"].as<string>();
        if (clip != "auto" && clip != "none") fail(where, "'Clip' must be auto or none");
        _listClip.push_back(clip);
        cout << "    Grain          = " << gk << "  (Clip: " << clip << ")\n";
      } else {
        if (dl["Grain"] || dl["Clip"]) fail(where, "'Grain' and 'Clip' need a polycrystal block");
        _listGrain.push_back(-1);
        _listClip.push_back("none");
      }
      _listM.push_back(m);
      _listN.push_back(n);
      _listShapes.push_back(shape);
      _listRadius.push_back(radius);
      _listHalfM.push_back(half_m);
      _listHalfP.push_back(half_p);
      _listCornerRadius.push_back(corner_r);
      _listSide.push_back(side);
      _listDir2.push_back(dir2);
      _listSides.push_back(sides);
      _listVertices.push_back(vertices);
    }
  } catch (const YAML::Exception& e) {
    cerr << "ERROR in input file " << _name << ": " << e.what() << "\n";
    abort();
  }
}

double input::getA0()         { return _a0; }
Lattice input::getLattice()   { return _lattice; }
string input::getElement()    { return _element; }
double input::getMass()       { return _mass; }
string input::getOutputFormat() { return _outputFormat; }
double input::getLooponSize() { return _looponSize; }
string input::getMethod()     { return _method; }
double input::getSegmentTolerance() { return _segmentTolerance; }

vector<int> input::getNumberOfCells() { return { _nx, _ny, _nz }; }
int         input::getNumberOfLoops() { return _numberOfLoops; }
int         input::getNumberOfAtoms()
{
  if (!_orient.active) return atomsPerCell(_lattice) * _nx * _ny * _nz;
  const vector<double> L = getBoxLengths();
  return (int)lround(atomsPerCell(_lattice) * L[0] * L[1] * L[2] / (_a0 * _a0 * _a0));
}

vector<double> input::getBoxLengths()
{
  if (_poly.active) return { _poly.box[0], _poly.box[1], _poly.box[2] };
  const int n[3] = { _nx, _ny, _nz };
  vector<double> L(3);
  for (int a = 0; a < 3; a++)
    L[a] = n[a] * (_orient.active ? latticePeriod(_lattice, _orient.dir[a], _a0) : _a0);
  return L;
}

string input::getOutputFileName()
{
  if (_poly.active) {
    const int G = _poly.list.empty() ? _poly.nGrains : (int)_poly.list.size();
    return "Poly_" + to_string(G) + "g_" + to_string(lround(_poly.box[0])) + "x"
           + to_string(lround(_poly.box[1])) + "x" + to_string(lround(_poly.box[2])) + ".xyz";
  }
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
string         input::getMatter(int d)       { return _listMatter[d]; }
int            input::getGrain(int d)        { return _listGrain[d]; }
string         input::getClip(int d)         { return _listClip[d]; }
int            input::getSides(int d)        { return _listSides[d]; }
vector< vector<double> > input::getVertices(int d) { return _listVertices[d]; }

#endif
