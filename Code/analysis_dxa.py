from ovito.io import import_file, export_file
from ovito.modifiers import DislocationAnalysisModifier
from ovito.data import DislocationNetwork
import ovito
ovito.enable_logging()
import sys

pipeline = import_file(f"dump.*.dump")

# Extract dislocation lines from a crystal with diamond structure:
modifier = DislocationAnalysisModifier()
modifier.input_crystal_structure = DislocationAnalysisModifier.Lattice.BCC
pipeline.modifiers.append(modifier)
data = pipeline.compute()

total_line_length = data.attributes['DislocationAnalysis.total_line_length']
cell_volume = data.attributes['DislocationAnalysis.cell_volume']
print("Dislocation density: %f" % (total_line_length / cell_volume))

# Export dislocation lines to a CA file:
export_file(pipeline, f"dislocations.*.ca", "ca", multiple_frames=True)

# Or export dislocations to a ParaView VTK file:
export_file(pipeline, f"dislocations.*.vtk", "vtk/disloc", multiple_frames=True)
