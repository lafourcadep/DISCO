# ==============================================================================
# Build the Kroupa dislocation-loop code
#
#   make            -> exec_cuda  GPU build (needs nvcc); `exec_cuda --cpu` or a
#                                 machine without GPU uses the OpenMP backend
#   make exec_cpu   -> exec_cpu   CPU-only build (g++ + OpenMP, no CUDA toolkit)
#   make all        -> both
#   make exec_mvcc  -> exec_mvcc  Apple-silicon GPU build through mvcc (CUDA on Metal,
#                                 https://github.com/doximity/mvcc): put mvcc's toolkit/bin
#                                 first in PATH, or make exec_mvcc MVCC_NVCC=/path/to/toolkit/bin/nvcc.
#                                 Apple GPUs emulate fp64: run it with --fp32.
#
# Adjust GPU_ARCH to match your hardware (make GPU_ARCH=sm_80):
#   sm_61  — GTX 1080 / 1070 (Pascal)
#   sm_70  — V100 (Volta)
#   sm_75  — RTX 2080 / Quadro RTX (Turing)
#   sm_80  — A100 (Ampere)
#   sm_86  — RTX 3090 / 3080 (Ampere)
#   sm_89  — RTX 4090 / 4080 (Ada Lovelace)
#   sm_90  — H100 (Hopper)
#   sm_120 — RTX PRO Blackwell
#
# Run  nvcc --list-gpu-arch  to see all options on your installation.
# ==============================================================================

NVCC      = nvcc
CXX       = g++
GPU_ARCH  = sm_120   # <-- change this to match your GPU

PARENT    = .
INCLUDES  = -I$(PARENT)

# yaml-cpp (input parser), first match wins:
#   1. make YAML_PREFIX=/path/to/yaml-cpp      (include/ and lib64/, lib/ or the prefix itself)
#      or make YAML_INCDIR=... YAML_LIBDIR=...
#   2. an environment module that sets YAMLCPP_INCDIR / YAMLCPP_LIBDIR
#      (e.g. `module load yaml-cpp` on the TGCC machines)
#   3. the spack install pinned below (so the build doesn't depend on whatever
#      version is `spack load`-ed: 0.6.3 vs 0.8.0)
#   4. pkg-config yaml-cpp
YAML_VERSION := 0.8.0
ifneq ($(strip $(YAML_PREFIX)),)
  YAML_INCDIR ?= $(YAML_PREFIX)/include
  YAML_LIBDIR ?= $(patsubst %/,%,$(dir $(firstword $(wildcard \
                   $(YAML_PREFIX)/lib64/libyaml-cpp.* $(YAML_PREFIX)/lib/libyaml-cpp.* $(YAML_PREFIX)/libyaml-cpp.*))))
else ifneq ($(strip $(YAMLCPP_INCDIR)),)
  YAML_INCDIR ?= $(YAMLCPP_INCDIR)
  YAML_LIBDIR ?= $(YAMLCPP_LIBDIR)
else
  SPACK_YAML := $(shell spack location -i yaml-cpp@$(YAML_VERSION) 2>/dev/null)
  ifneq ($(strip $(SPACK_YAML)),)
    YAML_INCDIR ?= $(SPACK_YAML)/include
    YAML_LIBDIR ?= $(firstword $(wildcard $(SPACK_YAML)/lib64 $(SPACK_YAML)/lib))
  endif
endif
ifneq ($(strip $(YAML_INCDIR)),)
  YAML_CFLAGS := -I$(YAML_INCDIR)
  YAML_LIBS   := -L$(YAML_LIBDIR) -Wl,-rpath,$(YAML_LIBDIR) -lyaml-cpp
  YAML_NVLIBS := -L$(YAML_LIBDIR) -Xlinker -rpath -Xlinker $(YAML_LIBDIR) -lyaml-cpp
else
  YAML_CFLAGS := $(shell pkg-config --cflags yaml-cpp)
  YAML_LIBS   := $(shell pkg-config --libs yaml-cpp)
  YAML_NVLIBS := $(YAML_LIBS)
endif

CXX_FLAGS  = -O3 -march=native -std=c++14 -fopenmp
NVCC_FLAGS = -O3 -arch=$(GPU_ARCH) -std=c++14 -DWITH_CUDA -Xcompiler -fopenmp,-O3,-march=native

SRC_HOST  = main.cpp input.cpp loop.cpp loop_deform.cpp prismatic.cpp backend_cpu.cpp
HEADERS   = input.hpp lattice.hpp loop.hpp loop_deform.hpp prismatic.hpp backend.hpp kroupa.hpp barnett.hpp

# ------------------------------------------------------------------------------
exec_cuda: $(SRC_HOST) backend_cuda.cu $(HEADERS)
	$(NVCC) $(NVCC_FLAGS) $(INCLUDES) $(YAML_CFLAGS) $(SRC_HOST) backend_cuda.cu -o $@ $(YAML_NVLIBS)

exec_cpu: $(SRC_HOST) $(HEADERS)
	$(CXX) $(CXX_FLAGS) $(INCLUDES) $(YAML_CFLAGS) $(SRC_HOST) -o $@ $(YAML_LIBS)

# mvcc (macOS, Apple silicon): host code is compiled by Apple clang, so OpenMP
# comes from Homebrew libomp (-Xpreprocessor -fopenmp)
MVCC_NVCC  ?= nvcc
OMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
MVCC_FLAGS  = -O3 -std=c++14 -DWITH_CUDA -Xcompiler -Xpreprocessor,-fopenmp,-I$(OMP_PREFIX)/include,-O3

exec_mvcc: $(SRC_HOST) backend_cuda.cu $(HEADERS)
	$(MVCC_NVCC) $(MVCC_FLAGS) $(INCLUDES) $(YAML_CFLAGS) $(SRC_HOST) backend_cuda.cu -o $@ \
	  $(YAML_NVLIBS) -L$(OMP_PREFIX)/lib -lomp

all: exec_cuda exec_cpu

clean:
	rm -f *.o *~ exec_cuda exec_cpu exec_mvcc

run: exec_cuda
	./exec_cuda ../Examples/quadrupole_circle.yaml

TEST_BIN  = tests/test_loop_deform
TEST_SRC  = tests/test_loop_deform.cpp

$(TEST_BIN): $(TEST_SRC) $(PARENT)/loop_deform.cpp $(PARENT)/loop_deform.hpp
	$(CXX) -O2 -std=c++14 -I$(PARENT) -o $(TEST_BIN) $(TEST_SRC) $(PARENT)/loop_deform.cpp

SNAP_BIN  = tests/test_snap

$(SNAP_BIN): tests/test_snap.cpp $(PARENT)/loop.cpp $(PARENT)/loop.hpp $(PARENT)/lattice.hpp
	$(CXX) -O2 -std=c++14 -I$(PARENT) -o $(SNAP_BIN) tests/test_snap.cpp $(PARENT)/loop.cpp

PRISM_BIN = tests/test_prismatic

$(PRISM_BIN): tests/test_prismatic.cpp $(PARENT)/loop.cpp $(PARENT)/loop.hpp $(PARENT)/prismatic.cpp $(PARENT)/prismatic.hpp $(PARENT)/kroupa.hpp $(PARENT)/lattice.hpp
	$(CXX) -O2 -std=c++14 -I$(PARENT) -o $(PRISM_BIN) tests/test_prismatic.cpp $(PARENT)/loop.cpp $(PARENT)/prismatic.cpp

BARN_BIN  = tests/test_barnett

$(BARN_BIN): tests/test_barnett.cpp $(PARENT)/loop.cpp $(PARENT)/loop.hpp $(PARENT)/barnett.hpp $(PARENT)/kroupa.hpp
	$(CXX) -O2 -std=c++14 -I$(PARENT) -o $(BARN_BIN) tests/test_barnett.cpp $(PARENT)/loop.cpp

test: $(TEST_BIN) $(SNAP_BIN) $(PRISM_BIN) $(BARN_BIN)
	./$(TEST_BIN)
	./$(SNAP_BIN)
	./$(PRISM_BIN)
	./$(BARN_BIN)

.PHONY: all clean run test
