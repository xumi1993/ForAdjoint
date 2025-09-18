# ForAdjoint

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language](https://img.shields.io/badge/Language-Fortran-blue.svg)](https://fortran-lang.org)
[![CMake](https://img.shields.io/badge/Build-CMake-green.svg)](https://cmake.org)

A high-performance Fortran library for calculating adjoint sources in seismic waveform inversion, providing efficient implementations of various misfit functions and their corresponding adjoint sources.

## 🚀 Features

- **Multiple Misfit Functions**: Comprehensive implementations of various misfit measurements
  - Exponentiated Phase Misfit
  - Multitaper Travel-time Misfit  
  - Cross-correlation Travel-time Misfit
  - Waveform Misfit (L2 norm)
  - Waveform Convolution Misfit
  - Receiver Function Misfit

- **High Performance**: Optimized Fortran implementations with efficient memory management
- **CMake Build System**: Cross-platform compilation with multiple compiler support
- **SAC Format Support**: Native support for SAC seismogram format
- **Python Integration**: Compatible with pyadjoint for validation and comparison

## 📦 Installation

### Prerequisites

- **Fortran Compiler**: Modern Fortran compiler (gfortran 9+, ifort 19+, or similar)
- **CMake**: Version 3.20 or higher  
- **LAPACK**: Linear algebra package

### Building from Source

```bash
# Clone the repository
git clone https://github.com/xumi1993/ForAdjoint.git
cd ForAdjoint

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build the library and examples
make

# Optional: Run tests
make test
```

### Compiler-Specific Flags

The build system automatically configures optimal flags for different compilers:

- **GNU Fortran**: `-std=f2008 -O3 -fimplicit-none -pedantic`
- **Intel Fortran**: `-std08 -O3 -xHost -fpe0 -assume buffered_io`
- **NVIDIA Fortran**: `-std=f2008 -O3 -Mcpu=auto -implicitnone`

## 📖 Usage

### Basic Example: Exponentiated Phase Misfit

```fortran
program example
  use exponentiated_phase_misfit
  use sacio
  implicit none
  
  type(ExponentiatedPhaseMisfit) :: epm
  real(kind=dp), dimension(:), allocatable :: observed, synthetic
  real(kind=dp), dimension(1,2) :: time_window
  real(kind=dp) :: dt
  
  ! Load seismogram data (SAC format)
  call read_sac("observed.sac", observed, dt)
  call read_sac("synthetic.sac", synthetic, dt)
  
  ! Define time window [start_time, end_time]
  time_window(1,:) = [50.0, 100.0]
  
  ! Calculate adjoint source
  call epm%calc_adjoint_source(observed, synthetic, dt, time_window)
  
  ! Access results
  print *, "Total misfit:", epm%total_misfit
  print *, "Adjoint source length:", size(epm%adj_src)
  
end program
```

### Command Line Tools

After building, executable tools are available in the `bin/` directory:

```bash
# Exponentiated phase misfit
./bin/xex_exp_ph_misfit observed.sac synthetic.sac 5.0 50.0 42.8 91.1 output_dir/

# Cross-correlation travel-time misfit  
./bin/xex_cc_misfit observed.sac synthetic.sac 5.0 50.0 42.8 91.1 output_dir/
```

**Parameters:**

- `observed.sac`: Path to observed seismogram
- `synthetic.sac`: Path to synthetic seismogram  
- `5.0 50.0`: Filter periods (short_period long_period)
- `42.8 91.1`: Time window (start_time end_time)
- `output_dir/`: Output directory for adjoint source

## 🧮 Supported Misfit Functions

### 1. Exponentiated Phase Misfit

Measures phase and amplitude differences simultaneously using complex-valued approach:

```math
χ(d,s) = 1/2 ∫ |d/|d|_w - s/|s|_w|² dt + 1/2 ∫ |H[d]/|d|_w - H[s]/|s|_w|² dt
```

**Features:**

- Robust to amplitude variations
- Suitable for regional and teleseismic waveforms
- Water-level stabilization for envelope normalization

### 2. Multitaper Travel-time Misfit

Uses multitaper method for stable frequency-domain measurements:

**Features:**

- Multiple orthogonal tapers for spectral estimation
- Travel-time and amplitude anomaly measurements
- Statistical error estimation
- Cycle-skipping detection

### 3. Cross-correlation Travel-time Misfit

Classic cross-correlation based travel-time measurement:

**Features:**

- Maximum cross-correlation lag measurement
- Amplitude ratio estimation
- Robust time-shift detection

### 4. Waveform Misfit

Simple L2 norm waveform difference:

```math
χ(d,s) = 1/2 ∫ |d(t) - s(t)|² dt
```

### 5. Receiver Function Misfit

Specialized for receiver function analysis:

**Features:**

- Deconvolution-based receiver function calculation
- Water-level deconvolution
- Radial/vertical component processing

## 🔧 Configuration Parameters

Key parameters can be configured through the adjoint configuration module:

```fortran
! Filter parameters
min_period = 5.0        ! Minimum period (s)
max_period = 50.0       ! Maximum period (s)

! Taper parameters  
taper_percentage = 0.3  ! Taper percentage (0.0-1.0)
itaper_type = 1        ! 1=Hanning, 2=Hamming, 3=Cosine, 4=Cosine^10

! Water level thresholds
water_threshold = 0.02  ! Deconvolution water level
wtr_env = 0.2          ! Envelope water level

! Multitaper parameters
mt_nw = 4.0           ! Time-bandwidth product
phase_step = 1.5      ! Phase step for cycle skipping
```

## 🧪 Testing and Validation

The library includes comprehensive test cases and validation against [pyadjoint](https://github.com/adjtomo/pyadjoint):

```bash
# Run Python validation tests (requires pyadjoint)
cd tests/
python test_comparison.py
```

Example validation for exponentiated phase misfit:

```python
import pyadjoint
import obspy

# Load test data
observed = obspy.read("example_data/TA.A38.BXZ.sac.obs")[0]
synthetic = obspy.read("example_data/TA.A38.BXZ.sac.syn")[0]

# Calculate with pyadjoint
config = pyadjoint.get_config(adjsrc_type="exponentiated_phase")
adj_py = pyadjoint.calculate_adjoint_source(config=config, 
                                           observed=observed, 
                                           synthetic=synthetic,
                                           windows=[(42.8, 91.1)])

# Compare with ForAdjoint results
# Results should match within numerical precision
```

## 🏗️ Project Structure

```text
ForAdjoint/
├── CMakeLists.txt              # Main CMake configuration
├── LICENSE                     # MIT License
├── README.md                   # This file
├── bin/                        # Compiled executables
├── build/                      # Build directory (generated)
├── src/                        # Source code
│   ├── adjoint_source_types/   # Misfit function implementations
│   │   ├── exponentiated_phase_misfit.f90
│   │   ├── mt_tt_misfit.f90
│   │   ├── cc_tt_misfit.f90
│   │   ├── waveform_misfit.f90
│   │   ├── waveform_conv_misfit.f90
│   │   └── rf_misfit.f90
│   └── utils/                  # Utility modules
│       ├── adj_config.f90      # Configuration parameters
│       ├── sacio.f90           # SAC I/O operations
│       ├── signal.f90          # Signal processing
│       ├── fftpack.f90         # FFT routines
│       └── dpss.f90            # Multitaper sequences
├── examples/                   # Example programs
│   ├── ex_exp_ph_misfit.f90
│   └── ex_cc_misfit.f90
├── tests/                      # Test cases and validation
├── example_data/               # Sample seismograms
└── docs/                       # Documentation
```

## 🔬 Scientific Applications

ForAdjoint is designed for:

- **Full Waveform Inversion (FWI)**: High-resolution subsurface imaging
- **Adjoint Tomography**: Large-scale seismic structure inversion  
- **Source Inversion**: Earthquake source parameter estimation
- **Receiver Function Analysis**: Crustal and mantle discontinuity imaging

## 🤝 Contributing

We welcome contributions! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Development Guidelines

- Follow Fortran 2008 standard
- Include comprehensive comments
- Add test cases for new features
- Ensure memory management is handled properly
- Validate against pyadjoint when applicable

## 📚 References

1. **Adjoint Methods**: Tromp, J., Tape, C., & Liu, Q. (2005). Seismic tomography, adjoint methods, time reversal and banana-doughnut kernels. *Geophysical Journal International*.

2. **Exponentiated Phase**: Yuan, Y. O., et al. (2016). Multiscale adjoint waveform tomography for surface and body waves. *Geophysics*.

3. **Multitaper Method**: Park, J., et al. (1987). Multiple-taper spectral analysis of high-frequency seismograms. *Journal of Geophysical Research*.

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Authors

- **Mijian Xu** - *Main Developer* - [xumi1993](https://github.com/xumi1993)

## 🙏 Acknowledgments

- [pyadjoint](https://github.com/adjtomo/pyadjoint) team for reference implementations
- ObsPy community for seismological data handling standards
- The broader seismological community for theoretical foundations

## 📧 Contact

For questions, suggestions, or collaboration opportunities:

- **Issues**: [GitHub Issues](https://github.com/xumi1993/ForAdjoint/issues)
