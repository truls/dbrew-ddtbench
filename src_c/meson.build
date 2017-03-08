# -*- mode: python -*-

sources = ['ddtbench.c', 'timing.c', 'utilities.c', 'timing_basic.c', 'wrapper.c']

if fft
  sources += 'timing_fft2d.c'
endif
if lammps
  sources += 'timing_lammps.c'
endif
if milc
  sources += 'timing_milc.c'
endif
if nas
  sources += 'timing_nas.c'
endif
if specfem3d
  sources += 'timing_specfem3d.c'
endif
if wrf
  sources += 'timing_wrf.c'
endif

ddtbench = executable('ddtbench', sources, dependencies: [mpi_dep, libdbrew, libpapi], include_directories: inc,
                      c_args: '-DHRT_ARCH=6')
