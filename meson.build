# -*- mode: python -*-

project('DDTBench', 'c', version: '1.2', license: 'NCSA',
        default_options: [ 'c_std=gnu11' ])

inc = include_directories('.')

conf = configuration_data()

mpi_impl = get_option('mpi_impl')
mpi_ver = mpi_impl == 'MPMPICH' ? 1 : 2

message(mpi_ver)

conf.set('MY_MPI_VERSION', mpi_ver)

enable_dbrew = get_option('enable_dbrew')
conf.set('HAVE_DBREW', enable_dbrew)

enable_papi = get_option('enable_papi')
conf.set('HAVE_PAPI', enable_papi)
if enable_papi
  conf.set('TEST_TYPE', 2)
else
  conf.set('TEST_TYPE', 1)
endif

verify_buffers = get_option('verify_buffers')
conf.set('VERIFY_BUFFERS', verify_buffers)
if (not enable_dbrew) and verify_buffers
  error('Buffer verification doesn\'t make sense without enabling DBrew')
endif

enable_all = get_option('enable_all')
milc = get_option('enable_milc') or enable_all
nas = get_option('enable_nas') or enable_all
specfem3d = get_option('enable_specfem3d') or enable_all
lammps = get_option('enable_lammps') or enable_all
fft = get_option('enable_fft')
wrf = get_option('enable_wrf')

# Fail if unsupported benchmarks are explicitly enabled on MPI2
if mpi_ver == 1 and fft
  error('FFT benchmark requires MPI2 compliant implementation')
elif mpi_ver == 1 and wrf
  error('WRF benchmark requires MPI2 compliant implementation')
elif mpi_ver == 2
  fft = fft or enable_all
  wrf = wrf or enable_all
endif

conf.set('WRF', wrf)
conf.set('MILC', milc)
conf.set('NAS', nas)
conf.set('FFT', fft)
conf.set('SPECFEM3D', specfem3d)
conf.set('LAMMPS', lammps)

# Configure MPI compile flags

path = ''
if mpi_impl == 'MPMPICH'
  path = get_option('MPMPICH_path')
elif mpi_impl == 'openmpi'
  path = get_option('openmpi_path')
elif mpi_impl == 'mpich'
  path = get_option('mpich_path')
endif
path = join_paths(path, 'bin')

# Append MPI
#env = environment()
#env.append('PATH', path, separator: ':')

mpicc = find_program(join_paths(path, 'mpicc'))


# TODO: Compiler may be something else than those listed
filtered = ['gcc', 'icc', 'clang', '-c', '']
link_flags = []
compile_flags = []
if mpi_impl == 'MPMPICH' or mpi_impl == 'mpich'
  cflags = run_command(mpicc, '-compile_info').stdout().strip().split(' ')
  foreach f : cflags
    if not filtered.contains(f)
      compile_flags += f
    endif
  endforeach
  ldflags = run_command(mpicc, '-link_info').stdout().strip().split(' ')
  foreach f : ldflags
    if not filtered.contains(f)
      link_flags += f
    endif
  endforeach
elif mpi_impl == 'openmpi'
  cflags = run_command(mpicc, '-showme').stdout().strip().split(' ')
  _flags = []
  foreach f : cflags
    if not filtered.contains(f)
      _flags += f
    endif
  endforeach
  foreach f : _flags
    if f.startswith('-I')
      compile_flags += f
    else
      link_flags += f
    endif
  endforeach
endif

libpapi = []
if enable_papi
  libpapi = dependency('papi', rewiured: true)
endif

# TODO: Do some sanity checks, e.g., check that mpi.h can actually be found

mpi_dep = declare_dependency(link_args: link_flags,
                             compile_args: compile_flags)

libdbrew = []
if enable_dbrew
  libdbrew = dependency('dbrew', required: true)
  conf.set('HAVE_DBREW', true)
endif

configure_file(input: 'config.h.in',
               output: 'config.h',
               configuration: conf)


subdir('src_c')

mpirun = find_program(join_paths(path, 'mpirun'))

test('Run compiled tets', mpirun, args: ['-n', '2', ddtbench.full_path()],
     timeout: 60*60)
