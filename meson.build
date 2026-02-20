project('pyenkf', 'c',
  version : '0.0.1',
  meson_version: '>=0.64.0',
  default_options: [
    'buildtype=debugoptimized']
 )

add_languages('fortran')
add_project_arguments(['-ffree-line-length-none','-fopenmp','-O2','-fdefault-real-8','-march=znver3'],language : 'fortran')
blas_dep=dependency('blas','openblas')

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

 incdir_numpy = run_command(py,
      [
       '-c',
       '''import os; os.chdir(".."); import numpy; print(numpy.get_include()) ''' ],
     check: true
   ).stdout().strip() 


incdir_f2py =run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

pyenkf_source = custom_target('pyenkfmodule.c',
  input : ['src/pyenkf/enksrf.F90'],  # .f so no F90 wrappers
  output : ['pyenkfmodule.c','pyenkf-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'pyenkf', '--lower']
)

inc_np = include_directories(incdir_numpy, incdir_f2py)

py.extension_module('pyenkf',
  ['src/pyenkf/enksrf.F90', pyenkf_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : [py_dep, blas_dep],
  install : true
)
