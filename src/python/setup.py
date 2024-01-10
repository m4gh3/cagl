''' This file is part of cagl.
 *
 * cagl is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * cagl is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cagl.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Lorenzo Magherini (m4gh3) '''

#NOTE: this file and this whole directory is mean to be copied to the build/ directory along with some other c files

from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np

#extensions = [ Extension("cagl", ["src/cython/cagl.pyx","src/fmpq_mpoly_matrix.c"], libraries=["flint", "gmp" ]) ]
extensions = [
    Extension(
        "cagl._cagl",
        ["cagl/_cagl.pyx", "fmpq_mpoly_matrix.c" ],
        libraries=["flint", "gmp", 'msolve', 'neogb' ],
        library_dirs=['cagl'],
        extra_link_args=["-Wl,-rpath=$ORIGIN"]
    )]

setup(
    name='cagl',
    version='0.0.2022a0',
    ext_modules=cythonize(extensions, compiler_directives={'language_level' : "3"} ),
    include_dirs=[np.get_include()]+[ "src/msolve/src/"+dir for dir in ["msolve", "usolve", "neogb", "fglm" ] ],
    packages=find_packages(), #['cagl',],
    include_package_data=True,
    package_data={'cagl': ['*.so']},
)

