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

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

#extensions = [ Extension("cagl", ["src/cython/cagl.pyx","src/fmpq_mpoly_matrix.c"], libraries=["flint", "gmp" ]) ]
extensions = [
    Extension(
        "cagl",
        ["src/cython/cagl.pyx", "src/fmpq_mpoly_matrix.c" ],
        libraries=["flint", "gmp", 'msolve', 'neogb' ],
        library_dirs=['libs'],
        extra_link_args=["-Wl,-rpath=$ORIGIN/libs"]
    )]

setup(
    name='cagl',
    version='0.0.2022a0',
    ext_modules=cythonize(extensions, compiler_directives={'language_level' : "3"}, build_dir="build" ),
    include_dirs=[np.get_include()]+[ "src/msolve/src/"+dir for dir in ["msolve", "usolve", "neogb", "fglm" ] ],
    package_data={'.': ['libs/*.so']},
)

