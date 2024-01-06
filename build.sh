#!/bin/sh

MAKEJOBS=2

rm -r build/python
cp -r src/python build
mkdir build/python/include
cp -r src/msolve src/include/fmpq_mpoly_matrix.h src/fmpq_mpoly_matrix.c build/python
cp -r src/include build/python

cd build/python/msolve
./autogen.sh
./configure
make -j$MAKEJOBS
cd ..

cp msolve/src/*/.libs/*.so cagl
python setup.py bdist_wheel
