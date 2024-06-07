rm -rv CMakeCache.txt
rm -rv CMakeFiles

cmake -DCMAKE_PREFIX_PATH=/home/garofalo/analysis/analysis_program/build_sanitizer/install_dir \
      -DCMAKE_BUILD_TYPE=DEBUG\
      -DCMAKE_CXX_FLAGS="-fopenmp   -pedantic  -g    -lm"  \
      -DWITH_ARB=OFF \
      -DCMAKE_VERBOSE_MAKEFILE=ON \
      -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
       \
      ..

