# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/_deps/googletest-src"
  "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/_deps/googletest-build"
  "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/_deps/googletest-subbuild/googletest-populate-prefix"
  "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/_deps/googletest-subbuild/googletest-populate-prefix/tmp"
  "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
  "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/_deps/googletest-subbuild/googletest-populate-prefix/src"
  "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
