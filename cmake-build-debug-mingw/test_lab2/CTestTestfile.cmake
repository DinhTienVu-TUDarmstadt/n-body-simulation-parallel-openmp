# CMake generated Testfile for 
# Source directory: D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/test_lab2
# Build directory: D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/test_lab2
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(SerialTests "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/cmake-build-debug-mingw/bin/lab2_test.exe")
set_tests_properties(SerialTests PROPERTIES  _BACKTRACE_TRIPLES "D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/test_lab2/CMakeLists.txt;66;add_test;D:/TU Darmstadt/Semester3/Parallele Programmierung/lab_2/test_lab2/CMakeLists.txt;0;")
subdirs("../_deps/googletest-build")
