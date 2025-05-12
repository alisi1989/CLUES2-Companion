# CMake generated Testfile for 
# Source directory: /Users/alessandrolisi1989/desktop/Michael_Covid_Polaris/Relate-Clues/relate-master
# Build directory: /Users/alessandrolisi1989/desktop/Michael_Covid_Polaris/Relate-Clues/relate-master/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTest "/Users/alessandrolisi1989/desktop/Michael_Covid_Polaris/Relate-Clues/relate-master/bin/Tests")
set_tests_properties(UnitTest PROPERTIES  _BACKTRACE_TRIPLES "/Users/alessandrolisi1989/desktop/Michael_Covid_Polaris/Relate-Clues/relate-master/CMakeLists.txt;54;add_test;/Users/alessandrolisi1989/desktop/Michael_Covid_Polaris/Relate-Clues/relate-master/CMakeLists.txt;0;")
subdirs("include/src")
subdirs("include/src/gzstream")
subdirs("include/file_formats/tskit")
subdirs("include/test")
subdirs("include/pipeline")
subdirs("include/evaluate")
subdirs("include/treeview")
subdirs("include/file_formats")
subdirs("include/extract")
