# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build"

# Include any dependencies generated for this target.
include EigenTutorial/CMakeFiles/EigenTutorial-test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include EigenTutorial/CMakeFiles/EigenTutorial-test.dir/compiler_depend.make

# Include the progress variables for this target.
include EigenTutorial/CMakeFiles/EigenTutorial-test.dir/progress.make

# Include the compile flags for this target's objects.
include EigenTutorial/CMakeFiles/EigenTutorial-test.dir/flags.make

EigenTutorial/CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.o: EigenTutorial/CMakeFiles/EigenTutorial-test.dir/flags.make
EigenTutorial/CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.o: ../EigenTutorial/unit_tests.cc
EigenTutorial/CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.o: EigenTutorial/CMakeFiles/EigenTutorial-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object EigenTutorial/CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.o"
	cd "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build/EigenTutorial" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT EigenTutorial/CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.o -MF CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.o.d -o CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.o -c "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/EigenTutorial/unit_tests.cc"

EigenTutorial/CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.i"
	cd "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build/EigenTutorial" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/EigenTutorial/unit_tests.cc" > CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.i

EigenTutorial/CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.s"
	cd "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build/EigenTutorial" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/EigenTutorial/unit_tests.cc" -o CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.s

# Object files for target EigenTutorial-test
EigenTutorial__test_OBJECTS = \
"CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.o"

# External object files for target EigenTutorial-test
EigenTutorial__test_EXTERNAL_OBJECTS =

Build/bin/EigenTutorial-test: EigenTutorial/CMakeFiles/EigenTutorial-test.dir/unit_tests.cc.o
Build/bin/EigenTutorial-test: EigenTutorial/CMakeFiles/EigenTutorial-test.dir/build.make
Build/bin/EigenTutorial-test: lib/libgtest.a
Build/bin/EigenTutorial-test: lib/libgtest_main.a
Build/bin/EigenTutorial-test: lib/libgtest.a
Build/bin/EigenTutorial-test: EigenTutorial/CMakeFiles/EigenTutorial-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../Build/bin/EigenTutorial-test"
	cd "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build/EigenTutorial" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/EigenTutorial-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
EigenTutorial/CMakeFiles/EigenTutorial-test.dir/build: Build/bin/EigenTutorial-test
.PHONY : EigenTutorial/CMakeFiles/EigenTutorial-test.dir/build

EigenTutorial/CMakeFiles/EigenTutorial-test.dir/clean:
	cd "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build/EigenTutorial" && $(CMAKE_COMMAND) -P CMakeFiles/EigenTutorial-test.dir/cmake_clean.cmake
.PHONY : EigenTutorial/CMakeFiles/EigenTutorial-test.dir/clean

EigenTutorial/CMakeFiles/EigenTutorial-test.dir/depend:
	cd "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6" "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/EigenTutorial" "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build" "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build/EigenTutorial" "/home/gobi/Documents/MyMaster/HS2024/Applied Optimiz/assignements/aopt-exercise6/build/EigenTutorial/CMakeFiles/EigenTutorial-test.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : EigenTutorial/CMakeFiles/EigenTutorial-test.dir/depend

