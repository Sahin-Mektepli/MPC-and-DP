# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.27.7/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.27.7/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build"

# Include any dependencies generated for this target.
include dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/compiler_depend.make

# Include the progress variables for this target.
include dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/progress.make

# Include the compile flags for this target's objects.
include dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/flags.make

dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o: dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/flags.make
dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o: /Users/dulumrae/Desktop/Kübra\ Hoca\ veya\ Ağ\ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/eigen/blas/testing/zblat1.f
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o"
	cd "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/dependencies/eigen/blas/testing" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/eigen/blas/testing/zblat1.f" -o CMakeFiles/zblat1.dir/zblat1.f.o

dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/zblat1.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/zblat1.dir/zblat1.f.i"
	cd "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/dependencies/eigen/blas/testing" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/eigen/blas/testing/zblat1.f" > CMakeFiles/zblat1.dir/zblat1.f.i

dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/zblat1.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/zblat1.dir/zblat1.f.s"
	cd "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/dependencies/eigen/blas/testing" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/eigen/blas/testing/zblat1.f" -o CMakeFiles/zblat1.dir/zblat1.f.s

# Object files for target zblat1
zblat1_OBJECTS = \
"CMakeFiles/zblat1.dir/zblat1.f.o"

# External object files for target zblat1
zblat1_EXTERNAL_OBJECTS =

dependencies/eigen/blas/testing/zblat1: dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/zblat1.f.o
dependencies/eigen/blas/testing/zblat1: dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/build.make
dependencies/eigen/blas/testing/zblat1: dependencies/eigen/blas/libeigen_blas.dylib
dependencies/eigen/blas/testing/zblat1: dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable zblat1"
	cd "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/dependencies/eigen/blas/testing" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/zblat1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/build: dependencies/eigen/blas/testing/zblat1
.PHONY : dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/build

dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/clean:
	cd "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/dependencies/eigen/blas/testing" && $(CMAKE_COMMAND) -P CMakeFiles/zblat1.dir/cmake_clean.cmake
.PHONY : dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/clean

dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/depend:
	cd "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/eigen/blas/testing" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/dependencies/eigen/blas/testing" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : dependencies/eigen/blas/testing/CMakeFiles/zblat1.dir/depend

