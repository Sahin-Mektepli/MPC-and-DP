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
include CMakeFiles/helper.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/helper.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/helper.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/helper.dir/flags.make

CMakeFiles/helper.dir/apps/helper.cpp.o: CMakeFiles/helper.dir/flags.make
CMakeFiles/helper.dir/apps/helper.cpp.o: /Users/dulumrae/Desktop/Kübra\ Hoca\ veya\ Ağ\ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/apps/helper.cpp
CMakeFiles/helper.dir/apps/helper.cpp.o: CMakeFiles/helper.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/helper.dir/apps/helper.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/helper.dir/apps/helper.cpp.o -MF CMakeFiles/helper.dir/apps/helper.cpp.o.d -o CMakeFiles/helper.dir/apps/helper.cpp.o -c "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/apps/helper.cpp"

CMakeFiles/helper.dir/apps/helper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/helper.dir/apps/helper.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/apps/helper.cpp" > CMakeFiles/helper.dir/apps/helper.cpp.i

CMakeFiles/helper.dir/apps/helper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/helper.dir/apps/helper.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/apps/helper.cpp" -o CMakeFiles/helper.dir/apps/helper.cpp.s

# Object files for target helper
helper_OBJECTS = \
"CMakeFiles/helper.dir/apps/helper.cpp.o"

# External object files for target helper
helper_EXTERNAL_OBJECTS =

helper: CMakeFiles/helper.dir/apps/helper.cpp.o
helper: CMakeFiles/helper.dir/build.make
helper: dependencies/cryptopp-cmake/cryptopp/libcryptopp.a
helper: CMakeFiles/helper.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable helper"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/helper.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/helper.dir/build: helper
.PHONY : CMakeFiles/helper.dir/build

CMakeFiles/helper.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/helper.dir/cmake_clean.cmake
.PHONY : CMakeFiles/helper.dir/clean

CMakeFiles/helper.dir/depend:
	cd "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/build/CMakeFiles/helper.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : CMakeFiles/helper.dir/depend

