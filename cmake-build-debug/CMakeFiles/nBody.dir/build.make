# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = "/Users/linus/Library/Application Support/JetBrains/Toolbox/apps/CLion/ch-0/202.7660.37/CLion.app/Contents/bin/cmake/mac/bin/cmake"

# The command to remove a file.
RM = "/Users/linus/Library/Application Support/JetBrains/Toolbox/apps/CLion/ch-0/202.7660.37/CLion.app/Contents/bin/cmake/mac/bin/cmake" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/linus/CLionProjects/nBody

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/linus/CLionProjects/nBody/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/nBody.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/nBody.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/nBody.dir/flags.make

CMakeFiles/nBody.dir/main.cpp.o: CMakeFiles/nBody.dir/flags.make
CMakeFiles/nBody.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/linus/CLionProjects/nBody/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/nBody.dir/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nBody.dir/main.cpp.o -c /Users/linus/CLionProjects/nBody/main.cpp

CMakeFiles/nBody.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nBody.dir/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/linus/CLionProjects/nBody/main.cpp > CMakeFiles/nBody.dir/main.cpp.i

CMakeFiles/nBody.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nBody.dir/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/linus/CLionProjects/nBody/main.cpp -o CMakeFiles/nBody.dir/main.cpp.s

CMakeFiles/nBody.dir/body.cpp.o: CMakeFiles/nBody.dir/flags.make
CMakeFiles/nBody.dir/body.cpp.o: ../body.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/linus/CLionProjects/nBody/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/nBody.dir/body.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nBody.dir/body.cpp.o -c /Users/linus/CLionProjects/nBody/body.cpp

CMakeFiles/nBody.dir/body.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nBody.dir/body.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/linus/CLionProjects/nBody/body.cpp > CMakeFiles/nBody.dir/body.cpp.i

CMakeFiles/nBody.dir/body.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nBody.dir/body.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/linus/CLionProjects/nBody/body.cpp -o CMakeFiles/nBody.dir/body.cpp.s

CMakeFiles/nBody.dir/io.cpp.o: CMakeFiles/nBody.dir/flags.make
CMakeFiles/nBody.dir/io.cpp.o: ../io.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/linus/CLionProjects/nBody/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/nBody.dir/io.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nBody.dir/io.cpp.o -c /Users/linus/CLionProjects/nBody/io.cpp

CMakeFiles/nBody.dir/io.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nBody.dir/io.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/linus/CLionProjects/nBody/io.cpp > CMakeFiles/nBody.dir/io.cpp.i

CMakeFiles/nBody.dir/io.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nBody.dir/io.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/linus/CLionProjects/nBody/io.cpp -o CMakeFiles/nBody.dir/io.cpp.s

CMakeFiles/nBody.dir/math_utils.cpp.o: CMakeFiles/nBody.dir/flags.make
CMakeFiles/nBody.dir/math_utils.cpp.o: ../math_utils.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/linus/CLionProjects/nBody/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/nBody.dir/math_utils.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nBody.dir/math_utils.cpp.o -c /Users/linus/CLionProjects/nBody/math_utils.cpp

CMakeFiles/nBody.dir/math_utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nBody.dir/math_utils.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/linus/CLionProjects/nBody/math_utils.cpp > CMakeFiles/nBody.dir/math_utils.cpp.i

CMakeFiles/nBody.dir/math_utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nBody.dir/math_utils.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/linus/CLionProjects/nBody/math_utils.cpp -o CMakeFiles/nBody.dir/math_utils.cpp.s

# Object files for target nBody
nBody_OBJECTS = \
"CMakeFiles/nBody.dir/main.cpp.o" \
"CMakeFiles/nBody.dir/body.cpp.o" \
"CMakeFiles/nBody.dir/io.cpp.o" \
"CMakeFiles/nBody.dir/math_utils.cpp.o"

# External object files for target nBody
nBody_EXTERNAL_OBJECTS =

nBody: CMakeFiles/nBody.dir/main.cpp.o
nBody: CMakeFiles/nBody.dir/body.cpp.o
nBody: CMakeFiles/nBody.dir/io.cpp.o
nBody: CMakeFiles/nBody.dir/math_utils.cpp.o
nBody: CMakeFiles/nBody.dir/build.make
nBody: /usr/local/lib/libmpi.dylib
nBody: CMakeFiles/nBody.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/linus/CLionProjects/nBody/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable nBody"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nBody.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/nBody.dir/build: nBody

.PHONY : CMakeFiles/nBody.dir/build

CMakeFiles/nBody.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/nBody.dir/cmake_clean.cmake
.PHONY : CMakeFiles/nBody.dir/clean

CMakeFiles/nBody.dir/depend:
	cd /Users/linus/CLionProjects/nBody/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/linus/CLionProjects/nBody /Users/linus/CLionProjects/nBody /Users/linus/CLionProjects/nBody/cmake-build-debug /Users/linus/CLionProjects/nBody/cmake-build-debug /Users/linus/CLionProjects/nBody/cmake-build-debug/CMakeFiles/nBody.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/nBody.dir/depend

