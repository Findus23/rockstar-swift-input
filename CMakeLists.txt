cmake_minimum_required(VERSION 3.10)
set(CMAKE_BUILD_TYPE DEBUG)
# set the project name
project(rockstar_cmake)

add_definitions(-DENABLE_HDF5 -DH5_USE_16_API -DENABLE_HDF5 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D _DEFAULT_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L  -D_DARWIN_C_SOURCE)
# add the executable

find_package(HDF5 REQUIRED)


add_executable(rockstar_cmake main.c rockstar.c check_syscalls.c fof.c groupies.c subhalo_metric.c potential.c nfw.c jacobi.c fun_times.c interleaving.c universe_time.c hubble.c integrate.c distance.c config_vars.c config.c bounds.c inthash.c io/read_config.c client.c server.c merger.c inet/socket.c inet/rsocket.c inet/address.c io/meta_io.c io/io_internal.c io/io_ascii.c io/stringparse.c io/io_gadget.c io/io_generic.c io/io_art.c io/io_bgc2.c io/io_util.c io/io_arepo.c io/io_swift.c io/io_hdf5.c)

target_include_directories(rockstar_cmake PRIVATE ${HDF5_INCLUDE_DIRS})

target_link_libraries(rockstar_cmake PRIVATE ${HDF5_LIBRARIES})
#target_link_libraries(rockstar_cmake PRIVATE m)

target_compile_options(rockstar_cmake PRIVATE -Wall -Wextra)
