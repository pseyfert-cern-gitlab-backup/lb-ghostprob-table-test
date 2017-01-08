################################################################################
# Package: TrackKernel
################################################################################
gaudi_subdir(TrackKernel v3r1)

gaudi_depends_on_subdirs(Det/DetDesc
			 Event/TrackEvent
                         Kernel/LHCbKernel
                         Kernel/LHCbMath)

find_package(Boost)
find_package(ROOT)
find_package(GSL)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})

gaudi_add_library(TrackKernel
                  src/*.cpp
                  INCLUDE_DIRS GSL Boost
                  PUBLIC_HEADERS TrackKernel DetDesc
                  LINK_LIBRARIES GSL Boost TrackEvent LHCbKernel LHCbMathLib)

gaudi_add_dictionary(TrackKernel
                     src/TrackKernelDict.h
                     src/TrackKernelDict.xml
                     INCLUDE_DIRS GSL Boost
                     LINK_LIBRARIES GSL Boost TrackEvent LHCbKernel LHCbMathLib TrackKernel
                     OPTIONS "-U__MINGW32__")

