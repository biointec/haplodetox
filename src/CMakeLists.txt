add_executable(haplodetox main.cpp settings.cpp dbgraph.cpp dsnode.cpp kmernpp.cpp
               util.cpp debug.cpp coverage.cpp alignment.cpp graphaln.cpp
               refcomp.cpp  correctgraph.cpp crfmult.cpp readaln.cpp)

target_link_libraries(haplodetox readfile kmer pgm pthread libdai gmp)

if (ZLIB_FOUND)
        target_link_libraries(haplodetox ${ZLIB_LIBRARY})
endif (ZLIB_FOUND)

target_link_libraries(haplodetox ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS haplodetox RUNTIME DESTINATION bin)
add_subdirectory(readfile)
add_subdirectory(kmer)
add_subdirectory(pgm)
add_subdirectory(libdai)

