cmake_minimum_required (VERSION 2.6)

project(contig)

set (APP_NAME contig)

set (HDR_DIR ${PROJECT_SOURCE_DIR}/contig)
	
set (TST_DIR ${PROJECT_SOURCE_DIR}/tests)

set (CMN_DIR ${PROJECT_SOURCE_DIR}/utils)
	
# DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL
set(CMAKE_BUILD_TYPE RELEASE)
IF(UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSEIF(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qstd=c++11")
ENDIF(UNIX)
	
include_directories(
	${PROJECT_SOURCE_DIR}/utils
	${PROJECT_SOURCE_DIR}/contig
	${PROJECT_SOURCE_DIR}/../../share/zlib-1.2.8
)

SET( SOURCES
	${CMN_DIR}/fasta_reader.cpp
	main.cpp
)

SET ( COMMON
	${CMN_DIR}/fasta_reader.h
	${CMN_DIR}/kseq.h
	${CMN_DIR}/read.h
	unittest.h
	)

SET( ANNOTATION
	${HDR_DIR}/annotation/annotation.hpp
	${HDR_DIR}/annotation/annotation_node.hpp
	${HDR_DIR}/annotation/annotation_record.hpp
)

SET( TRIE
	${HDR_DIR}/trie/trie_node.hpp
	${HDR_DIR}/trie/trie_impl.hpp
	${HDR_DIR}/trie/trie.hpp
	${HDR_DIR}/trie/trie_const_iterator.hpp
	${HDR_DIR}/trie/trie_iterator.hpp
	${HDR_DIR}/trie/trie_cont.hpp
)

SET( KSTAT
	${HDR_DIR}/kstat/kstat.hpp
)

SET( STORAGE
	${HDR_DIR}/contig.hpp
	${HDR_DIR}/contig_iterator.hpp
	${HDR_DIR}/contig_const_iterator.hpp
)

SET( ALIGN
	${HDR_DIR}/alicont/alialgo.hpp
	${HDR_DIR}/alicont/alicont.hpp
	${HDR_DIR}/alicont/alicont_line.hpp
	${HDR_DIR}/alicont/alimatrix.hpp
	${HDR_DIR}/alicont/score_matrix.hpp
)

SET ( TESTS
	${TST_DIR}/annotation_test.cpp
	${TST_DIR}/contig_test.cpp
	${TST_DIR}/fasta_test.cpp
	${TST_DIR}/kstat_test.cpp
	${TST_DIR}/trie_test.cpp
	${TST_DIR}/alicont_test.cpp
	)

IF(UNIX)
SET ( REQUIRED_LIBRARIES
        z 
)
ELSEIF(WIN32)	
SET ( REQUIRED_LIBRARIES
	zlib
)
ENDIF(UNIX)

source_group ("Annotation Files" FILES ${ANNOTATION})
source_group ("Trie Files" FILES ${TRIE})
source_group ("Kstat Files" FILES ${KSTAT})
source_group ("Storage Files" FILES ${STORAGE})
source_group ("Source Files" FILES ${SOURCES})
source_group ("Common Files" FILES ${COMMON})
source_group ("Test Files" FILES ${TESTS})
source_group ("Align Files" FILES ${ALIGN})

add_executable (${APP_NAME}
	${SOURCES}
	${ANNOTATION}
	${TRIE}
	${KSTAT}
	${STORAGE}
	${ALGORITHM}
	${COMMON}
	${TESTS}
	${ALIGN}
	)

target_link_libraries(${APP_NAME} 
	${REQUIRED_LIBRARIES}
	)
