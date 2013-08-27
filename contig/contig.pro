TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    utils/fasta_reader.cpp \
    tests/trie_test.cpp \
    tests/kstat_test.cpp \
    tests/annotation_test.cpp \
    tests/fasta_test.cpp \
    tests/contig_test.cpp \
    tests/alicont_test.cpp

HEADERS += \
    unittest.h \
    utils/read.h \
    utils/kseq.h \
    utils/fasta_reader.h \
    contig/trie/* \
    contig/annotation/* \
    contig/kstat/* \
    contig/contig.hpp \
    contig/contig_const_iterator.hpp \
    contig/contig_iterator.hpp \
    contig/alicont/alicont.hpp \
    contig/alicont/alialgo.hpp \
    contig/alicont/alimatrix.hpp \
    contig/alicont/score_matrix.hpp

#QMAKE_CC = icc
#QMAKE_CXX = icpc

#QMAKE_CC = clang
#QMAKE_CXX = clang++

QMAKE_CXXFLAGS += -std=c++11 -Wall -o parallel
LIBS += -lz -lirc
