cmake_minimum_required(VERSION 2.6)

project(gene_combinator)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(gene_combinator
	main.cpp
	read.cpp
	generic_exception.cpp
	fasta_reader.cpp
	reverse_complement.cpp)

target_link_libraries(gene_combinator z)
