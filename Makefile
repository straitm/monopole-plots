all: coverage.pdf

coverage.pdf: coverage.cc
	root -b -l -n -q coverage.cc
