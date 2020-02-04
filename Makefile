all: coverage.pdf

coverage.pdf: coverage.cc
	root -b -l -n -q coverage.cc

limit_plot.pdf: limit.cc
	root -b -l -n -q limit.cc+
