all: coverage.pdf limit_plot.pdf scatter.pdf r2min-n-1.pdf fmax-n-1.pdf

coverage.pdf: coverage.cc
	root -b -l -n -q coverage.cc

limit_plot.pdf: limit
	./limit

limit: limit.o Event_Info.o
	g++ `root-config --libs` limit.o Event_Info.o -o limit

limit.o: limit.cc Constants.hh
	g++ -Wall -Wextra -Werror -O3 `root-config --cflags` -c limit.cc

Event_Info.o: Event_Info.cc Constants.hh
	g++ -Wall -Wextra -Werror `root-config --cflags` -c Event_Info.cc

Event_List.o: Event_List.cc Constants.hh
	g++ -Wall -Wextra -Werror `root-config --cflags` -c Event_List.cc

scatter.pdf: scatter
	./scatter

scatter: scatter.o Event_Info.o Event_List.o
	g++ `root-config --libs` scatter.o Event_List.o Event_Info.o -o scatter

scatter.o: scatter.cc Event_Info.hh Event_List.hh MFRoot.hh
	g++ -Wall -Wextra -Werror `root-config --cflags` -c scatter.cc

r2min-n-1.pdf fmax-n-1.pdf: compare_data_and_mc.cc
	root -b -l -n -q compare_data_and_mc.cc+O

clean:
	rm -f *_C* *.d *.so *.o scatter limit *.pdf
