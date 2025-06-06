all: coverage.pdf \
  limit_plot.pdf limit_sensitivity_plot.pdf limit_sensitivity_fast_plot.pdf \
  scatterr2.pdf scatterfmax.pdf r2min-n-1.pdf \
  fmax-n-1.pdf fig-thetax.pdf fmax-n-2.pdf

fig-thetax.pdf: fig-thetax.C
	root -l -b -n -q fig-thetax.C

coverage.pdf: coverage.cc
	root -b -l -n -q coverage.cc

limit_plot.pdf: limit icecube.txt antares.txt
	./limit

limit: limit.o Event_Info.o
	g++ `root-config --libs` limit.o Event_Info.o -o limit

limit.o: limit.cc Constants.hh
	g++ -Wall -Wextra -Werror `root-config --cflags` -c limit.cc




limit_sensitivity_fast_plot.pdf: limit_sensitivity_fast icecube.txt antares.txt
	./limit_sensitivity_fast

limit_sensitivity_slow_plot.pdf: limit_sensitivity_slow icecube.txt antares.txt
	./limit_sensitivity_slow

limit_sensitivity_slowfast_plot.pdf: limit_sensitivity_slowfast icecube.txt antares.txt
	./limit_sensitivity_slowfast


limit_sensitivity_fast: limit_sensitivity_fast.o Event_Info.o make_sens_fastonly.C limitsensitivityfastdata.txt
	g++ `root-config --libs` limit_sensitivity_fast.o Event_Info.o -o limit_sensitivity_fast

limit_sensitivity_slow: limit_sensitivity_slow.o Event_Info.o make_sens_slowonly.C limitsensitivityslowdata.txt
	g++ `root-config --libs` limit_sensitivity_slow.o Event_Info.o -o limit_sensitivity_slow

limit_sensitivity_slowfast: limit_sensitivity_slowfast.o Event_Info.o make_sens_slowonly.C make_sens_fastonly.C limitsensitivityslowdata.txt limitsensitivityfastdata.txt
	g++ `root-config --libs` limit_sensitivity_slowfast.o Event_Info.o -o limit_sensitivity_slowfast


limitsensitivityslowdata.txt: make_sens_slowonly.C
	root -b -q -l -n make_sens_slowonly.C | awk 'NF == 3' | tee limitsensitivityslowdata.txt

limitsensitivityfastdata.txt: make_sens_fastonly.C
	root -b -q -l -n make_sens_fastonly.C | awk 'NF == 3' | tee limitsensitivityfastdata.txt


limit_sensitivity_fast.o: limit_sensitivity.cc Constants.hh
	g++ -DDRAWFAST -Wall -Wextra -Werror `root-config --cflags` -c limit_sensitivity.cc -o limit_sensitivity_fast.o

limit_sensitivity_slow.o: limit_sensitivity.cc Constants.hh
	g++ -DDRAWSLOW -Wall -Wextra -Werror `root-config --cflags` -c limit_sensitivity.cc -o limit_sensitivity_slow.o

limit_sensitivity_slowfast.o: limit_sensitivity.cc Constants.hh
	g++ -DDRAWSLOW -DDRAWFAST -Wall -Wextra -Werror `root-config --cflags` -c limit_sensitivity.cc -o limit_sensitivity_slowfast.o


Event_Info.o: Event_Info.cc Constants.hh
	g++ -Wall -Wextra -Werror `root-config --cflags` -c Event_Info.cc

Event_List.o: Event_List.cc Constants.hh
	g++ -Wall -Wextra -Werror `root-config --cflags` -c Event_List.cc

scatterfmax.pdf: scatter
	./scatter fmax

scatterr2.pdf: scatter
	./scatter r2

scatter: scatter.o Event_Info.o Event_List.o
	g++ -O3 `root-config --libs` scatter.o Event_List.o Event_Info.o -o scatter

scatter.o: scatter.cc Event_Info.hh Event_List.hh MFRoot.hh
	g++ -O3 -Wall -Wextra -Werror `root-config --cflags` -c scatter.cc

r2min-n-1.pdf fmax-n-1.pdf fmax-n-2.pdf: compare_data_and_mc.cc Event_List.cc \
 Event_Info.cc Constants.hh Event_Info.hh Event_List.hh
	root -b -l -n -q compare_data_and_mc.cc+O

clean:
	rm -f *_C* *.d *.so *.o scatter limit *.pdf *.pcm
