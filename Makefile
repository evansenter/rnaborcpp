FFTbor: convert_Vienna.o energy_par.o fold.o fold_vars.o misc.o pair_mat.o params.o RNAbor.o utils.o McCaskill.o
	g++ -O3 -o FFTbor -lfftw3 -lm convert_Vienna.o energy_par.o fold.o fold_vars.o McCaskill.o misc.o pair_mat.o params.o RNAbor.o utils.o

convert_Vienna.o: convert_Vienna.cpp convert_Vienna.h
	g++ -O3 -c convert_Vienna.cpp

energy_par.o: energy_par.cpp energy_par.h
	g++ -O3 -c energy_par.cpp

fold.o: fold.cpp fold.h
	g++ -O3 -c fold.cpp

fold_vars.o: fold_vars.cpp fold_vars.h
	g++ -O3 -c fold_vars.cpp

McCaskill.o: McCaskill.cpp McCaskill.h
	g++ -O3 -c McCaskill.cpp

misc.o: misc.cpp misc.h
	g++ -O3 -c misc.cpp

pair_mat.o: pair_mat.cpp pair_mat.h
	g++ -O3 -c pair_mat.cpp

params.o: params.cpp params.h
	g++ -O3 -c params.cpp

RNAbor.o: RNAbor.cpp RNAbor.h
	g++ -O3 -c RNAbor.cpp

utils.o: utils.cpp utils.h
	g++ -O3 -c utils.cpp

install: RNAbor
	cp FFTbor /usr/local/bin
	cp FFTbor.man /usr/local/man/man1/FFTbor.1
	
clean:
	rm -fr *.o
	rm -f ./FFTbor