FFTbor: convert_Vienna.o energy_par.o fold.o fold_vars.o misc.o pair_mat.o params.o RNAbor.o utils.o McCaskill.o
	g++ -o FFTbor -lfftw3 -lm convert_Vienna.o energy_par.o fold.o fold_vars.o McCaskill.o misc.o pair_mat.o params.o RNAbor.o utils.o

convert_Vienna.o: convert_Vienna.cpp convert_Vienna.h
	g++ -c convert_Vienna.cpp

energy_par.o: energy_par.cpp energy_par.h
	g++ -c energy_par.cpp

fold.o: fold.cpp fold.h
	g++ -c fold.cpp

fold_vars.o: fold_vars.cpp fold_vars.h
	g++ -c fold_vars.cpp

McCaskill.o: McCaskill.cpp McCaskill.h
	g++ -c McCaskill.cpp

misc.o: misc.cpp misc.h
	g++ -c misc.cpp

pair_mat.o: pair_mat.cpp pair_mat.h
	g++ -c pair_mat.cpp

params.o: params.cpp params.h
	g++ -c params.cpp

RNAbor.o: RNAbor.cpp RNAbor.h
	g++ -c RNAbor.cpp

utils.o: utils.cpp utils.h
	g++ -c utils.cpp

install: RNAbor
	cp FFTbor /usr/local/bin
	cp FFTbor.man /usr/local/man/man1/FFTbor.1
	
clean:
	rm -fr *.o
	rm -f ./FFTbor