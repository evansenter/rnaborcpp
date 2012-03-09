RNAhairpin: convert_Vienna.o energy_par.o fold.o fold_vars.o misc.o pair_mat.o params.o ParCal.o RNAhairpin.o utils.o McCaskill.o
	g++ -m32 -o RNAhairpin -lm convert_Vienna.o energy_par.o fold.o fold_vars.o ParCal.o McCaskill.o misc.o pair_mat.o params.o RNAhairpin.o utils.o

convert_Vienna.o: convert_Vienna.cpp convert_Vienna.h
	g++ -m32 -c convert_Vienna.cpp

energy_par.o: energy_par.cpp energy_par.h
	g++ -m32 -c energy_par.cpp

fold.o: fold.cpp fold.h
	g++ -m32 -c fold.cpp

fold_vars.o: fold_vars.cpp fold_vars.h
	g++ -m32 -c fold_vars.cpp

ParCal.o: ParCal.cpp ParCal.h
	g++ -m32 -c ParCal.cpp

McCaskill.o: McCaskill.cpp McCaskill.h
	g++ -m32 -c McCaskill.cpp

misc.o: misc.cpp misc.h
	g++ -m32 -c misc.cpp

pair_mat.o: pair_mat.cpp pair_mat.h
	g++ -m32 -c pair_mat.cpp

params.o: params.cpp params.h
	g++ -m32 -c params.cpp

RNAhairpin.o: RNAhairpin.cpp RNAhairpin.h
	g++ -m32 -c RNAhairpin.cpp

utils.o: utils.cpp utils.h
	g++ -m32 -c utils.cpp

install: RNAhairpin
	cp RNAhairpin /usr/local/bin
	cp RNAhairpin.man /usr/local/man/man1/RNAhairpin.1
	
clean:
	rm -fr *.o
	rm -f ./RNAhairpin