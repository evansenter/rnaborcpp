RNAhairpin: convert_Vienna.o energy_par.o fold.o fold_vars.o misc.o pair_mat.o params.o ParCal.o RNAhairpin.o utils.o McCaskill.o
	g++ -L /usr/local/lib -llapackpp -o RNAhairpin -lm convert_Vienna.o energy_par.o fold.o fold_vars.o ParCal.o McCaskill.o misc.o pair_mat.o params.o RNAhairpin.o utils.o

convert_Vienna.o: convert_Vienna.cpp convert_Vienna.h
	g++ -c convert_Vienna.cpp

energy_par.o: energy_par.cpp energy_par.h
	g++ -c energy_par.cpp

fold.o: fold.cpp fold.h
	g++ -c fold.cpp

fold_vars.o: fold_vars.cpp fold_vars.h
	g++ -c fold_vars.cpp

ParCal.o: ParCal.cpp ParCal.h
	g++ -c ParCal.cpp

McCaskill.o: McCaskill.cpp McCaskill.h
	g++ -I /usr/local/include/lapackpp -c McCaskill.cpp

misc.o: misc.cpp misc.h
	g++ -c misc.cpp

pair_mat.o: pair_mat.cpp pair_mat.h
	g++ -c pair_mat.cpp

params.o: params.cpp params.h
	g++ -c params.cpp

RNAhairpin.o: RNAhairpin.cpp RNAhairpin.h
	g++ -c RNAhairpin.cpp

utils.o: utils.cpp utils.h
	g++ -c utils.cpp

install: RNAhairpin
	cp RNAhairpin /usr/local/bin
	cp RNAhairpin.man /usr/local/man/man1/RNAhairpin.1
	
clean:
	rm -fr *.o
	rm -f ./RNAhairpin