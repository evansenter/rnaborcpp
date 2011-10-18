

RNAhairpin: convert_Vienna.o energy_par.o fold.o fold_vars.o misc.o pair_mat.o params.o ParCal.o RNAhairpin.o utils.o McCaskill.o
	gcc -o RNAhairpin -lm convert_Vienna.o energy_par.o fold.o fold_vars.o ParCal.o McCaskill.o misc.o pair_mat.o params.o RNAhairpin.o utils.o

convert_Vienna.o: convert_Vienna.c convert_Vienna.h
	gcc -c convert_Vienna.c

energy_par.o: energy_par.c energy_par.h
	gcc -c energy_par.c

fold.o: fold.c fold.h
	gcc -c fold.c

fold_vars.o: fold_vars.c fold_vars.h
	gcc -c fold_vars.c

ParCal.o: ParCal.c ParCal.h
	gcc -c ParCal.c

McCaskill.o: McCaskill.c McCaskill.h
	gcc -c McCaskill.c

misc.o: misc.c misc.h
	gcc -c misc.c

pair_mat.o: pair_mat.c pair_mat.h
	gcc -c pair_mat.c

params.o: params.c params.h
	gcc -c params.c

RNAhairpin.o: RNAhairpin.c RNAhairpin.h
	gcc -c RNAhairpin.c

utils.o: utils.c utils.h
	gcc -c utils.c

install: RNAhairpin
	cp RNAhairpin /usr/local/bin
	cp RNAhairpin.man /usr/local/man/man1/RNAhairpin.1
	
clean:
	rm -fr *.o


