.c.o:
	gcc -c ${CFLAGS} $<

SedBerg: SedBerg.o iceberg.o fjord.o random.o
	gcc -o SedBerg SedBerg.o iceberg1007.o fjord1007.o random1007.o

SedBerg.o: SedBerg.c constants.h
	gcc -c SedBerg.c

iceberg.o: iceberg1007.c iceberg.h
	gcc -c iceberg1007.c

fjord.o: fjord1007.c fjord.h
	gcc -c fjord1007.c

random.o: random1007.c random.h
	gcc -c random1007.c

run: SedBerg
	./SedBerg > logfile &
 

