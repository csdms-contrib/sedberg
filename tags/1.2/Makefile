PROGRAM=SedBerg
SOURCES=SedBerg.c iceberg1007.c fjord1007.c random1007.c
OBJS=$(SOURCES:.c=.o)

CC=gcc
CFLAGS=-O3
LIBS=-lm

.c.o:
	$(CC) -c ${CFLAGS} $<

$(PROGRAM): $(OBJS)
	$(CC) -o $(PROGRAM) $(OBJS) $(LIBS)

run: SedBerg
	./SedBerg > logfile &
 
clean:
	@rm -f $(OBJS) $(PROGRAM)

