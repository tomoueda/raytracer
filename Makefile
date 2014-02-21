CC = g++

CFLAG = -o as2

INSTALLPATH = -I./include

LINKFLAG = -llibgflags

main: as2.cpp 
	$(CC) $(CFLAG) $(INSTALLPATH) $(LINKFLAG) as2.cpp

test.o: test.cpp
	$(CC) $(INSTALLPATH) $(LINKFLAG) test.cpp

clean: 
	rm -rf *o as2