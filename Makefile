CC = clang++

INSTALLPATH = -I./include/ -I/opt/X11/include -g -DGL_GLEXT_PROTOTYPES -I/usr/X11/include -DOSX

LINKPATH = -L./lib -L/opt/X11/lib -L/usr/X11/lib

LDFLAG = -lgflags -lglog -lgtest -lX11

DYLDPATHC = export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:./lib"

CHANGE = -stdlib=libc++


main: as2.o classes.o
	$(CC) $(INSTALLPATH) $(LINKPATH) $(LDFLAG) $(CHANGE) -o as2 as2.o classes.o

as2.o: as2.cpp classes.cpp
	$(DYLDPATHC)
	$(CC) $(INSTALLPATH) $(CHANGE)  -c as2.cpp -c classes.cpp

test: test.o classes.o
	$(CC) $(INSTALLPATH) $(LINKPATH) $(LDFLAG) $(CHANGE) -o test test.o classes.o

test.o: test.cpp
	$(DYLDPATHC)
	$(CC) $(INSTALLPATH) $(CHANGE) -c test.cpp

clean: 
	rm -rf *o as2 test

fulltest:
	GLOG_logtostderr=1 ./as2 -unittest=1

unittest:
	GLOG_logtostderr=1 ./test

all:
	make
	make test