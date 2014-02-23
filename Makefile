CC = g++

MAGIC = `bin/GraphicsMagick++-config --cppflags --cxxflags --ldflags --libs`

INSTALLPATH = -I./include -I/opt/X11/include

LINKPATH = -L./lib -L/opt/X11/lib

LDFLAG = -lgflags -lglog -lX11 -lgtest

DYLDPATHC = export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:./lib"

UNAME = :$(shell, name)

main: as2.o
	$(CC) $(INSTALLPATH) $(LINKPATH) $(LDFLAG) -o as2 as2.o


as2.o: as2.cpp classes.cpp
	$(DYLDPATHC)
	$(CC) $(INSTALLPATH)  -c as2.cpp -c classes.cpp

test: test.o
	$(CC) $(INSTALLPATH) $(LINKPATH) $(LDFLAG) -o test test.o

test.o: test.cpp
	$(DYLDPATHC)
	$(CC) $(INSTALLPATH) -c test.cpp

clean: 
	rm -rf *o as2 test

fulltest:
	GLOG_logtostderr=1 ./as2 -unittest=1

unittest:
	./as2 -unittest=1