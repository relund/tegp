all: tegp

LIBXMLPATH = /usr/include/libxml2/

random.o: ../code/random.cpp ../code/random.hpp
	g++ -O -c ../code/random.cpp

tegp: ../code/teg-peak.cpp ../code/teg-peak.hpp random.o
	g++ random.o ../code/teg-peak.cpp -O -o tegp -I $(LIBXMLPATH) -lxml2
