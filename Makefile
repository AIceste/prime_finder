findprimes: prime_finder.o main.o
	g++ -O3 prime_finder.o main.o -o findprimes -lpthread -lgmp -lgmpxx

main.o: main.cpp
	g++ -c -O3 -Wall -Wextra -Wpedantic -std=c++11 main.cpp

prime_finder.o: prime_finder.cpp
	g++ -c -O3 -Wall -Wextra -Wpedantic -std=c++11 prime_finder.cpp

clean:
	rm -rf *.o findprimes
