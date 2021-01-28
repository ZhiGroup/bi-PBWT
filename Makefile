all: PBWT rPBWT

MemoryMap:
	g++ -std=c++17 -Wshadow -Wall -o PBWT PBWT.cpp -O2 -Wno-unused-result -DMEMORY_MAP -pthread -lboost_iostreams
	g++ -std=c++17 -Wshadow -Wall -o rPBWT rPBWT.cpp -O2 -Wno-unused-result -DMEMORY_MAP -pthread -lboost_iostreams

PBWT: PBWT.cpp
	g++ -std=c++17 -Wshadow -Wall -o PBWT PBWT.cpp -O2 -Wno-unused-result

rPBWT: rPBWT.cpp
	g++ -std=c++17 -Wshadow -Wall -o rPBWT rPBWT.cpp -O2 -Wno-unused-result

clean:
	rm -f PBWT rPBWT dPBWT drPBWT

.PHONY: all clean MemoryMap
