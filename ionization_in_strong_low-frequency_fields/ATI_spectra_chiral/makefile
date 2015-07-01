CC=g++ 
CFLAGS=-W -Wall -ansi -pedantic -std=c++11 -fopenmp 
LDFLAGS=-fopenmp
EXEC=exec.out
SOURCE= $(wildcard *.cpp)
OBJET= $(SOURCE:.cpp=.o)
HEADER= $(wildcard *.h)

run: exec.out
	@ ./exec.out

exec.out: $(OBJET)
	@ $(CC) -o $@ $^ $(LDFLAGS)

main.o: $(HEADER)

%.o: %.cpp
	@ $(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper run 

clean:
	@ rm -rf *.o

mrproper: clean
	@ rm -rf $(EXEC)

		
