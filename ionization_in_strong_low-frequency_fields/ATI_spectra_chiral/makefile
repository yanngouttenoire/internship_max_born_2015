CC=g++ 
CFLAGS=-ansi -pedantic -fopenmp 
WFLAGS=-W -Wall 
LDFLAGS=-fopenmp
EXEC=exec.out
SOURCE= $(wildcard *.cpp)
OBJET= $(SOURCE:.cpp=.o)
HEADER= $(wildcard *.h)

run: exec.out

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

		
