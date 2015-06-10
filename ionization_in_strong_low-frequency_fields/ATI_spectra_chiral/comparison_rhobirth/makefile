CC=g++ -std=c++0x
CFLAGS=-W -Wall -ansi -pedantic
LDFLAGS=
EXEC=exec.out
SOURCE= $(wildcard *.cpp)
OBJET= $(SOURCE:.cpp=.o)
HEADER= $(wildcard *.h)

run: exec.out
	@ ./exec.out

exec.out: $(OBJET)
	@ $(CC) -o $@ $^ $(LDFLAGS)

main.o: $(HEADER)

%.o: %.c
	@ $(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper run 

clean:
	@ rm -rf *.o

mrproper: clean
	@ rm -rf $(EXEC)

		
