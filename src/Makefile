CC = g++
SOURCE = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCE))
EXECUTABLE = lbm
CFLAGS =-O3 -Wno-deprecated -Wall

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXECUTABLE) 
	@echo DONE!

-include $(OBJECTS:.o=.d)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $*.cpp -o $*.o
	@$(CC) -MM $(CFLAGS) $*.cpp > $*.d
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean:
	rm -rf *.o *.d *.dat $(EXECUTABLE) *~
	@echo ALL CLEANED UP!

rebuild:
	make clean
	make
