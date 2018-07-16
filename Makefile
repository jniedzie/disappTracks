CC  = g++

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0
CCFLAGS = `root-config --cflags` -g -c -std=c++1y -Wall -Wextra -O0 -I./include/

TMP_DIR = tmp

all: plotDeDx

plotDeDx: ${TMP_DIR}/plotDeDx.o ${TMP_DIR}/Event.o ${TMP_DIR}/Track.o
	$(CC) $^ -o $@ $(LDFLAGS)

${TMP_DIR}/plotDeDx.o: plotDeDx.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/%.o: src/%.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f ${TMP_DIR}/*.o plotDeDx
