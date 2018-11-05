CC  = g++

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0 -lEve
CCFLAGS = `root-config --cflags` -g -c -Wall -Wextra -O0 -I./include/

TMP_DIR = tmp

all: plotDeDx display getFfactor scanCuts singleCutDetails

plotDeDx: ${TMP_DIR}/plotDeDx.o ${TMP_DIR}/Event.o ${TMP_DIR}/EventCut.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/JetCut.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o
	$(CC) $^ -o $@ $(LDFLAGS)

display: ${TMP_DIR}/display.o ${TMP_DIR}/Event.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o
	$(CC) $^ -o $@ $(LDFLAGS)

getFfactor: ${TMP_DIR}/getFfactor.o ${TMP_DIR}/Event.o ${TMP_DIR}/EventCut.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/JetCut.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o
	$(CC) $^ -o $@ $(LDFLAGS)

scanCuts: ${TMP_DIR}/scanCuts.o ${TMP_DIR}/Event.o ${TMP_DIR}/EventCut.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/JetCut.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o
	$(CC) $^ -o $@ $(LDFLAGS)

singleCutDetails: ${TMP_DIR}/singleCutDetails.o ${TMP_DIR}/Event.o ${TMP_DIR}/EventCut.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/JetCut.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o
	$(CC) $^ -o $@ $(LDFLAGS)

${TMP_DIR}/plotDeDx.o: plotDeDx.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/display.o: display.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/getFfactor.o: getFfactor.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/scanCuts.o: scanCuts.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/singleCutDetails.o: singleCutDetails.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/%.o: src/%.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f ${TMP_DIR}/*.o plotDeDx display getFfactor scanCuts singleCutDetails
