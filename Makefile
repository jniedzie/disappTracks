CC  = g++

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0 -lEve -lGeom
CCFLAGS = `root-config --cflags` -g -c -Wall -Wextra -O0 -I./include/ --std=c++17

TMP_DIR = tmp

all: runAnalysis eventDisplay helixFitter getFfactor scanCuts singleCutDetails

runAnalysis: ${TMP_DIR}/runAnalysis.o ${TMP_DIR}/Event.o ${TMP_DIR}/EventSet.o ${TMP_DIR}/EventCut.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/JetCut.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o
	$(CC) $^ -o $@ $(LDFLAGS)

eventDisplay: ${TMP_DIR}/eventDisplay.o ${TMP_DIR}/Event.o ${TMP_DIR}/EventSet.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o ${TMP_DIR}/Display.o ${TMP_DIR}/Fitter.o ${TMP_DIR}/Helix.o ${TMP_DIR}/Circle.o ${TMP_DIR}/Point.o ${TMP_DIR}/FitterConfig.o ${TMP_DIR}/PointsProcessor.o
	$(CC) $^ -o $@ $(LDFLAGS)

helixFitter: ${TMP_DIR}/helixFitter.o ${TMP_DIR}/Event.o ${TMP_DIR}/EventSet.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o ${TMP_DIR}/Display.o ${TMP_DIR}/Fitter.o ${TMP_DIR}/Helix.o ${TMP_DIR}/Circle.o ${TMP_DIR}/Point.o ${TMP_DIR}/FitterConfig.o ${TMP_DIR}/PointsProcessor.o ${TMP_DIR}/MonitorsManager.o
	$(CC) $^ -o $@ $(LDFLAGS)

getFfactor: ${TMP_DIR}/getFfactor.o ${TMP_DIR}/Event.o ${TMP_DIR}/EventSet.o ${TMP_DIR}/EventCut.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/JetCut.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o
	$(CC) $^ -o $@ $(LDFLAGS)

scanCuts: ${TMP_DIR}/scanCuts.o ${TMP_DIR}/Event.o ${TMP_DIR}/EventSet.o ${TMP_DIR}/EventCut.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/JetCut.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o
	$(CC) $^ -o $@ $(LDFLAGS)

singleCutDetails: ${TMP_DIR}/singleCutDetails.o ${TMP_DIR}/Event.o ${TMP_DIR}/EventSet.o ${TMP_DIR}/EventCut.o ${TMP_DIR}/Track.o ${TMP_DIR}/TrackCut.o ${TMP_DIR}/Jet.o ${TMP_DIR}/JetCut.o ${TMP_DIR}/HistSet.o ${TMP_DIR}/Lepton.o ${TMP_DIR}/LeptonCut.o
	$(CC) $^ -o $@ $(LDFLAGS)

${TMP_DIR}/runAnalysis.o: runAnalysis.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/eventDisplay.o: eventDisplay.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/helixFitter.o: helixFitter.cpp
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
	rm -f ${TMP_DIR}/*.o runAnalysis eventDisplay helixFitter getFfactor scanCuts singleCutDetails
