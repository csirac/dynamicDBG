LIBS=-l boost_system -l boost_log -l pthread
bds: BuildDataStructure.cpp formatutil.cpp TestUtil.cpp FDBG.cpp
	g++ -O3  -std=c++11 BuildDataStructure.cpp -o bds $(LIBS)
format: test_format.cpp
	g++ -O3  -std=c++11 test_format.cpp -o tf $(LIBS) -l bf
test: BuildDataStructure.cpp formatutil.cpp TestUtil.cpp FDBG.cpp
	g++ -O0 -g  -std=c++11 BuildDataStructure.cpp -o bds_test $(LIBS)
readsplit: ReadSplit.cpp
	g++ -O3 -std=c++11 ReadSplit.cpp -o readsplit $(LIBS)
countreads: CountReads.cpp
	g++ -O3 -std=c++11 CountReads.cpp -o countreads $(LIBS)
clean:
	rm bds *.bin
