bds: BuildDataStructure.cpp formatutil.cpp TestUtil.cpp FDBG.cpp
	g++ -O3  -std=c++11 BuildDataStructure.cpp -o bds -l boost_log -l pthread
format: test_format.cpp
	g++ -O3  -std=c++11 test_format.cpp -o tf -l boost_log -l pthread -l bf
test: BuildDataStructure.cpp formatutil.cpp TestUtil.cpp FDBG.cpp
	g++ -O0 -g  -std=c++11 BuildDataStructure.cpp -o bds_test -l boost_log -l pthread
readsplit: ReadSplit.cpp
	g++ -O3 -std=c++11 ReadSplit.cpp -o readsplit -l boost_log -l pthread
countreads: CountReads.cpp
	g++ -O3 -std=c++11 CountReads.cpp -o countreads -l boost_log -l pthread
clean:
	rm bds *.bin
