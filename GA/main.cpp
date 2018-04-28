#include "GA.h"
#include <iostream>
#include <cstdint>	//int16_t etc

//Testing
template<typename T> class TestGA : public GA<T> {
public:
	TestGA(int popsize, int genecount, int genesize) : GA<T>(popsize, genecount, genesize) {}

	double phenotype(int ind) { //individual
		double val = (double)GA<T>::pop.get(ind, 0);

		return 96.6 - pow(val/10 - 24.5, 2);
	}

};

TestGA<uint16_t> a = TestGA<uint16_t >(10000, 1, 16);

int main()
{

	triple top_pair = { 0, 0, 0.0 };

	for (int i = 0; i<50000; i++) {
		top_pair = a.evolve();
		//std::cout << top_pair.first << "," << top_pair.second << "," << top_pair.third << "," << endl;
		cout << a.phenotype(top_pair.first) << "(" << a.pop.get(top_pair.first,0) << ")" << endl;
	}



	/*
	GenePop<int16_t> g;
	GA<int16_t unsigned> a(10,10,8);
	g.set_all(38);
	g.serialise();
	g.deserialiseInit();

	GenePop<int8_t> g1("GA1.pop");
	g1.archive_file = "GA1.pop";
	g1.serialise();
	*/

	getchar();
	return 0;
}
