#include "GA.h"
#include <iostream>
#include <cstdint>	//int16_t etc

//Testing

template<typename T> class TestPhenotype: public Phenotype<T> {
public:
	using Phenotype::Phenotype;
	double calc(GenePop<T>& pop, int ind) {
		double val = (double)pop.get(ind, 0);
		return 96.6 - pow(val / 10 - 24.5, 2);
	}
};

int main()
{
	//TODO: is passing pointer to object on the stack a bad idea?

	TestPhenotype<uint16_t>* tp = new TestPhenotype<uint16_t>(1);
	GA<uint16_t> a = GA<uint16_t>(tp, 10000);

	triple top_pair = { 0, 0, 0.0 };

	for (int i = 0; i<500000; i++) {
		top_pair = a.evolve();
		//std::cout << top_pair.first << "," << top_pair.second << "," << top_pair.third << "," << endl;
		cout << a.calc(top_pair.first) << "(" << a.get(top_pair.first,0) << ")" << endl;
	}

	delete tp;
	

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
