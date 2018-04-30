#pragma once

/*
	Innject the model to optimise
*/

template<typename T> class Phenotype {

public:
	string SIGNATURE = "BASE";

	int genecount;
	int genesize;

	Phenotype() {
		//remember to init!
	}

	Phenotype(int genecount) {
		init(genecount);
	}

	virtual double calc(GenePop<T>& pop, int ind) {
		return 1;
	}

	void init(int genecount) {
		this->genecount = genecount;
		this->genesize = sizeof(T) * 8;
	}
};
