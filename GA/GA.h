#pragma once

/*
	GA will maximise the phenotype function


*/

#define _GA_VERSION_ "2.1.2013.6.19.13.00"

#include "GenePop.h"
#include "Phenotype.h"
#include <algorithm>//sort
#include <vector>//moving individuals about

#define MAX_double 1e38
#define MIN_SUM_ABV_MIN 1e-10           //limit phenotype resolution so that rounding errors don't change pop size!
//investigate improving this
/*
Bitmem
v2.0.2013.5.29.21:57    //removed tight packing genes
v2.1.2013.5.31.21:28    //working
v2.2.2013.5.31.21:49    //working - removed Bitmem (all in GenePop class)
V2.3.2013.6.13.21.55
V3.0.2013.6.19.12.58

NOTES:
In reality haploid m/f join. Cross overs + muts occur in parent only!
Adapted here so that cross overs occur between haploid parents.
*/

/*
== DOCS ==

evolve	main loop
	calcFitness
	marry()
	breed
	return top couple

*/

using namespace std;


//Phenotype declaration
//template<typename T> double phenotype(T* pop, int ind, int genecount);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Genetic Algorithm Class
//////////////////////////////////////////////////////////////////////////////////////////////////////////

//use with couples (16bytes)
struct triple {
	int first;      //parent1
	int second;     //parent2
	double third;   //combined fitness
};

template<typename T>
class GA {
	Phenotype<T>* phenotype;
	GenePop<T> pop;
	//12bytes
	pair<int, double>* fitness; //collect individual fitnesses 1st 8Bytes fitness, 2nd 4B num
	int numcouples;
	triple* couples;            //marriage lists

	void init() {//init rest of variables
		double memMB = (pop.pop_size / 2 + 28 * pop.pop_size) / 1048576;

		cout << "************************************************" << endl;
		cout << "*                                              *" << endl;
		cout << "* GENETIC ALGORITHM                            *" << endl;
		cout << "* Vers. " << _GA_VERSION_ << "                    *" << endl;
		cout << "*                                              *" << endl;
		cout << "************************************************" << endl;
		cout << "Memory footprint: " << to_string(memMB) << "MB" << endl;

		fitness = new pair<int, double>[pop.pop_size];
		numcouples = static_cast<const int> (ceil(pop.pop_size / 2));
		couples = new triple[numcouples];
	}

	bool calcFitness(double& pop_min_f, double& f_unit) {	//
		//returns false if all same
		double sum = 0;
		double mn = MAX_double;
		double f;
		double popsize = static_cast<double>(pop.pop_size);

		for (int i = 0; i < pop.pop_size; i++)
		{
			////////////////////////////////////////////////////////
			// fitness run here
			f = phenotype->calc(pop, i);
			//
			////////////////////////////////////////////////////////

			fitness[i].first = i;   //id pop posn
			fitness[i].second = f;
			//calc pop fitness params
			sum += f;                   //sum
			if (f < mn)                 //min
			{
				mn = f;
				pop_min_f = f;
			}
		}
		//this is num children you get per couples fitness above min
		//can be near zero so handle!
		double denom = (sum - popsize * pop_min_f);

		if (denom < MIN_SUM_ABV_MIN) return false;

		f_unit = popsize / denom;

		return true;
	}

	void sortFitness() {
		std::sort(fitness, fitness + pop.pop_size, GA::paircmp);
	}

	void marry() {
		//sort fitness first
		int i2;
		for (int i = 0; i<numcouples; i++)   //marry sorted neighbours 0...half, 1...half+1
		{
			i2 = 2 * i;
			couples[i].first = fitness[i2].first;
			couples[i].second = fitness[i2 + 1].first;
			couples[i].third = fitness[i2].second + fitness[i2 + 1].second;
		}
	}

	void breed(double pop_min_f, double f_unit, bool not_same) {
		//Get positions of all cross overs
		pop.populateCrossOverMap();

		int child_count = 0;        //Num children so far
		double total_children = 0;  //Total children to aim for
		int num_children = 0;       //rounded Total children (don't round before else rounding artefacts)
		for (int i = 0; i < numcouples; i++) {
			int parent1 = couples[i].first;
			int parent2 = couples[i].second;
			double sum_fitness = couples[i].third;
			double couple_f_abv_min = sum_fitness - 2 * pop_min_f;
			//fitnesses too similar
			if (not_same)
				total_children += couple_f_abv_min * f_unit;
			else
				total_children += 2.0;
			/*
			add before rounding to avoid rounding artifacts
			instead of looping num_children which would amplify the error.
			*/
			num_children = floor(total_children);   //total_children -> popsize
			while (child_count < (const int)num_children)
			{
				pop.child(child_count++, parent1, parent2);
			}
		}
		//Children made, swap populations around
		pop.swapPopulations();
		//irradiate the main population cause mutations
		pop.irradiate();
	}

public:

	static bool paircmp(pair<int, double> a, pair<int, double> b) {
		return a.second > b.second;
	}

	GA(string archive) : pop(archive) {
		init();
	}

	GA(Phenotype<T>* phenotype, int popsize) : phenotype(phenotype), pop(popsize, phenotype->genecount) {
		init();
		cout << pop.pretties();
	}

	~GA() {
		delete[] couples;
		delete[] fitness;
	}

	//funcs control underlying population parameters
	void setProbCO(double _prob_cross) { pop.setProbCO(_prob_cross); }
	void setProbMut(double _prob_mut) { pop.setProbMut(_prob_mut); }
	void resizePop(int newSize) { pop.resizePop(newSize); }

	void serialise(string archive_name) {
		string archive_header = phenotype->SIGNATURE;
		pop.serialise(archive_name, archive_header);
	}

	triple evolve() {
		double pop_min_f = 0;   //sum, min
		double f_unit = 0;      //converts fitness abv min -> children
		bool not_same = calcFitness(pop_min_f, f_unit);		//TODO: ignoring all same
		sortFitness();
		marry();
		pop.gen++;
		breed(pop_min_f, f_unit, not_same);		
		return couples[0];
		//Swapped at end of breed... move into evolve?
	}

	double calc(int ind) { 
		return phenotype->calc(pop, ind);
	}

	//TODO: check Casting big to avoid using T
	int get(int ind, int gene) {
		return pop.get(ind, gene);
	}

	int gen() {
		return pop.gen;
	}

	int getGeneCount() {
		return pop.gene_count;
	}

	enum POSITION {TOP,BOTTOM};

	bool copyDNA(vector<T> dna, POSITION pos) {
		if (pos == TOP) {
			return pop.copyDNA(dna, 0);
		}
	}

	void replaceDNA(vector<T> dna, POSITION pos) {
		if (pos == BOTTOM) {
			pop.replaceDNA(dna, pop.gene_count);
		}
	}

};