#pragma once

#define _GA_VERSION_ "2.1.2013.6.19.13.00"

#include "GenePop.h"
#include <algorithm>//sort

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

template<typename T> class Phenotype {
public:
	int genecount;
	int genesize;

	Phenotype(int genecount) {
		this->genecount = genecount;
		this->genesize = sizeof(T) * 8;
	}

	virtual double calc(GenePop<T>& pop, int ind) {
		return 1;
	}
};


//Phenotype declaration
//template<typename T> double phenotype(T* pop, int ind, int genecount);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Genetic Algorithm Class
//////////////////////////////////////////////////////////////////////////////////////////////////////////

//use with couples
struct triple {
	int first;      //parent1
	int second;     //parent2
	double third;   //combined fitness
};

template<typename T> class GA {
	Phenotype<T>* phenotype;
	GenePop<T> pop;

	pair<int, double>* fitness; //collect individual fitnesses 1st 6Bytes fitness, 2nd 2B num
	int numcouples;
	triple* couples;            //marriage lists

	void init() {//init rest of variables
		cout << "************************************************" << endl;
		cout << "*                                              *" << endl;
		cout << "* GENETIC ALGORITHM                            *" << endl;
		cout << "* Vers. " << _GA_VERSION_ << "                    *" << endl;
		cout << "*                                              *" << endl;
		cout << "************************************************" << endl << endl;

		fitness = new pair<int, double>[pop.pop_size];
		numcouples = static_cast<const int> (ceil(pop.pop_size / 2));
		couples = new triple[numcouples];
	}

	void serialise() {
		pop.serialise();
	}

	bool calcFitness(double& pop_min_f, double& f_unit) {//returns false if all same
		double sum = 0;
		double mn = MAX_double;
		double f;
		double popsize = static_cast<double>(pop.pop_size);

		for (int i = 0; i<pop.pop_size; i++)
		{
			f = phenotype->calc(pop, i);
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

	void marry() {
		//sort fitness<pop id, phenotype> by phenotype
		std::sort(fitness, fitness + pop.pop_size, GA::paircmp);
		int i2;
		for (int i = 0; i<numcouples; i++)   //marry sorted neighbours
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
			*/
			num_children = floor(total_children);   //total_children -> popsize
			while (child_count < (const int)num_children)
			{
				child(child_count++, parent1, parent2);
			}
		}
		//Children made, swap populations around
		pop.swapPopulations();
		//irradiate the main population cause mutations
		pop.irradiate();
	}

	inline T compl(T byte) {
		return (T)~byte;
	}

	//Give birth
	void child(int ind, int parent1, int parent2) {

		bool s = gsl_rng_uniform_int(pop.r, 2);					//random parent to start (false = parent1)
		/*
		cross-over mask
		00011100 means bits from parent 0 or 1
		then swapped by random parent swap

		*/
		T comask;
		
		for (int g = 0; g<pop.gene_count; g++)					//go thru each gene position
		{
			comask = pop.get_mask_inc(ind, g);

			//TODO need to maintain state of last bit (flip s by end bit)
			if (comask == 0)
			{
				//0 = No crossover i.e. just copy parent
				pop.setChild(
					ind,
					g,
					((s) ? pop.get(parent2, g) : pop.get(parent1, g))
				);
			}
			else
			{
				T mask1, mask2;
				if (s)
				{
					mask1 = comask;           //get mask orig. starting 0
					mask2 = compl(comask);    //complement for parent2
				}
				else
				{
					mask1 = compl(comask);    //get mask orig. starting 1
					mask2 = comask;
				}

				//get parental genes
				T g1 = pop.get(parent1, g);
				T g2 = pop.get(parent2, g);

				//Apply cross over mask
				g1 &= mask1;
				g2 &= mask2;

				//OR combine into child
				T child = g1 | g2;
				//Save child

				//pop.setChild(ind, g, pop.get(parent1,g));
				pop.setChild(ind, g, child);

				//which parent uppermost after crossovers
				//Look at last bit (should match start next mask)
				//s = ((mask1 & pop.bmask[pop.sizeofT8-1]) == pop.bmask[pop.sizeofT8-1]);
				s = (mask1 & 1);
			}
		}
	}

public:

	static bool paircmp(pair<int, double> a, pair<int, double> b) {
		return a.second > b.second;
	}

	GA(string archive) : pop(archive) {
		init();
	}

	GA(Phenotype<T>* phenotype, int popsize) : pop(popsize, phenotype->genecount) {
		this->phenotype = phenotype;
		init();
	}

	~GA() {
		delete[] couples;
		delete[] fitness;
	}

	//funcs control underlying population parameters
	void setProbCO(double _prob_cross) { pop.setProbCO(_prob_cross); }
	void setProbMut(double _prob_mut) { pop.setProbMut(_prob_mut); }
	void setGeneSize(int gs_) { pop.setGeneSize(gs_); }
	void setArchiveFile(string af_) { pop.setArchiveFile(af_); }
	void setArchiveHdr(string hdr) { pop.setArchiveHdr(hdr); }
	void resizePop(int newSize) { pop.resizePop(newSize); }

	triple evolve() {
		double pop_min_f = 0;   //sum, min
		double f_unit = 0;      //converts fitness abv min -> children
		bool not_same = calcFitness(pop_min_f, f_unit);
		marry();
		pop.gen++;
		breed(pop_min_f, f_unit, not_same);
		return couples[0];
		//Swapped at end of breed... move into evolve?
	}

	double calc(int ind) { 
		return phenotype->calc(pop, ind);
	}

	T get(int ind, int gene) {
		return pop.get(ind, gene);
	}
};