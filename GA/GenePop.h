#pragma once
/*

Holds population and encapsulates gene ops

*/
#define _GP_VERSION_ "30.4.2018.0.0.0.1"		//resurrected. Removed bit resolution. Only evolves bytes now
//#define _GP_VERSION_ "3.1.2013.10.17.18.30"   //serialise bug fixed

#include <map>
#include <time.h>   //to seed rand
#include <fstream>
#include "gsl/gsl_randist.h"
//#include "gsl/gsl_cdf.h"
#include "myfilefuncs.h"

using namespace std;

struct GAException : public exception
{
	string s;
	GAException(string s_) : s(s_) {}
	~GAException() throw () {}
	const char* what() const throw() { return s.c_str(); }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gene/Population Management Class
//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T> class GenePop {
	//Population parameters
	double prob_cross;
	double prob_mut;
	string archive_file;    //stores the serialised pop
	string archive_hdr;     //info about pop training (env files/stamps etc).
							//derivative parameters
	int num_units;      //num storage units (n)
	int sizeofT;        //bytes in T
	int sizeofT8;       //bits in T
	int population_bits;       //Total bits for global functions
	float RANDEXPDISTSCALE;         //multiple for randExpDist
	float RANDEXPDISTCORRECTION;    //corrects distribution to unit probability
	T* units;           //the mem
	T* units2;          //mem for children
	T* bmask;           //2^x ( to get bit 0,1,2,3: & with bmask[0,1,2,3] ) :bmask[5] = 100000 : used (i=0;i++<sizeof(T);)
	T* mask;            //2^x - 1 ( to limit to x bits)                     :mask[5] = 11111 : used to mask sizeof(T)
	/*
		CROSS-OVERS
		
	
	*/
	map<int, T> comap;  //map for filing gene cross overs
	typename map<int, T>::iterator it_comap;

	void init() {//once pop_size, gene_count, prob_cross, prob_mut loaded then init
		sizeofT = sizeof(T);
		sizeofT8 = sizeofT * 8;
			initMem();
		//set dynamic masks for bit manipulation
		//mask to extract bits
		bmask = new T[sizeofT8+1];		//0000 0001,0000 0010,0000 0100,..,1000 0000,0000 0000 (i=T then = 0)
		//mask to remove left bits
		mask = new T[sizeofT8+1];		//0000 0000,0000 0001,0000 0011,..,0111 1111,1111 1111

		for (int i = 0; i<sizeofT8+1; i++)
		{
			bmask[i] = pow(2, i);           
			mask[i] = bmask[i] - 1;
		}

		//Get ptr to random number generator
		GenePop<T>::r = gsl_rng_alloc(gsl_rng_taus);
		//seed generator
		gsl_rng_set(r, time(NULL));
	}
	void memParameters()
	{
		num_units = pop_size * gene_count;
		population_bits = num_units * sizeofT8;
	}
	void initMem()  //separated so can call separately when mem resized
	{
		memParameters();
		//Create units
		units = new T[num_units];           //create memory
		units2 = new T[num_units];           //create child memory
	}

	T* delete_old_mem()
	{
		//Returns a *copy old memory
		T* units3 = new T[num_units];
		//Copy old mem into units3 storage
		memcpy(units3, units, num_units);
		//delete existing memory
		delete[] units;
		delete[] units2;
		return units3;
	}

	T getRandGene()
	{
		//Simply find a random number between 0 and 111...111 (size of T)
		static unsigned long int geneinfosize = (unsigned long int)(mask[sizeofT8]);
		return gsl_rng_uniform_int(r, geneinfosize);  //limit size (bits)
	}

	void set_rand() {       //byte by byte
		for (int i = 0; i<num_units; i++)
			units[i] = getRandGene();
	}

	/////////////////////////////////////

	void deserialiseInit() {
		//relative to current directory
		cout << "\n\nInitialising & Loading archive...\n";
		cout << archive_file << "\n";

		ifstream fp(archive_file.c_str(), ios::binary);
		if (!fp) {
			throw GAException("Invalid file.");
		}
		//Environment files/ id stamps
		fp >> archive_hdr;
		cout << "\nHeader : " << archive_hdr << endl;

		int unitsize;
		fp >> unitsize;
		if (sizeof(T) < unitsize)
		{
			throw GAException("Non matching unitsizes! " + patch::to_string(sizeofT) + "!=" + patch::to_string(unitsize));
		}
		fp >> gen;
		fp >> pop_size;
		fp >> gene_count;
		fp >> sizeofT8;
		fp >> prob_cross;
		fp >> prob_mut;

		if (sizeofT8 != sizeof(T) * 8) {
			//TODO: actually this can just fit to new size with a warning
			throw GAException("Failed to initialize. Popstore genesize doesn't fit object T size.")
		}

		init();         //initialise rest of class

						//strip the final eol as binary read pointer moves to end of line?
						//fp.ignore() dont use as wish to see which character is ignored
		char eol;
		fp.get(eol);
		cout << "stripped: ASCII eol: " << (int)eol << ". Deserialising Binary...";

		// deserialise binary file
		char buffer[unitsize];
		for (int i = 0; i<num_units; i++) {
			fp.read(buffer, unitsize);
			units[i] = *((T*)buffer);
		}

		fp.close();
		cout << "done.\n\n";
	}

	//gets the index of the gene in the memory
	inline int gene_pos(int ind, int gene) {
		return gene_count * ind + gene;
	}

	void set(int ind, int gene, int value) {
		int offset = gene_pos(ind, gene);
		units[offset] = static_cast<T>(value);
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///
	///
	///		PUBLIC
	///
	///
	///////////////////////////////////////////////////////////////////////////////////////////

public:
	static gsl_rng* r;         //ptr random number generator

	int gen;            //generation
	int pop_size;       //num individuals
	int gene_count;     //num genes per individual

	//Init from archive
	GenePop(string archive) :archive_file(archive) {
		deserialiseInit();
		pretties();
	}

	//Create new population
	GenePop(int _pop_size = 10, int _gene_count = 1, double _prob_cross = 0.002, double _prob_mut = 0.00002, string af_ = "GA.pop", string ah_ = "") :
		gen(0), pop_size(_pop_size), gene_count(_gene_count), prob_cross(_prob_cross), prob_mut(_prob_mut), archive_file(af_), archive_hdr(ah_)
	{
		init();
		set_rand();
		pretties();
	}

	~GenePop() {
		//Release resources used by random number generator
		gsl_rng_free(r);
		//Clear memory
		delete[] units;
		delete[] units2;
		delete[] bmask;
		delete[] mask;
		cout << "Memory freed." << endl;
	}

	//Set the population parameters
	void setProbCO(double _prob_cross) { prob_cross = _prob_cross; cout << "Prob cross: " << prob_cross << endl; }
	void setProbMut(double _prob_mut) { prob_mut = _prob_mut; cout << "Prob mut: " << prob_mut << endl; }
	void setArchiveFile(string af) { archive_file = af; }
	void setArchiveHdr(string hdr) { archive_hdr = hdr; cout << "Header: " << archive_hdr << endl; }

	//Set population to fixed value
	void set_all(T val) {       //byte by byte
		for (int i = 0; i<num_units; i++)
			units[i] = static_cast<T>(val);
	}

	void setChild(int ind, int gene, int value) {
		int offset = gene_pos(ind, gene);
		units2[offset] = static_cast<T>(value);
	}

	//Return a gene
	T get(int ind, int gene) {
		int offset = gene_pos(ind, gene);
		return static_cast<T>(units[offset]);
	}

	void serialise() {  //allow upto 8 bytes read in
						//relative to current directory
		ofstream fp(archive_file.c_str(), ios::binary);
		if (archive_hdr.length() == 0) archive_hdr = "(none)";   //ensure not empty so that >> can read back
		fp << archive_hdr << endl;
		fp << sizeofT << endl;
		fp << gen << endl;
		fp << pop_size << endl;
		fp << gene_count << endl;
		fp << sizeofT8 << endl;
		fp << prob_cross << endl;
		fp << prob_mut << endl;             //seems to add too many 0As

		for (int i = 0; i<num_units; i++)
		{
			fp.write((char*)&(units[i]), sizeofT);
		}
		fp.close();

	}

	void irradiate() {
		/*
		Using Binomial to find out how many mutations in whole population.
		Then, evenly distributing over dna.
		Essentially a dna_bits size sequence of n sided dice throws.

		Is it for binomial generator to work with smaller n? To `
		*/

		int num_muts = static_cast<int>(gsl_ran_binomial(r, prob_mut, population_bits));
		//Irradiate entire population
		int c = 0;
		int gene = 0;
		int bit = 0;

		//No collisions (ignoring for now)
		while (c++ < num_muts)
		{
			//Get the gene and bit position of the next mutated bit
			gene = gsl_rng_uniform_int(r, num_units);          //0 <= x < num_units
			bit = gsl_rng_uniform_int(r, sizeofT8);           //This is the actual limit to number size not store!!!
															   //If T = bit then = 0
			if (bmask[bit]) units[gene] ^= bmask[bit];
		}
	}

	void swapPopulations() {
		//make children then new pop
		T* temp = units;
		units = units2;
		units2 = temp;
	}

	///////////////////////////////////
	// CROSS_OVERS
	///////////////////////////////////

	//CROSS OVER MAP FUNCTIONS
	void populateCrossOverMap() {
		//Foretell crossovers masks for all children in advance of breeding!
		comap.clear();
		int num_co = static_cast<int>(gsl_ran_binomial(r, prob_cross, population_bits - 1));
		//crosses occur at boundaries so n-1 boundaries
		//look into co masking find out whether before boundary or after

		int gene = 0;
		while (num_co-- != 0)
		{
			gene = gsl_rng_uniform_int(r, num_units);        //0 <= x < num_units
															 //XOR with 00011111 type mask
															 // causes 0 masked bits stay same & 1 masked bits to flip
															 //comap 11111111(start swith co so all flipped)...1 (co b4 last one)
															 // != 0 since decided ther is a co here
			
			if (comap[gene] == 0) //new crossover in gene
				comap[gene] = mask[1 + gsl_rng_uniform_int(r, sizeofT8 - 1)];       //not 0 nor gene_size since must be a crossover! (note: mask[gene_size] -> all bits set)
			else
				comap[gene] ^= mask[1 + gsl_rng_uniform_int(r, sizeofT8 - 1)];       //not 0 nor gene_size since must be a crossover!
		}

		//Init iterator for comaps
		it_comap = comap.begin();
	}

	//cycle thorough crossover points map 
	T get_mask_inc(int ind, int gene)
	{
		//TODO: should maintain state of last byte!

		if (gene_pos(ind, gene) < it_comap->first)
			return 0;
		else
		{
			//testing for current mask, grab it
			T mask = it_comap->second;
			//increment iterator;
			it_comap++;
			//return it
			return mask;
		}
	}
	//temporary for testing (delete)
	void resetCOInterator()
	{
		it_comap = comap.begin();
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Resize Model
	///////////////////////////////////////////////////////////////////////////////////////////////

	void resizePop(int newSize)     //change population size
	{
		if (newSize == pop_size) return;
		//Keep old pop size
		int oldSize = pop_size;
		//Resize memory+init memory and get *old memory
		T* units3 = delete_old_mem();
		//Set new popsize
		pop_size = newSize;
		//recalc and allocate memory
		initMem();
		//Copy old units into new memory repeating as needed (using new num_units)
		for (int i = 0; i != num_units; i++)
		{
			units[i] = units3[i % oldSize];
		}
		delete[] units3;
		cout << "Pop resize: " << pop_size << endl;
	}

	void resizeGeneCount(int gc_)   //change DNA size
	{
		if (gc_ == gene_count) return;
		//Keep old gene_count
		int oldGC = gene_count;
		//Set new gene count
		gene_count = gc_;
		//Delete old mem + return *old memory
		T* units3 = delete_old_mem();
		//recalc memory parameters and allocate memory
		initMem();
		//Copy old individuals into new individuals, looping if necessary
		for (int i = 0; i != pop_size; i++)
			for (int j = 0; j != gene_count; j++)
			{
				units[i*gene_count + j] = units3[i * oldGC + j % oldGC];
			}
		delete[] units3;
		cout << "Gene count resize: " << gene_count << endl;
	}

	void resizeGeneSize(int gs_)    //change number bits
	{
		//Mask data to new gene_size
		if (gs_ < sizeofT8)
		{
			sizeofT8 = gs_;
			//Mask genes to new gene_size
			for (int i = 0; i != num_units; ++i)
			{
				units[i] &= mask[sizeofT8];
			}
		}
		//Set population bits for irradiation
		memParameters();
	}
	/*
	//Init from command line
	void cmd(int argc, char* argv[])
	{
	vector< vector<string> > q;
	for(int i=1;i<argc;i++)
	{
	q.push_back( strSplit(argv[i], "=") );
	}

	for_each(
	q.begin(),
	q.end(),
	[](vector<string> vs)
	{
	switch (vs[0][1])
	{
	case 97:    //a archive
	setArchiveFile(vs[1]);
	break;
	case 99:    //c crossovers prob
	setProbCO(vs[1]);
	break;
	case 101:    //e environment
	setEnvironmentFile(vs[1]);
	break;
	case 103:    //g environment
	setGeneSize(vs[1]);
	break;
	case 109:    //m mutation prob
	setProbMut(vs[1]);
	break;
	case 112:    //p population
	resizePop(vs[1]);
	break;
	default:
	cout <<
	"\t-a=archive_file\n\t-e=environment_file\n" <<
	"\t-m=prob_mutation\n\t-c=prob_crossover\n" <<
	"\t-g=gene_size\n\t-p=population size\n";
	break;
	}
	}
	);
	}
	*/

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Reports
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	string memout() {
		string str = "";
		for (int i = 0; i<num_units; i++)
			str += displayByte(units[i], i);
		return str;
	}

	string crossoversout() {
		string str = "CrossOver Map (" + to_string(comap.size()) + ")\n";
		for (it_comap = comap.begin(); it_comap != comap.end(); ++it_comap) {
			str += displayByte(it_comap->second, it_comap->first);
		}
		return str;
	}

	inline string displayByte(T byte, int i = 0) {
		return to_string(i) + "\t" + byteToString(byte) + "\t(" + to_string((int)byte) + ")\n";
	}

	string pretties()
	{
		string str = "";
		str += "************************************************\n";
		str += "*                                              *\n";
		str += "* GENEPOP                                      *\n";
		str += "* Vers. " + string(_GP_VERSION_) + "                   *\n";
		str += "*                                              *\n";
		str += "************************************************\n\n";
		str += "Pop size:\t" + to_string(pop_size) + "\nGene Count:\t" + to_string(gene_count) + "\nGene size\t" + to_string(sizeofT8) + "\nProb cross\t" + to_string(prob_cross) + "\nProb Mutate\t" + to_string(prob_mut) + "\n\n";
		return str;
	}

	static string byteToString(T byte) {
		string str = "";
		for (int i = sizeof(T) * 8 - 1; i>-1; i--)
		{
			str += (byte & (T)pow(2, i)) ? "1" : "0";
		}
		return str;
	}
};

//Declaration of ptr to random number generator
template<typename T> gsl_rng* GenePop<T>::r;