#include "CppUnitTest.h"
#include "CppUnitTestLogger.h"

#include "windows.h"
#include "../GA/GenePop.h"
#include <iostream>
#include <string>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace GenePopTest
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		Assert a;
		//GenePop(int _pop_size = 10, int _gene_count = 1, int _gene_size = 8, double _prob_cross = 0.002, double _prob_mut = 0.00002, string af_ = "GA.pop", string ah_ = "") :

		TEST_METHOD(TestMethod1)
		{
			Logger::WriteMessage("#######################################################");

			//10 individuals with 2x 8bit genes
			GenePop<uint8_t> pop(10,2,8);

			a.IsTrue(pop.pop_size == 10);
			a.IsTrue(pop.gene_count == 2);

			a.IsTrue(pop.byteToString(17) == "00010001");

			Logger::WriteMessage(("\n" + pop.memout()).c_str() );

			Logger::WriteMessage("1st gene by get()\n");
			uint8_t gene1 = pop.get(0, 0);
			Logger::WriteMessage(to_string(gene1).c_str());

			pop.swapPopulations();
			uint8_t gene2 = pop.get(0, 0);
			a.IsFalse(gene1 == gene2);

			pop.swapPopulations();
			uint8_t gene3 = pop.get(0, 0);
			a.IsTrue(gene3 == gene1);

			pop.serialise("test.txt");

		}

		TEST_METHOD(TestCrossovers)
		{
			//10 individuals with 2x 8bit genes
			GenePop<uint8_t> pop(10, 2, 8);

			Logger::WriteMessage("#######################################################");
			//Test Reproduction

			pop.setProbCO(0);

			Logger::WriteMessage("Crossover Map\n");
			pop.populateCrossOverMap();

			Logger::WriteMessage(pop.crossoversout().c_str());

			pop.resetCOInterator();
			string str = "\n";
			for (int i = 0; i < pop.pop_size; i++) {
				for (int j = 0; j < pop.gene_count; j++) {
					str += (to_string(i) + "," + to_string(j) + "\t" + pop.byteToString(pop.get_mask_inc(i, j)) + "\n");
				}
			}
			Logger::WriteMessage(str.c_str());

		}

		TEST_METHOD(TestBreeding)
		{
			int child_count = 0;        //Num children so far

			double size = 11;

			double fitness[11] = { 0.781007, 0.987095, 0.816114, 0.333005, 0.838099, 0.319576, 0.476979,
				0.384771, 0.425035, 0.263178, 0.530664 };

			double pop_min_f = 0.263178;
			double sum = 6.15552 - size * pop_min_f;
			double f_unit = size / sum;

			double total_children = -1;		//Total children couple allowed to aim for (double)
											//must get to 0 before get first child!
			for (int i = 0; i < size; i++) {
				double couple_f_abv_min = fitness[i] - pop_min_f;		//couple f above min
				total_children += couple_f_abv_min * f_unit;

				Logger::WriteMessage(std::to_string(total_children).c_str());

				while (child_count <= total_children)
				{
					Logger::WriteMessage(std::to_string(child_count).c_str());

					child_count++;
				}
			}

			a.IsTrue(child_count == size);
			
		}

	};
}