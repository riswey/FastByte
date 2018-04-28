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

			pop.setArchiveFile("test.txt");
				
			pop.serialise();

		}

		TEST_METHOD(TestBreeding)
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

		TEST_METHOD(GA)
		{
		}


	};
}