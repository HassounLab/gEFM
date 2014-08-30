#include <stdio.h>
#include <stdlib.h>

#include "EFMGenerator.h"
#include "BitVector.h"
#include "Rank.h"

/**
 * Network for which EFMs are to be computed.
 */
Network network;

/**
 * Computes elementary flux modes for the loaded network. The function is
 * templated for BitVector used to store reactions.
 */
template<class BitVector>
void execute() {
   reversibleTreeIndices = NULL;

   reversibleTreeIndices = (int*) malloc(MAX_PATHWAY_POOL_SIZE * sizeof (int));
   if (reversibleTreeIndices == NULL) {
      printf("Error allocating memory for reversible tree indices\n");
      return;
   }

   Pathway<BitVector>* pathways = (Pathway<BitVector>*) malloc(numReactions * sizeof (Pathway<BitVector>));
   for (int r = 0; r < numReactions; r++) {
      reactions.push_back(network.reactions[r]);
      pathways[r] = Pathway<BitVector>(r);
      reversible.push_back(network.reversible[r]);
      for (int m = 0; m < numMetabolites; m++) {
         pathways[r].setMetaboliteCoefficient(m, network.s[m][r]);
      }
   }

   reversiblePairs = malloc(reversiblePairCount * sizeof (BitVector));
   BitVector* revPairs = (BitVector*) reversiblePairs;
   for (unsigned int r = 0; r < network.reversiblePairs.size(); r++) {
      revPairs[r] = BitVector();
      revPairs[r].setBit(network.reversiblePairs[r], true);
      revPairs[r].setBit(network.reversiblePairs[r] + 1, true);
   }

   EFMGenerator<BitVector> efmgenerator(pathways);
   efmgenerator.genenrateEFMs();
   efmgenerator.printEFMs();
   free(pathways);
}

/**
 * Frees allocated memory resources.
 */
void freeResources() {
   if (reversibleTreeIndices != NULL) {
      free(reversibleTreeIndices);
   }
   if (reversiblePairs != NULL) {
      free(reversiblePairs);
   }
}

/**
 * Computes maximum cardinality for the rank test.
 */
void computeMaxCardinality() {
   double* matrix = (double*) malloc(sizeof (double) * numMetabolites * numReactions);
   int rows = 0;
   for (int m = 0, i = 0; m < numMetabolites; m++) {
      if (!network.external[m]) {
         rows++;
         for (int r = 0; r < numReactions; r++, i++) {
            matrix[i] = network.s[m][r];
         }
      }
   }
   maxCardinality = computeRank(matrix, numReactions, rows, NULL, NULL) + 1;
   free(matrix);
}

/**
 * Computes elementary flux modes for the given network.
 * @param file Network file.
 */
void execute(const char* file) {
   if (!network.readNetworkFile(file)) {
      printf("Error loading network file\n");
      return;
   }
   numReactions = network.reactions.size();
   numMetabolitesRemaining = numMetabolites = network.metabolites.size();
   for (int m = 0; m < numMetabolites; m++) {
      metabolites.push_back(network.metabolites[m]);
      externalMetabolites.push_back(network.external[m]);
   }
   reversiblePairCount = network.reversiblePairs.size();
   computeMaxCardinality();
   if (numReactions <= 32) {
      execute<BitVector32>();
   } else if (numReactions <= 64) {
      execute<BitVector64>();
   } else if (numReactions <= 96) {
      execute<BitVector96>();
   } else if (numReactions <= 128) {
      execute<BitVector128>();
   } else if (numReactions <= 160) {
      execute<BitVector160>();
   } else if (numReactions <= 192) {
      execute<BitVector192>();
   } else if (numReactions <= 448) {
      execute<BitVector448>();
   } else {
      printf("Maximum number of reactions supported is 448. Number of reactions in the network is %d", numReactions);
   }
}

/**
 * Main entry point for the code.
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code.
 */
int main(int argc, char** argv) {
   if (argc == 2) {
      execute(argv[1]);
      freeResources();
   } else {
      printf("Please specify the network file.");
   }
   return (EXIT_SUCCESS);
}
