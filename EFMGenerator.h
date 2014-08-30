#ifndef EFMGENERATOR_H
#define	EFMGENERATOR_H

#include <iostream>
#include <climits>
#include "BitPatternTree.h"
#include "PathwayPool.h"
#include "ReversibleTree.h"
#include "Stopwatch.h"
#include "ReactionCoefficients.h"

using namespace std;

template <class BitVector>
class EFMGenerator {
private:
   //List of partial pathways
   PathwayPool<BitVector> pathways;
   int combinationCount;
   //Total number of combinations generated
   long combinationsGenerated;
   //Total number combinations stored in the bitpattern tree
   //This is equal to total number of combinations - combinations rejected on the fly
   long combinationsStored;
   //Timer for removal of a metabolite
   Stopwatch timeMetaboliteRemvoed;
   BitPatternTree<BitVector> bptInput, bptOutput, bptNonpart;

   //Format time

   float format(double time) {
      return time;
   }

   //This function removes all the pathways that are inputs to the metabolite at metaboliteIndex

   void removeUnusedInputs(int metaboliteIndex) {
      for (unsigned int i = 0; i < pathways.size(); i++) {
         if (pathways[i].isInput(metaboliteIndex)) {
            pathways.remove(i);
            i--;
         }
      }
   }

   //This function removes all the pathways that are outputs to the metabolite at metaboliteIndex

   void removeUnusedOutputs(int metaboliteIndex) {
      for (unsigned int i = 0; i < pathways.size(); i++) {
         if (pathways[i].isOutput(metaboliteIndex)) {
            pathways.remove(i);
            i--;
         }
      }
   }

   //This function finds the index of metabolite to be removed. The next metabolite is selected
   //on the basis of minimum number of combinations of pathways generated if a metabolite is removed

   void findNextMetaboliteToRemove(int& metaboliteIndex, int& inputCount, int& outputCount) {
      metaboliteIndex = 0;
      inputCount = 0;
      outputCount = 0;
      long prod, minProd = LONG_MAX;
      int inputs, outputs;
      for (int m = 0; m < numMetabolitesRemaining; m++) {
         inputs = outputs = 0;
         for (int p = pathways.size() - 1; p >= 0; p--) {
            if (pathways[p].isInput(m)) {
               inputs++;
            } else if (pathways[p].isOutput(m)) {
               outputs++;
            }
         }
         prod = inputs;
         prod *= outputs;
         if (prod < minProd) {
            metaboliteIndex = m;
            inputCount = inputs;
            outputCount = outputs;
            minProd = prod;
         }
      }
   }


   //Process external metabolites
   //Processed metabolites are placed before unprocessed
   //This method reorders metabolites by moving external metabolites after internal metabolites

   void processExternals() {
      //Start from end of the
      for (int i = numMetabolitesRemaining - 1; i >= 0; i--) {
         //If i is external ignore. Otherwise swap it with an external that is at lower index
         if (externalMetabolites[i]) {
            numMetabolitesRemaining--;
            continue;
         }
         for (int j = i - 1; j >= 0; j--) {
            if (externalMetabolites[j]) {
               reorderMetabolites(i, j);
               numMetabolitesRemaining--;
               break;
            }
         }
      }
   }

   //Eliminates an unprocessed metabolite

   void removeNextMetabolite() {
      int metaboliteIndex, inputCount, outputCount;
      //Find next metabolite to remove
      findNextMetaboliteToRemove(metaboliteIndex, inputCount, outputCount);
      //Decrement number of remaining metabolites
      numMetabolitesRemaining--;
      //Move the metabolite to the end of unprocessed metabolites
      reorderMetabolites(metaboliteIndex, numMetabolitesRemaining);
      metaboliteIndex = numMetabolitesRemaining;
      long totalCombinations = inputCount;
      totalCombinations *= outputCount;
      cout << metabolites[metaboliteIndex] << "\t" << inputCount << "\t" << outputCount << "\t" << (pathways.size() - inputCount - outputCount) << "\t" << totalCombinations << "\t";
      if (inputCount == 0) {
         removeUnusedOutputs(metaboliteIndex);
      } else if (outputCount == 0) {
         removeUnusedInputs(metaboliteIndex);
      } else {
         generateCombinations(metaboliteIndex, inputCount, outputCount);
      }
   }

   //Prints data after a metabolite is removed

   void printMetaboliteData() {
      cout << combinationsGenerated << "\t" << combinationsStored << "\t" << pathways.size() << "\t" << format(timeMetaboliteRemvoed.elapsedTime) << endl;
   }

   //This function reorders metabolites by swapping metabolite at index1 and index2

   void reorderMetabolites(int index1, int index2) {
      if (index1 == index2) {
         return;
      }
      string name = metabolites[index1];
      metabolites[index1] = metabolites[index2];
      metabolites[index2] = name;

      bool external = externalMetabolites[index1];
      externalMetabolites[index1] = externalMetabolites[index2];
      externalMetabolites[index2] = external;

      double coeff;

      for (int i = pathways.size() - 1; i >= 0; i--) {
         coeff = pathways[i].metaboliteCoefficients[index1];
         pathways[i].metaboliteCoefficients[index1] = pathways[i].metaboliteCoefficients[index2];
         pathways[i].metaboliteCoefficients[index2] = coeff;
      }
   }

   //This function generates combinations for the metabolite at given index

   void generateCombinations(int metaboliteIndex, int inputCount, int outputCount) {
      bptInput.init();
      bptOutput.init();
      bptNonpart.init();
      //Identify input, output and non-participating pathways for the given metabolite
      int prevSize = pathways.size();
      //Sorting input pathways
      for (int i = 0, j = prevSize - 1; i < inputCount;) {
         if (pathways[i].isInput(metaboliteIndex)) {
            i++;
         } else if (pathways[j].isInput(metaboliteIndex)) {
            pathways.swap(i, j);
         } else {
            j--;
         }
      }
      //Sorting output pathways
      for (int i = inputCount, j = prevSize - 1; i < inputCount + outputCount;) {
         if (pathways[i].isOutput(metaboliteIndex)) {
            i++;
         } else if (pathways[j].isOutput(metaboliteIndex)) {
            pathways.swap(i, j);
            i++;
            j--;
         } else {
            j--;
         }
      }
      //Create reversible trees for inputs and outputs
      ReversibleTree<BitVector> inputTree(pathways, 0, inputCount);
      ReversibleTree<BitVector> outputTree(pathways, inputCount, inputCount + outputCount);
      for (int i = prevSize - 1, in = 0, out = inputCount; i >= 0; i--) {
         if (pathways[i].isInput(metaboliteIndex)) {
            bptInput.addPathway(&pathways[i]);
            in++;
         } else if (pathways[i].isOutput(metaboliteIndex)) {
            bptOutput.addPathway(&pathways[i]);
            out++;
         } else {
            bptNonpart.addPathway(&pathways[i]);
         }
      }
      //Generate combinations
      generateCombinations(inputTree.getRoot(), outputTree.getRoot());
      //Update list of pathways in the network
      for (int i = pathways.size() - 1; i >= prevSize; i--) {
         pathways[i].updateMetaboliteCoefficients(metaboliteIndex);
      }
      for (int i = prevSize - 1; i >= 0; i--) {
         if (pathways[i].isInput(metaboliteIndex) || pathways[i].isOutput(metaboliteIndex)) {
            pathways.remove(i);
         }
      }
      clearBPTNodePool();
      clearRevNodePool();
   }

   //This function recursively generates combinations from reversible tree nodes

   void generateCombinations(ReversibleTreeNode<BitVector>* input, ReversibleTreeNode<BitVector>* output) {
      BitVector comboLabel = input->getBitsUsed();
      comboLabel.bitwiseOr(output->getBitsUsed());
      if (!isValidPathway<BitVector>(comboLabel) || (comboLabel.getCardinality() > maxCardinality)) {
         return;
      }
      if (input->isLeaf()) {
         if (output->isLeaf()) {
            if (bptInput.isSuperSet(comboLabel, input->getBitsUsed()) || bptOutput.isSuperSet(comboLabel, output->getBitsUsed()) || bptNonpart.isSuperSet(comboLabel)) {
               return;
            }
            int iStart = input->getStart(), iEnd = input->getEnd();
            int oStart = output->getStart(), oEnd = output->getEnd();
            Pathway<BitVector> *combo = ++pathways;
            Pathway<BitVector> *in;
            Pathway<BitVector> *out;
            for (int i = iStart; i < iEnd; i++) {
               for (int o = oStart; o < oEnd; o++) {
                  in = &pathways[i];
                  out = &pathways[o];
                  combo->setParents(in, out);
                  if (combo->getCardinality() <= maxCardinality) {
                     combinationsGenerated++;
                     if (!(bptInput.isSuperSet(combo, in) || bptOutput.isSuperSet(combo, out) || bptNonpart.isSuperSet(combo))) {
                        combinationsStored++;
                        combo = ++pathways;
                     }
                  }
               }
            }
            --pathways;
         } else {
            generateCombinations(input, output->getNode0());
            generateCombinations(input, output->getNode1());
         }
      } else {
         if (output->isLeaf()) {
            generateCombinations(input->getNode0(), output);
            generateCombinations(input->getNode1(), output);
         } else {
            generateCombinations(input->getNode0(), output->getNode0());
            generateCombinations(input->getNode0(), output->getNode1());
            generateCombinations(input->getNode1(), output->getNode0());
            generateCombinations(input->getNode1(), output->getNode1());
         }
      }
   }


public:
   bool debug_met;

   EFMGenerator(Pathway<BitVector>* pths) {
      for (int i = 0; i < numReactions; i++) {
         pathways << pths[i];
      }
   }

   /**
    * Generates EFMs
    */
   void genenrateEFMs() {
      processExternals();
      cout << "Reactions\t" << numReactions - reversiblePairCount << endl;
      cout << "Reversible\t" << reversiblePairCount << endl;
      cout << "Metabolites\t" << numMetabolites << endl;
      cout << "Internal\t" << numMetabolitesRemaining << endl;
      cout << "Serial\tMetabolite Removed\tInputs\tOutputs\tNon-participating\tTotal Combinations\tCombinations Generated\tIndependent Pathways\tTotal Pathways\tTime" << endl;
      int serial = 1;
      double totalTime = 0;
      initBPTNodePool(sizeof (BitPatternTreeNode<BitVector>));
      initRevNodePool(sizeof (ReversibleTreeNode<BitVector>));
      while (numMetabolitesRemaining > 0) {
         debug_met = (serial == 1);
         combinationsGenerated = 0;
         combinationsStored = 0;
         timeMetaboliteRemvoed.restart();
         cout << serial << "\t";
         removeNextMetabolite();
         timeMetaboliteRemvoed.stop();
         totalTime += timeMetaboliteRemvoed.elapsedTime;
         printMetaboliteData();
         serial++;
      }
      freeBPTNodePool();
      freeRevNodePool();
      printf("Total Time : %g\n", totalTime);
   }

   /**
    * Prints EFMs after computing reaction coefficients.
    */
   void printEFMs() {
      ReactionCoefficients<BitVector> coeff;
      cout << "*** EFMs ***" << endl;
      coeff.printHeader();
      for (int i = 0; i < pathways.size(); i++) {
         cout << (i + 1);
         coeff.computeCoefficients(pathways[i].reactionBitData);
      }
   }

};

#endif	/* EFMGENERATOR_H */

