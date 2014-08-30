#ifndef REACTIONCOEFFICIENTS_H
#define	REACTIONCOEFFICIENTS_H

#include "BitVector.h"
#include "GlobalData.h"
#include "Rank.h"

template <class BitVector>
class ReactionCoefficients {
private:

   double* base;
   double* matrix;
   double* coefficients;
   int* colId;
   int rows;
public:

   ReactionCoefficients() {
      int m, r, i;
      rows = 1;
      for (m = 0; m < numMetabolites; m++) {
         if (!network.external[m]) {
            rows++;
         }
      }
      base = (double*) malloc(rows * numReactions * sizeof (double));
      matrix = (double*) malloc(rows * numReactions * sizeof (double));
      coefficients = (double*) malloc(rows * sizeof (double));
      colId = (int*) malloc(numReactions * sizeof (int));
      for (m = 0, i = 0; m < numMetabolites; m++) {
         if (!network.external[m]) {
            for (r = 0; r < numReactions; r++, i++) {
               base[i] = network.s[m][r];
            }
         }
      }
      for (r = 0; r < numReactions; r++, i++) {
         base[i] = 0;
      }
   }

   ~ReactionCoefficients() {
      if (base != NULL) {
         free(base);
      }
      if (matrix != NULL) {
         free(matrix);
      }
      if (coefficients != NULL) {
         free(coefficients);
      }
   }

   void printHeader() {
      cout << "EFM";
      for (int i = 0; i < numReactions; i++) {
         cout << "\t" << reactions[i];
         if (reversible[i]) {
            i++;
         }
      }
      cout << endl;
   }

   void computeCoefficients(BitVector pathway) {
      int index, m, r;
      //Initialize column space
      for (r = 0; r < numReactions; r++) {
         colId[r] = r;
      }
      //Copy base matrix
      for (m = 0, index = 0; m < rows; m++) {
         coefficients[m] = 0;
         for (r = 0; r < numReactions; r++, index++) {
            matrix[index] = pathway[r] ? base[index] : 0;
         }
      }
      //Normalization vector
      coefficients[rows - 1] = 1;
      for (r = 0; r < numReactions; r++) {
         if (pathway[r]) {
            matrix[(rows - 1) * numReactions + r] = 1;
            break;
         }
      }
      //Row echelon form
      int rank = computeRank(matrix, numReactions, rows, colId, coefficients);
      //Row reduced echelon form
      for (r = rank - 1; r >= 0; r--) {
         for (m = 0; m < r; m++) {
            index = m * numReactions + r;
            coefficients[m] -= coefficients[r] * matrix[index];
            matrix[index] = 0;
         }
      }
      //Rearrange column space
      for (int i = 0; i < numReactions; i++) {
         matrix[i] = 0;
      }
      for (int i = 0; i < rank; i++) {
         matrix[colId[i]] = coefficients[i];
      }
      //Printing mode
      double coeff;
      for (int i = 0; i < numReactions; i++) {
         coeff = matrix[i];
         if (reversible[i]) {
            if (!pathway[i]) {
               coeff = -matrix[i + 1];
            }
            i++;
         }
         if (coeff == -0.0) {
            coeff = 0;
         }
         cout << "\t" << coeff;
      }
      cout << endl;
   }
};

#endif	/* REACTIONCOEFFICIENTS_H */

