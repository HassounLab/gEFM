#ifndef REVERSIBLETREE_H
#define	REVERSIBLETREE_H

#include "PathwayPool.h"
#include "ReversibleTreeNode.h"

using namespace std;

template <class BitVector>
class ReversibleTree {
public:

   ReversibleTree(PathwayPool<BitVector>& list, int start, int end) {
      BitVector bitsCannotbeUsed(true);
      getBitsCannotBeUsed(bitsCannotbeUsed, list, start, end);
      BitVector bitsUsed;
      root = (ReversibleTreeNode<BitVector>*) nextRevNode(sizeof (ReversibleTreeNode<BitVector>));
      root->init(list, start, end, bitsCannotbeUsed, bitsUsed);
   }

   ~ReversibleTree() {
   }

   ReversibleTreeNode<BitVector>* getRoot() {
      return root;
   }


private:
   //Root of the tree
   ReversibleTreeNode<BitVector>* root;

   void getBitsCannotBeUsed(BitVector& bitsCannotbeUsed, PathwayPool<BitVector>& list, int start, int end) {
      unsigned int count;
      for (int r = 0; r < numReactions; r++) {
         bitsCannotbeUsed.setBit(r, !reversible[r]);
         if (!bitsCannotbeUsed[r]) {
            count = 0;
            for (int p = start; p < end; p++) {
               if (list[p].reactionBitData[r]) {
                  count++;
               } else {
                  count--;
               }
            }
            if (count < 0) {
               count = -count;
            }
         }
      }
   }
};

#endif	/* REVERSIBLETREE_H */

