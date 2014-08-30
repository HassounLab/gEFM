#include <cassert>

#include "BitPatternTreeNode.h"
#include "ReversibleTreeNode.h"

//Total number of reactions in the network
int numReactions;
//Total number of metabolites in the network
int numMetabolites;
//Boolean flag for external metabolites (true if a metabolite is external)
vector<bool> externalMetabolites;
//Boolean flag for reversible reactions (true if a reaction is reversible)
vector<bool> reversible;
//Number of metabolites to be processed 
int numMetabolitesRemaining;
//Reaction names
vector<string> reactions;
//Metabolites names
vector<string> metabolites;
//Maximum cardinality for EFM computed from rank condition
int maxCardinality;

//Array for reversible tree indices
int* reversibleTreeIndices;

//Reversible reaction pairs
void* reversiblePairs;
//Number of reversible reaction pairs
int reversiblePairCount;

//Maximum number of Bitpattern Tree Nodes
int maxPoolSizeBPT;
//Pointer array of Bitpattern Tree Nodes
void* poolBPT;
//Number of Bitpattern Tree Nodes
int poolSizeBPT;

//Maximum number of Reversible Tree Nodes
int maxPoolSizeRev;
//Pointer array of Reversible Tree Nodes
void* poolRev;
//Number of Reversible Tree Nodes
int poolSizeRev;

//Lookup table to compute cardinality of bit vectors
int BIT_COUNT_LOOKUP[] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};

/**
 * Initializes Bitpattern Tree Node pool
 * @param size Size of Bitpattern Tree Node pool
 */
void initBPTNodePool(int size) {
   maxPoolSizeBPT = MAX_BPT_TREE_NODES;
   poolBPT = malloc(maxPoolSizeBPT * size);
   assert(poolBPT != NULL);
   poolSizeBPT = 0;
}

/**
 * Free memory resources allocated for Bitpattern Tree Node pool
 */
void freeBPTNodePool() {
   free(poolBPT);
   poolBPT = NULL;
}

/**
 * Clears the Bitpattern Tree Node pool
 */
void clearBPTNodePool() {
   poolSizeBPT = 0;
}

/**
 * Returns pointer to next free Bitpattern Tree Node in the pool
 * @param size Node size in bytes
 * @return Pointer to the node
 */
void* nextBPTNode(int size) {
   poolSizeBPT++;
   assert(poolSizeBPT < maxPoolSizeBPT);
   return ((char*) poolBPT) + ((poolSizeBPT - 1) * size);
}

/**
 * Initializes Reversible Tree Node pool
 * @param size Size of Reversible Tree Node pool
 */
void initRevNodePool(int size) {
   maxPoolSizeRev = MAX_REVERSIBLE_TREE_NODES;
   poolRev = malloc(maxPoolSizeRev * size);
   assert(poolRev != NULL);
   poolSizeRev = 0;
}

/**
 * Free memory resources allocated for Reversible Tree Node pool
 */
void freeRevNodePool() {
   free(poolRev);
   poolRev = NULL;
}

/**
 * Clears the Reversible Tree Node pool
 */
void clearRevNodePool() {
   poolSizeRev = 0;
}

/**
 * Returns pointer to next free Reversible Tree Node in the pool
 * @param size Node size in bytes
 * @return Pointer to the node
 */
void* nextRevNode(int size) {
   poolSizeRev++;
   assert(poolSizeRev < maxPoolSizeRev);
   return ((char*) poolRev) + ((poolSizeRev - 1) * size);
}
