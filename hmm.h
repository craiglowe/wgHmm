/*

hmm.h
Written by Craig Lowe

*/

#ifndef HMM_H
#define HMM_H

double addLog(double x, double y );

double subtractLog(double x, double y);

struct bed *genomeWideViterbiPath(struct bed *noGapList, unsigned int numStates, struct hash *emissionProbsHash, double **transitionProbs, double *startingProbs);

#endif

