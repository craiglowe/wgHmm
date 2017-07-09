/*

hmm.c
Written by Craig Lowe

*/

#include "common.h"
#include "hash.h"
#include "memalloc.h"
#include "obscure.h"
#include "sqlNum.h"
#include "hmm.h"
#include "bed.h"

double dTwoMax(double a, double b, unsigned short *passBack)
{
	if(a >= b)
	{
		*passBack=0;
		return(a);
	}
	else
	{
		*passBack=1;
		return(b);
	}
}

double dThreeMax(double a, double b, double c, unsigned short *passBack)
{
	if(a >= b && a >= c)
	{
		*passBack=0;
		return(a);
	}
	else if(b >= c)
	{
		*passBack=1;
		return(b);
	}
	else
	{
		*passBack=2;
		return(c);
	}
}

double dNMax(double *nums, unsigned int len, unsigned short *passBack)
{
	unsigned int x = 0;
	unsigned short bestIndex = 0;

	for(x=1, bestIndex=0; x<len; x++)
	{
		if(nums[x] > nums[bestIndex]){bestIndex = x;}
	}
	
	*passBack = bestIndex;
	return(nums[bestIndex]);
}

unsigned int dNMaxIndex(double *nums, unsigned int len)
{
	unsigned int x = 0, bestIndex = 0;

	for(x=1; x<len; x++)
	{
		if(nums[x] > nums[bestIndex]){bestIndex = x;}
	}

	return(bestIndex);
}


struct bed *pathToBed(char *chrom, unsigned int start, unsigned int end, unsigned short **path, unsigned int endState)
{
	unsigned int tempEnd = 0;
	unsigned int x = 0, xGenome = 0, prevRow = 0, currRow = 0;
	struct bed *answer = NULL, *currBed = NULL;
	char buffer[256];

	tempEnd = end;
	prevRow = endState;
	currRow = endState;
	for(x = end - start - 1; x > 1; x--)
	{
		xGenome = x + start;
		if(currRow != prevRow)
		{
			verbose(4, "change from %u to %u at %u\n", currRow, prevRow, xGenome);
			AllocVar(currBed);
			currBed->next = answer;
			currBed->chrom = cloneString(chrom);
			currBed->chromStart = xGenome;
			currBed->chromEnd = tempEnd;
			sprintf(buffer, "%u", prevRow);
			currBed->name = cloneString(buffer);
			answer = currBed;
			tempEnd = xGenome;
		}
		
		prevRow = currRow;
		currRow = path[currRow][x];
	}

	AllocVar(currBed);
	currBed->next = answer;
	currBed->chrom = cloneString(chrom);
	currBed->chromStart = start;
	currBed->chromEnd = tempEnd;
	sprintf(buffer, "%u", currRow);
	currBed->name = cloneString(buffer);
	answer = currBed;

	verbose(4, "done with path finding.  There were %d states.\n", slCount(answer));
	return(answer);
}

double addLog(double x, double y )
{
	if(x == -INFINITY){return(y);}
	if(y == -INFINITY){return(x);}
	if(x >= y){return(x + log(1 + exp(y-x)));}
	else{return(y + log(1 + exp(x-y)));}
}

double multiplyLog(double x, double y)
{
	if(x == -INFINITY || y == -INFINITY){return(-INFINITY);}
	else{return(x + y);}
}

double subtractLog(double x, double y)
{
	if(x < y){errAbort("Error: Taking log of negative number\n");}
	if(x == y){return(-INFINITY);}
	if(y == -INFINITY){return(x);}
	return(x + log(1 - exp(y-x)));
}


struct bed *viterbiPath(struct bed *roi, unsigned int numStates, struct hash *emissionProbsHash, double **transitionProbs, double *startingProbs)
{
	unsigned int prevState = 0, state = 0, x = 0, xGenome = 0, endState = 0, start = 0, end = 0;
	unsigned short **path = NULL;
	double **stateProb = NULL,  **emissionProbs = NULL;
	struct bed *regions = NULL;
	char *chrom = NULL;
	double probSum = 0;
	double *possibleOrigin = NULL;

	chrom = roi->chrom;
	start = roi->chromStart;
	end = roi->chromEnd;
	emissionProbs = hashMustFindVal(emissionProbsHash, roi->chrom);
	AllocArray(path, numStates);
	AllocArray(stateProb, numStates);
	AllocArray(possibleOrigin, numStates);
	for(state = 0; state < numStates; state++)
	{
		AllocArray(stateProb[state], end - start + 1);
		AllocArray(path[state], end - start + 1);
		path[state][0] = -1;
		stateProb[state][0] = multiplyLog(startingProbs[state], emissionProbs[state][start]);
	}

	verbose(4, "starting calculations\n");
	for(x = 1; x < end - start + 1; x++)
	{
		xGenome = x + start;

		for(state=0, probSum=-INFINITY; state<numStates; state++)
		{
			for(prevState=0; prevState<numStates; prevState++)
			{
				possibleOrigin[prevState] = stateProb[prevState][x-1] + transitionProbs[prevState][state] + emissionProbs[state][xGenome];
			}
			stateProb[state][x] = dNMax(possibleOrigin, numStates, &(path[state][x]));
			//stateProb[state][x] = dThreeMax(stateProb[0][x-1] + transitionProbs[0][state] + emissionProbs[state][xGenome],
			//		stateProb[1][x-1] + transitionProbs[1][state] + emissionProbs[state][xGenome],
			//		stateProb[2][x-1] + transitionProbs[2][state] + emissionProbs[state][xGenome],
			//		&(path[state][x]));
			probSum = addLog(probSum, stateProb[state][x]);
		}
		for(state=0; state<numStates; state++)
		{
			//stateProb[state][x] = stateProb[state][x] - probSum;
			verbose(5, "%u %u %e\n", state, x, stateProb[state][x]);
		}

		//uglyf("%d %d : %E %E %E : %E %E %E : %d %d %d\n", x, xGenome, stateProb[0][x], stateProb[1][x], stateProb[2][x] , emissionProbs[0][xGenome], emissionProbs[1][xGenome], emissionProbs[2][xGenome], (int)path[0][x], (int)path[1][x], (int)path[2][x]);
		/*for(state=0; state<numStates; state++)
		{
			if(stateProb[state][x] == -INFINITY){errAbort("Error: a probability has gone to zero\n");}
		}*/
	}

	for(state=0; state<numStates; state++)
	{
		possibleOrigin[state] = stateProb[state][end-start];
		verbose(4, "end prob for state %u is %e %e\n", state, stateProb[state][end-start], stateProb[state][end-start-1]);
	}

	endState = dNMaxIndex(possibleOrigin, numStates);
	//if(stateProb[0][end-start] >= stateProb[1][end-start] && stateProb[0][end-start] >= stateProb[2][end-start]){endState = 0;}
	//else if(stateProb[1][end-start] >= stateProb[2][end-start]){endState = 1;}
	//else{endState = 2;}

	verbose(4, "trying to figure out path\n");
	regions = pathToBed(chrom, start, end, path, endState);

	/* freeing memory */
	for(state = 0; state < numStates; state++)
	{
		freeMem(stateProb[state]);
		freeMem(path[state]);
	}
	freeMem(path);
	freeMem(stateProb);
	freeMem(possibleOrigin);

	return(regions);
}

void addLodScore(struct hash *emissionProbsHash, struct bed *regionList, int defaultState)
{
	struct bed *region = NULL;
	unsigned int x = 0;
	double sum = 0;
	double **emissionProbs = NULL;

	for(region=regionList; region != NULL; region=region->next)
	{
		emissionProbs = hashMustFindVal(emissionProbsHash, region->chrom);
		sum = 0;
		for(x=region->chromStart; x<region->chromEnd; x++)
		{
			sum += emissionProbs[sqlUnsigned(region->name)][x] - emissionProbs[defaultState][x];
		}
		region->score = (int)(sum/log(10));
	}
}


struct bed *bedCat(struct bed *a, struct bed *b)
{
	struct bed *curr = NULL;
	if(a == NULL){return(b);}
	if(b == NULL){return(a);}
	for(curr=a;curr->next!=NULL;curr=curr->next)
	{
	}
	curr->next=b;
	return(a);
}


struct bed *genomeWideViterbiPath(struct bed *noGapList, unsigned int numStates, struct hash *emissionProbsHash, double **transitionProbs, double *startingProbs, int defaultState)
{
	struct bed *roi = NULL, *altRegions = NULL, *temp = NULL;

	for(roi = noGapList; roi != NULL; roi = roi->next)
	{
		verbose(4, "starting region: %s %u %u\n", roi->chrom, roi->chromStart, roi->chromEnd);
		temp = viterbiPath(roi, numStates, emissionProbsHash, transitionProbs, startingProbs);
		if(defaultState >= 0){addLodScore(emissionProbsHash, temp, defaultState);}
		altRegions = bedCat(altRegions, temp);
	}
	return(altRegions);
}


