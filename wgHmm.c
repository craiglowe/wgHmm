/*

wgHmm.c
Written by Craig Lowe

*/

#include "common.h"
#include "linefile.h"
#include "options.h"
#include "memalloc.h"
#include "bed.h"
#include "hmm.h"
#include "sqlNum.h"

/*---------------------------------------------------------------------------*/

static struct optionSpec optionSpecs[] =
{
	{"defaultState", OPTION_INT},
	{"inputIsLog", OPTION_BOOLEAN},
	{NULL, 0}
};

boolean optInputIsLog = FALSE;
int optDefaultState = -1;

/*---------------------------------------------------------------------------*/

void usage()
{
errAbort(
	"wgHmm - a simple and general implementation of HMMs for genomics\n"
	"usage:\n"
	"   wgHmm chrom.sizes noGaps.bed numberOfStates transition.matrix emissionProbs.wig output.bed\n"
	"options:\n"
	" -inputIsLog   (false)  This flag says that the transition matrix and the\n"
	"                        emisison probabilities are both in log-space.  This\n"
	"                        is probably a good idea in most cases.\n"
	" -defaultState (none)   This will report of Lod-Score comparing the selected\n"
	"                         state to the default state\n"
	"notes:\n"
	" chrom.sizes - tab separated file where the first column is the name\n"
	"                of the chromosome and the second column is the length\n"
	" noGaps.bed - bed file of regions no overlapping an assembly gap\n"
	" numberOfStates - number of states in the HMM/transducer\n"
	" transition.matrix - is a tab separated file where column X\n"
	"                      row Y gives the transition probabiliy\n"
	"                      from state X to state Y\n"
	" emissionProbs.wig - based on a fixed-step wiggle file, but with multiple columns\n"
	"                      one for each state.\n"
	);
}

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/

struct slChrom
{
	struct slChrom *next;
	char *name;
	unsigned int length;
};


struct slChrom *createChrom(char *name, unsigned int length)
{
	struct slChrom *ret;
	AllocVar(ret);
	ret->name = cloneString(name);
	ret->length = length;
	return(ret);
}


struct slChrom *bedListSortedToChromMax(struct bed *bedList)
{
	char *currName = "";
	struct bed *currBed = NULL;
	struct slChrom *head = NULL, *temp = NULL;
	for(currBed=bedList; currBed!=NULL; currBed=currBed->next)
	{
		if(differentString(currName,currBed->chrom))
		{
			currName = currBed->chrom;
			temp = createChrom(currBed->chrom, currBed->chromEnd);
			temp->next = head;
			head = temp;
		}
		else
		{
			temp->length = currBed->chromEnd;
		}
	}
	return(head);
}

struct slChrom *readChromInfoFile(char *chromSizeFilename)
{
	struct lineFile *lf = lineFileOpen(chromSizeFilename, TRUE);
	char *row[2];
	unsigned int chromLength = 0;
	struct slChrom *chromList = NULL, *currChrom = NULL;
	char *chrom = NULL;

	while(lineFileRow(lf, row))
	{
		chrom = row[0];
		chromLength = sqlUnsigned(row[1]);
		currChrom = createChrom(row[0], chromLength);
		currChrom->next = chromList;
		chromList = currChrom;
	}
	lineFileClose(&lf);
	return(chromList);
}

void printBedList(char *outputFilename, struct bed *head)
{
	FILE *outFile = mustOpen(outputFilename, "w");
	struct bed *curr = NULL;

	for(curr=head; curr != NULL; curr=curr->next)
	{
		if(optDefaultState >= 0){bedTabOutN(curr, 5, outFile);}
		else{bedTabOutN(curr, 4, outFile);}
	}
	carefulClose(&outFile);
}

struct hash *allocateEmissionHash(unsigned int numStates, struct slChrom *chromList)
{
	struct hash *answer = newHash(6);
	struct slChrom *currChrom = NULL;
	unsigned int x = 0;
	double **emissionProbs = NULL;

	for(currChrom=chromList; currChrom!=NULL; currChrom=currChrom->next)
	{
		verbose(3, "  Allocating %u for %s\n", currChrom->length, currChrom->name);
		emissionProbs = NULL;
		AllocArray(emissionProbs, numStates);
		for(x=0; x<numStates; x++)
		{
			AllocArray(emissionProbs[x], currChrom->length);
		}
		hashAdd(answer, currChrom->name, emissionProbs);
	}
	return(answer);
}


void setBaseDataFromWig(unsigned int numStates, struct hash *emissionHash, char *wigFilename)
{
	struct lineFile *lf = NULL;
	unsigned int rows = 0;
	unsigned int base = 0, x = 0;
	char *chrom = NULL;
	double **emissionByBase = NULL;

	if(numStates <= 4){rows=4;}
	else{rows = numStates;}
	char *row[rows];

	lf = lineFileOpen(wigFilename, TRUE);
	while(lineFileChop(lf, row))
	{
		if(sameString(row[0], "fixedStep"))
		{
			chrom = row[1] + 6;
			base = sqlUnsigned(row[2] + 6) - 1;
			emissionByBase = ((double **)hashMustFindVal(emissionHash, chrom));
		}
		else
		{
			if(emissionByBase == NULL){errAbort("wig file not in right format, don't know the chrom\n");}
			for(x=0; x<numStates; x++)
			{
				if(optInputIsLog){emissionByBase[x][base] = sqlDouble(row[x]);}
				else{emissionByBase[x][base] = log(sqlDouble(row[x]));}
			}
			base++;
		}
	}
	lineFileClose(&lf);
}


struct hash *createEmissionHash(unsigned int numStates, struct slChrom *chromList, char *inputWigFilename)
{
	struct hash *emissionHash = NULL;

	verbose(2, " Allocating memory for emission probs\n");
	emissionHash = allocateEmissionHash(numStates, chromList);
	verbose(2, " Reading file: %s\n", inputWigFilename);
	setBaseDataFromWig(numStates, emissionHash, inputWigFilename);
	return(emissionHash);
}


void setTransitionProbs(unsigned int numStates, char *transitionFilename, double *startingProbs, double **transitionProbs)
{
	struct lineFile *lf = NULL;
	char *row[numStates];
	unsigned int x = 0, line = 0;

	lf = lineFileOpen(transitionFilename, TRUE);
	while(lineFileChop(lf, row))
	{
		for(x=0; x<numStates; x++)
		{
			if(line == 0)
			{
				if(optInputIsLog){startingProbs[x] = sqlDouble(row[x]);}
				else{startingProbs[x] = log(sqlDouble(row[x]));}
			}
			else
			{
				if(optInputIsLog){transitionProbs[line-1][x] = sqlDouble(row[x]);}
				else{transitionProbs[line-1][x] = log(sqlDouble(row[x]));}
			}
		}
		line++;
	}
	lineFileClose(&lf);
}


/*---------------------------------------------------------------------------*/

void wgHmm(unsigned int numStates, char *chromInfoFilename, char *noGapFilename, char *transitionMatrixFilename, char *inputWigFilename, char *outputBedFilename)
{
	struct bed *noGapBedList = NULL, *regions = NULL;
	struct hash *emission = NULL;
	struct slChrom *chromList = NULL;
	double *startingProbs = NULL;
	double **transitionProbs = NULL;
	unsigned int x = 0;

	verbose(2, "Reading chromInfo file\n");
	chromList = readChromInfoFile(chromInfoFilename);
	verbose(2, "Reading noGap file\n");
	noGapBedList = bedLoadNAll(noGapFilename, 3);
	verbose(2, "Sorting noGap regions\n");
	slSort(&noGapBedList, bedCmp);	

	verbose(2, "Making the table of emission probabilities\n");
	emission = createEmissionHash(numStates, chromList, inputWigFilename);

	verbose(2, "Making the table of transition probabilities\n");
	AllocArray(startingProbs, numStates);
	AllocArray(transitionProbs, numStates);
	for(x=0;x<numStates;x++){AllocArray(transitionProbs[x], numStates);}
	setTransitionProbs(numStates, transitionMatrixFilename, startingProbs, transitionProbs);

	verbose(2, "Finding the viterbi path\n");
	regions = genomeWideViterbiPath(noGapBedList, numStates, emission, transitionProbs, startingProbs, optDefaultState);

	verbose(2, "Printing results\n");
	printBedList(outputBedFilename, regions);

	verbose(2, "Done\n");
}

/*---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
/* Process command line. */
{
	optionInit(&argc, argv, optionSpecs);
	if (argc != 7){usage();}

	optInputIsLog = optionExists("inputIsLog");
	optDefaultState = optionInt("defaultState", optDefaultState);

	char *chromInfoFilename = argv[1];
	char *noGapFilename = argv[2];
	unsigned int numStates = sqlUnsigned(argv[3]);
	char *transitionMatrixFilename = argv[4];
	char *inputWigFilename = argv[5];
	char *outputBedFilename = argv[6];

	wgHmm(numStates, chromInfoFilename, noGapFilename, transitionMatrixFilename, inputWigFilename, outputBedFilename);

	return(0);
}

