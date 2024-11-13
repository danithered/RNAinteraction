#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ViennaRNA/part_func_up.h>
#include <ViennaRNA/fold.h>

struct RNA{
	char* seq;
	char* str;
	double mfe;
	unsigned int length;
	unsigned int max_length;
} RNA_def = {NULL, NULL, 0.0, 0, 0};

int initRNA(const unsigned int len, struct RNA *rna){
	if(rna->max_length < len){
		if(rna->seq || rna->str){
			free(rna->seq);
			free(rna->str);
		}
		
		// init seq
		rna->seq = calloc(len+1, sizeof(char));
		if(!rna->seq){ 
			printf("ERROR: initRNA: could not init!\n");
			return(0);
		}
		memset(rna->seq, '\0', len+1);

		// init str
		rna->str = calloc(len+1, sizeof(char));
		if(!rna->str){ 
			free(rna->seq);
			printf("ERROR: initRNA: could not init!\n");
			return(0);
		}
		memset(rna->seq, '\0', len+1);

		// set all the other values
		rna->length = 0;
		rna->mfe=0.0;
		rna->max_length = len;
		return(1);
	}
	return (2);
}

void freeRNA(struct RNA *rna){
	if(rna->max_length){
		if(rna->seq) free(rna->seq);
		if(rna->str) free(rna->str);
	}
}

int addRNA(char* seq, struct RNA *rna){
	const unsigned int len = strlen(seq);
	if(!initRNA(len, rna)){ // make sure it has enough space
		printf("ERROR: addRNA: could not init rna!\n");
		return(0);
	}

	strcpy(rna->seq, seq);
	rna->mfe = (double) vrna_fold(seq, rna->str);
	rna->length = len;

	return(1);
}

void printRNA(struct RNA *rna){
	printf("%s %s [%d] (%g)\n", rna->seq, rna->str, rna->length, rna->mfe);
}

int main(){
	struct RNA rna1 = RNA_def, rna2 = RNA_def;
	addRNA("GCCCGAAAAUUUUUU\0", &rna1);
	addRNA("GCCCGAAAAUUUUUU\0", &rna2);

	printRNA(&rna1);
	printRNA(&rna2);

	// function part
	struct RNA *left= &rna1, *right= $rna2;

	// init char storing structural constraint
	unsigned int constraint_length = rna1->length+rna2->length+1;
	char *constraint = (char*) calloc(constraint_length, sizeof(char));
	if(!constraint){
		printf("ERROR: Could not initialise array for constraint in size %d\n", constraint_length);
		return(0);
	}
	memset(constraint, 'x', --constraint_length); // init all to have no intermoleclar interaction
	constraint[constraint_length] = '\0'; // do not forget about terminator character 

	// have to find the start og left 5' dangling end
	unsigned int left_start= left->length; 
	while(left_start && left->str[left_start-1] == '.'){--left_start;}

	// have to find the end of right 3' dangling end
	int right_end= -1; 
	for(const unsigned int end = right->length-1; right_end < end && right->str[right_end+1] == '.'; ++right_end){}

	// which one is longer
	struct RNA *longer, *shorter;
	if(left->length < right->length){
		longer = right;
		shorter = left;
	} else {
		longer = left;
		shorter = right;
	}

	// find binding energy
	if(right_end >= 0 && left_start != left->length ) { // if there are dangling ends
		interact *interaction =pf_interact(
				longer->seq, // *s1 - the longer sequence
				shorter->seq, // *s2, - the shorter sequence
				NULL, // pu_contrib *p_c, - probabiltity of the bases of the first sequence to be unpaired
				NULL, // pu_contrib *p_c2, - unpaired probs for the second sequence
				, // int max_w, - maximal length of interaction
				constraint, // char *cstruc, - structural constraint: 'x' for no intermolecular interaction, '|' intermolecular binding, '.' no constraint
				0 , // int incr3, - inclusion of unpaired residues right of the region of interaction in ‘s1’
				0 // int incr5 - inclusion of unpaired residues left of the region of interaction in ‘s1’
				);

		// free data
		free(interact);
	}

	free(constraint);

	// end of function part

	// free
	freeRNA(&rna1);
	freeRNA(&rna2);
	
	return(0);
}

