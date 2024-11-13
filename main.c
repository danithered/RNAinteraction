#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ViennaRNA/part_func_up.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/constraints/basic.h>
#include <ViennaRNA/partfunc/global.h>

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
	// binding rna1 5' dangling end to rna2 3' dangling end
	struct RNA *left= &rna1, *right= &rna2;

	// init char storing structural constraint
	unsigned int constraint_length = left->length+right->length+2;
	char *constraint = (char*) calloc(constraint_length, sizeof(char));
	if(!constraint){
		printf("ERROR: Could not initialise array for constraint in size %d\n", constraint_length);
		return(0);
	}
	memset(constraint, 'x', --constraint_length); // init all to have no intermoleclar interaction (decrement constraint_length)
	constraint[constraint_length] = '\0'; // do not forget about terminator character 

	// have to find the start og left 5' dangling end
	unsigned int left_start= left->length; 
	while(left_start && left->str[left_start-1] == '.'){--left_start;}

	// have to find the end of right 3' dangling end
	int right_end= -1; 
	for(const unsigned int end = right->length-1; right_end < end && right->str[right_end+1] == '.'; ++right_end){}

	// which one is longer
	struct RNA *longer, *shorter;
	unsigned int length5 = 0, length3 = 0; 

	if(left->length < right->length){
		longer = right;
		shorter = left;

		// construct structural binding: starting with right, and then left

		// have to find the end of right 3' dangling end
		for(unsigned int i = 0; i != right->length && right->str[i] == '.'; ++i){
			constraint[i] = '.';
			++length3;
		}

		// put separator character '&'. Left will start at constraint[right->length+1]
		constraint[right->length] = '&';
		
		// have to find the start og left 5' dangling end (do not forget: constraint_length = len1 + len2)
		for(int i = left->length-1, start = right->length +1; i != -1 && left->str[i] == '.'; --i){
			constraint[start + i] = '.';
			++length5;
		}

	} else {
		longer = left;
		shorter = right;
		// construct structural binding: starting with left, and then right
		
		// left 5' dangling end
		for(int i = left->length-1; i != -1 && left->str[i] == '.'; --i){
			constraint[i] = '.';
			++length5;
		}

		// put separator character '&'. Left will start at constraint[right->length+1]
		constraint[left->length] = '&';
		
		// right 3' dangling end 
		for(unsigned int i = 0, start = left->length+1; i != right->length && right->str[i] == '.'; ++i){
			constraint[start+i] = '.';
			++length3;
		}
	}
	printf("constraint computed: %s\n", constraint);
	
	// maximal window length
	const int w = (length3<length5)?length3:length5;

	// calculating real probabilities
	vrna_fold_compound_t *fc = vrna_fold_compound(longer->seq, NULL, VRNA_OPTION_DEFAULT);	
	vrna_constraints_add(fc, longer->str, VRNA_CONSTRAINT_DB_DEFAULT);
	//vrna_pf(fc, NULL); 
	//vrna_fold_compound_free(fc);
	
	//double temp = pf_fold_par(longer->seq, longer->str, NULL, 0, 1, 0);
	pf_fold(longer->seq, NULL);

	pu_contrib *pu = pf_unstru(longer->seq, w);
	free_pf_arrays();

	// calculating fake probabilities
	/*pu_contrib *pu = get_pu_contrib_struct(longer->length, w);
	for (unsigned int i = 1; i <= longer->length; ++i){
		for (unsigned int j = 0; j <= w; ++j) {
			pu->H[i][j] = pu->I[i][j] = pu->M[i][j] = pu->E[i][j] = (constraint[i-1]=="x")?0.0:2.5;
		}
	}*/


	// find binding energy
	if(length3 | length5) { // if there are dangling ends
		interact *interaction =pf_interact(
				longer->seq, // *s1 - the longer sequence
				shorter->seq, // *s2, - the shorter sequence
				pu, // pu_contrib *p_c, - probabiltity of the bases of the first sequence to be unpaired
				NULL, // pu_contrib *p_c2, - unpaired probs for the second sequence
				w, // int max_w, - maximal length of interaction
				constraint, // char *cstruc, - structural constraint: 'x' for no intermolecular interaction, '|' intermolecular binding, '.' no constraint
				0 , // int incr3, - inclusion of unpaired residues right of the region of interaction in ‘s1’
				0 // int incr5 - inclusion of unpaired residues left of the region of interaction in ‘s1’
				);

		// free data
		free_interact(interaction);
	}
	
	free_pu_contrib(pu);
	free(constraint);

	// end of function part

	// free
	freeRNA(&rna1);
	freeRNA(&rna2);
	
	return(0);
}
