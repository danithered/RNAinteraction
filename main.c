#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ViennaRNA/part_func_up.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/constraints/basic.h>
#include <ViennaRNA/partfunc/global.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/utils/strings.h>

// RNA struct - you do not need this, just for myself
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

// functions - need this!
double fn(struct RNA *left, struct RNA *right){
	// binding rna1 5' dangling end to rna2 3' dangling end

	// init char storing structural constraint
	unsigned int constraint_length = left->length+right->length+2;
	char *constraint = (char*) calloc(constraint_length, sizeof(char));
	if(!constraint){
		printf("ERROR: Could not initialise array for constraint in size %d\n", constraint_length);
		return(0);
	}
	memset(constraint, 'x', constraint_length-1); // init all to have no intermoleclar interaction 
	constraint[constraint_length-1] = '\0'; // do not forget about terminator character 

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
		
		// have to find the start og left 5' dangling end 
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
	// vrna_fold_compound_t *fc = vrna_fold_compound(longer->seq, NULL, VRNA_OPTION_DEFAULT);	
	// vrna_constraints_add(fc, longer->str, VRNA_CONSTRAINT_DB_DEFAULT);
	//vrna_pf(fc, NULL); 
	//vrna_fold_compound_free(fc);
	
	//double temp = pf_fold_par(longer->seq, longer->str, NULL, 0, 1, 0);
	pf_fold(longer->seq, NULL);

	pu_contrib *pu = pf_unstru(longer->seq, w);
	free_pf_arrays();

	// calculating fake probabilities
	//pu_contrib *pu = get_pu_contrib_struct(longer->length, w);
	//for (unsigned int i = 1; i <= longer->length; ++i){
	//	for (unsigned int j = 0; j <= w; ++j) {
	//		pu->H[i][j] = pu->I[i][j] = pu->M[i][j] = pu->E[i][j] = (constraint[i-1]=="x")?0.0:2.5;
	//	}
	//}


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

	return(0.0);
}

///  binding left 5' dangling end to right 3' dangling end
/**
 * Function that computes the MFE of two RNA-s binding with their 5' and 3' ends. The first RNA binds with its 5' end to the 3' end of the second. Binding energy can be calculated with the substraction of structure MFEs from the MFE returned from this function. Please note, if you request the output of complex sequence and structure, those will have to be freed at the end, as they are allocated during the run of the function!
 *
 * @param[in] left_seq Pointer to the sequence of the RNA, whose 5' dangling end assotiaties
 * @param[in] left_str Pointer to the dot-bracket 2D structure of the RNA, whose 5' dangling end assotiaties.
 * @param[in] right_seq Pointer to the sequence of the RNA, whose 3' dangling end assotiaties
 * @param[in] right_str Pointer to the dot-bracket 2D structure of the RNA, whose 3' dangling end assotiaties.
 * @param[out] compl_seq Sequence of the complex. If NULL no output will be written.
 * @param[out] compl_str 2D stucture of the complex. If NULL no output will be written.
 *
 * @return MFE of the composit. If it is positive, some error happened.
 */
double fn2(
		char *left_seq, char *left_str,
		char *right_seq, char *right_str,
		char **compl_seq, char **compl_str)
{
	const unsigned int left_length = strlen(left_seq), right_length = strlen(right_seq);
	unsigned int length5 = 0, length3 = 0; 

	// init dynamic data
	unsigned int constraint_length = left_length + right_length + 2;
	char *constraint = (char*) calloc(constraint_length, sizeof(char));
	if( !constraint ){
		printf("ERROR: Could not initialise arrays in size %d\n", constraint_length);
		return(1.0);
	}
	char *concat     = (char*) calloc(constraint_length, sizeof(char));
	if( !concat ){
		printf("ERROR: Could not initialise arrays in size %d\n", constraint_length);
		free(constraint);
		return(1.0);
	}	
	char *concatstr = NULL;
	if(compl_str){
		concatstr  = (char*) calloc(constraint_length, sizeof(char));
		if( !concatstr ){
			printf("ERROR: Could not initialise arrays in size %d\n", constraint_length);
			free(constraint);
			free(concat);
			return(1.0);
		}
		memset(concatstr, '\0', constraint_length);
	}

	// create concatenated string
	strcpy(concat, left_seq);
	concat[left_length] = '&';
	strcpy(concat + left_length + 1, right_seq);
	concat[constraint_length-1] = '\0';
	if(compl_seq){
		*compl_seq = (char*) calloc(constraint_length, sizeof(char));
		if(*compl_seq) strcpy(*compl_seq, concat);
	}

	// create constraint
	{ // left 5' dangling end
		int i = left_length-1;
		for(; i != -1 && left_str[i] == '.'; --i){
			constraint[i] = 'e';
			++length5;
		}
		for(; i != -1; --i){
			if(left_str[i] == '.') constraint[i] = 'x';
			else constraint[i] = left_str[i];
		}
	}

	constraint[left_length] = '&'; // separator

	{ // right 3' dangling end 
		unsigned int i = 0;
		const unsigned int start = left_length+1;
		for(; i != right_length && right_str[i] == '.'; ++i){
			constraint[start+i] = 'e';
			++length3;
		}
		for(; i != right_length; ++i){
			if(right_str[i] == '.') constraint[start+i] = 'x';
			else constraint[start+i] = right_str[i];
		}
	}

	constraint[constraint_length-1] = '\0'; // do not forget about terminator character 
	
//	printf("constraint: %s\n", constraint);


	/* create a new model details structure to store the Model Settings */
	vrna_md_t md;
	vrna_md_set_default(&md); 
	vrna_fold_compound_t *fc = vrna_fold_compound(concat,
                                                &md,//&(opt->md),
                                                VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);

	// add hard constraint
	vrna_constraints_add(fc, constraint,
			//VRNA_CONSTRAINT_DB | 
			VRNA_CONSTRAINT_DB_X | 
			VRNA_CONSTRAINT_DB_INTERMOL |
			//VRNA_CONSTRAINT_DB_INTRAMOL |
			VRNA_CONSTRAINT_DB_DEFAULT |
			//VRNA_CONSTRAINT_DB_ENFORCE_BP |
			VRNA_CONSTRAINT_DB_PIPE
			);
	
	// compute dimer structure
	float mfe = vrna_mfe_dimer(fc, concatstr);
      	
	// write out structure
	if(compl_str){
		// get string
		char * concatstr2;
		for(unsigned int sep = 1; sep < fc->strands; ++sep){
			concatstr2 = (char *) vrna_cut_point_insert(concatstr, (int)fc->strand_start[sep] + (sep-1) );
			free(concatstr);
			concatstr = concatstr2;
		}
		// alloc memory for external usage
		*compl_str = (char*) calloc(constraint_length, sizeof(char));

		// write it out
		if(*compl_str) strcpy(*compl_str, concatstr);

		// free
		free(concatstr); // concatstr2 is always the same as this one
	}

	// free stuff
	vrna_fold_compound_free(fc);
	free(constraint);
	free(concat);

	return(mfe);
}

	// compute binding energy
	//mfe -= left_mfe + right_mfe;

int main(int argc, char** argv){
	// read in or load rna-s + also compute str and mfe
	struct RNA rna1 = RNA_def, rna2 = RNA_def, rna3 = RNA_def;
	if(argc < 4){
		addRNA("AUAUAAUUUGGGGGAUAUACCCCCCGGGGGGG\0", &rna1);
		addRNA("CCCCCCCCCGGGGGAUAUACCCCCCUUUUUU\0", &rna2);
		addRNA("AAAAAAAAAGGGGGAUAUACCCCCCU\0", &rna3);
	} else {
		addRNA(argv[1], &rna1);
		addRNA(argv[2], &rna2);
		addRNA(argv[3], &rna3);
	}

	// print RNAs
	printRNA(&rna1);
	printRNA(&rna2);
	printRNA(&rna3);

	
	// step one: bind rna2 to rna1
	char *outstr, *outseq; // here will be stored the structure of the complex
	const double mfe_comp = fn2(rna1.seq, rna1.str, rna2.seq, rna2.str, &outseq, &outstr); // run calculation
	const double be = mfe_comp - rna1.mfe - rna2.mfe; // compute binding energy	

	printf("complex (1+2):\n%s\n%s [%f]\nbinding energy: %f\n", outseq, outstr, mfe_comp, be);
	printf("Based on energies the binding%s occour.\n", (be<0.0)?"":" do not");

	// step two: bind rna3 to the 1st complex
	char *outstr2; // here will be stored the structure of the complex
	const double mfe_comp2 = fn2(outseq, outstr, rna3.seq, rna3.str, NULL, &outstr2);	
	const double be2 = mfe_comp2 - rna3.mfe - mfe_comp; // compute binding energy	

	printf("complex ((1+2)+3):\n%s [%f]\nbinding energy: %f\n", outstr2, mfe_comp2, be2);
	printf("Based on energies the binding%s occour.\n", (be2<0.0)?"":" do not");

	// step three: bind rna3 to rna2
	char *outstr3, *outseq3; // here will be stored the structure of the complex
	const double mfe_comp3 = fn2(rna2.seq, rna2.str, rna3.seq, rna3.str, &outseq3, &outstr3); // run calculation
	const double be3 = mfe_comp3 - rna3.mfe - rna2.mfe; // compute binding energy	

	printf("complex (2+3):\n%s\n%s [%f]\nbinding energy: %f\n", outseq3, outstr3, mfe_comp3, be3);
	printf("Based on energies the binding%s occour.\n", (be3<0.0)?"":" do not");

	// free
	freeRNA(&rna1);
	freeRNA(&rna2);
	free(outstr);
	free(outstr2);
	free(outseq);
	
	return(0);
}

