#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/subopt/wuchty.h>  
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>

vrna_subopt_solution_t* fn3(char *rna1, char *rna2, const double temperature){
	// dynamic allocation for concatenated string
	const unsigned int length1=strlen(rna1), length2=strlen(rna2);
	const unsigned int length=length1+length2+1;
	
	char *concatenated = (char*) calloc(length+1, sizeof(char));
	if(!concatenated){
		printf("could not allocate memory is size sizeof(char) * %d\n", length+1);
		return(NULL);
	}
	char *cstr = (char*) calloc(length+1, sizeof(char)); // while the doc says that it does not need it, vrna_dimer still uses this memory
	if(!cstr){
		printf("could not allocate memory is size sizeof(char) * %d\n", length+1);
		free(concatenated);
		return(NULL);
	}
	strcpy(concatenated, rna1);
	concatenated[length1] = '&';
	strcat(concatenated + length1 + 1, rna2);	

	/* initialize random number generator */
	vrna_init_rand();

	char* constraint = (char*) calloc(length+1, sizeof(char));
	if(!constraint){
		printf("could not allocate memory is size sizeof(char) * %d\n", length+1);
		free(concatenated);
		free(cstr);
		return(NULL);
	}
	
	// create constraint
	memset(constraint, 'e', length);
	constraint[length1]='&';

	//printf("%s\n", constraint);
	printf("%s\n", concatenated);

	/* create a new model details structure to store the Model Settings */
	vrna_md_t md;
	vrna_md_set_default(&md); 
	md.uniq_ML = 1; // keep stuff for suboptim
	
	// set temperature in celsius degrees. The default is 37 Celsius
	md.temperature = temperature;
	
	// create fold compound
	vrna_fold_compound_t *fc = vrna_fold_compound(concatenated,
                                                &md,//&(opt->md),
                                                VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);

	// tell them binding has to be external binding
	vrna_constraints_add(fc, constraint,
			//VRNA_CONSTRAINT_DB | 
			VRNA_CONSTRAINT_DB_X | 
			VRNA_CONSTRAINT_DB_INTERMOL |
			//VRNA_CONSTRAINT_DB_INTRAMOL |
			VRNA_CONSTRAINT_DB_DEFAULT |
			VRNA_CONSTRAINT_DB_ENFORCE_BP |
			VRNA_CONSTRAINT_DB_PIPE
			);

	// compute dimer structure
	float mfe = vrna_mfe_dimer(fc, cstr);
	if(mfe >= 0.0) {
		vrna_fold_compound_free(fc);
		free(constraint);
		free(concatenated);
		free(cstr);
		return(NULL);
	}

	// collect suboptimal structures
	int delta = (int)(-mfe*100.0)+1; // the range of subopt structures calculated around the optimal: delta * 0.01 kcal/mol 
	vrna_subopt_solution_t *subopts = vrna_subopt(fc, delta, 1, NULL); 

	// free
	vrna_fold_compound_free(fc);
	free(constraint);
	free(concatenated);
	free(cstr);

	return( subopts );
}

void freeSubopt(vrna_subopt_solution_t *l){
      	for (unsigned int i = 0; l[i].structure; i++) free(l[i].structure);
	free(l);
}

unsigned int countLength(vrna_subopt_solution_t* x){
	if(!x) return(0);
	unsigned int length=0;
	while(x[length].structure && x[length].energy<0.0) {++length;}
	return(length);
}

void makeStickyEnds(char* begin, char* end){
	if(!begin || !end || begin >= end) {
		fprintf(stderr, "ERROR: non valid begin and end pointers supplied for makeStickyEnds!\n");
		return;
	}

	char *c = begin; // iterator from first char
	// 3' dangling end 
	for(; c != end && *c == '.'; ++c){
		*c = 'e';
	}
	
	if(c != end){
		char *e = end-1; // iterator from last char
		// 5' dangling end
		for(; *e == '.'; --e){ // no need to guard begin, as there is structure in the beginning
			*e = 'e';
		}

		// inner part
		for(; c != e; ++c){
			if(*c == '.') *c = 'x';
		}
	}

}

///  binding a single sequence to the 4 possible dangling ends of a cofolded duplex
/**
 * Function that computes the MFE of a duplex and a simplex RNA-s binding with their 5' and 3' ends. Binding energy can be calculated with the substraction of structure MFEs from the MFE returned from this function. Please note, if you request the output of complex sequence and structure, those will have to be freed at the end, as they are allocated during the run of the function!
 *
 * @param[in] duplex1_seq Pointer to the sequence of the first member of the cofolded duplex RNA
 * @param[in] duplex2_seq Pointer to the sequence of the second member of the cofolded duplex RNA
 * @param[in] duplex_str Pointer to the dot-bracket 2D structure of the duplex RNA
 * @param[in] single_seq Pointer to the sequence of the unfolded, single RNA
 * @param[in] temperature Temperature of binding in Celsius degrees
 *
 * @return MFE of the composit. If it is positive, some error happened.
 */
vrna_subopt_solution_t* connect3(
		char *duplex1_seq, 
		char *duplex2_seq, 
		char *duplex_str,
		char *single_seq,
		const double temperature
)
{
	const unsigned int duplex1_length = strlen(duplex1_seq);
	const unsigned int duplex_length = strlen(duplex_str);
	const unsigned int single_length = strlen(single_seq);

	// init dynamic data
	const unsigned int constraint_length = duplex_length + single_length + 2;
	char *constraint = (char*) calloc(constraint_length, sizeof(char));
	if( !constraint ){
		fprintf(stderr, "ERROR: Could not initialise arrays in size %d\n", constraint_length);
		return(NULL);
	}
	char *concat     = (char*) calloc(constraint_length, sizeof(char));
	if( !concat ){
		fprintf(stderr, "ERROR: Could not initialise arrays in size %d\n", constraint_length);
		free(constraint);
		return(NULL);
	}	

	// create concatenated string
	strcpy(concat, duplex1_seq);
	concat[duplex1_length] = '&';
	strcpy(concat + duplex1_length + 1, duplex2_seq);
	concat[duplex_length] = '&';
	strcpy(concat + duplex_length + 1, single_seq);
	concat[constraint_length-1] = '\0';

	// create constraint
	// { // duplex 5' dangling end
		// int i = duplex_length-1;
		// for(; i != -1 && duplex_str[i] == '.'; --i){
			// constraint[i] = 'e';
			// ++length5;
		// }
		// for(; i != -1; --i){
			// if(duplex_str[i] == '.') constraint[i] = 'x';
			// else constraint[i] = duplex_str[i];
		// }
	// }
// 
	// { // duplex 3' dangling end 
		// for(unsigned int i = 0; i != duplex_length && duplex_str[i] == '.'; ++i){
			// constraint[i] = 'e';
		// }
	// }
	// constraint[duplex_length] = '&'; // separator
// 
	// { // single 3' dangling end 
		// unsigned int i = 0;
		// const unsigned int start = duplex_length+1;
		// for(; i != single_length && single_str[i] == '.'; ++i){
			// constraint[start+i] = 'e';
		// }
		// for(; i != single_length; ++i){
			// if(single_str[i] == '.') constraint[start+i] = 'x';
			// else constraint[start+i] = single_str[i];
		// }
	// }
// 
	// { // single 5' dangling end
		// unsigned int i = duplex_length-1;
		// const unsigned int start = ; 
		// for(; i != -1 && duplex_str[i] == '.'; --i){
			// constraint[i] = 'e';
			// ++length5;
		// }
	// }
	// constraint[constraint_length-1] = '\0'; // do not forget about terminator character 
	
	strcpy(constraint, duplex_str);
	constraint[duplex_length] = '&';
	// strcpy(constraint + duplex_length + 1, single_str);
	memset(constraint + duplex_length + 1, 'e', single_length);
	constraint[constraint_length-1] = '\0';

	// find & separating first and second part of duplex - if it has none it is undefined behaviour!
	// char *end = constraint + 1;
	// while(*end != '&'){++end;}
	
	// make ends sticky
	makeStickyEnds(constraint, constraint + duplex1_length); // make sticky ends for the first part of the duplex
	makeStickyEnds(constraint + duplex1_length + 1, constraint + duplex_length); // make sticky ends for the second part of the duplex
							      //
	// makeStickyEnds(constraint, end-1); // make sticky ends for the first part of the duplex
	// makeStickyEnds(end+1, constraint + duplex_length -1); // make sticky ends for the second part of the duplex
	// makeStickyEnds(constraint + duplex_length + 1, constraint + constraint_length - 1); // make single ends sticky

	// printf("constraint: %s\n", constraint);


	/* create a new model details structure to store the Model Settings */
	vrna_md_t md;
	vrna_md_set_default(&md); 
	
	// set temperature in celsius degrees
	md.temperature = temperature;
	
	// create fold compound
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
			VRNA_CONSTRAINT_DB_ENFORCE_BP |
			VRNA_CONSTRAINT_DB_PIPE
			);
	
	// compute dimer structure
	float mfe = vrna_mfe_dimer(fc, NULL);
	if(mfe >= 0.0) {
		vrna_fold_compound_free(fc);
		free(constraint);
		free(concat);
		return(NULL);
	}
      	
	// collect suboptimal structures
	int delta = (int)(-mfe*100.0)+1; // the range of subopt structures calculated around the optimal: delta * 0.01 kcal/mol 
	vrna_subopt_solution_t *subopts = vrna_subopt(fc, delta, 1, NULL); 

	// free stuff
	vrna_fold_compound_free(fc);
	free(constraint);
	free(concat);

	return( subopts );
}


int main(int argc, char** argv){
	char rna1[] = "ACCCCC\0";
	char rna2[] = "UGGGCC\0";
	char rna3[] = "AGGG\0";

	printf("Connecting RNA sequences %s and %s\n", rna1, rna2);

	vrna_subopt_solution_t *subopts = fn3(rna1, rna2, VRNA_MODEL_DEFAULT_TEMPERATURE );

	// print
	if(!subopts){ // if there was no beneficial structure print it
		printf("Could not find optimal solutions.\n");
		return(1);
	}

	// iterating tru list
	for(vrna_subopt_solution_t *subopt = subopts; subopt->structure; ++subopt){
		printf("%s\t%f\n", subopt->structure, subopt->energy); // printing it
	}

	const unsigned int no_subopt = countLength(subopts);
	if(no_subopt){
		printf("Connecting first RNA duplex of %d suboptimal structures and rna %s\n", no_subopt, rna3);
		vrna_subopt_solution_t *subopts2 = connect3(
				rna1, rna2, // the seq of the 1st and 2nd member of the cofolded complex
				subopts[0].structure, // structure of the cofolded complex
				rna3, // the single rna to bind to the complex
				VRNA_MODEL_DEFAULT_TEMPERATURE); // temperature
		// print
		if(!subopts2){ // if there was no beneficial structure print it
			printf("Could not find optimal solutions.\n");
			return(1);
		}

		// iterating tru list
		for(vrna_subopt_solution_t *subopt = subopts2; subopt->structure; ++subopt){
			printf("%s\t%f\n", subopt->structure, subopt->energy); // printing it
		}
	}

	// free and return
	freeSubopt(subopts);
	return 0;
}

