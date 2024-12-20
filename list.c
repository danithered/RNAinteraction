#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/subopt/wuchty.h>  
#include <ViennaRNA/utils/basic.h>

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
			//VRNA_CONSTRAINT_DB_ENFORCE_BP |
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
	while(x[length].structure && x[length].energy>0.0) {++length;}
	return(length);
}

int main(int argc, char** argv){
	char rna1[] = "AC\0";
	char rna2[] = "UG\0";

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

	// free and return
	freeSubopt(subopts);
	return 0;
}

