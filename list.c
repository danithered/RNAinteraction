#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/subopt/wuchty.h>  
#include <ViennaRNA/utils/basic.h>

int main(int argc, char** argv){
	char rna1[] = "GGGG\0";
	char rna2[] = "CCCC\0";

	// dynamic allocation for concatenated string
	const unsigned int length1=strlen(rna1), length2=strlen(rna2);
	const unsigned int length=length1+length2+1;
	
	char *concatenated = (char*) calloc(length+1, sizeof(char));
	if(!concatenated){
		printf("could not allocate memory is size sizeof(char) * %d\n", length+1);
		return(2);
	}
	strcpy(concatenated, rna1);
	concatenated[length1] = '&';
	strcat(concatenated + length1 + 1, rna2);	

	/* initialize random number generator */
	vrna_init_rand();

	char* constraint = (char*) calloc(length+1, sizeof(char));
	if(!constraint){
		printf("could not allocate memory is size sizeof(char) * %d\n", length+1);
		return(2);
	}
	
	// create constraint
	memset(constraint, 'e', length);
	char *pos = strchr(argv[1], '&');
	if(pos){
		constraint[pos-argv[1]] = '&';
		//printf("%s\n", constraint);
	}

	/* create a new model details structure to store the Model Settings */
	vrna_md_t md;
	vrna_md_set_default(&md); 
	md.uniq_ML = 1; // keep stuff for suboptim
	
	// set temperature in celsius degrees. The default is 37 Celsius
	md.temperature = (argc < 3)?VRNA_MODEL_DEFAULT_TEMPERATURE:atof(argv[2]);
	
	// create fold compound
	vrna_fold_compound_t *fc = vrna_fold_compound(argv[1],
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
	float mfe = vrna_mfe_dimer(fc, NULL);

	// collect suboptimal structures
	int delta = (int)(-mfe*100.0)+1; // the range of subopt structures calculated around the optimal: delta * 0.01 kcal/mol 
	vrna_subopt_solution_t *subopts = vrna_subopt(fc, delta, 1, NULL);

	// print
	for(vrna_subopt_solution_t *subopt = subopts; subopt->structure; ++subopt){
		printf("%s\t%f\n", subopt->structure, subopt->energy);
	}

	// free
	vrna_fold_compound_free(fc);
	free(constraint);
	free(concatenated);
	return 0;
}

