#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/subopt/wuchty.h>  
#include <ViennaRNA/utils/basic.h>

int main(int argc, char** argv){
	/* initialize random number generator */
	vrna_init_rand();

	if(argc < 2){
		printf("please supply sequence as argument!\n");
		return(1);
	}

	//const unsigned int length=strlen(argv[1]);

	// reserve memory for structure
	//char* str = (char*) calloc(length+1, sizeof(char));
	//if(!str){
	//	printf("could not allocate memory is size sizeof(char) * %d\n", length+1);
	//	return(2);
	//}
	
	/* create a new model details structure to store the Model Settings */
	vrna_md_t md;
	vrna_md_set_default(&md); 
	md.uniq_ML = 1; // keep stuff for suboptim
	
	// set temperature in celsius degrees. The default is 37 Celsius
	md.temperature = (argc < 3)?VRNA_MODEL_DEFAULT_TEMPERATURE:atof(argv[2]);
	
	// create fold compound
	vrna_fold_compound_t *fc = vrna_fold_compound(argv[1],
                                                &md,//&(opt->md),
                                                VRNA_OPTION_DEFAULT);

	// compute dimer structure
	//float mfe = vrna_mfe_dimer(fc, str);
	float mfe = vrna_mfe(fc, NULL);

	// collect suboptimal structures
	int delta = (int)(-mfe*100.0)+1; // the range of subopt structures calculated around the optimal: delta * 0.01 kcal/mol 
	vrna_subopt_solution_t *subopts = vrna_subopt(fc, delta, 1, NULL);

	// print
	for(vrna_subopt_solution_t *subopt = subopts; subopt->structure; ++subopt){
		printf("%s\t%f\n", subopt->structure, subopt->energy);
	}

	// free
	vrna_fold_compound_free(fc);
	//free(str);
	return 0;
}

