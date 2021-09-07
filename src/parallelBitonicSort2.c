/*
   This file is part of mpiMarkDup
   
   Copyright Institut Curie 2020
   
   This software is a computer program whose purpose is to sort SAM file and mark duplicates.
   
   You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
   The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
   The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
   Module:
     parallelBitonicSort2.c

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* parallel_bitonic.c -- parallel bitonic sort of randomly generated list
 *     of integers
 *
 * Input:
 *     n: the global length of the list -- must be a power of 2.
 *
 * Output:
 *     The sorted list.
 *
 * Notes:
 *     1.  Assumes the number of processes p = 2^d and p divides n.
 *     2.  The lists are statically allocated -- size specified in MAX.
 *     3.  Keys are in the range 0 -- KEY_MAX-1.
 *     4.  Implementation can be made much more efficient by using
 *         pointers and avoiding re-copying lists in merges.
 *
 * See Chap 14, pp. 320 & ff. in PPMPI.
 *
 * https://www.eecis.udel.edu/~saunders/courses/372/01f/ppmpi_c/chap08/
 * https://www.eecis.udel.edu/~saunders/courses/372/01f/ppmpi_c/chap14a/
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <malloc.h>

#include "compat.h"

#include "parallelBitonicSort2.h"

// we limit to 1gb per proc

#define LOW 0
#define HIGH 1

// the number of key max is 1gb
#define key_mpi_t MPI_LONG_LONG_INT
#define index_mpi_t MPI_LONG_LONG_INT

static MPI_Comm COMM_WORLD;

/********************************************************************/

size_t *base_arr2;


static int compare_size_t_V2(const void *a, const void *b){

	if (*(const size_t *)a > *(const size_t *)b)
		return 1;
	else if (*(const size_t *)a < *(const size_t *)b)
		return -1;
	else
		return 0;

}

static int compare_size_t(const void *a, const void *b){

	 size_t aa = *(size_t *)a, bb = *(size_t *)b;
	 //fprintf(stderr, "rank %d :::::[MPIBITONIC2] we call parallel bitonic sort \n", split_rank);

	 if (base_arr2[aa] > base_arr2[bb])
		return 1;
	else if (base_arr2[aa] < base_arr2[bb])
		return -1;
	else
		return 0;
}


void ParallelBitonicSort2(
		MPI_Comm split_comm,
		int my_rank,
		int dimension,
		size_t *local_list, 	//stand for coordinates
		int *local_list1, 		//stand for read sizes
		int *local_list2, 		//stand for read rank
		size_t *local_list3,	//stand for read offset source
		int *local_list4,    	//stand for original rank offset source
		size_t list_size) {

	// In this version of bitonic there is no more index.
	// Everything is sorted according
	// the local_list vector.

    int       proc_set_size;
    unsigned  and_bit;
    size_t k = 0;
    COMM_WORLD = split_comm;

    Local_sort2(list_size, local_list, local_list1, local_list2, local_list3, local_list4);

    //we check the local_list is sorted
    for (k = 0; k < (list_size - 1); k++){
    	assert(local_list[k] <= local_list[k + 1]);
    }

    /* and_bit is a bitmask that, when "anded" with  */
    /* my_rank, tells us whether we're working on an */
    /* increasing or decreasing list                 */
    for (proc_set_size = 2, and_bit = 2; proc_set_size <= dimension;
    		proc_set_size = proc_set_size*2, and_bit = and_bit << 1){

        if ((my_rank & and_bit) == 0){

            Par_bitonic_sort_incr2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		proc_set_size,
            		my_rank);
        }
        else{

            Par_bitonic_sort_decr2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		proc_set_size,
            		my_rank);
        }
    }

    for (k = 0; k < (list_size - 1); k++){
      	assert(local_list[k] <= local_list[k + 1]);
    }
}

void ParallelBitonicSortAll(
	MPI_Comm split_comm,
	int my_rank,
	int dimension,
	size_t *local_list, 	          //stand for coordinates
	int *local_list1, 		          //stand for read sizes
	int *local_list2, 		          //stand for read rank
	size_t *local_list3,	          //stand for read offset source
	int *local_list4,    	          //stand for original rank offset source
	int *local_list5,                 //stand for qname keys 
	unsigned int *local_list6,        //stand for flags
	unsigned int *local_list7,        //stand for pair nums
	unsigned int *local_list8,        //stand for orientations
	int *local_list9,                 //stand for mate scores
	int *local_list10,                //stand for read LB 
	int *local_list11,                //stand for read chr names
	int *local_list12,                //stand for mate chr names
	int *local_list13,                //stand for mate ranks 
	int *local_list14,                //stand for physical location x
	int *local_list15, 
	int *local_list16,               //stand for physical location y 
	size_t *local_list17,             //stand for mate coordinates
	size_t *local_list18,
	size_t *local_list19,             //stand for unclipped_coordinates
	size_t list_size
){

	int       proc_set_size;
    unsigned  and_bit;
    size_t k = 0;
    COMM_WORLD = split_comm;

	Local_sortAll(list_size, local_list, local_list1, local_list2, local_list3, local_list4, local_list5, local_list6, local_list7, local_list8, local_list9, local_list10, local_list11, local_list12, local_list13, local_list14, local_list15, local_list16, local_list17, local_list18, local_list19);


    //we check the local_list is sorted
    for (k = 0; k < (list_size - 1); k++){
    	assert(local_list[k] <= local_list[k + 1]);
    }

	/* and_bit is a bitmask that, when "anded" with  */
    /* my_rank, tells us whether we're working on an */
    /* increasing or decreasing list                 */
    for (proc_set_size = 2, and_bit = 2; proc_set_size <= dimension;
    		proc_set_size = proc_set_size*2, and_bit = and_bit << 1){

        if ((my_rank & and_bit) == 0){

            Par_bitonic_sort_incrAll(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
            		proc_set_size,
            		my_rank);
        }
        else{

            Par_bitonic_sort_decrAll(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
            		proc_set_size,
            		my_rank);
        }
    }

    for (k = 0; k < (list_size - 1); k++){
      	assert(local_list[k] <= local_list[k + 1]);
	}

}

void ParallelBitonicSortAll2(
	MPI_Comm split_comm,
	int my_rank,
	int dimension,
	size_t *local_list, 	          //stand for coordinates
	int *local_list1, 		          //stand for read sizes
	int *local_list2, 		          //stand for read rank
	size_t *local_list3,	          //stand for read offset source
	size_t *local_list4,	          //stand for read offset source
	int *local_list5,    	          //stand for original rank offset source
	int *local_list6,                 //stand for qname keys 
	unsigned int *local_list7,        //stand for flags
	unsigned int *local_list8,        //stand for pair nums
	unsigned int *local_list9,        //stand for orientations
	int *local_list10,                 //stand for mate scores
	int *local_list11,                //stand for read LB 
	int *local_list12,                //stand for read chr names
	int *local_list13,                //stand for mate chr names
	int *local_list14,                //stand for mate ranks 
	int *local_list15,                //stand for physical location x
	int *local_list16,                //stand for physical location y 
	int *local_list17,
	size_t *local_list18,             //stand for mate coordinates
	size_t *local_list19,             //stand for unclipped_coordinates
	size_t *local_list20,
	size_t list_size
){
	int       proc_set_size;
    unsigned  and_bit;
    size_t k = 0;
    COMM_WORLD = split_comm;

	Local_sortAll2(list_size, local_list, local_list1, local_list2, local_list3, local_list4, local_list5, local_list6, local_list7, local_list8, local_list9, local_list10, local_list11, local_list12, local_list13, local_list14, local_list15, local_list16, local_list17, local_list18, local_list19, local_list20);

    //we check the local_list is sorted
    for (k = 0; k < (list_size - 1); k++){
    	assert(local_list[k] <= local_list[k + 1]);
    }

	/* and_bit is a bitmask that, when "anded" with  */
    /* my_rank, tells us whether we're working on an */
    /* increasing or decreasing list                 */
    for (proc_set_size = 2, and_bit = 2; proc_set_size <= dimension;
    		proc_set_size = proc_set_size*2, and_bit = and_bit << 1){

        if ((my_rank & and_bit) == 0){

            Par_bitonic_sort_incrAll2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
					local_list20,
            		proc_set_size,
            		my_rank);

        }
        else{
            Par_bitonic_sort_decrAll2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
					local_list20,
            		proc_set_size,
            		my_rank);
        }

    }

    for (k = 0; k < (list_size - 1); k++){
      	assert(local_list[k] <= local_list[k + 1]);
	}

}


/*********************************************************************/
void Local_sort2(
         size_t   list_size     /* in     */,
         size_t  *local_keys    /* in/out */,
         int     *local_keys1   /* in/out */,
         int     *local_keys2   /* in/out */,
         size_t  *local_keys3   /* in/out */,
         int     *local_keys4   /* in/out */
         ) {

	//we create an index vector
	size_t *local_keys_temp   	 = malloc(sizeof(size_t)*list_size);
	int    *local_keys_temp1     = malloc(sizeof(int)*list_size);
	int    *local_keys_temp2     = malloc(sizeof(int)*list_size);
	size_t *local_keys_temp3  	 = malloc(sizeof(size_t)*list_size);
	int    *local_keys_temp4  	 = malloc(sizeof(int)*list_size);

	assert(local_keys_temp);
	assert(local_keys_temp1);
	assert(local_keys_temp2);
	assert(local_keys_temp3);
	assert(local_keys_temp4);

	size_t *index_vector = (size_t *)malloc(sizeof(size_t)*list_size);

	size_t j = 0;

	for(j = 0; j < list_size; j++){
			index_vector[j] = j;
	}

	base_arr2 = local_keys;
	bitonic_qksort2(index_vector, list_size, 0, list_size - 1, compare_size_t);

	//then we apply loac index to local_keys
	for(j = 0; j < list_size; j++){
		local_keys_temp[j]  = local_keys[index_vector[j]];
		local_keys_temp1[j] = local_keys1[index_vector[j]];
		local_keys_temp2[j] = local_keys2[index_vector[j]];
		local_keys_temp3[j] = local_keys3[index_vector[j]];
		local_keys_temp4[j] = local_keys4[index_vector[j]];
	}

	for(j = 0; j < list_size; j++){
		local_keys[j]  = local_keys_temp[j];
		local_keys1[j] = local_keys_temp1[j];
		local_keys2[j] = local_keys_temp2[j];
		local_keys3[j] = local_keys_temp3[j];
		local_keys4[j] = local_keys_temp4[j];
	}

	free(index_vector);
	free(local_keys_temp);
	free(local_keys_temp1);
	free(local_keys_temp2);
	free(local_keys_temp3);
	free(local_keys_temp4);
	malloc_trim(0);
}

void Local_sortAll(
	size_t   list_size     /* in     */,
	size_t  *local_keys    /* in/out */,
	int     *local_keys1   /* in/out */,
	int     *local_keys2   /* in/out */,
	size_t  *local_keys3   /* in/out */,
	int     *local_keys4   /* in/out */,
	int     *local_keys5,
	unsigned int *local_keys6,
	unsigned int *local_keys7,
	unsigned int *local_keys8,
	int *local_keys9,
	int *local_keys10,
	int *local_keys11,
	int *local_keys12,
	int *local_keys13,
	int *local_keys14,
	int *local_keys15,
	int *local_keys16,
	size_t *local_keys17,
	size_t *local_keys18,
	size_t *local_keys19){

		//we create an index vector
	size_t *local_keys_temp   	 = malloc(sizeof(size_t)*list_size);
	int    *local_keys_temp1     = malloc(sizeof(int)*list_size);
	int    *local_keys_temp2     = malloc(sizeof(int)*list_size);
	size_t *local_keys_temp3  	 = malloc(sizeof(size_t)*list_size);
	int    *local_keys_temp4  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp5  	 = malloc(sizeof(int)*list_size);
	unsigned int    *local_keys_temp6  	 = malloc(sizeof(unsigned int)*list_size);
	unsigned int    *local_keys_temp7  	 = malloc(sizeof(unsigned int)*list_size);
	unsigned int    *local_keys_temp8  	 = malloc(sizeof(unsigned int)*list_size);
	int    *local_keys_temp9  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp10  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp11  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp12  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp13 	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp14  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp15  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp16  	 = malloc(sizeof(int)*list_size);
	size_t    *local_keys_temp17  	 = malloc(sizeof(size_t)*list_size);
	size_t    *local_keys_temp18 	 = malloc(sizeof(size_t)*list_size);
	size_t    *local_keys_temp19 	 = malloc(sizeof(size_t)*list_size);

	
	assert(local_keys_temp);
	assert(local_keys_temp1);
	assert(local_keys_temp2);
	assert(local_keys_temp3);
	assert(local_keys_temp4);
	assert(local_keys_temp5);
	assert(local_keys_temp6);
	assert(local_keys_temp7);
	assert(local_keys_temp8);
	assert(local_keys_temp9);
	assert(local_keys_temp10);
	assert(local_keys_temp11);
	assert(local_keys_temp12);
	assert(local_keys_temp13);
	assert(local_keys_temp14);
	assert(local_keys_temp15);
	assert(local_keys_temp16);
	assert(local_keys_temp17);
	assert(local_keys_temp18);
	assert(local_keys_temp19);

	size_t *index_vector = (size_t *)malloc(sizeof(size_t)*list_size);

	size_t j = 0;

	for(j = 0; j < list_size; j++){
			index_vector[j] = j;
	}

	base_arr2 = local_keys;
	bitonic_qksort2(index_vector, list_size, 0, list_size - 1, compare_size_t);

	//then we apply loac index to local_keys
	for(j = 0; j < list_size; j++){
		local_keys_temp[j]  = local_keys[index_vector[j]];
		local_keys_temp1[j] = local_keys1[index_vector[j]];
		local_keys_temp2[j] = local_keys2[index_vector[j]];
		local_keys_temp3[j] = local_keys3[index_vector[j]];
		local_keys_temp4[j] = local_keys4[index_vector[j]];
		local_keys_temp5[j] = local_keys5[index_vector[j]];
		local_keys_temp6[j] = local_keys6[index_vector[j]];
		local_keys_temp7[j] = local_keys7[index_vector[j]];
		local_keys_temp8[j] = local_keys8[index_vector[j]];
		local_keys_temp9[j] = local_keys9[index_vector[j]];
		local_keys_temp10[j] = local_keys10[index_vector[j]];
		local_keys_temp11[j] = local_keys11[index_vector[j]];
		local_keys_temp12[j] = local_keys12[index_vector[j]];
		local_keys_temp13[j] = local_keys13[index_vector[j]];
		local_keys_temp14[j] = local_keys14[index_vector[j]];
		local_keys_temp15[j] = local_keys15[index_vector[j]];
		local_keys_temp16[j] = local_keys16[index_vector[j]];
		local_keys_temp17[j] = local_keys17[index_vector[j]];
		local_keys_temp18[j] = local_keys18[index_vector[j]];
		local_keys_temp19[j] = local_keys19[index_vector[j]];


	}

	for(j = 0; j < list_size; j++){
		local_keys[j]  = local_keys_temp[j];
		local_keys1[j] = local_keys_temp1[j];
		local_keys2[j] = local_keys_temp2[j];
		local_keys3[j] = local_keys_temp3[j];
		local_keys4[j] = local_keys_temp4[j];
		local_keys5[j] = local_keys_temp5[j];
		local_keys6[j] = local_keys_temp6[j];
		local_keys7[j] = local_keys_temp7[j];
		local_keys8[j] = local_keys_temp8[j];
		local_keys9[j] = local_keys_temp9[j];
		local_keys10[j] = local_keys_temp10[j];
		local_keys11[j] = local_keys_temp11[j];
		local_keys12[j] = local_keys_temp12[j];
		local_keys13[j] = local_keys_temp13[j];
		local_keys14[j] = local_keys_temp14[j];
		local_keys15[j] = local_keys_temp15[j];
		local_keys16[j] = local_keys_temp16[j];
		local_keys17[j] = local_keys_temp17[j];
		local_keys18[j] = local_keys_temp18[j];
		local_keys19[j] = local_keys_temp19[j];


	}

	free(index_vector);
	free(local_keys_temp);
	free(local_keys_temp1);
	free(local_keys_temp2);
	free(local_keys_temp3);
	free(local_keys_temp4);
	free(local_keys_temp5);
	free(local_keys_temp6);
	free(local_keys_temp7);
	free(local_keys_temp8);
	free(local_keys_temp9);
	free(local_keys_temp10);
	free(local_keys_temp11);
	free(local_keys_temp12);
	free(local_keys_temp13);
	free(local_keys_temp14);
	free(local_keys_temp15);
	free(local_keys_temp16);
	free(local_keys_temp17);
	free(local_keys_temp18);
	free(local_keys_temp19);

	malloc_trim(0);
}

void Local_sortAll2(
	size_t   list_size     /* in     */,
	size_t  *local_keys    /* in/out */,
	int     *local_keys1   /* in/out */,
	int     *local_keys2   /* in/out */,
	size_t  *local_keys3   /* in/out */,
	size_t  *local_keys4   /* in/out */,
	int     *local_keys5   /* in/out */,
	int     *local_keys6,
	unsigned int *local_keys7,
	unsigned int *local_keys8,
	unsigned int *local_keys9,
	int *local_keys10,
	int *local_keys11,
	int *local_keys12,
	int *local_keys13,
	int *local_keys14,
	int *local_keys15,
	int *local_keys16,
	int *local_keys17,
	size_t *local_keys18,
	size_t *local_keys19,
	size_t *local_keys20){

		//we create an index vector
	size_t *local_keys_temp   	 = malloc(sizeof(size_t)*list_size);
	int    *local_keys_temp1     = malloc(sizeof(int)*list_size);
	int    *local_keys_temp2     = malloc(sizeof(int)*list_size);
	size_t *local_keys_temp3  	 = malloc(sizeof(size_t)*list_size);
	size_t *local_keys_temp4  	 = malloc(sizeof(size_t)*list_size);
	int    *local_keys_temp5  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp6  	 = malloc(sizeof(int)*list_size);
	unsigned int    *local_keys_temp7  	 = malloc(sizeof(unsigned int)*list_size);
	unsigned int    *local_keys_temp8  	 = malloc(sizeof(unsigned int)*list_size);
	unsigned int    *local_keys_temp9  	 = malloc(sizeof(unsigned int)*list_size);
	int    *local_keys_temp10  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp11  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp12  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp13  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp14 	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp15  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp16  	 = malloc(sizeof(int)*list_size);
	int    *local_keys_temp17  	 = malloc(sizeof(int)*list_size);
	size_t    *local_keys_temp18  	 = malloc(sizeof(size_t)*list_size);
	size_t    *local_keys_temp19 	 = malloc(sizeof(size_t)*list_size);
	size_t    *local_keys_temp20 	 = malloc(sizeof(size_t)*list_size);

	
	assert(local_keys_temp);
	assert(local_keys_temp1);
	assert(local_keys_temp2);
	assert(local_keys_temp3);
	assert(local_keys_temp4);
	assert(local_keys_temp5);
	assert(local_keys_temp6);
	assert(local_keys_temp7);
	assert(local_keys_temp8);
	assert(local_keys_temp9);
	assert(local_keys_temp10);
	assert(local_keys_temp11);
	assert(local_keys_temp12);
	assert(local_keys_temp13);
	assert(local_keys_temp14);
	assert(local_keys_temp15);
	assert(local_keys_temp16);
	assert(local_keys_temp17);
	assert(local_keys_temp18);
	assert(local_keys_temp19);
	assert(local_keys_temp20);

	size_t *index_vector = (size_t *)malloc(sizeof(size_t)*list_size);

	size_t j = 0;

	for(j = 0; j < list_size; j++){
			index_vector[j] = j;
	}

	base_arr2 = local_keys;
	bitonic_qksort2(index_vector, list_size, 0, list_size - 1, compare_size_t);

	//then we apply loac index to local_keys
	for(j = 0; j < list_size; j++){
		local_keys_temp[j]  = local_keys[index_vector[j]];
		local_keys_temp1[j] = local_keys1[index_vector[j]];
		local_keys_temp2[j] = local_keys2[index_vector[j]];
		local_keys_temp3[j] = local_keys3[index_vector[j]];
		local_keys_temp4[j] = local_keys4[index_vector[j]];
		local_keys_temp5[j] = local_keys5[index_vector[j]];
		local_keys_temp6[j] = local_keys6[index_vector[j]];
		local_keys_temp7[j] = local_keys7[index_vector[j]];
		local_keys_temp8[j] = local_keys8[index_vector[j]];
		local_keys_temp9[j] = local_keys9[index_vector[j]];
		local_keys_temp10[j] = local_keys10[index_vector[j]];
		local_keys_temp11[j] = local_keys11[index_vector[j]];
		local_keys_temp12[j] = local_keys12[index_vector[j]];
		local_keys_temp13[j] = local_keys13[index_vector[j]];
		local_keys_temp14[j] = local_keys14[index_vector[j]];
		local_keys_temp15[j] = local_keys15[index_vector[j]];
		local_keys_temp16[j] = local_keys16[index_vector[j]];
		local_keys_temp17[j] = local_keys17[index_vector[j]];
		local_keys_temp18[j] = local_keys18[index_vector[j]];
		local_keys_temp19[j] = local_keys19[index_vector[j]];
		local_keys_temp20[j] = local_keys20[index_vector[j]];


	}

	for(j = 0; j < list_size; j++){
		local_keys[j]  = local_keys_temp[j];
		local_keys1[j] = local_keys_temp1[j];
		local_keys2[j] = local_keys_temp2[j];
		local_keys3[j] = local_keys_temp3[j];
		local_keys4[j] = local_keys_temp4[j];
		local_keys5[j] = local_keys_temp5[j];
		local_keys6[j] = local_keys_temp6[j];
		local_keys7[j] = local_keys_temp7[j];
		local_keys8[j] = local_keys_temp8[j];
		local_keys9[j] = local_keys_temp9[j];
		local_keys10[j] = local_keys_temp10[j];
		local_keys11[j] = local_keys_temp11[j];
		local_keys12[j] = local_keys_temp12[j];
		local_keys13[j] = local_keys_temp13[j];
		local_keys14[j] = local_keys_temp14[j];
		local_keys15[j] = local_keys_temp15[j];
		local_keys16[j] = local_keys_temp16[j];
		local_keys17[j] = local_keys_temp17[j];
		local_keys18[j] = local_keys_temp18[j];
		local_keys19[j] = local_keys_temp19[j];
		local_keys20[j] = local_keys_temp20[j];


	}

	free(index_vector);
	free(local_keys_temp);
	free(local_keys_temp1);
	free(local_keys_temp2);
	free(local_keys_temp3);
	free(local_keys_temp4);
	free(local_keys_temp5);
	free(local_keys_temp6);
	free(local_keys_temp7);
	free(local_keys_temp8);
	free(local_keys_temp9);
	free(local_keys_temp10);
	free(local_keys_temp11);
	free(local_keys_temp12);
	free(local_keys_temp13);
	free(local_keys_temp14);
	free(local_keys_temp15);
	free(local_keys_temp16);
	free(local_keys_temp17);
	free(local_keys_temp18);
	free(local_keys_temp19);
	free(local_keys_temp20);

	malloc_trim(0);
}
	


/*********************************************************************/
int Key_compare2(const size_t* p, const size_t* q) {

    if (*p < *q)
        return -1;
    else if (*p == *q)
        return 0;
    else /* *p > *q */
        return 1;

}  /* Key_compare */


/********************************************************************/
void Par_bitonic_sort_incr2(
        size_t      list_size      /* in     */,
        size_t*    	local_list     /* in/out */,
        int*    	local_list1    /* in/out */,
        int*    	local_list2    /* in/out */,
        size_t*    	local_list3    /* in/out */,
        int*    	local_list4    /* in/out */,
        int       	proc_set_size  /* in     */,
        int 	  	my_rank
        ) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    //int       my_rank;

    MPI_Comm_rank(COMM_WORLD, &my_rank);

    proc_set_dim = log2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);

    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;
        if (my_rank < partner){

            Merge_split2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		LOW,
            		partner
            		);
        }
        else{

            Merge_split2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		HIGH,
            		partner
            		);
        }
        eor_bit = eor_bit >> 1;
    }
}  /* Par_bitonic_sort_incr */


void Par_bitonic_sort_incrAll(
        size_t      list_size      /* in     */,
        size_t*    	local_list     /* in/out */,
        int*    	local_list1    /* in/out */,
        int*    	local_list2    /* in/out */,
        size_t*    	local_list3    /* in/out */,
        int*    	local_list4    /* in/out */,
		int* local_list5,
		unsigned int* local_list6,
		unsigned int* local_list7,
		unsigned int* local_list8,
		int* local_list9,
		int* local_list10,
		int* local_list11,
		int* local_list12,
		int* local_list13,
		int* local_list14,
		int* local_list15,
		int* local_list16,
		size_t* local_list17,
		size_t* local_list18,
		size_t* local_list19,
        int       	proc_set_size  /* in     */,
        int 	  	my_rank
        ) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    //int       my_rank;

    MPI_Comm_rank(COMM_WORLD, &my_rank);

    proc_set_dim = log2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);

    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;
        if (my_rank < partner){

            Merge_splitAll(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
            		LOW,
            		partner
            		);
        }
        else{

            Merge_splitAll(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
            		HIGH,
            		partner
            		);
        }
        eor_bit = eor_bit >> 1;
    }
}  /* Par_bitonic_sort_incr */

void Par_bitonic_sort_incrAll2(
        size_t      list_size      /* in     */,
        size_t*    	local_list     /* in/out */,
        int*    	local_list1    /* in/out */,
        int*    	local_list2    /* in/out */,
        size_t*    	local_list3    /* in/out */,
        size_t*    	local_list4    /* in/out */,
        int*    	local_list5    /* in/out */,
		int         *local_list6,
		unsigned int *local_list7,
		unsigned int *local_list8,
		unsigned int *local_list9,
		int *local_list10,
		int *local_list11,
		int *local_list12,
		int *local_list13,
		int *local_list14,
		int *local_list15,
		int *local_list16,
		int* local_list17,
		size_t* local_list18,
		size_t* local_list19,
		size_t* local_list20,
        int       	proc_set_size  /* in     */,
        int 	  	my_rank
        ) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    //int       my_rank;

    MPI_Comm_rank(COMM_WORLD, &my_rank);

    proc_set_dim = log2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);

    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;
        if (my_rank < partner){

            Merge_splitAll2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
					local_list20,
            		LOW,
            		partner
            		);
        }
        else{

            Merge_splitAll2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
					local_list20,
            		HIGH,
            		partner
            		);
        }
        eor_bit = eor_bit >> 1;
    }
}  /* Par_bitonic_sort_incr */



/********************************************************************/
void Par_bitonic_sort_decr2(
        size_t	list_size      /* in     */,
        size_t* local_list     /* in/out */,
        int*    local_list1    /* in/out */,
        int*    local_list2    /* in/out */,
        size_t* local_list3    /* in/out */,
        int*    local_list4    /* in/out */,
        int     proc_set_size  /* in     */,
        int 	my_rank
        ) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    //int       my_rank;

    MPI_Comm_rank(COMM_WORLD, &my_rank);

    proc_set_dim = log2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);
    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;

        if (my_rank > partner){
            Merge_split2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		LOW,
            		partner
            );
        }
        else{
            Merge_split2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		HIGH,
            		partner
            );
        }
        eor_bit = eor_bit >> 1;
    }

} /* Par_bitonic_sort_decr */


void Par_bitonic_sort_decrAll(
        size_t	list_size      /* in     */,
        size_t* local_list     /* in/out */,
        int*    local_list1    /* in/out */,
        int*    local_list2    /* in/out */,
        size_t* local_list3    /* in/out */,
        int*    local_list4    /* in/out */,
		int* local_list5,
		unsigned int* local_list6,
		unsigned int* local_list7,
		unsigned int* local_list8,
		int* local_list9,
		int* local_list10,
		int* local_list11,
		int* local_list12,
		int* local_list13,
		int* local_list14,
		int* local_list15,
		int* local_list16,
		size_t* local_list17,
		size_t* local_list18,
		size_t* local_list19,
        int     proc_set_size  /* in     */,
        int 	my_rank
        ) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    //int       my_rank;

    MPI_Comm_rank(COMM_WORLD, &my_rank);

    proc_set_dim = log2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);
    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;

        if (my_rank > partner){
            Merge_splitAll(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
            		LOW,
            		partner
            );
        }
        else{
            Merge_splitAll(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
            		HIGH,
            		partner
            );
        }
        eor_bit = eor_bit >> 1;
    }

} /* Par_bitonic_sort_decr */

void Par_bitonic_sort_decrAll2(
        size_t      list_size      /* in     */,
        size_t*    	local_list     /* in/out */,
        int*    	local_list1    /* in/out */,
        int*    	local_list2    /* in/out */,
        size_t*    	local_list3    /* in/out */,
        size_t*    	local_list4    /* in/out */,
        int*    	local_list5    /* in/out */,
		int         *local_list6,
		unsigned int *local_list7,
		unsigned int *local_list8,
		unsigned int *local_list9,
		int *local_list10,
		int *local_list11,
		int *local_list12,
		int *local_list13,
		int *local_list14,
		int *local_list15,
		int *local_list16,
		int* local_list17,
		size_t* local_list18,
		size_t* local_list19,
		size_t* local_list20,
        int       	proc_set_size  /* in     */,
        int 	  	my_rank
        ) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    //int       my_rank;

    MPI_Comm_rank(COMM_WORLD, &my_rank);

    proc_set_dim = log2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);
    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;

        if (my_rank > partner){

            Merge_splitAll2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
					local_list20,
            		LOW,
            		partner
            		);
        }
        else{

            Merge_splitAll2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
					local_list5,
					local_list6,
					local_list7,
					local_list8,
					local_list9,
					local_list10,
					local_list11,
					local_list12,
					local_list13,
					local_list14,
					local_list15,
					local_list16,
					local_list17,
					local_list18,
					local_list19,
					local_list20,
            		HIGH,
            		partner
            		);
        }
        eor_bit = eor_bit >> 1;
    }
} 


/********************************************************************/
void Merge_split2(
        size_t    list_size     /* in     */,
        size_t    *local_list1   /* in/out */,
        int	      *local_list2  /* in/out */,
        int       *local_list3  /* in/out */,
        size_t    *local_list4  /* in/out */,
        int       *local_list5  /* in/out */,
        int       which_keys    /* in     */,
        int       partner       /* in     */) {

	 int number_amount;
	 size_t k=0;

	 size_t *temp_key_list1 	= malloc(list_size*sizeof(size_t));
	 size_t *temp_key_list4   	= malloc(list_size*sizeof(size_t));
	 int *temp_key_list2   		= malloc(list_size*sizeof(int));
	 int *temp_key_list3 		= malloc(list_size*sizeof(int));
	 int *temp_key_list5   		= malloc(list_size*sizeof(int));

	 assert(temp_key_list1 != 0);
	 assert(temp_key_list2 != 0);
	 assert(temp_key_list3 != 0);
	 assert(temp_key_list4 != 0);
	 assert(temp_key_list5 != 0);

	 temp_key_list1 = memset(temp_key_list1, 0, sizeof(size_t)*list_size);
	 temp_key_list2 = memset(temp_key_list2, 0, sizeof(int)*list_size);
	 temp_key_list3 = memset(temp_key_list3, 0, sizeof(int)*list_size);
	 temp_key_list4 = memset(temp_key_list4, 0, sizeof(size_t)*list_size);
	 temp_key_list5 = memset(temp_key_list5, 0, sizeof(int)*list_size);


	 /*
	  * we pack the data into 1 vector called interbuff
	  *
	  */

	 MPI_Status status;
	 int res;
	 size_t *interbuff = NULL;
	 res = MPI_Alloc_mem((5*list_size)*sizeof(size_t), MPI_INFO_NULL, &interbuff);
	 assert(res == MPI_SUCCESS);
	 
	 for ( k = 0 ; k < list_size; k++ ){
		 interbuff[k] 				=  local_list1[k];
		 interbuff[k + list_size] 	=  (size_t)local_list2[k];
		 interbuff[k + 2*list_size] =  (size_t)local_list3[k];
		 interbuff[k + 3*list_size] =  local_list4[k];
		 interbuff[k + 4*list_size] =  (size_t)local_list5[k];
	 }

	 size_t *interbuff2;
	 res = MPI_Alloc_mem((5*list_size)*sizeof(size_t), MPI_INFO_NULL, &interbuff2);
	 assert(res == MPI_SUCCESS);


	 MPI_Sendrecv(interbuff,
			 	  5*list_size,
			 	  MPI_LONG_LONG_INT,
			 	  partner,
			 	  0,
			 	  interbuff2,
			 	  5*list_size,
			 	  MPI_LONG_LONG_INT,
			 	  partner,
			 	  0,
			 	  COMM_WORLD,
			 	  &status);


	 MPI_Get_count(&status, MPI_PACKED, &number_amount);

	 for ( k = 0 ; k < list_size; k++ ){
		 temp_key_list1[k] = (size_t)interbuff2[k];
		 temp_key_list2[k] = (int)   interbuff2[k + list_size];
		 temp_key_list3[k] = (int)   interbuff2[k + 2*list_size];
		 temp_key_list4[k] = (size_t)interbuff2[k + 3*list_size];
		 temp_key_list5[k] = (int)   interbuff2[k + 4*list_size];
	 }


    if (which_keys == HIGH){
    	Merge_list_high2(
    			 list_size,
    			 local_list1,
    			 local_list2,
    			 local_list3,
    			 local_list4,
    			 local_list5,
    			 temp_key_list1,
    			 temp_key_list2,
    			 temp_key_list3,
    			 temp_key_list4,
    			 temp_key_list5
    			 );
    }
    else{
        Merge_list_low2(
        		list_size,
        		local_list1,
        		local_list2,
        		local_list3,
        		local_list4,
        		local_list5,
        		temp_key_list1,
        		temp_key_list2,
        		temp_key_list3,
        		temp_key_list4,
        		temp_key_list5
        		);

    }

    if (temp_key_list1) free(temp_key_list1);
    if (temp_key_list2) free(temp_key_list2);
    if (temp_key_list3) free(temp_key_list3);
    if (temp_key_list4) free(temp_key_list4);
    if (temp_key_list5) free(temp_key_list5);

	MPI_Free_mem(interbuff);
	MPI_Free_mem(interbuff2);
    //malloc_trim(0);

} /* Merge_split */

void Merge_splitAll(
        size_t    list_size     /* in     */,
        size_t    *local_list1   /* in/out */,
        int	      *local_list2  /* in/out */,
        int       *local_list3  /* in/out */,
        size_t    *local_list4  /* in/out */,
        int       *local_list5  /* in/out */,
		int         *local_list6,
		unsigned int *local_list7,
		unsigned int *local_list8,
		unsigned int *local_list9,
		int *local_list10,
		int *local_list11,
		int *local_list12,
		int *local_list13,
		int *local_list14,
		int *local_list15,
		int *local_list16,
		int *local_list17,
		size_t *local_list18,
		size_t *local_list19,
		size_t *local_list20,
        int       which_keys    /* in     */,
        int       partner       /* in     */) {

	 int number_amount;
	 size_t k=0;

	 size_t *temp_key_list1 	= malloc(list_size*sizeof(size_t));
	 size_t *temp_key_list4   	= malloc(list_size*sizeof(size_t));
	 int *temp_key_list2   		= malloc(list_size*sizeof(int));
	 int *temp_key_list3 		= malloc(list_size*sizeof(int));
	 int *temp_key_list5   		= malloc(list_size*sizeof(int));
	 int *temp_key_list6   		= malloc(list_size*sizeof(int));
	 unsigned int *temp_key_list7   		= malloc(list_size*sizeof(unsigned int));
	 unsigned int *temp_key_list8   		= malloc(list_size*sizeof(unsigned int));
	 unsigned int *temp_key_list9   		= malloc(list_size*sizeof(unsigned int));
	 int *temp_key_list10   		= malloc(list_size*sizeof(int));
	 int *temp_key_list11   		= malloc(list_size*sizeof(int));
	 int *temp_key_list12  		= malloc(list_size*sizeof(int));
	 int *temp_key_list13   		= malloc(list_size*sizeof(int));
	 int *temp_key_list14   		= malloc(list_size*sizeof(int));
	 int *temp_key_list15  		= malloc(list_size*sizeof(int));
	 int *temp_key_list16  		= malloc(list_size*sizeof(int));
	 int *temp_key_list17  		= malloc(list_size*sizeof(int));

	 size_t *temp_key_list18   		= malloc(list_size*sizeof(size_t));
	 size_t *temp_key_list19   		= malloc(list_size*sizeof(size_t));	 
	 size_t *temp_key_list20   		= malloc(list_size*sizeof(size_t));	 

	 assert(temp_key_list1 != 0);
	 assert(temp_key_list2 != 0);
	 assert(temp_key_list3 != 0);
	 assert(temp_key_list4 != 0);
	 assert(temp_key_list5 != 0);
	 assert(temp_key_list6 != 0);
	 assert(temp_key_list7 != 0);
	 assert(temp_key_list8 != 0);
	 assert(temp_key_list9 != 0);
	 assert(temp_key_list10 != 0);
	 assert(temp_key_list11 != 0);
	 assert(temp_key_list12 != 0);
	 assert(temp_key_list13 != 0);
	 assert(temp_key_list14 != 0);
	 assert(temp_key_list15 != 0);
	 assert(temp_key_list16 != 0);
	 assert(temp_key_list17 != 0);
	 assert(temp_key_list18 != 0);
	 assert(temp_key_list19 != 0);
	 assert(temp_key_list20 != 0);

	 temp_key_list1 = memset(temp_key_list1, 0, sizeof(size_t)*list_size);
	 temp_key_list2 = memset(temp_key_list2, 0, sizeof(int)*list_size);
	 temp_key_list3 = memset(temp_key_list3, 0, sizeof(int)*list_size);
	 temp_key_list4 = memset(temp_key_list4, 0, sizeof(size_t)*list_size);
	 temp_key_list5 = memset(temp_key_list5, 0, sizeof(int)*list_size);
	 temp_key_list6 = memset(temp_key_list6, 0, sizeof(int)*list_size);
	 temp_key_list7 = memset(temp_key_list7, 0, sizeof(unsigned int)*list_size);
	 temp_key_list8 = memset(temp_key_list8, 0, sizeof(unsigned int)*list_size);
	 temp_key_list9 = memset(temp_key_list9, 0, sizeof(unsigned int)*list_size);
	 temp_key_list10 = memset(temp_key_list10, 0, sizeof(int)*list_size);
	 temp_key_list11 = memset(temp_key_list11, 0, sizeof(int)*list_size);
	 temp_key_list12 = memset(temp_key_list12, 0, sizeof(int)*list_size);
	 temp_key_list13 = memset(temp_key_list13, 0, sizeof(int)*list_size);
	 temp_key_list14 = memset(temp_key_list14, 0, sizeof(int)*list_size);
	 temp_key_list15 = memset(temp_key_list15, 0, sizeof(int)*list_size);
	 temp_key_list16 = memset(temp_key_list16, 0, sizeof(int)*list_size);
	 temp_key_list17 = memset(temp_key_list17, 0, sizeof(int)*list_size);
	 temp_key_list18 = memset(temp_key_list18, 0, sizeof(size_t)*list_size);
	 temp_key_list19 = memset(temp_key_list19, 0, sizeof(size_t)*list_size);
	 temp_key_list20 = memset(temp_key_list20, 0, sizeof(size_t)*list_size);



	 /*
	  * we pack the data into 1 vector called interbuff
	  *
	  */

	 MPI_Status status;
	 int res;
	 size_t *interbuff = NULL;
	 res = MPI_Alloc_mem((20*list_size)*sizeof(size_t), MPI_INFO_NULL, &interbuff);
	 assert(res == MPI_SUCCESS);
	 
	 for ( k = 0 ; k < list_size; k++ ){
		 interbuff[k] 				=  local_list1[k];
		 interbuff[k + list_size] 	=  (size_t)local_list2[k];
		 interbuff[k + 2*list_size] =  (size_t)local_list3[k];
		 interbuff[k + 3*list_size] =  local_list4[k];
		 interbuff[k + 4*list_size] =  (size_t)local_list5[k];
		 interbuff[k + 5*list_size] =  (size_t)local_list6[k];
		 interbuff[k + 6*list_size] =  (size_t)local_list7[k];
		 interbuff[k + 7*list_size] =  (size_t)local_list8[k];
		 interbuff[k + 8*list_size] =  (size_t)local_list9[k];
		 interbuff[k + 9*list_size] =  (size_t)local_list10[k];
		 interbuff[k + 10*list_size] =  (size_t)local_list11[k];
		 interbuff[k + 11*list_size] =  (size_t)local_list12[k];
		 interbuff[k + 12*list_size] =  (size_t)local_list13[k];
		 interbuff[k + 13*list_size] =  (size_t)local_list14[k];
		 interbuff[k + 14*list_size] =  (size_t)local_list15[k];
		 interbuff[k + 15*list_size] =  (size_t)local_list16[k];
		 interbuff[k + 16*list_size] =  (size_t)local_list17[k];
		 interbuff[k + 17*list_size] =  local_list18[k];
		 interbuff[k + 18*list_size] =  local_list19[k];
		 interbuff[k + 19*list_size] =  local_list20[k];

	 }

	 size_t *interbuff2;
	 res = MPI_Alloc_mem((20*list_size)*sizeof(size_t), MPI_INFO_NULL, &interbuff2);
	 assert(res == MPI_SUCCESS);


	 MPI_Sendrecv(interbuff,
			 	  20*list_size,
			 	  MPI_LONG_LONG_INT,
			 	  partner,
			 	  0,
			 	  interbuff2,
			 	  20*list_size,
			 	  MPI_LONG_LONG_INT,
			 	  partner,
			 	  0,
			 	  COMM_WORLD,
			 	  &status);


	 MPI_Get_count(&status, MPI_PACKED, &number_amount);

	 for ( k = 0 ; k < list_size; k++ ){
		 temp_key_list1[k] = (size_t)interbuff2[k];
		 temp_key_list2[k] = (int)   interbuff2[k + list_size];
		 temp_key_list3[k] = (int)   interbuff2[k + 2*list_size];
		 temp_key_list4[k] = (size_t)interbuff2[k + 3*list_size];
		 temp_key_list5[k] = (int)   interbuff2[k + 4*list_size];
		 temp_key_list6[k] = (int)   interbuff2[k + 5*list_size];
		 temp_key_list7[k] = (unsigned int)   interbuff2[k + 6*list_size];
		 temp_key_list8[k] = (unsigned int)   interbuff2[k + 7*list_size];
		 temp_key_list9[k] = (unsigned int)   interbuff2[k + 8*list_size];
		 temp_key_list10[k] = (int)   interbuff2[k + 9*list_size];
		 temp_key_list11[k] = (int)   interbuff2[k + 10*list_size];
		 temp_key_list12[k] = (int)   interbuff2[k + 11*list_size];
		 temp_key_list13[k] = (int)   interbuff2[k + 12*list_size];
		 temp_key_list14[k] = (int)   interbuff2[k + 13*list_size];
		 temp_key_list15[k] = (int)   interbuff2[k + 14*list_size];
		 temp_key_list16[k] = (int)   interbuff2[k + 15*list_size];
		 temp_key_list17[k] = (int)   interbuff2[k + 16*list_size];
		 temp_key_list18[k] = (size_t)   interbuff2[k + 17*list_size];
		 temp_key_list19[k] = (size_t)   interbuff2[k + 18*list_size];
		 temp_key_list20[k] = (size_t)   interbuff2[k + 19*list_size];

	 }


    if (which_keys == HIGH){
    	Merge_list_highAll(
    			 list_size,
    			 local_list1,
    			 local_list2,
    			 local_list3,
    			 local_list4,
    			 local_list5,
				 local_list6,
				 local_list7,
				 local_list8,
				 local_list9,
				 local_list10,
				 local_list11,
				 local_list12,
				 local_list13,
				 local_list14, 
				 local_list15,
				 local_list16,
				 local_list17,
				 local_list18,
				 local_list19,
				 local_list20,
    			 temp_key_list1,
    			 temp_key_list2,
    			 temp_key_list3,
    			 temp_key_list4,
    			 temp_key_list5,
				 temp_key_list6,
				 temp_key_list7,
				 temp_key_list8,
				 temp_key_list9, 
				 temp_key_list10,
				 temp_key_list11,
				 temp_key_list12,
				 temp_key_list13,
				 temp_key_list14,
				 temp_key_list15,
				 temp_key_list16,
				 temp_key_list17,
				 temp_key_list18,
				 temp_key_list19,
				 temp_key_list20
    			 );
    }
    else{
        Merge_list_lowAll(
        		list_size,
    			 local_list1,
    			 local_list2,
    			 local_list3,
    			 local_list4,
    			 local_list5,
				 local_list6,
				 local_list7,
				 local_list8,
				 local_list9,
				 local_list10,
				 local_list11,
				 local_list12,
				 local_list13,
				 local_list14, 
				 local_list15,
				 local_list16,
				 local_list17,
				 local_list18,
				 local_list19,
				 local_list20,
    			 temp_key_list1,
    			 temp_key_list2,
    			 temp_key_list3,
    			 temp_key_list4,
    			 temp_key_list5,
				 temp_key_list6,
				 temp_key_list7,
				 temp_key_list8,
				 temp_key_list9, 
				 temp_key_list10,
				 temp_key_list11,
				 temp_key_list12,
				 temp_key_list13,
				 temp_key_list14,
				 temp_key_list15,
				 temp_key_list16,
				 temp_key_list17,
				 temp_key_list18,
				 temp_key_list19,
				 temp_key_list20
        		);

    }

    if (temp_key_list1) free(temp_key_list1);
    if (temp_key_list2) free(temp_key_list2);
    if (temp_key_list3) free(temp_key_list3);
    if (temp_key_list4) free(temp_key_list4);
    if (temp_key_list5) free(temp_key_list5);
    if (temp_key_list6) free(temp_key_list6);
    if (temp_key_list7) free(temp_key_list7);
    if (temp_key_list8) free(temp_key_list8);
    if (temp_key_list9) free(temp_key_list9);
    if (temp_key_list10) free(temp_key_list10);
    if (temp_key_list11) free(temp_key_list11);
    if (temp_key_list12) free(temp_key_list12);
    if (temp_key_list13) free(temp_key_list13);
    if (temp_key_list14) free(temp_key_list14);
    if (temp_key_list15) free(temp_key_list15);
    if (temp_key_list16) free(temp_key_list16);
    if (temp_key_list17) free(temp_key_list17);
    if (temp_key_list18) free(temp_key_list18);
	if (temp_key_list19) free(temp_key_list19);
    if (temp_key_list20) free(temp_key_list20);


	MPI_Free_mem(interbuff);
	MPI_Free_mem(interbuff2);
    //malloc_trim(0);

} /* Merge_split */


void Merge_splitAll2(
        size_t    list_size     /* in     */,
        size_t    *local_list1   /* in/out */,
        int	      *local_list2  /* in/out */,
        int       *local_list3  /* in/out */,
        size_t    *local_list4  /* in/out */,
		size_t    *local_list5  /* in/out */,
        int       *local_list6  /* in/out */,
		int         *local_list7,
		unsigned int *local_list8,
		unsigned int *local_list9,
		unsigned int *local_list10,
		int *local_list11,
		int *local_list12,
		int *local_list13,
		int *local_list14,
		int *local_list15,
		int *local_list16,
		int *local_list17,
		int *local_list18,
		size_t *local_list19,
		size_t *local_list20,
		size_t *local_list21,
        int       which_keys    /* in     */,
        int       partner       /* in     */) {

	 int number_amount;
	 size_t k=0;

	 size_t *temp_key_list1 	= malloc(list_size*sizeof(size_t));
	 size_t *temp_key_list4   	= malloc(list_size*sizeof(size_t));
	 size_t *temp_key_list5   	= malloc(list_size*sizeof(size_t));
	 int *temp_key_list2   		= malloc(list_size*sizeof(int));
	 int *temp_key_list3 		= malloc(list_size*sizeof(int));
	 int *temp_key_list6   		= malloc(list_size*sizeof(int));
	 int *temp_key_list7   		= malloc(list_size*sizeof(int));
	 unsigned int *temp_key_list8   		= malloc(list_size*sizeof(unsigned int));
	 unsigned int *temp_key_list9   		= malloc(list_size*sizeof(unsigned int));
	 unsigned int *temp_key_list10   		= malloc(list_size*sizeof(unsigned int));
	 int *temp_key_list11   		= malloc(list_size*sizeof(int));
	 int *temp_key_list12   		= malloc(list_size*sizeof(int));
	 int *temp_key_list13  		= malloc(list_size*sizeof(int));
	 int *temp_key_list14   		= malloc(list_size*sizeof(int));
	 int *temp_key_list15   		= malloc(list_size*sizeof(int));
	 int *temp_key_list16  		= malloc(list_size*sizeof(int));
	 int *temp_key_list17  		= malloc(list_size*sizeof(int));
	 int *temp_key_list18  		= malloc(list_size*sizeof(int));

	 size_t *temp_key_list19   		= malloc(list_size*sizeof(size_t));
	 size_t *temp_key_list20   		= malloc(list_size*sizeof(size_t));	 
	 size_t *temp_key_list21   		= malloc(list_size*sizeof(size_t));	 

	 assert(temp_key_list1 != 0);
	 assert(temp_key_list2 != 0);
	 assert(temp_key_list3 != 0);
	 assert(temp_key_list4 != 0);
	 assert(temp_key_list5 != 0);
	 assert(temp_key_list6 != 0);
	 assert(temp_key_list7 != 0);
	 assert(temp_key_list8 != 0);
	 assert(temp_key_list9 != 0);
	 assert(temp_key_list10 != 0);
	 assert(temp_key_list11 != 0);
	 assert(temp_key_list12 != 0);
	 assert(temp_key_list13 != 0);
	 assert(temp_key_list14 != 0);
	 assert(temp_key_list15 != 0);
	 assert(temp_key_list16 != 0);
	 assert(temp_key_list17 != 0);
	 assert(temp_key_list18 != 0);
	 assert(temp_key_list19 != 0);
	 assert(temp_key_list20 != 0);
	 assert(temp_key_list21 != 0);

	 temp_key_list1 = memset(temp_key_list1, 0, sizeof(size_t)*list_size);
	 temp_key_list2 = memset(temp_key_list2, 0, sizeof(int)*list_size);
	 temp_key_list3 = memset(temp_key_list3, 0, sizeof(int)*list_size);
	 temp_key_list4 = memset(temp_key_list4, 0, sizeof(size_t)*list_size);
	 temp_key_list5 = memset(temp_key_list5, 0, sizeof(size_t)*list_size);
	 temp_key_list6 = memset(temp_key_list6, 0, sizeof(int)*list_size);
	 temp_key_list7 = memset(temp_key_list7, 0, sizeof(int)*list_size);
	 temp_key_list8 = memset(temp_key_list8, 0, sizeof(unsigned int)*list_size);
	 temp_key_list9 = memset(temp_key_list9, 0, sizeof(unsigned int)*list_size);
	 temp_key_list10 = memset(temp_key_list10, 0, sizeof(unsigned int)*list_size);
	 temp_key_list11 = memset(temp_key_list11, 0, sizeof(int)*list_size);
	 temp_key_list12 = memset(temp_key_list12, 0, sizeof(int)*list_size);
	 temp_key_list13 = memset(temp_key_list13, 0, sizeof(int)*list_size);
	 temp_key_list14 = memset(temp_key_list14, 0, sizeof(int)*list_size);
	 temp_key_list15 = memset(temp_key_list15, 0, sizeof(int)*list_size);
	 temp_key_list16 = memset(temp_key_list16, 0, sizeof(int)*list_size);
	 temp_key_list17 = memset(temp_key_list17, 0, sizeof(int)*list_size);
	 temp_key_list18 = memset(temp_key_list18, 0, sizeof(int)*list_size);
	 temp_key_list19 = memset(temp_key_list19, 0, sizeof(size_t)*list_size);
	 temp_key_list20 = memset(temp_key_list20, 0, sizeof(size_t)*list_size);
	 temp_key_list21 = memset(temp_key_list21, 0, sizeof(size_t)*list_size);



	 /*
	  * we pack the data into 1 vector called interbuff
	  *
	  */

	 MPI_Status status;
	 int res;
	 size_t *interbuff = NULL;
	 res = MPI_Alloc_mem((21*list_size)*sizeof(size_t), MPI_INFO_NULL, &interbuff);
	 assert(res == MPI_SUCCESS);
	 
	 for ( k = 0 ; k < list_size; k++ ){
		 interbuff[k] 				=  local_list1[k];
		 interbuff[k + list_size] 	=  (size_t)local_list2[k];
		 interbuff[k + 2*list_size] =  (size_t)local_list3[k];
		 interbuff[k + 3*list_size] =  local_list4[k];
		 interbuff[k + 4*list_size] =  local_list5[k];
		 interbuff[k + 5*list_size] =  (size_t)local_list6[k];
		 interbuff[k + 6*list_size] =  (size_t)local_list7[k];
		 interbuff[k + 7*list_size] =  (size_t)local_list8[k];
		 interbuff[k + 8*list_size] =  (size_t)local_list9[k];
		 interbuff[k + 9*list_size] =  (size_t)local_list10[k];
		 interbuff[k + 10*list_size] =  (size_t)local_list11[k];
		 interbuff[k + 11*list_size] =  (size_t)local_list12[k];
		 interbuff[k + 12*list_size] =  (size_t)local_list13[k];
		 interbuff[k + 13*list_size] =  (size_t)local_list14[k];
		 interbuff[k + 14*list_size] =  (size_t)local_list15[k];
		 interbuff[k + 15*list_size] =  (size_t)local_list16[k];
		 interbuff[k + 16*list_size] =  (size_t)local_list17[k];
		 interbuff[k + 17*list_size] =  (size_t)local_list18[k];
		 interbuff[k + 18*list_size] =  local_list19[k];
		 interbuff[k + 19*list_size] =  local_list20[k];
		 interbuff[k + 20*list_size] =  local_list21[k];

	 }

	 size_t *interbuff2;
	 res = MPI_Alloc_mem((21*list_size)*sizeof(size_t), MPI_INFO_NULL, &interbuff2);
	 assert(res == MPI_SUCCESS);


	 MPI_Sendrecv(interbuff,
			 	  21*list_size,
			 	  MPI_LONG_LONG_INT,
			 	  partner,
			 	  0,
			 	  interbuff2,
			 	  21*list_size,
			 	  MPI_LONG_LONG_INT,
			 	  partner,
			 	  0,
			 	  COMM_WORLD,
			 	  &status);


	 MPI_Get_count(&status, MPI_PACKED, &number_amount);

	 for ( k = 0 ; k < list_size; k++ ){
		 temp_key_list1[k] = (size_t)interbuff2[k];
		 temp_key_list2[k] = (int)   interbuff2[k + list_size];
		 temp_key_list3[k] = (int)   interbuff2[k + 2*list_size];
		 temp_key_list4[k] = (size_t)interbuff2[k + 3*list_size];
		 temp_key_list5[k] = (size_t)interbuff2[k + 4*list_size];
		 temp_key_list6[k] = (int)   interbuff2[k + 5*list_size];
		 temp_key_list7[k] = (int)   interbuff2[k + 6*list_size];
		 temp_key_list8[k] = (unsigned int)   interbuff2[k + 7*list_size];
		 temp_key_list9[k] = (unsigned int)   interbuff2[k + 8*list_size];
		 temp_key_list10[k] = (unsigned int)   interbuff2[k + 9*list_size];
		 temp_key_list11[k] = (int)   interbuff2[k + 10*list_size];
		 temp_key_list12[k] = (int)   interbuff2[k + 11*list_size];
		 temp_key_list13[k] = (int)   interbuff2[k + 12*list_size];
		 temp_key_list14[k] = (int)   interbuff2[k + 13*list_size];
		 temp_key_list15[k] = (int)   interbuff2[k + 14*list_size];
		 temp_key_list16[k] = (int)   interbuff2[k + 15*list_size];
		 temp_key_list17[k] = (int)   interbuff2[k + 16*list_size];
		 temp_key_list18[k] = (int)   interbuff2[k + 17*list_size];
		 temp_key_list19[k] = (size_t)   interbuff2[k + 18*list_size];
		 temp_key_list20[k] = (size_t)   interbuff2[k + 19*list_size];
		 temp_key_list21[k] = (size_t)   interbuff2[k + 20*list_size];

	 }


    if (which_keys == HIGH){
    	Merge_list_highAll2(
    			 list_size,
    			 local_list1,
    			 local_list2,
    			 local_list3,
    			 local_list4,
    			 local_list5,
				 local_list6,
				 local_list7,
				 local_list8,
				 local_list9,
				 local_list10,
				 local_list11,
				 local_list12,
				 local_list13,
				 local_list14, 
				 local_list15,
				 local_list16,
				 local_list17,
				 local_list18,
				 local_list19,
				 local_list20,
				 local_list21,
    			 temp_key_list1,
    			 temp_key_list2,
    			 temp_key_list3,
    			 temp_key_list4,
    			 temp_key_list5,
				 temp_key_list6,
				 temp_key_list7,
				 temp_key_list8,
				 temp_key_list9, 
				 temp_key_list10,
				 temp_key_list11,
				 temp_key_list12,
				 temp_key_list13,
				 temp_key_list14,
				 temp_key_list15,
				 temp_key_list16,
				 temp_key_list17,
				 temp_key_list18,
				 temp_key_list19,
				 temp_key_list20,
				 temp_key_list21
    			 );
    }
    else{
        Merge_list_lowAll2(
        		list_size,
    			 local_list1,
    			 local_list2,
    			 local_list3,
    			 local_list4,
    			 local_list5,
				 local_list6,
				 local_list7,
				 local_list8,
				 local_list9,
				 local_list10,
				 local_list11,
				 local_list12,
				 local_list13,
				 local_list14, 
				 local_list15,
				 local_list16,
				 local_list17,
				 local_list18,
				 local_list19,
				 local_list20,
				 local_list21,
    			 temp_key_list1,
    			 temp_key_list2,
    			 temp_key_list3,
    			 temp_key_list4,
    			 temp_key_list5,
				 temp_key_list6,
				 temp_key_list7,
				 temp_key_list8,
				 temp_key_list9, 
				 temp_key_list10,
				 temp_key_list11,
				 temp_key_list12,
				 temp_key_list13,
				 temp_key_list14,
				 temp_key_list15,
				 temp_key_list16,
				 temp_key_list17,
				 temp_key_list18,
				 temp_key_list19,
				 temp_key_list20,
				 temp_key_list21
        		);

    }


    if (temp_key_list1) free(temp_key_list1);
    if (temp_key_list2) free(temp_key_list2);
    if (temp_key_list3) free(temp_key_list3);
    if (temp_key_list4) free(temp_key_list4);
    if (temp_key_list5) free(temp_key_list5);
    if (temp_key_list6) free(temp_key_list6);
    if (temp_key_list7) free(temp_key_list7);
    if (temp_key_list8) free(temp_key_list8);
    if (temp_key_list9) free(temp_key_list9);
    if (temp_key_list10) free(temp_key_list10);
    if (temp_key_list11) free(temp_key_list11);
    if (temp_key_list12) free(temp_key_list12);
    if (temp_key_list13) free(temp_key_list13);
    if (temp_key_list14) free(temp_key_list14);
    if (temp_key_list15) free(temp_key_list15);
    if (temp_key_list16) free(temp_key_list16);
    if (temp_key_list17) free(temp_key_list17);
    if (temp_key_list18) free(temp_key_list18);
	if (temp_key_list19) free(temp_key_list19);
    if (temp_key_list20) free(temp_key_list20);
    if (temp_key_list21) free(temp_key_list21);

	fprintf(stderr,"TEST TEST\n");


	MPI_Free_mem(interbuff);
	MPI_Free_mem(interbuff2);
    malloc_trim(0);

} /* Merge_split */



/********************************************************************/
/* Merges the contents of the two lists. */
/* Returns the smaller keys in list1     */
void Merge_list_lowAll(
        size_t   list_size  	/* in     */,
        size_t  *list_key    	/* in/out */,
        int		*list_key1    	/* in/out */,
        int		*list_key2    	/* in/out */,
        size_t  *list_key3    	/* in/out */,
        int     *list_key4    	/* in/out */,
        int     *list_key5    	/* in/out */,
        unsigned int     *list_key6    	/* in/out */,
        unsigned int     *list_key7    	/* in/out */,
        unsigned int     *list_key8    	/* in/out */,
        int     *list_key9    	/* in/out */,
        int     *list_key10    	/* in/out */,
        int     *list_key11    	/* in/out */,
        int     *list_key12    	/* in/out */,
        int     *list_key13    	/* in/out */,
        int     *list_key14    	/* in/out */,
        int     *list_key15    	/* in/out */,
        int     *list_key16    	/* in/out */,
        size_t     *list_key17    	/* in/out */,
        size_t     *list_key18    	/* in/out */,
        size_t     *list_key19    	/* in/out */,
        size_t  *list_tmp_key   /* in     */,
        int		*list_tmp_key1   /* in     */,
        int		*list_tmp_key2   /* in     */,
        size_t  *list_tmp_key3   /* in     */,
        int     *list_tmp_key4   /* in     */,
        int     *list_tmp_key5   /* in     */,
        unsigned int     *list_tmp_key6   /* in     */,
        unsigned int     *list_tmp_key7   /* in     */,
        unsigned int     *list_tmp_key8   /* in     */,
        int     *list_tmp_key9   /* in     */,
        int     *list_tmp_key10   /* in     */,
        int     *list_tmp_key11   /* in     */,
        int     *list_tmp_key12  /* in     */,
        int     *list_tmp_key13   /* in     */,
        int     *list_tmp_key14   /* in     */,
        int     *list_tmp_key15   /* in     */,
        int     *list_tmp_key16   /* in     */,
        size_t     *list_tmp_key17   /* in     */,
        size_t    *list_tmp_key18   /* in     */,
        size_t     *list_tmp_key19   /* in     */
        ) {

	size_t  i;
    size_t  index1 = 0;
    size_t  index2 = 0;

    size_t *scratch_list_key  = malloc(list_size*sizeof(size_t));
    int	   *scratch_list_key1 = malloc(list_size*sizeof(size_t));
    int    *scratch_list_key2 = malloc(list_size*sizeof(size_t));
    size_t *scratch_list_key3 = malloc(list_size*sizeof(size_t));
    int    *scratch_list_key4 = malloc(list_size*sizeof(int));
    int    *scratch_list_key5 = malloc(list_size*sizeof(int));
    unsigned int    *scratch_list_key6 = malloc(list_size*sizeof(unsigned int));
    unsigned int    *scratch_list_key7 = malloc(list_size*sizeof(unsigned int));
    unsigned int    *scratch_list_key8 = malloc(list_size*sizeof(unsigned int));
    int    *scratch_list_key9 = malloc(list_size*sizeof(int));
    int    *scratch_list_key10 = malloc(list_size*sizeof(int));
    int    *scratch_list_key11 = malloc(list_size*sizeof(int));
    int    *scratch_list_key12 = malloc(list_size*sizeof(int));
    int    *scratch_list_key13 = malloc(list_size*sizeof(int));
    int    *scratch_list_key14 = malloc(list_size*sizeof(int));
    int    *scratch_list_key15 = malloc(list_size*sizeof(int));
    int    *scratch_list_key16 = malloc(list_size*sizeof(int));
    size_t    *scratch_list_key17 = malloc(list_size*sizeof(size_t));
    size_t    *scratch_list_key18 = malloc(list_size*sizeof(size_t));
    size_t    *scratch_list_key19 = malloc(list_size*sizeof(size_t));



    scratch_list_key[0]  = 0;
    scratch_list_key1[0] = 0;
    scratch_list_key2[0] = 0;
    scratch_list_key3[0] = 0;
    scratch_list_key4[0] = 0;
    scratch_list_key5[0] = 0;
    scratch_list_key6[0] = 0;
    scratch_list_key7[0] = 0;
    scratch_list_key8[0] = 0;
    scratch_list_key9[0] = 0;
    scratch_list_key10[0] = 0;
    scratch_list_key11[0] = 0;
    scratch_list_key12[0] = 0;
    scratch_list_key13[0] = 0;
    scratch_list_key14[0] = 0;
    scratch_list_key15[0] = 0;
    scratch_list_key16[0] = 0;
    scratch_list_key17[0] = 0;
    scratch_list_key18[0] = 0;
    scratch_list_key19[0] = 0;


    for (i = 0; i < list_size; i++){
        if (list_key[index1] <= list_tmp_key[index2]) {

        	scratch_list_key[i]  = list_key[index1];
        	scratch_list_key1[i] = list_key1[index1];
        	scratch_list_key2[i] = list_key2[index1];
        	scratch_list_key3[i] = list_key3[index1];
        	scratch_list_key4[i] = list_key4[index1];
        	scratch_list_key5[i] = list_key5[index1];
        	scratch_list_key6[i] = list_key6[index1];
        	scratch_list_key7[i] = list_key7[index1];
        	scratch_list_key8[i] = list_key8[index1];
        	scratch_list_key9[i] = list_key9[index1];
        	scratch_list_key10[i] = list_key10[index1];
        	scratch_list_key11[i] = list_key11[index1];
        	scratch_list_key12[i] = list_key12[index1];
        	scratch_list_key13[i] = list_key13[index1];
        	scratch_list_key14[i] = list_key14[index1];
        	scratch_list_key15[i] = list_key15[index1];
        	scratch_list_key16[i] = list_key16[index1];
        	scratch_list_key17[i] = list_key17[index1];
        	scratch_list_key18[i] = list_key18[index1];
        	scratch_list_key19[i] = list_key19[index1];
            index1++;

        } else {

        	scratch_list_key[i]  = list_tmp_key[index2];
        	scratch_list_key1[i] = list_tmp_key1[index2];
        	scratch_list_key2[i] = list_tmp_key2[index2];
        	scratch_list_key3[i] = list_tmp_key3[index2];
        	scratch_list_key4[i] = list_tmp_key4[index2];
        	scratch_list_key5[i] = list_tmp_key5[index2];
        	scratch_list_key6[i] = list_tmp_key6[index2];
        	scratch_list_key7[i] = list_tmp_key7[index2];
        	scratch_list_key8[i] = list_tmp_key8[index2];
        	scratch_list_key9[i] = list_tmp_key9[index2];
        	scratch_list_key10[i] = list_tmp_key10[index2];
        	scratch_list_key11[i] = list_tmp_key11[index2];
        	scratch_list_key12[i] = list_tmp_key12[index2];
        	scratch_list_key13[i] = list_tmp_key13[index2];
        	scratch_list_key14[i] = list_tmp_key14[index2];
        	scratch_list_key15[i] = list_tmp_key15[index2];
        	scratch_list_key16[i] = list_tmp_key16[index2];
        	scratch_list_key17[i] = list_tmp_key17[index2];
        	scratch_list_key18[i] = list_tmp_key18[index2];
        	scratch_list_key19[i] = list_tmp_key19[index2];
            index2++;
        }
    }
    for (i = 0; i < list_size; i++){
    	list_key[i] = scratch_list_key[i];
    	list_key1[i] = scratch_list_key1[i];
    	list_key2[i] = scratch_list_key2[i];
    	list_key3[i] = scratch_list_key3[i];
    	list_key4[i] = scratch_list_key4[i];
    	list_key5[i] = scratch_list_key5[i];
    	list_key6[i] = scratch_list_key6[i];
    	list_key7[i] = scratch_list_key7[i];
    	list_key8[i] = scratch_list_key8[i];
    	list_key9[i] = scratch_list_key9[i];
    	list_key10[i] = scratch_list_key10[i];
    	list_key11[i] = scratch_list_key11[i];
    	list_key12[i] = scratch_list_key12[i];
    	list_key13[i] = scratch_list_key13[i];
    	list_key14[i] = scratch_list_key14[i];
    	list_key15[i] = scratch_list_key15[i];
    	list_key16[i] = scratch_list_key16[i];
    	list_key17[i] = scratch_list_key17[i];
    	list_key18[i] = scratch_list_key18[i];
    	list_key19[i] = scratch_list_key19[i];

    }

    free(scratch_list_key);
    free(scratch_list_key1);
    free(scratch_list_key2);
    free(scratch_list_key3);
    free(scratch_list_key4);
    free(scratch_list_key5);
    free(scratch_list_key6);
    free(scratch_list_key7);
    free(scratch_list_key8);
    free(scratch_list_key9);
    free(scratch_list_key10);
    free(scratch_list_key11);
    free(scratch_list_key12);
    free(scratch_list_key13);
    free(scratch_list_key14);
    free(scratch_list_key15);
    free(scratch_list_key16);
    free(scratch_list_key17);
    free(scratch_list_key18);
    free(scratch_list_key19);

    //malloc_trim(0);
}  /* Merge_list_low */

void Merge_list_lowAll2(
        size_t   list_size  	/* in     */,
        size_t  *list_key    	/* in/out */,
        int		*list_key1    	/* in/out */,
        int		*list_key2    	/* in/out */,
        size_t  *list_key3    	/* in/out */,
        size_t  *list_key4    	/* in/out */,
        int     *list_key5    	/* in/out */,
        int     *list_key6    	/* in/out */,
        unsigned int     *list_key7    	/* in/out */,
        unsigned int     *list_key8    	/* in/out */,
        unsigned int     *list_key9    	/* in/out */,
        int     *list_key10    	/* in/out */,
        int     *list_key11    	/* in/out */,
        int     *list_key12    	/* in/out */,
        int     *list_key13    	/* in/out */,
        int     *list_key14    	/* in/out */,
        int     *list_key15    	/* in/out */,
        int     *list_key16    	/* in/out */,
        int     *list_key17    	/* in/out */,
        size_t     *list_key18    	/* in/out */,
        size_t     *list_key19    	/* in/out */,
        size_t     *list_key20    	/* in/out */,
        size_t  *list_tmp_key   /* in     */,
        int		*list_tmp_key1   /* in     */,
        int		*list_tmp_key2   /* in     */,
        size_t  *list_tmp_key3   /* in     */,
        size_t  *list_tmp_key4   /* in     */,
        int     *list_tmp_key5   /* in     */,
        int     *list_tmp_key6   /* in     */,
        unsigned int     *list_tmp_key7   /* in     */,
        unsigned int     *list_tmp_key8   /* in     */,
        unsigned int     *list_tmp_key9   /* in     */,
        int     *list_tmp_key10   /* in     */,
        int     *list_tmp_key11   /* in     */,
        int     *list_tmp_key12   /* in     */,
        int     *list_tmp_key13  /* in     */,
        int     *list_tmp_key14   /* in     */,
        int     *list_tmp_key15   /* in     */,
        int     *list_tmp_key16   /* in     */,
        int     *list_tmp_key17   /* in     */,
        size_t     *list_tmp_key18   /* in     */,
        size_t    *list_tmp_key19   /* in     */,
        size_t     *list_tmp_key20   /* in     */
        ) {

	size_t  i;
    size_t  index1 = 0;
    size_t  index2 = 0;

    size_t *scratch_list_key  = malloc(list_size*sizeof(size_t));
    int	   *scratch_list_key1 = malloc(list_size*sizeof(size_t));
    int    *scratch_list_key2 = malloc(list_size*sizeof(size_t));
    size_t *scratch_list_key3 = malloc(list_size*sizeof(size_t));
    size_t *scratch_list_key4 = malloc(list_size*sizeof(size_t));
    int    *scratch_list_key5 = malloc(list_size*sizeof(int));
    int    *scratch_list_key6 = malloc(list_size*sizeof(int));
    unsigned int    *scratch_list_key7 = malloc(list_size*sizeof(unsigned int));
    unsigned int    *scratch_list_key8 = malloc(list_size*sizeof(unsigned int));
    unsigned int    *scratch_list_key9 = malloc(list_size*sizeof(unsigned int));
    int    *scratch_list_key10 = malloc(list_size*sizeof(int));
    int    *scratch_list_key11 = malloc(list_size*sizeof(int));
    int    *scratch_list_key12 = malloc(list_size*sizeof(int));
    int    *scratch_list_key13 = malloc(list_size*sizeof(int));
    int    *scratch_list_key14 = malloc(list_size*sizeof(int));
    int    *scratch_list_key15 = malloc(list_size*sizeof(int));
    int    *scratch_list_key16 = malloc(list_size*sizeof(int));
    int    *scratch_list_key17 = malloc(list_size*sizeof(int));
    size_t    *scratch_list_key18 = malloc(list_size*sizeof(size_t));
    size_t    *scratch_list_key19 = malloc(list_size*sizeof(size_t));
    size_t    *scratch_list_key20 = malloc(list_size*sizeof(size_t));



    scratch_list_key[0]  = 0;
    scratch_list_key1[0] = 0;
    scratch_list_key2[0] = 0;
    scratch_list_key3[0] = 0;
    scratch_list_key4[0] = 0;
    scratch_list_key5[0] = 0;
    scratch_list_key6[0] = 0;
    scratch_list_key7[0] = 0;
    scratch_list_key8[0] = 0;
    scratch_list_key9[0] = 0;
    scratch_list_key10[0] = 0;
    scratch_list_key11[0] = 0;
    scratch_list_key12[0] = 0;
    scratch_list_key13[0] = 0;
    scratch_list_key14[0] = 0;
    scratch_list_key15[0] = 0;
    scratch_list_key16[0] = 0;
    scratch_list_key17[0] = 0;
    scratch_list_key18[0] = 0;
    scratch_list_key19[0] = 0;
    scratch_list_key20[0] = 0;


    for (i = 0; i < list_size; i++){
        if (list_key[index1] <= list_tmp_key[index2]) {

        	scratch_list_key[i]  = list_key[index1];
        	scratch_list_key1[i] = list_key1[index1];
        	scratch_list_key2[i] = list_key2[index1];
        	scratch_list_key3[i] = list_key3[index1];
        	scratch_list_key4[i] = list_key4[index1];
        	scratch_list_key5[i] = list_key5[index1];
        	scratch_list_key6[i] = list_key6[index1];
        	scratch_list_key7[i] = list_key7[index1];
        	scratch_list_key8[i] = list_key8[index1];
        	scratch_list_key9[i] = list_key9[index1];
        	scratch_list_key10[i] = list_key10[index1];
        	scratch_list_key11[i] = list_key11[index1];
        	scratch_list_key12[i] = list_key12[index1];
        	scratch_list_key13[i] = list_key13[index1];
        	scratch_list_key14[i] = list_key14[index1];
        	scratch_list_key15[i] = list_key15[index1];
        	scratch_list_key16[i] = list_key16[index1];
        	scratch_list_key17[i] = list_key17[index1];
        	scratch_list_key18[i] = list_key18[index1];
        	scratch_list_key19[i] = list_key19[index1];
        	scratch_list_key20[i] = list_key20[index1];
            index1++;

        } else {

        	scratch_list_key[i]  = list_tmp_key[index2];
        	scratch_list_key1[i] = list_tmp_key1[index2];
        	scratch_list_key2[i] = list_tmp_key2[index2];
        	scratch_list_key3[i] = list_tmp_key3[index2];
        	scratch_list_key4[i] = list_tmp_key4[index2];
        	scratch_list_key5[i] = list_tmp_key5[index2];
        	scratch_list_key6[i] = list_tmp_key6[index2];
        	scratch_list_key7[i] = list_tmp_key7[index2];
        	scratch_list_key8[i] = list_tmp_key8[index2];
        	scratch_list_key9[i] = list_tmp_key9[index2];
        	scratch_list_key10[i] = list_tmp_key10[index2];
        	scratch_list_key11[i] = list_tmp_key11[index2];
        	scratch_list_key12[i] = list_tmp_key12[index2];
        	scratch_list_key13[i] = list_tmp_key13[index2];
        	scratch_list_key14[i] = list_tmp_key14[index2];
        	scratch_list_key15[i] = list_tmp_key15[index2];
        	scratch_list_key16[i] = list_tmp_key16[index2];
        	scratch_list_key17[i] = list_tmp_key17[index2];
        	scratch_list_key18[i] = list_tmp_key18[index2];
        	scratch_list_key19[i] = list_tmp_key19[index2];
        	scratch_list_key20[i] = list_tmp_key20[index2];
            index2++;
        }
    }
    for (i = 0; i < list_size; i++){
    	list_key[i] = scratch_list_key[i];
    	list_key1[i] = scratch_list_key1[i];
    	list_key2[i] = scratch_list_key2[i];
    	list_key3[i] = scratch_list_key3[i];
    	list_key4[i] = scratch_list_key4[i];
    	list_key5[i] = scratch_list_key5[i];
    	list_key6[i] = scratch_list_key6[i];
    	list_key7[i] = scratch_list_key7[i];
    	list_key8[i] = scratch_list_key8[i];
    	list_key9[i] = scratch_list_key9[i];
    	list_key10[i] = scratch_list_key10[i];
    	list_key11[i] = scratch_list_key11[i];
    	list_key12[i] = scratch_list_key12[i];
    	list_key13[i] = scratch_list_key13[i];
    	list_key14[i] = scratch_list_key14[i];
    	list_key15[i] = scratch_list_key15[i];
    	list_key16[i] = scratch_list_key16[i];
    	list_key17[i] = scratch_list_key17[i];
    	list_key18[i] = scratch_list_key18[i];
    	list_key19[i] = scratch_list_key19[i];
    	list_key20[i] = scratch_list_key20[i];

    }
	
    free(scratch_list_key);
    free(scratch_list_key1);
    free(scratch_list_key2);
    free(scratch_list_key3);
    free(scratch_list_key4);
    free(scratch_list_key5);
    free(scratch_list_key6);
    free(scratch_list_key7);
    free(scratch_list_key8);
    free(scratch_list_key9);
    free(scratch_list_key10);
    free(scratch_list_key11);
    free(scratch_list_key12);
    free(scratch_list_key13);
    free(scratch_list_key14);
    free(scratch_list_key15);
    free(scratch_list_key16);
    free(scratch_list_key17);
    free(scratch_list_key18);
    free(scratch_list_key19);
	free(scratch_list_key20);

    //malloc_trim(0);
}

void Merge_list_low2(
        size_t   list_size  	/* in     */,
        size_t  *list_key    	/* in/out */,
        int		*list_key1    	/* in/out */,
        int		*list_key2    	/* in/out */,
        size_t  *list_key3    	/* in/out */,
        int     *list_key4    	/* in/out */,
        size_t  *list_tmp_key   /* in     */,
        int		*list_tmp_key1   /* in     */,
        int		*list_tmp_key2   /* in     */,
        size_t  *list_tmp_key3   /* in     */,
        int     *list_tmp_key4   /* in     */
        ) {

	size_t  i;
    size_t  index1 = 0;
    size_t  index2 = 0;

    size_t *scratch_list_key  = malloc(list_size*sizeof(size_t));
    int	   *scratch_list_key1 = malloc(list_size*sizeof(size_t));
    int    *scratch_list_key2 = malloc(list_size*sizeof(size_t));
    size_t *scratch_list_key3 = malloc(list_size*sizeof(size_t));
    int    *scratch_list_key4 = malloc(list_size*sizeof(int));

    scratch_list_key[0]  = 0;
    scratch_list_key1[0] = 0;
    scratch_list_key2[0] = 0;
    scratch_list_key3[0] = 0;
    scratch_list_key4[0] = 0;

    for (i = 0; i < list_size; i++){
        if (list_key[index1] <= list_tmp_key[index2]) {

        	scratch_list_key[i]  = list_key[index1];
        	scratch_list_key1[i] = list_key1[index1];
        	scratch_list_key2[i] = list_key2[index1];
        	scratch_list_key3[i] = list_key3[index1];
        	scratch_list_key4[i] = list_key4[index1];
            index1++;

        } else {

        	scratch_list_key[i]  = list_tmp_key[index2];
        	scratch_list_key1[i] = list_tmp_key1[index2];
        	scratch_list_key2[i] = list_tmp_key2[index2];
        	scratch_list_key3[i] = list_tmp_key3[index2];
        	scratch_list_key4[i] = list_tmp_key4[index2];
            index2++;
        }
    }
    for (i = 0; i < list_size; i++){
    	list_key[i] = scratch_list_key[i];
    	list_key1[i] = scratch_list_key1[i];
    	list_key2[i] = scratch_list_key2[i];
    	list_key3[i] = scratch_list_key3[i];
    	list_key4[i] = scratch_list_key4[i];

    }

    free(scratch_list_key);
    free(scratch_list_key1);
    free(scratch_list_key2);
    free(scratch_list_key3);
    free(scratch_list_key4);

    //malloc_trim(0);
}  /* Merge_list_low */




/********************************************************************/
/* Returns the larger keys in list 1.    */
void Merge_list_high2(
		 size_t   list_size  	/* in     */,
		 size_t  *list_key    	/* in/out */,
		 int	 *list_key1    	/* in/out */,
		 int	 *list_key2    	/* in/out */,
		 size_t  *list_key3    	/* in/out */,
		 int     *list_key4    	/* in/out */,
		 size_t  *list_tmp_key   /* in     */,
		 int	 *list_tmp_key1   /* in     */,
		 int     *list_tmp_key2   /* in     */,
		 size_t  *list_tmp_key3   /* in     */,
		 int     *list_tmp_key4   /* in     */
		 ) {

    size_t  i;
    size_t  index1 = list_size - 1;
    size_t  index2 = list_size - 1;

    size_t  *scratch_list_key  = malloc(list_size*sizeof(size_t));
    int 	*scratch_list_key1 = malloc(list_size*sizeof(int));
    int 	*scratch_list_key2 = malloc(list_size*sizeof(int));
    size_t  *scratch_list_key3 = malloc(list_size*sizeof(size_t));
    int 	*scratch_list_key4 = malloc(list_size*sizeof(int));

    scratch_list_key[0] =0;
    scratch_list_key1[0]=0;
    scratch_list_key2[0]=0;
    scratch_list_key3[0]=0;
    scratch_list_key4[0]=0;

    size_t counter =0;
    int  rank;
    MPI_Comm_rank(COMM_WORLD, &rank);
    for (i = list_size - 1;; i--){

        if (list_key[index1] >= list_tmp_key[index2]) {

        	scratch_list_key[i]  = list_key[index1];
        	scratch_list_key1[i] = list_key1[index1];
        	scratch_list_key2[i] = list_key2[index1];
        	scratch_list_key3[i] = list_key3[index1];
        	scratch_list_key4[i] = list_key4[index1];
        	index1--;
        	counter++;

        } else {

        	scratch_list_key[i]  = list_tmp_key[index2];
        	scratch_list_key1[i] = list_tmp_key1[index2];
        	scratch_list_key2[i] = list_tmp_key2[index2];
        	scratch_list_key3[i] = list_tmp_key3[index2];
        	scratch_list_key4[i] = list_tmp_key4[index2];
            index2--;
            counter++;
        }

        if (counter >= list_size)
            break;

        if (i == 0) break;
    }

    for (i = 0; i < list_size; i++){

        	list_key[i]   = scratch_list_key[i];
        	list_key1[i]  = scratch_list_key1[i];
        	list_key2[i]  = scratch_list_key2[i];
        	list_key3[i]  = scratch_list_key3[i];
        	list_key4[i]  = scratch_list_key4[i];
    }

    free(scratch_list_key);
    free(scratch_list_key1);
    free(scratch_list_key2);
    free(scratch_list_key3);
    free(scratch_list_key4);
}  /* Merge_list _high */

void Merge_list_highAll(size_t   list_size  	/* in     */,
        size_t  *list_key    	/* in/out */,
        int		*list_key1    	/* in/out */,
        int		*list_key2    	/* in/out */,
        size_t  *list_key3    	/* in/out */,
        int     *list_key4    	/* in/out */,
        int     *list_key5    	/* in/out */,
        unsigned int     *list_key6    	/* in/out */,
        unsigned int     *list_key7    	/* in/out */,
        unsigned int     *list_key8    	/* in/out */,
        int     *list_key9    	/* in/out */,
        int     *list_key10    	/* in/out */,
        int     *list_key11    	/* in/out */,
        int     *list_key12    	/* in/out */,
        int     *list_key13    	/* in/out */,
        int     *list_key14    	/* in/out */,
        int     *list_key15    	/* in/out */,
        int     *list_key16    	/* in/out */,
        size_t     *list_key17    	/* in/out */,
        size_t     *list_key18    	/* in/out */,
        size_t     *list_key19    	/* in/out */,
        size_t  *list_tmp_key   /* in     */,
        int		*list_tmp_key1   /* in     */,
        int		*list_tmp_key2   /* in     */,
        size_t  *list_tmp_key3   /* in     */,
        int     *list_tmp_key4   /* in     */,
        int     *list_tmp_key5   /* in     */,
        unsigned int     *list_tmp_key6   /* in     */,
        unsigned int     *list_tmp_key7   /* in     */,
        unsigned int     *list_tmp_key8   /* in     */,
        int     *list_tmp_key9   /* in     */,
        int     *list_tmp_key10   /* in     */,
        int     *list_tmp_key11   /* in     */,
        int     *list_tmp_key12  /* in     */,
        int     *list_tmp_key13   /* in     */,
        int     *list_tmp_key14   /* in     */,
        int     *list_tmp_key15   /* in     */,
        int     *list_tmp_key16   /* in     */,
        size_t     *list_tmp_key17   /* in     */,
        size_t    *list_tmp_key18   /* in     */,
        size_t     *list_tmp_key19   /* in     */
) {

    size_t  i;
    size_t  index1 = list_size - 1;
    size_t  index2 = list_size - 1;

    size_t  *scratch_list_key  = malloc(list_size*sizeof(size_t));
    int 	*scratch_list_key1 = malloc(list_size*sizeof(int));
    int 	*scratch_list_key2 = malloc(list_size*sizeof(int));
    size_t  *scratch_list_key3 = malloc(list_size*sizeof(size_t));
    int 	*scratch_list_key4 = malloc(list_size*sizeof(int));
	int    *scratch_list_key5 = malloc(list_size*sizeof(int));
    unsigned int    *scratch_list_key6 = malloc(list_size*sizeof(unsigned int));
    unsigned int    *scratch_list_key7 = malloc(list_size*sizeof(unsigned int));
    unsigned int    *scratch_list_key8 = malloc(list_size*sizeof(unsigned int));
    int    *scratch_list_key9 = malloc(list_size*sizeof(int));
    int    *scratch_list_key10 = malloc(list_size*sizeof(int));
    int    *scratch_list_key11 = malloc(list_size*sizeof(int));
    int    *scratch_list_key12 = malloc(list_size*sizeof(int));
    int    *scratch_list_key13 = malloc(list_size*sizeof(int));
    int    *scratch_list_key14 = malloc(list_size*sizeof(int));
    int    *scratch_list_key15 = malloc(list_size*sizeof(int));
    int    *scratch_list_key16 = malloc(list_size*sizeof(int));
    size_t    *scratch_list_key17 = malloc(list_size*sizeof(size_t));
    size_t    *scratch_list_key18 = malloc(list_size*sizeof(size_t));
    size_t    *scratch_list_key19 = malloc(list_size*sizeof(size_t));


    scratch_list_key[0] =0;
    scratch_list_key1[0]=0;
    scratch_list_key2[0]=0;
    scratch_list_key3[0]=0;
    scratch_list_key4[0]=0;
	scratch_list_key5[0] = 0;
    scratch_list_key6[0] = 0;
    scratch_list_key7[0] = 0;
    scratch_list_key8[0] = 0;
    scratch_list_key9[0] = 0;
    scratch_list_key10[0] = 0;
    scratch_list_key11[0] = 0;
    scratch_list_key12[0] = 0;
    scratch_list_key13[0] = 0;
    scratch_list_key14[0] = 0;
    scratch_list_key15[0] = 0;
    scratch_list_key16[0] = 0;
    scratch_list_key17[0] = 0;
    scratch_list_key18[0] = 0;
    scratch_list_key19[0] = 0;

    size_t counter =0;
    int  rank;
    MPI_Comm_rank(COMM_WORLD, &rank);
    for (i = list_size - 1;; i--){

        if (list_key[index1] >= list_tmp_key[index2]) {

        	scratch_list_key[i]  = list_key[index1];
        	scratch_list_key1[i] = list_key1[index1];
        	scratch_list_key2[i] = list_key2[index1];
        	scratch_list_key3[i] = list_key3[index1];
        	scratch_list_key4[i] = list_key4[index1];
			scratch_list_key5[i] = list_key5[index1];
        	scratch_list_key6[i] = list_key6[index1];
        	scratch_list_key7[i] = list_key7[index1];
        	scratch_list_key8[i] = list_key8[index1];
        	scratch_list_key9[i] = list_key9[index1];
        	scratch_list_key10[i] = list_key10[index1];
        	scratch_list_key11[i] = list_key11[index1];
        	scratch_list_key12[i] = list_key12[index1];
        	scratch_list_key13[i] = list_key13[index1];
        	scratch_list_key14[i] = list_key14[index1];
        	scratch_list_key15[i] = list_key15[index1];
        	scratch_list_key16[i] = list_key16[index1];
        	scratch_list_key17[i] = list_key17[index1];
        	scratch_list_key18[i] = list_key18[index1];
        	scratch_list_key19[i] = list_key19[index1];
        	index1--;
        	counter++;

        } else {

        	scratch_list_key[i]  = list_tmp_key[index2];
        	scratch_list_key1[i] = list_tmp_key1[index2];
        	scratch_list_key2[i] = list_tmp_key2[index2];
        	scratch_list_key3[i] = list_tmp_key3[index2];
        	scratch_list_key4[i] = list_tmp_key4[index2];
			scratch_list_key5[i] = list_tmp_key5[index2];
        	scratch_list_key6[i] = list_tmp_key6[index2];
        	scratch_list_key7[i] = list_tmp_key7[index2];
        	scratch_list_key8[i] = list_tmp_key8[index2];
        	scratch_list_key9[i] = list_tmp_key9[index2];
        	scratch_list_key10[i] = list_tmp_key10[index2];
        	scratch_list_key11[i] = list_tmp_key11[index2];
        	scratch_list_key12[i] = list_tmp_key12[index2];
        	scratch_list_key13[i] = list_tmp_key13[index2];
        	scratch_list_key14[i] = list_tmp_key14[index2];
        	scratch_list_key15[i] = list_tmp_key15[index2];
        	scratch_list_key16[i] = list_tmp_key16[index2];
        	scratch_list_key17[i] = list_tmp_key17[index2];
        	scratch_list_key18[i] = list_tmp_key18[index2];
        	scratch_list_key19[i] = list_tmp_key19[index2];
            index2--;
            counter++;
        }

        if (counter >= list_size)
            break;

        if (i == 0) break;
    }

    for (i = 0; i < list_size; i++){

        	list_key[i]   = scratch_list_key[i];
        	list_key1[i]  = scratch_list_key1[i];
        	list_key2[i]  = scratch_list_key2[i];
        	list_key3[i]  = scratch_list_key3[i];
        	list_key4[i]  = scratch_list_key4[i];
			list_key5[i] = scratch_list_key5[i];
			list_key6[i] = scratch_list_key6[i];
			list_key7[i] = scratch_list_key7[i];
			list_key8[i] = scratch_list_key8[i];
			list_key9[i] = scratch_list_key9[i];
			list_key10[i] = scratch_list_key10[i];
			list_key11[i] = scratch_list_key11[i];
			list_key12[i] = scratch_list_key12[i];
			list_key13[i] = scratch_list_key13[i];
			list_key14[i] = scratch_list_key14[i];
			list_key15[i] = scratch_list_key15[i];
			list_key16[i] = scratch_list_key16[i];
			list_key17[i] = scratch_list_key17[i];	
			list_key18[i] = scratch_list_key18[i];	
			list_key19[i] = scratch_list_key19[i];	
    }

    free(scratch_list_key);
    free(scratch_list_key1);
    free(scratch_list_key2);
    free(scratch_list_key3);
    free(scratch_list_key4);
	free(scratch_list_key5);
    free(scratch_list_key6);
    free(scratch_list_key7);
    free(scratch_list_key8);
    free(scratch_list_key9);
    free(scratch_list_key10);
    free(scratch_list_key11);
    free(scratch_list_key12);
    free(scratch_list_key13);
    free(scratch_list_key14);
    free(scratch_list_key15);
    free(scratch_list_key16);
    free(scratch_list_key17);
    free(scratch_list_key18);
    free(scratch_list_key19);
}  /* Merge_list _high */


void Merge_list_highAll2(size_t   list_size  	/* in     */,
        size_t  *list_key    	/* in/out */,
        int		*list_key1    	/* in/out */,
        int		*list_key2    	/* in/out */,
        size_t  *list_key3    	/* in/out */,
        size_t  *list_key4    	/* in/out */,
        int     *list_key5    	/* in/out */,
        int     *list_key6    	/* in/out */,
        unsigned int     *list_key7    	/* in/out */,
        unsigned int     *list_key8    	/* in/out */,
        unsigned int     *list_key9    	/* in/out */,
        int     *list_key10    	/* in/out */,
        int     *list_key11    	/* in/out */,
        int     *list_key12    	/* in/out */,
        int     *list_key13    	/* in/out */,
        int     *list_key14    	/* in/out */,
        int     *list_key15    	/* in/out */,
        int     *list_key16    	/* in/out */,
        int     *list_key17    	/* in/out */,
        size_t     *list_key18    	/* in/out */,
        size_t     *list_key19    	/* in/out */,
        size_t     *list_key20    	/* in/out */,
        size_t  *list_tmp_key   /* in     */,
        int		*list_tmp_key1   /* in     */,
        int		*list_tmp_key2   /* in     */,
        size_t  *list_tmp_key3   /* in     */,
        size_t  *list_tmp_key4   /* in     */,
        int     *list_tmp_key5   /* in     */,
        int     *list_tmp_key6   /* in     */,
        unsigned int     *list_tmp_key7   /* in     */,
        unsigned int     *list_tmp_key8   /* in     */,
        unsigned int     *list_tmp_key9   /* in     */,
        int     *list_tmp_key10   /* in     */,
        int     *list_tmp_key11   /* in     */,
        int     *list_tmp_key12   /* in     */,
        int     *list_tmp_key13  /* in     */,
        int     *list_tmp_key14   /* in     */,
        int     *list_tmp_key15   /* in     */,
        int     *list_tmp_key16   /* in     */,
        int     *list_tmp_key17   /* in     */,
        size_t     *list_tmp_key18   /* in     */,
        size_t    *list_tmp_key19   /* in     */,
        size_t     *list_tmp_key20   /* in     */
) {

    size_t  i;
    size_t  index1 = list_size - 1;
    size_t  index2 = list_size - 1;

    size_t  *scratch_list_key  = malloc(list_size*sizeof(size_t));
    int 	*scratch_list_key1 = malloc(list_size*sizeof(int));
    int 	*scratch_list_key2 = malloc(list_size*sizeof(int));
    size_t  *scratch_list_key3 = malloc(list_size*sizeof(size_t));
    size_t  *scratch_list_key4 = malloc(list_size*sizeof(size_t));
    int 	*scratch_list_key5 = malloc(list_size*sizeof(int));
	int    *scratch_list_key6 = malloc(list_size*sizeof(int));
    unsigned int    *scratch_list_key7 = malloc(list_size*sizeof(unsigned int));
    unsigned int    *scratch_list_key8 = malloc(list_size*sizeof(unsigned int));
    unsigned int    *scratch_list_key9 = malloc(list_size*sizeof(unsigned int));
    int    *scratch_list_key10 = malloc(list_size*sizeof(int));
    int    *scratch_list_key11 = malloc(list_size*sizeof(int));
    int    *scratch_list_key12 = malloc(list_size*sizeof(int));
    int    *scratch_list_key13 = malloc(list_size*sizeof(int));
    int    *scratch_list_key14 = malloc(list_size*sizeof(int));
    int    *scratch_list_key15 = malloc(list_size*sizeof(int));
    int    *scratch_list_key16 = malloc(list_size*sizeof(int));
    int    *scratch_list_key17 = malloc(list_size*sizeof(int));
    size_t    *scratch_list_key18 = malloc(list_size*sizeof(size_t));
    size_t    *scratch_list_key19 = malloc(list_size*sizeof(size_t));
    size_t    *scratch_list_key20 = malloc(list_size*sizeof(size_t));


    scratch_list_key[0] =0;
    scratch_list_key1[0]=0;
    scratch_list_key2[0]=0;
    scratch_list_key3[0]=0;
    scratch_list_key4[0]=0;
	scratch_list_key5[0] = 0;
    scratch_list_key6[0] = 0;
    scratch_list_key7[0] = 0;
    scratch_list_key8[0] = 0;
    scratch_list_key9[0] = 0;
    scratch_list_key10[0] = 0;
    scratch_list_key11[0] = 0;
    scratch_list_key12[0] = 0;
    scratch_list_key13[0] = 0;
    scratch_list_key14[0] = 0;
    scratch_list_key15[0] = 0;
    scratch_list_key16[0] = 0;
    scratch_list_key17[0] = 0;
    scratch_list_key18[0] = 0;
    scratch_list_key19[0] = 0;
    scratch_list_key20[0] = 0;

    size_t counter =0;
    int  rank;
    MPI_Comm_rank(COMM_WORLD, &rank);
    for (i = list_size - 1;; i--){

        if (list_key[index1] >= list_tmp_key[index2]) {

        	scratch_list_key[i]  = list_key[index1];
        	scratch_list_key1[i] = list_key1[index1];
        	scratch_list_key2[i] = list_key2[index1];
        	scratch_list_key3[i] = list_key3[index1];
        	scratch_list_key4[i] = list_key4[index1];
			scratch_list_key5[i] = list_key5[index1];
        	scratch_list_key6[i] = list_key6[index1];
        	scratch_list_key7[i] = list_key7[index1];
        	scratch_list_key8[i] = list_key8[index1];
        	scratch_list_key9[i] = list_key9[index1];
        	scratch_list_key10[i] = list_key10[index1];
        	scratch_list_key11[i] = list_key11[index1];
        	scratch_list_key12[i] = list_key12[index1];
        	scratch_list_key13[i] = list_key13[index1];
        	scratch_list_key14[i] = list_key14[index1];
        	scratch_list_key15[i] = list_key15[index1];
        	scratch_list_key16[i] = list_key16[index1];
        	scratch_list_key17[i] = list_key17[index1];
        	scratch_list_key18[i] = list_key18[index1];
        	scratch_list_key19[i] = list_key19[index1];
        	scratch_list_key20[i] = list_key20[index1];
        	index1--;
        	counter++;

        } else {

        	scratch_list_key[i]  = list_tmp_key[index2];
        	scratch_list_key1[i] = list_tmp_key1[index2];
        	scratch_list_key2[i] = list_tmp_key2[index2];
        	scratch_list_key3[i] = list_tmp_key3[index2];
        	scratch_list_key4[i] = list_tmp_key4[index2];
			scratch_list_key5[i] = list_tmp_key5[index2];
        	scratch_list_key6[i] = list_tmp_key6[index2];
        	scratch_list_key7[i] = list_tmp_key7[index2];
        	scratch_list_key8[i] = list_tmp_key8[index2];
        	scratch_list_key9[i] = list_tmp_key9[index2];
        	scratch_list_key10[i] = list_tmp_key10[index2];
        	scratch_list_key11[i] = list_tmp_key11[index2];
        	scratch_list_key12[i] = list_tmp_key12[index2];
        	scratch_list_key13[i] = list_tmp_key13[index2];
        	scratch_list_key14[i] = list_tmp_key14[index2];
        	scratch_list_key15[i] = list_tmp_key15[index2];
        	scratch_list_key16[i] = list_tmp_key16[index2];
        	scratch_list_key17[i] = list_tmp_key17[index2];
        	scratch_list_key18[i] = list_tmp_key18[index2];
        	scratch_list_key19[i] = list_tmp_key19[index2];
        	scratch_list_key20[i] = list_tmp_key20[index2];
            index2--;
            counter++;
        }

        if (counter >= list_size)
            break;

        if (i == 0) break;
    }

    for (i = 0; i < list_size; i++){

        	list_key[i]   = scratch_list_key[i];
        	list_key1[i]  = scratch_list_key1[i];
        	list_key2[i]  = scratch_list_key2[i];
        	list_key3[i]  = scratch_list_key3[i];
        	list_key4[i]  = scratch_list_key4[i];
			list_key5[i] = scratch_list_key5[i];
			list_key6[i] = scratch_list_key6[i];
			list_key7[i] = scratch_list_key7[i];
			list_key8[i] = scratch_list_key8[i];
			list_key9[i] = scratch_list_key9[i];
			list_key10[i] = scratch_list_key10[i];
			list_key11[i] = scratch_list_key11[i];
			list_key12[i] = scratch_list_key12[i];
			list_key13[i] = scratch_list_key13[i];
			list_key14[i] = scratch_list_key14[i];
			list_key15[i] = scratch_list_key15[i];
			list_key16[i] = scratch_list_key16[i];
			list_key17[i] = scratch_list_key17[i];	
			list_key18[i] = scratch_list_key18[i];	
			list_key19[i] = scratch_list_key19[i];	
			list_key20[i] = scratch_list_key20[i];	
    }

    free(scratch_list_key);
    free(scratch_list_key1);
    free(scratch_list_key2);
    free(scratch_list_key3);
    free(scratch_list_key4);
	free(scratch_list_key5);
    free(scratch_list_key6);
    free(scratch_list_key7);
    free(scratch_list_key8);
    free(scratch_list_key9);
    free(scratch_list_key10);
    free(scratch_list_key11);
    free(scratch_list_key12);
    free(scratch_list_key13);
    free(scratch_list_key14);
    free(scratch_list_key15);
    free(scratch_list_key16);
    free(scratch_list_key17);
    free(scratch_list_key18);
    free(scratch_list_key19);
	free(scratch_list_key20);

}  /* Merge_list _high */

/*
 * -------------------------            qksort          ------------------------------
 */

int bitonic_qksort2(void *data, size_t size, size_t i, size_t k, int (*compare)(const void *key1, const void *key2)){

	int j;

	/*
	 * stop recursion when it is not possible to partition further
	 * when calling qksort:
	 * i = 0
	 * k = size-1
	 */

	while (i < k){

		/*
		 * find where to partition the elements
		 */

		if ((j = bitonic_partition2(data, i, k, compare)) < 0){
			return -1;
		}

		/*
		 * recursively sort the left partition
		 */
		if (bitonic_qksort2(data, size, i, j, compare) < 0)
			return -1;

		/*
		 * iterate and sort the right partition
		 */
		i = j + 1;
	}
	return 0;
}


int bitonic_issort2(void *data, size_t size, int (*compare)(const void *key, const void *key2)){

	size_t *a = data;
	size_t *key;
	size_t i,j;

	if ((key = malloc(sizeof(size_t))) == NULL)
		return -1;

	for ( j = 1; j < size; j++){
		memcpy(key, &a[j], sizeof(size_t));
		i = j - 1;
		while ( compare(&a[i], key) > 0 ){
			memcpy(&a[(i + 1)], &a[i], sizeof(size_t));
			if (i == 0) break;
            i--;
		}
		memcpy(&a[(i + 1)], key, sizeof(size_t));
	}

	free(key);
	return 0;

}

int bitonic_partition2(void *data, size_t i, size_t k, int (*compare)(const void *key1, const void *key2)){

	size_t *a = data;
	size_t *pval, *temp;
	size_t r[3];
	/*
	 * allocate value for the partition value and swapping
	 */
	if ((pval = malloc(sizeof(size_t))) == NULL)
		return -1;
	if ((temp = malloc(sizeof(size_t))) == NULL){
		free(pval);
		return -1;
	}


	/*
	 * Use the median-of-three method to find partition value
	 */

	int part = (k - i + 1);
	if(part <= 0)
		part = 1;

	r[0] = (rand() % part + i);
	r[1] = (rand() % part + i);
	r[2] = (rand() % part + i);

	/*
	 * TODO: replace the qsort with issort
	 *
	 * issort(r, 3, sizeof(size_t), compare_size_t_V2);
	 */
	qsort(r, 3, sizeof(size_t),compare_size_t_V2);
	memcpy(pval, &a[r[1]], sizeof(size_t));

	/*
	 * Create 2 partitions around the partition value
	 */
	i--;
	k++;
	while(1) {

		/*
		 * move left until an element is found in the wrong partition
		 */

		do {
			k--;
		} while (compare(&a[k], pval) > 0);

		/*
		 * move right until an element is found in the wrong partition
		 */

		do {
			i++;
		} while (compare(&a[i], pval) < 0);

		if (i >= k){
			/*
			 * break when left and right counter cross
			 */
			break;
		}

		else{
			// swap element under the left and right counters
			memcpy(temp, &a[i], sizeof(size_t));
			memcpy(&a[i], &a[k], sizeof(size_t));
			memcpy(&a[k], temp, sizeof(size_t));
		}
	}

	/*
	 * free the storage allocated for partitioning
	 */
	free(pval);
	free(temp);

	/*
	 * return position dividing the two partition
	 */
	return k;

}
