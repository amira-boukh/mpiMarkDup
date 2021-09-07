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
     parallelBitonicSort2.h

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#include <mpi.h>

//typedef size_t size_t;
void Generate_local_list2(
		size_t list_size,
		size_t local_list[]
		);
void Local_sort2(
		size_t list_size,
		size_t local_keys[],
		int    local_keys1[],
		int    local_keys2[],
		size_t local_keys3[],
		int local_keys4[]
		);

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
	size_t *local_keys19); 

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
	size_t *local_keys20); 


int Key_compare2(const size_t* p, const size_t* q);
void Par_bitonic_sort_incr2(
		size_t list_size,
		size_t local_list[],
		int    local_list1[],
		int    local_list2[],
		size_t local_list3[],
		int local_list4[],
		int    proc_set_size,
        int    rank
        );

void Par_bitonic_sort_incrAll(
        size_t      list_size      /* in     */,
        size_t*    	local_list     /* in/out */,
        int*    	local_list1    /* in/out */,
        int*    	local_list2    /* in/out */,
        size_t*    	local_list3    /* in/out */,
        int*    	local_list4    /* in/out */,
		int         *local_list5,
		unsigned int *local_list6,
		unsigned int *local_list7,
		unsigned int *local_list8,
		int *local_list9,
		int *local_list10,
		int *local_list11,
		int *local_list12,
		int *local_list13,
		int *local_list14,
		int *local_list15,
		int* local_list16,
		size_t* local_list17,
		size_t* local_list18,
		size_t* local_list19,
        int       	proc_set_size  /* in     */,
        int 	  	my_rank
        );

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
        );


void Par_bitonic_sort_decr2(
		size_t list_size,
		size_t local_list[],
		int    local_list1[],
		int    local_list2[],
		size_t local_list3[],
		int local_list4[],
		int    proc_set_size,
        int    rank
        );

void Par_bitonic_sort_decrAll(
        size_t	list_size      /* in     */,
        size_t* local_list     /* in/out */,
        int*    local_list1    /* in/out */,
        int*    local_list2    /* in/out */,
        size_t* local_list3    /* in/out */,
        int*    local_list4    /* in/out */,
		int         *local_list5,
		unsigned int *local_list6,
		unsigned int *local_list7,
		unsigned int *local_list8,
		int *local_list9,
		int *local_list10,
		int *local_list11,
		int *local_list12,
		int *local_list13,
		int *local_list14,
		int *local_list15,
		int* local_list16,
		size_t* local_list17,
		size_t* local_list18,
		size_t* local_list19,
        int     proc_set_size  /* in     */,
        int 	my_rank
        );

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
        );

void Merge_split2(
		size_t list_size,
		size_t local_list[],
		int    local_list1[],
		int    local_list2[],
		size_t local_list3[],
		int local_list4[],
		int    which_keys,
		int    partner
		);

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
        int       partner       /* in     */);

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
        int       partner       /* in     */);


void Merge_list_low2(
		size_t  list_size,
		size_t  list_key[],
		int     list_key1[],
		int     list_key2[],
		size_t  list_key3[],
		int  list_key4[],
		size_t  list_tmp_key[],
		int     list_tmp_key1[],
		int     list_tmp_key2[],
		size_t  list_tmp_key3[],
		int  list_tmp_key4[]
		);

void Merge_list_lowAll(
        size_t   list_size,
        size_t  list_key[],
        int		list_key1[],
        int		list_key2[],
        size_t  list_key3[],
        int     list_key4[],
        int     list_key5[],
        unsigned int     list_key6[],
        unsigned int     list_key7[],
        unsigned int     list_key8[],
        int     list_key9[],
        int     list_key10[],
        int     list_key11[],
        int     list_key12[],
        int     list_key13[],
        int     list_key14[],
        int     list_key15[],
        int     list_key16[],
        size_t     list_key17[],
        size_t     list_key18[],
        size_t     list_key19[],
        size_t  list_tmp_key[],
        int		list_tmp_key1[],
        int		list_tmp_key2[],
        size_t  list_tmp_key3[],
        int     list_tmp_key4[],
        int     list_tmp_key5[],
        unsigned int     list_tmp_key6[],
        unsigned int     list_tmp_key7[],
        unsigned int     list_tmp_key8[],
        int     list_tmp_key9[],
        int     list_tmp_key10[],
        int     list_tmp_key11[],
        int     list_tmp_key12[],
        int     list_tmp_key13[],
        int     list_tmp_key14[],
        int     list_tmp_key15[],
        int     list_tmp_key16[],
        size_t     list_tmp_key17[],
        size_t    list_tmp_key18[],
        size_t    list_tmp_key19[]
        );

void Merge_list_lowAll2(
        size_t   list_size,
        size_t  list_key[],
        int		list_key1[],
        int		list_key2[],
        size_t  list_key3[],
        size_t  list_key4[],
        int     list_key5[],
        int     list_key6[],
        unsigned int     list_key7[],
        unsigned int     list_key8[],
        unsigned int     list_key9[],
        int     list_key10[],
        int     list_key11[],
        int     list_key12[],
        int     list_key13[],
        int     list_key14[],
        int     list_key15[],
        int     list_key16[],
        int     list_key17[],
        size_t     list_key18[],
        size_t     list_key19[],
        size_t     list_key20[],
        size_t  list_tmp_key[],
        int		list_tmp_key1[],
        int		list_tmp_key2[],
        size_t  list_tmp_key3[],
        size_t  list_tmp_key4[],
        int     list_tmp_key5[],
        int     list_tmp_key6[],
        unsigned int     list_tmp_key7[],
        unsigned int     list_tmp_key8[],
        unsigned int     list_tmp_key9[],
        int     list_tmp_key10[],
        int     list_tmp_key11[],
        int     list_tmp_key12[],
        int     list_tmp_key13[],
        int     list_tmp_key14[],
        int     list_tmp_key15[],
        int     list_tmp_key16[],
        int     list_tmp_key17[],
        size_t     list_tmp_key18[],
        size_t    list_tmp_key19[],
        size_t    list_tmp_key20[]
        );


void Merge_list_high2(
		size_t  list_size,
		size_t  list_key[],
		int     list_key1[],
		int     list_key2[],
		size_t  list_key3[],
		int     list_key4[],
		size_t  list_tmp_key[],
		int     list_tmp_key1[],
		int     list_tmp_key2[],
		size_t  list_tmp_key3[],
		int     list_tmp_key4[]
		);
void Merge_list_highAll(size_t   list_size,
        size_t  list_key[],
        int		list_key1[],
        int		list_key2[],
        size_t  list_key3[],
        int     list_key4[],
        int     list_key5[],
        unsigned int     list_key6[],
        unsigned int     list_key7[],
        unsigned int     list_key8[],
        int     list_key9[],
        int     list_key10[],
        int     list_key11[],
        int     list_key12[],
        int     list_key13[],
        int     list_key14[],
        int     list_key15[],
        int     list_key16[],
        size_t     list_key17[],
        size_t     list_key18[],
        size_t     list_key19[],
        size_t  list_tmp_key[],
        int		list_tmp_key1[],
        int		list_tmp_key2[],
        size_t  list_tmp_key3[],
        int     list_tmp_key4[],
        int     list_tmp_key5[],
        unsigned int     list_tmp_key6[],
        unsigned int     list_tmp_key7[],
        unsigned int     list_tmp_key8[],
        int     list_tmp_key9[],
        int     list_tmp_key10[],
        int     list_tmp_key11[],
        int     list_tmp_key12[],
        int     list_tmp_key13[],
        int     list_tmp_key14[],
        int     list_tmp_key15[],
        int     list_tmp_key16[],
        size_t  list_tmp_key17[],
        size_t    list_tmp_key18[],
		size_t    list_tmp_key19[]

);

void Merge_list_highAll2(
        size_t   list_size,
        size_t  list_key[],
        int		list_key1[],
        int		list_key2[],
        size_t  list_key3[],
        size_t  list_key4[],
        int     list_key5[],
        int     list_key6[],
        unsigned int     list_key7[],
        unsigned int     list_key8[],
        unsigned int     list_key9[],
        int     list_key10[],
        int     list_key11[],
        int     list_key12[],
        int     list_key13[],
        int     list_key14[],
        int     list_key15[],
        int     list_key16[],
        int     list_key17[],
        size_t     list_key18[],
        size_t     list_key19[],
        size_t     list_key20[],
        size_t  list_tmp_key[],
        int		list_tmp_key1[],
        int		list_tmp_key2[],
        size_t  list_tmp_key3[],
        size_t  list_tmp_key4[],
        int     list_tmp_key5[],
        int     list_tmp_key6[],
        unsigned int     list_tmp_key7[],
        unsigned int     list_tmp_key8[],
        unsigned int     list_tmp_key9[],
        int     list_tmp_key10[],
        int     list_tmp_key11[],
        int     list_tmp_key12[],
        int     list_tmp_key13[],
        int     list_tmp_key14[],
        int     list_tmp_key15[],
        int     list_tmp_key16[],
        int     list_tmp_key17[],
        size_t     list_tmp_key18[],
        size_t    list_tmp_key19[],
        size_t    list_tmp_key20[]
        );


		 
int bitonic_qksort2(
		void   *data,
		size_t size,
		size_t i,
		size_t k,
		int    (*compare)(const void *key1, const void *key2)
		);
int bitonic_partition2(
		void   *data,
		size_t i,
		size_t k,
		int    (*compare)(const void *key1, const void *key2)
		);
void ParallelBitonicSort2(
		MPI_Comm split_comm,
		int      my_rank,
		int      dimension,
		size_t   *local_list,
		int      *local_list1,
		int      *local_list2,
		size_t   *local_list3,
		int      *local_list4,
		size_t   list_size
		);

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
	int *local_list15,                //stand for physical location y 
	int *local_list16,
	size_t *local_list17,             //stand for mate coordinates
	size_t *local_list18,             //stand for unclipped_coordinates
	size_t *local_lits19,
	size_t list_size
);

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
	size_t *local_lits20,
	size_t list_size
);
