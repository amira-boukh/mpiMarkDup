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
     llist.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/


#ifndef LLIST_H
#define LLIST_H 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "reads.h"
#include "log.h"

/**
 * @file llist.h
 * @brief Circular doubly linked list definition and basic operations
 * @version 0.1
 * @date 01/01/2018
 */

typedef struct lnode_t {
    struct lnode_t* prev; /**< previous node */
    readInfo* read;       /**< data pointer */
    struct lnode_t* next; /**< next node */
} lnode_t;

typedef struct llist_t {
    lnode_t* nil;  /**< sentinel (dummy node) */
    int size;      /**< count element amount  */
    lnode_t* head; /**< head node */
    lnode_t* tail; /**< tail node */
} llist_t;

int compareRead(readInfo *read1, readInfo *read2, const int isFragment) ;
llist_t*  llist_create();
void llist_clear(llist_t *l) ;
lnode_t *llist_append(llist_t *l, readInfo *e) ;
void  llist_destruct(llist_t* l);
void llist_merge_sort(llist_t* llist, const int isFragment) ;
lnode_t *llist_insert_node(llist_t *l, lnode_t* lnode, readInfo *e) ;
lnode_t*  llist_insert_index(llist_t* l, int n, readInfo* e);
readInfo* llist_delete_ptr(llist_t *l, lnode_t* p) ;
readInfo* llist_delete(llist_t* l, int n);
void  llist_readInfo_print(llist_t* l, MPI_Comm comm);
lnode_t*  llist_lsearch(llist_t* l, int n);
#endif /* ifndef LLIST_H */
