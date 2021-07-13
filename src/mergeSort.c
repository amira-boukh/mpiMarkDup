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
     mergeSort.c

   Authors:
    Frederic Jarlier,   Institut Curie
    Nicolas Joly,       Institut Pasteur
    Nicolas Fedy,       Institut Curie
    Leonor Sirotti,     Institut Curie
    Thomas Magalhaes,   Institut Curie
    Paul Paganiban,     Institut Curie
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mergeSort.h"
#include "reads.h"

readInfo *mergeSort(readInfo *c, size_t n) {
    size_t q, p;
    readInfo *d;

    q = n / 2;
    p = n - q;

    if (p > 1) {
        d = mergeSort(c, p);

        if (q > 1) {
            mergeSort(d, q);
        }

    } else {
        d = c->next;
    }

    d = structMerge(c, p, d, q);

    return d;
}

readInfo *structMerge(readInfo *c, size_t p, readInfo *d, size_t q) {

    readInfo *t;

    while (1) {

        if (c->next->coordPos > d->next->coordPos) {
            t = d->next;
            d->next = t->next;
            t->next = c->next;
            c->next = t;

            if (q == 1) {
                break;
            }

            --q;

        } else {
            if (p == 1) {
                while (q > 0) {
                    d = d->next;
                    --q;
                }

                break;
            }

            --p;
        }

        c = c->next;
    }

    return d;
}
