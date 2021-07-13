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
     parser.c

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

#include "parser.h"
#include "log.h"
#include "reads.h"
#include "mpiMarkDuplicatesUtils.h"
#include "markDuplicates.h"
#include "createLBList.h"

char parse_mode;

size_t hash_name(char *line, size_t max) {
    
    size_t i, j, result = 0;

    for (i = 0, j = 0; i < strlen(line); i++) {
        if ( (line[i] >= 'A' && line[i] <= 'Z') || (line[i] >= 'a' && line[i] <= 'z') ) {
            j = i;
        }
    }

    i = j + 1;

    for (j = 0; i < strlen(line) && j < max; i++) {
        if (line[i] >= '0' && line[i] <= '9') {
            char tmp;
            tmp = line[i] - '0';
            result <<= 4;
            result += tmp;
            j++;
        }
    }

    return result;
}

void init_goff(MPI_File mpi_filed, unsigned int headerSize, size_t fsize, int numproc, size_t *goff) {


    char *current_line = NULL;
    MPI_Status status;
    int i = 0;
    size_t j = 0;

    size_t lsize = fsize / numproc;
    goff[0] = headerSize;

    for (i = 1; i < numproc; i++) {
        goff[i] = lsize * i + headerSize;
    }

    goff[numproc] = fsize;

    for (i = 1; i < numproc; i++) {
        current_line = (char *)calloc(1000, sizeof(char));
        MPI_File_read_at(mpi_filed, (MPI_Offset)goff[i], current_line, 1000, MPI_CHAR, &status);
        assert(strlen(current_line) != 0);
        j = 0;

        while (j < fsize && current_line[j] != '\n') {
            j++;
        }

        goff[i] += (j + 1);
        free(current_line);
    }


}

void parser_paired(char *localData, int rank, size_t start_offset, unsigned char threshold,
                   int nbchrom, size_t **preadNumberByChr, char **chrNames, readInfo ***preads, int optical_distance, char* header) {

    char *currentCarac;
    char *cigar; 
    char *cigar_mate;
    char currentLine[MAX_LINE_SIZE];
    unsigned char quality;
    int i, chr, mchr, nbchr = 0;
    int lastChr = -1;
    int next, mate_score, key;
    size_t lineSize, offset_read_in_source_file;
    size_t coord, coord_mate,flag;
    size_t *readNumberByChr;
    size_t counter = 0;
    readInfo **reads = *preads;
    char *qname;  
    unsigned int firstInPair;
    //lbInfo *lb; 

    chrInfo chromosome;
    lbInfo lb;
    initChrInfo(&chromosome, header);
    initLbInfo(&lb, header);



    for (i = 0; i < MAX_LINE_SIZE; i++) {
        currentLine[i] = 0;
    }

    //we take the first line *
    //before calling parsepaired, we know that localdata is at the begining of a read
    next = tokenizer(localData, '\n', currentLine);
    offset_read_in_source_file = start_offset;

    nbchr = nbchrom;
    readNumberByChr = (size_t *)calloc(nbchr, sizeof(size_t));

    while (next) {

        lineSize = strlen(currentLine) + 1;

        //we update the offset in the
        //source file
        currentLine[lineSize - 1] = '\n';
        currentLine[lineSize] = '\0';

        //GO TO FLAG
        currentCarac = strstr(currentLine, "\t");
        //Qname 
        char *tokenCar;
        currentCarac = currentLine;
        int count = getTokenTab(&currentCarac, &tokenCar);
        qname = tokenCar;
        
        //free(tokenCar);


        //*currentCarac = '\0';
        //currentCarac++;
        //currentCarac = strstr(currentCarac, "\t");
        flag = strtoull(currentCarac, &currentCarac,10);
        currentCarac = strstr(currentCarac, "\t");
        //fprintf(stdout,"%d\n",flag);
        

        if (lastChr == (nbchr - 1)) {
            chr = (nbchr - 1);

        } else {
            chr = getChr(currentCarac, chrNames, nbchr);
        }

        //GO TO COORD
        currentCarac = strstr(currentCarac + 1, "\t");


    
        key = strtoull(currentLine, NULL, strlen(currentLine));
        key = hash_name(currentLine, 16);

        //strtoull(currentCarac, &currentCarac, 10);

    
        //TAKE COORD AND GO TO MAPQ
        coord = strtoull(currentCarac, &currentCarac, 10);

        
       //fprintf(stderr,"qname : %s, key : %d, coord : %zu\n",qname, key,coord);
        
    

        //TAKE MAPQ AND GO TO CIGAR
        quality = strtoull(currentCarac, &currentCarac, 10);
        currentCarac = strstr(currentCarac, "\t");
        currentCarac = currentCarac +1;
        count = getTokenTab(&currentCarac, &tokenCar);
        cigar = tokenCar;
        //free(tokenCar);


        //GO TO RNEXT
        currentCarac = strstr(currentCarac-1, "\t");

        if (currentCarac[1] == '=') {
            mchr = chr;

        } else if (currentCarac[1] == '*') {
            mchr = (nbchr - 1);

        } else {
            mchr = getChr(currentCarac, chrNames, nbchr);
        }
    
        currentCarac = strstr(currentCarac+1, "\t");
        coord_mate = strtoull(currentCarac,&currentCarac, 10); 

        currentCarac = strstr(currentCarac+1, "\t");
        currentCarac = strstr(currentCarac+1, "\t");
        currentCarac = strstr(currentCarac+1, "\t");
        currentCarac = currentCarac +1;


        
        while (getTokenTab(&currentCarac, &tokenCar)) {
            if (strncmp(tokenCar, "MC:Z:", strlen("MC:Z:")) == 0) {
                cigar_mate = getReadTagValue(tokenCar, "MC:Z:");

            }
            if (strncmp(tokenCar, "ms:i:", strlen("ms:i:")) == 0) {
                char *score = getReadTagValue(tokenCar, "ms:i:");
                mate_score = strtoull(score,&score, 10); 

            }
            if (strncmp(tokenCar, "LB:Z:", strlen("LB:Z:")) == 0) {
            //fillReadLBValue(tokenCar, read);
            char *lbName = getReadTagValue(tokenCar, "LB:Z:");

            if (strcmp((&lb)->lbList[(&lb)->lastLb], lbName) == 0) {
                reads[chr]->next->readLb = (&lb)->lastLb;

            } else {
                for (i = 0; i < (&lb)->lbNum; ++i) {
                    if (strcmp((&lb)->lbList[i], lbName) == 0) {
                        reads[chr]->next->readLb = i;
                    }
                }
            }

            free(lbName);
            free(tokenCar);
            break;
        }

        free(tokenCar);
        }

        //fprintf(stderr,"qname : %s, mate score : %d cigar mate : %s\n",qname, mate_score, cigar_mate);


        //first we check if reads mapped on the same chromosome

        if ((chr < nbchr - 2) && (chr == mchr) ) {
            //then we found concordant reads
            if (quality >= threshold) {

                reads[chr]->next = malloc(sizeof(readInfo));
                reads[chr]->next->Qname = qname;
                reads[chr]->next->valueFlag = flag;
                reads[chr]->next->coordPos = coord;
                reads[chr]->next->coordMatePos = coord_mate;
                reads[chr]->next->readChromosome = chr;
                reads[chr]->next->mateChromosome = mchr;
                reads[chr]->next->quality = quality;
                reads[chr]->next->cigar = cigar;
                reads[chr]->next->mate_cigar = cigar_mate;
                reads[chr]->next->offset_source_file = offset_read_in_source_file;
                reads[chr]->next->offset = lineSize;
                fillUnclippedCoord(reads[chr]->next);
                reads[chr]->next->orientation = getOrientation(reads[chr]->next, 0);
                if (optical_distance > 0) {
                    reads[chr]->next->physicalLocation = computePhysicalLocation(reads[chr]->next);
                }
                firstInPair = readBits((unsigned int)reads[chr]->next->valueFlag, 6);
                if (firstInPair == 1)
                    reads[chr]->next->pair_num = 1;
                else
                    reads[chr]->next->pair_num = 2;

                reads[chr]->next->fingerprint = read2Fingerprint(reads[chr]->next);
                reads[chr]->next->mate_fingerprint = read2mateFP(reads[chr]->next);
                reads[chr] = reads[chr]->next;
                readNumberByChr[chr]++;
                //fprintf(stdout,"flag %ld rmane %d quality %d\n",flag,chr,quality);
            }

        } else if ((chr < (nbchr - 2)) && ( mchr < (nbchr - 2))) {

            
            if (quality >= threshold) {
                
                //we found discordant reads
                reads[nbchr - 1]->next = malloc(sizeof(readInfo));
                reads[nbchr - 1]->next->Qname = qname;
                reads[nbchr - 1]->next->valueFlag = flag;
                reads[nbchr - 1]->next->coordPos = coord;
                reads[nbchr - 1]->next->coordMatePos = coord_mate;
                reads[nbchr - 1]->next->readChromosome = chr;
                reads[nbchr - 1]->next->mateChromosome = mchr;
                reads[nbchr - 1]->next->cigar = cigar;
                reads[nbchr - 1]->next->mate_cigar = cigar_mate;
                reads[nbchr - 1]->next->quality = quality;
                reads[nbchr - 1]->next->offset_source_file = offset_read_in_source_file;
                reads[nbchr - 1]->next->offset = lineSize;
                reads[nbchr - 1]->next->orientation = getOrientation(reads[nbchr - 1]->next, 0);
                fillUnclippedCoord(reads[nbchr - 1]->next);
                if (optical_distance > 0) {
                    reads[nbchr - 1]->next->physicalLocation = computePhysicalLocation(reads[nbchr-1]->next);
                }
               firstInPair = readBits((unsigned int)reads[nbchr - 1]->next->valueFlag, 6);
                if (firstInPair == 1)
                    reads[nbchr - 1]->next->pair_num = 1;
                else
                    reads[nbchr - 1]->next->pair_num = 2;
                reads[nbchr - 1]->next->fingerprint = read2Fingerprint(reads[nbchr - 1]->next);
                reads[nbchr - 1]->next->mate_fingerprint = read2mateFP(reads[nbchr - 1]->next);
                reads[nbchr - 1] = reads[nbchr - 1]->next;
                readNumberByChr[nbchr - 1]++;

                //we add it in reads[chr] too
                //as we want to keep it in the bam
                reads[chr]->next = malloc(sizeof(readInfo));
                reads[chr]->next->Qname = qname;
                reads[chr]->next->valueFlag = flag;
                reads[chr]->next->coordPos = coord;
                reads[chr]->next->coordMatePos = coord_mate;
                reads[chr]->next->quality = quality;
                reads[chr]->next->readChromosome = chr;
                reads[chr]->next->mateChromosome = mchr;
                reads[chr]->next->cigar = cigar;
                reads[chr]->next->mate_cigar = cigar_mate;
                reads[chr]->next->offset_source_file = offset_read_in_source_file;
                reads[chr]->next->offset = lineSize;
                reads[chr]->next->orientation = getOrientation(reads[chr]->next, 0);
                fillUnclippedCoord(reads[chr]->next);
                if (optical_distance > 0) {
                    reads[chr]->next->physicalLocation = computePhysicalLocation(reads[chr]->next);
                }
                firstInPair = readBits((unsigned int)reads[chr]->next->valueFlag, 6);
                if (firstInPair == 1)
                    reads[chr]->next->pair_num = 1;
                else
                    reads[chr]->next->pair_num = 2;
                reads[chr]->next->fingerprint = read2Fingerprint(reads[chr]->next);
                reads[chr]->next->mate_fingerprint = read2mateFP(reads[chr]->next);
                reads[chr] = reads[chr]->next;
                readNumberByChr[chr]++;
      
            }

        } else if ((chr == '*') && ( mchr < (nbchr - 2))) {

            //we found discordant reads with one pair unmapped
            reads[nbchr - 2]->next = malloc(sizeof(readInfo));
            reads[nbchr - 2]->next->Qname = qname;
            reads[nbchr - 2]->next->valueFlag = flag;
            reads[nbchr - 2]->next->offset_source_file = offset_read_in_source_file;
            reads[nbchr - 2]->next->offset = lineSize;
            reads[nbchr - 2] = reads[nbchr - 2]->next;
            reads[nbchr - 2]->next->readChromosome = -1;
            readNumberByChr[nbchr - 2]++;

        } else if ((mchr == '*') && ( chr < (nbchr - 2))) {

            //we found discordant reads with one pair unmapped
            reads[nbchr - 2]->next = malloc(sizeof(readInfo));
            reads[nbchr - 2]->next->Qname = qname;
            reads[nbchr - 2]->next->valueFlag = flag;
            reads[nbchr - 2]->next->offset_source_file = offset_read_in_source_file;
            reads[nbchr - 2]->next->offset = lineSize;
            reads[nbchr - 2]->next->readChromosome = -1;
            reads[nbchr - 2] = reads[nbchr - 2]->next;
            readNumberByChr[nbchr - 2]++;
        }

        else {
            //we found unmapped pairs reads
            reads[nbchr - 2]->next = malloc(sizeof(readInfo));
            reads[nbchr - 2]->next->Qname = qname;
            reads[nbchr - 2]->next->valueFlag = flag;
            reads[nbchr - 2]->next->offset_source_file = offset_read_in_source_file;
            reads[nbchr - 2]->next->offset = lineSize;
            reads[nbchr - 2]->next->readChromosome = -1;
            reads[nbchr - 2] = reads[nbchr - 2]->next;
            readNumberByChr[nbchr - 2]++;
        }



        //we update the offset_read_in_source_file
        offset_read_in_source_file += lineSize;
        //we read the next line

        for (i = 0; i < MAX_LINE_SIZE; i++) {
            currentLine[i] = 0;
        }

        next = tokenizer(NULL, '\n', currentLine);

        counter++;
    }

    //fprintf(stderr, "rank %d ::: counter = %zu \n", rank, counter);
    md_log_rank_trace(rank, "counter = %zu\n", counter);

    for (i = 0; i < nbchr; i++) {
        preadNumberByChr[0][i] += readNumberByChr[i];
    }

    free(readNumberByChr);
}

int getChr(char *str, char **chrNames, int nbchr) {
    int i = 0, found = 0, size;
    char *str1 = str, *str2;

    str2 = str1 + 1;

    for (; *str2 != '\t'; str2++);

    size = strlen(str1) - strlen(str2);

    char *tmp_chr = (char *)malloc(sizeof(char) * size + 1);
    tmp_chr[0] = 0;

    for (i = 0; i < size; i++) {
        tmp_chr[i] = str1[i + 1];
    }

    tmp_chr[size - 1] = 0;

    assert(strlen(tmp_chr) != 0);

    for (i = 0, found = 0; i < nbchr && !found; i++) {
        found = !strcmp(tmp_chr, chrNames[i]);
    }

    free(tmp_chr);
    return i - 1;
}
