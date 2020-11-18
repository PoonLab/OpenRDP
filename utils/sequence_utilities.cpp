#include "omp.h"
#include <stdexcept>
#include <malloc.h>

using namespace std;

namespace MathFuncs {
    /*
    Categorizes sites based on number of nucleotide polymorphisms and the presence and absence of markers

    @param <int> NextNo - the number of sequences
    @param <int> LSeq - the length of the strain sequence
    @param <int> UBSN1 - the length of the sequence number
    @param <int> UBSCC1 - the length of the sequence category count
    @param <int> StartPosInAlign - the start position in the alignment
    @param <int> EndPosInAlign - the end position in the alignment
    @param <int array> SeqCatCount
    @param <int array> AA
    @param <int array> SeqNum - matrix of sequences
    @param <string> NucMat
    @param <string> SeqSpace - stores if there are gaps at a site (0 if no gaps, 1 if gaps)
    @param <string> NucMatB
    @param <string> NucMat2
    @param <string> flp
    @param <string> ml - stores if there are gaps at a site (0 if no gaps, 1 if gaps)
    @param <string> nl - stores number of variants at a site
    */
    int MyMathFuncs::MakeSeqCatCount2P(int Nextno, int LSeq, int UBSN1, int UBSCC1, int StartPosInAlign,
                                       int EndPosInAlign, int * SeqCatCount, int * AA, short int * SeqNum,
                                       unsigned char * NucMat, unsigned char * SeqSpace, unsigned char * NucMatB,
                                       unsigned char * NucMat2, unsigned char * flp, unsigned char * ml,
                                       unsigned char * nl) {

        int X, Y, os1, os2, ANum, NumNucs, off3, off2, aapos;

        os1 = UBSN1 + 1;
        os2 = UBSCC1 + 1;
        ANum = 0;
        SeqCatCount[0] = 0;
        SeqCatCount[1] = 0;
        SeqCatCount[2] = 0;
        SeqCatCount[3] = 0;
        SeqCatCount[4] = 0;
        SeqCatCount[5] = 0;
        SeqCatCount[6] = 0;
        SeqCatCount[7] = 0;
        SeqCatCount[8] = 0;
        SeqCatCount[9] = 0;
        off3 = Nextno + 1;

        if (StartPosInAlign <= EndPosInAlign) {
            // Loop for each site in the aligned  sequences
            for (X = StartPosInAlign; X <= EndPosInAlign; X++) {
                off2 = off3*X;   // offset2 = curr_pos + 1 * site in alignment
                AA[0] = 0;
                AA[1] = 0;
                AA[2] = 0;
                AA[3] = 0;
                AA[4] = 0;

                // Stores whether there was a gap, A, C, G, T in each column
                #pragma omp parallel
                for private(Y, aapos)   // not shared between threads
                for (Y = 0; Y <= Nextno; Y++) {
                    // Get sequence number Y at site X
                    // Look in NucMat to see the site is a 1 (A), 2 (C), 3 (G), 4 (T)
                    aapos = NucMat[SeqNum[X + Y*os1]];

                    // If the position corresponds to a gap (0) in AA, set AA[pos] to 1
                    if (AA[aapos] == 0)
                        AA[aapos] = 1;
                }

                // Number of different variants (A, C, G, T)
                ANum = (int)(AA[1] + AA[2] + AA[3] + AA[4]);

                // Gap at site A
                SeqSpace[X] = AA[0];

                // 1 if there was a gap, 0 otherwise
                NumNucs = (int)(AA[0]);

                // Increment category counts (2 * 5 array)
                SeqCatCount[NumNucs + ANum*os2] = SeqCatCount[NumNucs + ANum*os2] + 1;

                // Set NucMat 2 to store 0-3 if gaps or 1-4 is no gaps
                for (Y = 1; Y <= 4; Y++) {
                    if (AA[Y] == 1) {
                        NucMat2[NucMatB[Y]] = NumNucs;
                        NumNucs = NumNucs + 1;
                    }
                }

                ml[X] = AA[0];      // Store if gap at site in alignment
                nl[X] = ANum;       // Store number of different variants

                //
                #pragma omp parallel
                for private(Y)
                for (Y = 0; Y <= Nextno; Y++)
                    flp[Y + off2] = NucMat2[(int)(SeqNum[X + Y*os1])];
            }
        }

        // Wrap around if the end position is before the start
        else if (StartPosInAlign > EndPosInAlign) {
            // Loop over sites from the start position to the end
            for (X = StartPosInAlign; X <= LSeq; X++) {
                off2 = off3*X;
                AA[0] = 0;
                AA[1] = 0;
                AA[2] = 0;
                AA[3] = 0;
                AA[4] = 0;

                #pragma omp parallel
                for private(Y, aapos)
                for (Y = 0; Y <= Nextno; Y++) {
                    aapos = NucMat[SeqNum[X + Y*os1]];
                    if (AA[aapos] == 0)
                        AA[aapos] = 1;
                }
                ANum = (int)(AA[1] + AA[2] + AA[3] + AA[4]);

                SeqSpace[X] = AA[0];
                NumNucs = (int)(AA[0]);
                SeqCatCount[NumNucs + ANum*os2] = SeqCatCount[NumNucs + ANum*os2] + 1;
                NumNucs = AA[0];

                for (Y = 1; Y <= 4; Y++) {
                    if (AA[Y] == 1) {
                        NucMat2[NucMatB[Y]] = NumNucs;
                        NumNucs = NumNucs + 1;
                    }
                }
                ml[X] = AA[0];
                nl[X] = ANum;

                #pragma omp parallel
                for private(Y)
                for (Y = 0; Y <= Nextno; Y++)
                    flp[Y + off2] = NucMat2[(int)(SeqNum[X + Y*os1])];
            }

            // Loop over sites from the beginning of the sequence to the end position
            for (X = 1; X <= EndPosInAlign; X++) {
                off2 = off3*X;
                AA[0] = 0;
                AA[1] = 0;
                AA[2] = 0;
                AA[3] = 0;
                AA[4] = 0;

                #pragma omp parallel
                for private(Y, aapos)
                for (Y = 0; Y <= Nextno; Y++) {
                    aapos = NucMat[SeqNum[X + Y*os1]];
                    if (AA[aapos] == 0)
                        AA[aapos] = 1;
                }

                ANum = (int)(AA[1] + AA[2] + AA[3] + AA[4]);

                SeqSpace[X] = AA[0];
                NumNucs = (int)(AA[0]);
                SeqCatCount[NumNucs + ANum*os2] = SeqCatCount[NumNucs + ANum*os2] + 1;
                NumNucs = AA[0];

                for (Y = 1; Y <= 4; Y++) {
                    if (AA[Y] == 1) {
                        NucMat2[NucMatB[Y]] = NumNucs;
                        NumNucs = NumNucs + 1;
                    }
                }
                ml[X] = AA[0];
                nl[X] = ANum;

                #pragma omp parallel
                for private(Y)
                for (Y = 0; Y <= Nextno; Y++)
                    flp[Y + off2] = NucMat2[(int)(SeqNum[X + Y*os1])];
            }

        }

        return(1);
}
