#include "omp.h"
#include <stdexcept>
#include <malloc.h>

using namespace std;

namespace MathFuncs {
    // Categorizes sites based on number of nucleotide polymorphisms and the presence and absence of markers
    // P: parallelized
    // **Called by FastDistanceCalcZ in Module 4 **
    int MyMathFuncs::MakeSeqCatCount2P(int Nextno, int LSeq, int UBSN1, int UBSCC1, int StartPosInAlign,
                                       int EndPosInAlign, int *SeqCatCount, int *AA, short int *SeqNum,
                                       unsigned char *NucMat, unsigned char *SeqSpace, unsigned char *NucMatB,
                                       unsigned char *NucMat2, unsigned char *flp, unsigned char *ml,
                                       unsigned char *nl) {

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

                for (X = StartPosInAlign; X <= EndPosInAlign; X++) {
                    off2 = off3*X;
                    AA[0] = 0;
                    AA[1] = 0;
                    AA[2] = 0;
                    AA[3] = 0;
                    AA[4] = 0;

    #pragma omp parallel for private(Y, aapos)
                    for (Y = 0; Y <= Nextno; Y++) {
                        aapos = NucMat[SeqNum[X + Y*os1]];
                        if (AA[aapos] == 0)
                            AA[aapos] = 1;
                    }
                    ANum = (int)(AA[1] + AA[2] + AA[3] + AA[4]);

                    SeqSpace[X] = AA[0];
                    NumNucs = (int)(AA[0]);
                    SeqCatCount[NumNucs + ANum*os2] = SeqCatCount[NumNucs + ANum*os2] + 1;

                    for (Y = 1; Y <= 4; Y++) {
                        if (AA[Y] == 1) {
                            NucMat2[NucMatB[Y]] = NumNucs;
                            NumNucs = NumNucs + 1;
                        }
                    }
                    ml[X] = AA[0];
                    nl[X] = ANum;

    #pragma omp parallel for private(Y)
                    for (Y = 0; Y <= Nextno; Y++)
                        flp[Y + off2] = NucMat2[(int)(SeqNum[X + Y*os1])];
                }
            }
            else if (StartPosInAlign > EndPosInAlign) {
                for (X = StartPosInAlign; X <= LSeq; X++) {
                    off2 = off3*X;
                    AA[0] = 0;
                    AA[1] = 0;
                    AA[2] = 0;
                    AA[3] = 0;
                    AA[4] = 0;

    #pragma omp parallel for private(Y, aapos)
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
    #pragma omp parallel for private(Y)
                    for (Y = 0; Y <= Nextno; Y++)
                        flp[Y + off2] = NucMat2[(int)(SeqNum[X + Y*os1])];
                }
                for (X = 1; X <= EndPosInAlign; X++) {
                    off2 = off3*X;
                    AA[0] = 0;
                    AA[1] = 0;
                    AA[2] = 0;
                    AA[3] = 0;
                    AA[4] = 0;
    #pragma omp parallel for private(Y, aapos)
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
    #pragma omp parallel for private(Y)
                    for (Y = 0; Y <= Nextno; Y++)
                        flp[Y + off2] = NucMat2[(int)(SeqNum[X + Y*os1])];
                }

            }

            return(1);
        }
