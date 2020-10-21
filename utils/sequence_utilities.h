#ifndef SEQUENCE_UTILITIES_H
#define SEQUENCE_UTILITIES_H


namespace MathFuncs {
    class MyMathFuncs {

    public:
        static MATHFUNCSDLL_API int _stdcall MakeSeqCatCount2P(int Nextno, int LSeq, int UBSN1, int UBSCC1,
                                                               int StartPosInAlign, int EndPosInAlign, int *SeqCatCount,
                                                               int *AA, short int *SeqNum, unsigned char *NucMat,
                                                               unsigned char *SeqSpace, unsigned char *NucMatB,
                                                               unsigned char *NucMat2, unsigned char *flp,
                                                               unsigned char *ml, unsigned char *nl);
    }
}

#endif