#ifndef PROTEIN_H
#define PROTEIN_H


#include "vector"
#include "string"

using namespace std;

class protein
{
    public:
        protein();
        virtual ~protein();
        void load_file( string filename );

        void gen_score1();
        void gen_score2();
        void gen_score3();
        void gen_score4();
        void gen_score5();

        vector<double *> score1;
        vector<double *> score2;
        vector<double *> score3;
        vector<double *> score4;
        vector<double *> score5;

        vector<string> all_name;
    private:

        string data;
        double *ten;
        vector<string> all_protein;

    private:
        void lq( double *score , int len , double *ten );

};

#endif // PROTEIN_H
