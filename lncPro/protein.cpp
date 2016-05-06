#include "protein.h"
#include "fstream"
#include "iostream"
#include "cmath"
#include "cstring"

#define PI 3.141592653589793

protein::protein()
{
    //ctor
}

protein::~protein()
{
    //dtor
}
void protein::load_file( string filename)
{
    data = filename;
    ifstream in( filename.c_str() );
    string tmp;
    string buff;
    getline( in , tmp );
    while( 1 )
    {
        buff = "";
        for( int i = 1 ; i < (int)tmp.size() ; ++i )
            buff += tmp[i];
        all_name.push_back( buff );
        buff = "";
        while( getline( in , tmp ) )
        {
            if( tmp[0] == '>' )
            {
                all_protein.push_back( buff );
                break;
            }
            else
            {
                buff += tmp;
            }
        }
        if( in.eof() )
        {
            all_protein.push_back( buff );
            break;
        }
    }
//    cout << all_protein.size() << endl;
//    cout << all_name.size() << endl;
//    for( int i = 0 ;i < all_protein.size() ; ++i )
//    {
//       cout << all_name[i] << "\t" << all_protein[i] << endl;
//    }
}

void protein::gen_score1()
{
    string structure_tmp;
    vector< string > structure;
    string cmd = "./lncPro/predator -a ";
    cmd += data;
    FILE *pp;
    if( ( pp = popen(cmd.c_str() , "r") ) == NULL )
    {
        cout << "popen error!\n";
    }
    char a[200];
    fgets( a , 200 , pp );
    while( !feof(pp) )
    {
        //printf( "%s" , a );
        if( a[0] == '>' )
        {
            structure_tmp.clear();
            fgets( a , 200 , pp );
            fgets( a , 200 , pp );
            while( a[0] >= '0' && a[0] <= '9' )
            {
                fgets( a , 200 , pp );
                for( int i = 0 ; i < (int)strlen( a ) ; ++i )
                {
                    if( a[i] == 'E') structure_tmp += "e";
                    else if( a[i] == 'H' ) structure_tmp += "h";
                    else if( a[i] == '_' ) structure_tmp += "c";
                }
                fgets( a, 200 , pp );
                fgets( a, 200 , pp );
                fgets( a, 200 , pp );
            }
            //cout << structure_tmp << endl;
            structure.push_back( structure_tmp );
        }
        else
        {
            fgets( a , 200 , pp );
        }
    }
    pclose(pp);
    double code[100000];
    for( int i = 0 ; i < (int)structure.size() ; ++i )
    {
       for( int j = 0 ; j < (int)structure[i].size() ; ++j )
        {
            ten = new double[10];
            if( structure[i][j] == 'h')
            {
                if( all_protein[i][j] == 'A' )
                    code[j] = 1.42;
                else if( all_protein[i][j] == 'R' )
                    code[j] = 0.98;
                else if( all_protein[i][j] == 'N' )
                    code[j] = 0.67;
                else if( all_protein[i][j] == 'D' )
                    code[j] = 1.01;
                else if( all_protein[i][j] == 'C' )
                    code[j] = 0.7;
                else if( all_protein[i][j] == 'Q' )
                    code[j] = 1.11;
                else if( all_protein[i][j] == 'E' )
                    code[j] = 1.51;
                else if( all_protein[i][j] == 'G' )
                    code[j] = 0.57;
                else if( all_protein[i][j] == 'H' )
                    code[j] = 1.0;
                else if( all_protein[i][j] == 'I' )
                    code[j] = 1.08;
                else if( all_protein[i][j] == 'L' )
                    code[j] = 1.21;
                else if( all_protein[i][j] == 'K' )
                    code[j] = 1.16;
                else if( all_protein[i][j] == 'M' )
                    code[j] = 1.45;
                else if( all_protein[i][j] == 'F' )
                    code[j] = 1.13;
                else if( all_protein[i][j] == 'P' )
                    code[j] = 0.57;
                else if( all_protein[i][j] == 'S' )
                    code[j] = 0.77;
                else if( all_protein[i][j] == 'T' )
                    code[j] = 0.83;
                else if( all_protein[i][j] == 'W' )
                    code[j] = 1.08;
                else if( all_protein[i][j] == 'Y' )
                    code[j] = 0.69;
                else if( all_protein[i][j] == 'V' )
                    code[j] = 1.06;
                else
                    code[j] = 0.0;
            }
            else if( structure[i][j] == 'c')
            {
                if( all_protein[i][j] == 'A' )
                    code[j] = 0.66;
                else if( all_protein[i][j] == 'R' )
                    code[j] = 0.95;
                else if( all_protein[i][j] == 'N' )
                    code[j] = 1.56;
                else if( all_protein[i][j] == 'D' )
                    code[j] = 1.46;
                else if( all_protein[i][j] == 'C' )
                    code[j] = 1.19;
                else if( all_protein[i][j] == 'Q' )
                    code[j] = 0.98;
                else if( all_protein[i][j] == 'E' )
                    code[j] = 0.74;
                else if( all_protein[i][j] == 'G' )
                    code[j] = 1.56;
                else if( all_protein[i][j] == 'H' )
                    code[j] = 0.95;
                else if( all_protein[i][j] == 'I' )
                    code[j] = 0.47;
                else if( all_protein[i][j] == 'L' )
                    code[j] = 0.59;
                else if( all_protein[i][j] == 'K' )
                    code[j] = 1.01;
                else if( all_protein[i][j] == 'M' )
                    code[j] = 0.6;
                else if( all_protein[i][j] == 'F' )
                    code[j] = 1.13;
                else if( all_protein[i][j] == 'P' )
                    code[j] = 1.52;
                else if( all_protein[i][j] == 'S' )
                    code[j] = 1.43;
                else if( all_protein[i][j] == 'T' )
                    code[j] = 0.96;
                else if( all_protein[i][j] == 'W' )
                    code[j] = 0.96;
                else if( all_protein[i][j] == 'Y' )
                    code[j] = 1.14;
                else if( all_protein[i][j] == 'V' )
                    code[j] = 0.5;
                else
                    code[j] = 0.0;
            }
            else if( structure[i][j] == 'e')
            {
                if( all_protein[i][j] == 'A' )
                    code[j] = 0.83;
                else if( all_protein[i][j] == 'R' )
                    code[j] = 0.93;
                else if( all_protein[i][j] == 'N' )
                    code[j] = 0.89;
                else if( all_protein[i][j] == 'D' )
                    code[j] = 0.54;
                else if( all_protein[i][j] == 'C' )
                    code[j] = 1.19;
                else if( all_protein[i][j] == 'Q' )
                    code[j] = 1.1;
                else if( all_protein[i][j] == 'E' )
                    code[j] = 0.37;
                else if( all_protein[i][j] == 'G' )
                    code[j] = 0.75;
                else if( all_protein[i][j] == 'H' )
                    code[j] = 0.87;
                else if( all_protein[i][j] == 'I' )
                    code[j] = 1.6;
                else if( all_protein[i][j] == 'L' )
                    code[j] = 1.3;
                else if( all_protein[i][j] == 'K' )
                    code[j] = 0.74;
                else if( all_protein[i][j] == 'M' )
                    code[j] = 1.05;
                else if( all_protein[i][j] == 'F' )
                    code[j] = 1.38;
                else if( all_protein[i][j] == 'P' )
                    code[j] = 0.55;
                else if( all_protein[i][j] == 'S' )
                    code[j] = 0.75;
                else if( all_protein[i][j] == 'T' )
                    code[j] = 1.19;
                else if( all_protein[i][j] == 'W' )
                    code[j] = 1.37;
                else if( all_protein[i][j] == 'Y' )
                    code[j] = 1.47;
                else if( all_protein[i][j] == 'V' )
                    code[j] = 1.7;
                else
                    code[j] = 0.0;
            }
        }
        lq( code , structure[i].size() , ten );
        score1.push_back( ten );
    }
    /* ofstream out1( "protein_fea_1.txt" );;
    for( int i = 0 ; i < 1 ; ++i )
    {
        for( int j = 0 ; j < 9 ; ++j )
            out1 << score1[i][j] << "\t";
        out1 << score1[i][9] << "\n";
    }
out1.close();
*/
}
void protein::gen_score2()
{
    double code[100000];
    for( int i = 0 ; i < (int)all_protein.size() ; ++i )
    {
        for( int j = 0 ; j < (int)all_protein[i].size() ; ++j )
        {
            ten = new double[10];
            if( all_protein[i][j] == 'A' )
                code[j] = 8.1;
            else if( all_protein[i][j] == 'R' )
                code[j] = 10.5;
            else if( all_protein[i][j] == 'N' )
                code[j] = 11.6;
            else if( all_protein[i][j] == 'D' )
                code[j] = 13.0;
            else if( all_protein[i][j] == 'C' )
                code[j] = 5.5;
            else if( all_protein[i][j] == 'Q' )
                code[j] = 10.5;
            else if( all_protein[i][j] == 'E' )
                code[j] = 12.3;
            else if( all_protein[i][j] == 'G' )
                code[j] = 9;
            else if( all_protein[i][j] == 'H' )
                code[j] = 10.4;
            else if( all_protein[i][j] == 'I' )
                code[j] = 5.2;
            else if( all_protein[i][j] == 'L' )
                code[j] = 4.9;
            else if( all_protein[i][j] == 'K' )
                code[j] = 11.3;
            else if( all_protein[i][j] == 'M' )
                code[j] = 5.7;
            else if( all_protein[i][j] == 'F' )
                code[j] = 5.2;
            else if( all_protein[i][j] == 'P' )
                code[j] = 8;
            else if( all_protein[i][j] == 'S' )
                code[j] = 9.2;
            else if( all_protein[i][j] == 'T' )
                code[j] = 8.6;
            else if( all_protein[i][j] == 'W' )
                code[j] = 5.4;
            else if( all_protein[i][j] == 'Y' )
                code[j] = 6.2;
            else if( all_protein[i][j] == 'V' )
                code[j] = 5.9;
            else
                cout << "error" << endl;
        }
        lq( code , all_protein[i].size() , ten );
        score2.push_back( ten );
    }
    /* ofstream out1( "protein_fea_2.txt" );;
    for( int i = 0 ; i < 1 ; ++i )
    {
        for( int j = 0 ; j < 9 ; ++j )
            out1 << score2[i][j] << "\t";
        out1 << score2[i][9] << "\n";
    }
out1.close();*/
//    cout << score2.size() << endl;
//    for( int i = 0 ; i < score2.size() ; ++i )
//    {
//        for( int j = 0 ; j < 10 ; ++ j)
//            cout << score2[i][j] << "\t";
//        cout << endl;
//    }
}

void protein::gen_score3()
{
    double code[100000];
    for( int i = 0 ; i < (int)all_protein.size() ; ++i )
    {
        for( int j = 0 ; j < (int)all_protein[i].size() ; ++j )
        {
            ten = new double[10];
            if( all_protein[i][j] == 'A' )
                code[j] = 0.0;
            else if( all_protein[i][j] == 'R' )
                code[j] = 52;
            else if( all_protein[i][j] == 'N' )
                code[j] = 3.38;
            else if( all_protein[i][j] == 'D' )
                code[j] = 49.7;
            else if( all_protein[i][j] == 'C' )
                code[j] = 1.48;
            else if( all_protein[i][j] == 'Q' )
                code[j] = 3.53;
            else if( all_protein[i][j] == 'E' )
                code[j] = 49.9;
            else if( all_protein[i][j] == 'G' )
                code[j] = 0.0;
            else if( all_protein[i][j] == 'H' )
                code[j] = 51.6;
            else if( all_protein[i][j] == 'I' )
                code[j] = 0.13;
            else if( all_protein[i][j] == 'L' )
                code[j] = 0.13;
            else if( all_protein[i][j] == 'K' )
                code[j] = 49.5;
            else if( all_protein[i][j] == 'M' )
                code[j] = 1.43;
            else if( all_protein[i][j] == 'F' )
                code[j] = 0.35;
            else if( all_protein[i][j] == 'P' )
                code[j] = 1.58;
            else if( all_protein[i][j] == 'S' )
                code[j] = 1.67;
            else if( all_protein[i][j] == 'T' )
                code[j] = 1.66;
            else if( all_protein[i][j] == 'W' )
                code[j] = 2.1;
            else if( all_protein[i][j] == 'Y' )
                code[j] = 1.61;
            else if( all_protein[i][j] == 'V' )
                code[j] = 0.13;
            else
                cout << "error" << endl;
        }
        lq( code , all_protein[i].size() , ten );
        score3.push_back( ten );
    }
    /* ofstream out1( "protein_fea_3.txt" );;
    for( int i = 0 ; i < 1 ; ++i )
    {
        for( int j = 0 ; j < 9 ; ++j )
            out1 << score3[i][j] << "\t";
        out1 << score3[i][9] << "\n";
    }
out1.close();*/
//    cout << score3.size() << endl;
//    for( int i = 0 ; i < score3.size() ; ++i )
//    {
//        for( int j = 0 ; j < 10 ; ++ j)
//            cout << score3[i][j] << "\t";
//        cout << endl;
//    }
}
void protein::gen_score4()
{
    double code[100000];
    for( int i = 0 ; i < (int)all_protein.size() ; ++i )
    {
        for( int j = 0 ; j < (int)all_protein[i].size() ; ++j )
        {
            ten = new double[10];
            if( all_protein[i][j] == 'A' )
                code[j] = 1.8;
            else if( all_protein[i][j] == 'R' )
                code[j] = -4.5;
            else if( all_protein[i][j] == 'N' )
                code[j] = -3.5;
            else if( all_protein[i][j] == 'D' )
                code[j] = -3.5;
            else if( all_protein[i][j] == 'C' )
                code[j] = 2.5;
            else if( all_protein[i][j] == 'Q' )
                code[j] = -3.5;
            else if( all_protein[i][j] == 'E' )
                code[j] = -3.5;
            else if( all_protein[i][j] == 'G' )
                code[j] = -0.4;
            else if( all_protein[i][j] == 'H' )
                code[j] = -3.2;
            else if( all_protein[i][j] == 'I' )
                code[j] = 4.5;
            else if( all_protein[i][j] == 'L' )
                code[j] = 3.8;
            else if( all_protein[i][j] == 'K' )
                code[j] = -3.9;
            else if( all_protein[i][j] == 'M' )
                code[j] = 1.9;
            else if( all_protein[i][j] == 'F' )
                code[j] = 2.8;
            else if( all_protein[i][j] == 'P' )
                code[j] = -1.6;
            else if( all_protein[i][j] == 'S' )
                code[j] = -0.8;
            else if( all_protein[i][j] == 'T' )
                code[j] = -0.7;
            else if( all_protein[i][j] == 'W' )
                code[j] = -0.9;
            else if( all_protein[i][j] == 'Y' )
                code[j] = -1.3;
            else if( all_protein[i][j] == 'V' )
                code[j] = 4.2;
            else
                cout << "error" << endl;
        }
        lq( code , all_protein[i].size() , ten );
        score4.push_back( ten );
    }
    /* ofstream out1( "protein_fea_4.txt" );;
    for( int i = 0 ; i < 1 ; ++i )
    {
        for( int j = 0 ; j < 9 ; ++j )
            out1 << score4[i][j] << "\t";
        out1 << score4[i][9] << "\n";
    }
out1.close();*/
//    cout << score4.size() << endl;
//    for( int i = 0 ; i < score4.size() ; ++i )
//    {
//        for( int j = 0 ; j < 10 ; ++ j)
//            cout << score4[i][j] << "\t";
//        cout << endl;
//    }
}

void protein::gen_score5()
{
    double code[100000];
    for( int i = 0 ; i < (int)all_protein.size() ; ++i )
    {
        for( int j = 0 ; j < (int)all_protein[i].size() ; ++j )
        {
            ten = new double[10];
            if( all_protein[i][j] == 'A' )
                code[j] = 0.61;
            else if( all_protein[i][j] == 'R' )
                code[j] = 0.69;
            else if( all_protein[i][j] == 'N' )
                code[j] = 0.89;
            else if( all_protein[i][j] == 'D' )
                code[j] = 0.61;
            else if( all_protein[i][j] == 'C' )
                code[j] = 0.36;
            else if( all_protein[i][j] == 'Q' )
                code[j] = 0.97;
            else if( all_protein[i][j] == 'E' )
                code[j] = 0.51;
            else if( all_protein[i][j] == 'G' )
                code[j] = 0.81;
            else if( all_protein[i][j] == 'H' )
                code[j] = 0.69;
            else if( all_protein[i][j] == 'I' )
                code[j] = -1.45;
            else if( all_protein[i][j] == 'L' )
                code[j] = -1.65;
            else if( all_protein[i][j] == 'K' )
                code[j] = 0.46;
            else if( all_protein[i][j] == 'M' )
                code[j] = -0.66;
            else if( all_protein[i][j] == 'F' )
                code[j] = -1.52;
            else if( all_protein[i][j] == 'P' )
                code[j] = -0.17;
            else if( all_protein[i][j] == 'S' )
                code[j] = 0.42;
            else if( all_protein[i][j] == 'T' )
                code[j] = 0.29;
            else if( all_protein[i][j] == 'W' )
                code[j] = -1.2;
            else if( all_protein[i][j] == 'Y' )
                code[j] = -1.43;
            else if( all_protein[i][j] == 'V' )
                code[j] = -0.75;
            else
                cout << "error" << endl;
        }
        lq( code , all_protein[i].size() , ten );
        score5.push_back( ten );
    }
     /*ofstream out1( "protein_fea.txt" );
    int protein_size = score5.size();
    for( int i = 0 ; i < protein_size ; ++i )
    {
	out1 <<all_name[i] << "\t";
        /*if (all_protein[i].length() < 30){
	    for( int j = 0 ; j < 49 ; ++j )
                out1 << 0 << "\t";
            out1 << 0 << "\n";   
	}*/
	//else{
       /* for( int j = 0 ; j < 10 ; ++j ){
            out1 << score1[i][j] << "\t";
	    out1 << score2[i][j] << "\t";
	    out1 << score3[i][j] << "\t";
            out1 << score4[i][j] << "\t";
	}
        for( int j = 0 ; j < 9 ; ++j )
            out1 << score5[i][j] << "\t";
        out1 << score5[i][9] << "\n";
       // }
    }*/
    /*for( int i = 0 ; i < protein_size ; ++i )
    {
        for( int j = 0 ; j < 10 ; ++j )
            out1 << score2[i][j] << "\t";
        //out1 << score5[i][9] << "\n";
    }
    for( int i = 0 ; i < 1 ; ++i )
    {
        for( int j = 0 ; j < 10 ; ++j )
            out1 << score3[i][j] << "\t";
        //out1 << score5[i][9] << "\n";
    }
    for( int i = 0 ; i < 1 ; ++i )
    {
        for( int j = 0 ; j < 10 ; ++j )
            out1 << score4[i][j] << "\t";
        //out1 << score5[i][9] << "\n";
    }
	for( int i = 0 ; i < 1 ; ++i )
    {
        for( int j = 0 ; j < 9 ; ++j )
            out1 << score5[i][j] << "\t";
        out1 << score5[i][9] << "\n";
    }*/
	//out1.close();
//    cout << score5.size() << endl;
//    for( int i = 0 ; i < score5.size() ; ++i )
//    {
//        for( int j = 0 ; j < 10 ; ++ j)
//            cout << score5[i][j] << "\t";
//        cout << endl;
//    }
}
void protein::lq( double *score , int len , double *ten )
{
    double tmp = 0 ;
    for( int i = 0 ; i < 10 ; ++i )
    {
        tmp = 0 ;
        for( int j = 0 ; j < len ; ++j )
        {
            tmp += double(score[j]) * cos( (double)PI * ( (double)j + 0.5 ) * ( (double)( i + 1 ) - 0.5 ) / double(len) );
        }
        ten[i] = sqrt( 2 / double(len) ) * tmp ;
    }
}
