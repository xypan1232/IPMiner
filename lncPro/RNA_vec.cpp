#include "RNA.h"
#include "fstream"
#include "iostream"
#include "cmath"

using namespace std;

#define PI 3.1415926

RNA::RNA()
{
    len_score1 = len_score2 = len_score3 = 0 ;
    len_w = 100;
}

void RNA::init( string filename )
{
    ifstream in( filename.c_str() );
    string tmp , buff ;
    getline( in , buff );
    while( getline( in , buff ) )
    {
//        cout << buff << endl;
        if( buff[0] == '>' )
            break;
        tmp += buff;
    }
    RNA_string = tmp;
//    cout << "RNA_String : " << RNA_string << endl;
    len_score1 = len_RNA_string = RNA_string.size();
    //cout << "RNA string is : " << RNA_string << endl << "length is : " << len_RNA_string << endl;
    in.close();
}

void RNA::stepFirst()
{
    //cout << "score1 is running........." << endl;
    encode1();
    //cout << "score2 is running........." << endl;
    encode2();
    //cout << "score3 is running........." << endl;
    encode3();
}

void RNA::stepSecond()
{
    lq( score1 , len_score1 , ten1 );
    lq( score2 , len_score2 , ten2 );
    lq( score3 , len_score3 , ten3 );
      /* ofstream out1( "RNA_fea.txt" );
        for( int i = 0 ; i < 10 ; ++i )
        {
            out1 << ten1[i] << "\t";

        }
        for( int i = 0 ; i < 10 ; ++i )
        {
            out1 << ten2[i] << "\t";

        }
        for( int i = 0 ; i < 9 ; ++i )
        {
            out1 << ten3[i] << "\t";

        }
	out1 << ten3[9] << "\n";
        out1.close();*/
/*if (score1)
{delete [] score1;
score1 = NULL;
}

if (score2)
{delete [] score2;
score2 = NULL;
}

if (score3)
{delete [] score3;
score3 = NULL;
}
*/
}

void RNA::stepThird()
{

    load_w( "w_1.txt" , w1 );
    load_w( "w_2.txt" , w2 );
    load_w( "w_3.txt" , w3 );
    load_w( "w_4.txt" , w4 );
    load_w( "w_5.txt" , w5 );

    load_m( "m_1_1.txt" , m11 );
    load_m( "m_1_2.txt" , m12 );
    load_m( "m_2_1.txt" , m21 );
    load_m( "m_2_2.txt" , m22 );
    load_m( "m_3_1.txt" , m31 );
    load_m( "m_3_2.txt" , m32 );
    load_m( "m_4_1.txt" , m41 );
    load_m( "m_4_2.txt" , m42 );
    load_m( "m_5_1.txt" , m51 );
    load_m( "m_5_2.txt" , m52 );

    fscore1 = compute_score( humanprotein_1 , ten1 , w1 , fscore1 , len_humanprotein_1 );
    fscore2 = compute_score( humanprotein_2 , ten2 , w2 , fscore2 , len_humanprotein_2 );
    fscore3 = compute_score( humanprotein_3 , ten2 , w3 , fscore3 , len_humanprotein_3 );
    fscore4 = compute_score( humanprotein_4 , ten3 , w4 , fscore4 , len_humanprotein_4 );
    fscore5 = compute_score( humanprotein_5 , ten3 , w5 , fscore5 , len_humanprotein_5 );



    transfer_score( fscore1 , m11 , m12 , w1 , len_humanprotein_1 );
    transfer_score( fscore2 , m21 , m22 , w2 , len_humanprotein_2 );
    transfer_score( fscore3 , m31 , m32 , w3 , len_humanprotein_3 );
    transfer_score( fscore4 , m41 , m42 , w4 , len_humanprotein_4 );
    transfer_score( fscore5 , m51 , m52 , w5 , len_humanprotein_5 );

    /*
        ofstream out1( "fscore1.txt" );
        for( int i = 0 ; i < len_humanprotein_1 ; ++i )
        {
            out1 << fscore1[i] << "\r\n";

        }
        out1.close();
        ofstream out2( "fscore2.txt" );
        for( int i = 0 ; i < len_humanprotein_1 ; ++i )
        {
            out2 << fscore2[i] << "\r\n";
        }
        out2.close();
        ofstream out3( "fscore3.txt" );
        for( int i = 0 ; i < len_humanprotein_1 ; ++i )
        {
            out3 << fscore3[i] << "\r\n";
        }
        out3.close();
        ofstream out4( "fscore4.txt" );
        for( int i = 0 ; i < len_humanprotein_1 ; ++i )
        {
            out4 << fscore4[i] << "\r\n";
        }
        out4.close();
        ofstream out5( "fscore5.txt" );
        for( int i = 0 ; i < len_humanprotein_5 ; ++i )
        {
            out5 << fscore5[i] << "\r\n";
        }
        out5.close();
    */
    cal_fscore();
}

void RNA::cal_fscore()
{
    double tmp = 0 ;
    fscore = new double[len_humanprotein_1];
    for( int i = 0 ; i < len_humanprotein_1 ; ++i )
    {
        tmp = fscore1[i] + fscore2[i] + fscore3[i] + fscore4[i] + fscore5[i];
        fscore[i] = tmp / 5;
    }
}


void RNA::transfer_score( double *score , double *m1 , double *m2 , double *w , int len )
{
    double cc1 , cc2 ;
    cc1 = cc2 = 0;
    for( int i = 0 ;  i < 100 ; ++i )
    {
        cc1 += m1[i] * w[i];
        cc2 += m2[i] * w[i];
    }
    double cc = ( cc1 + cc2 ) / 2;
    for( int i = 0 ; i < len ; ++i )
    {
        score[i] = ( 100 / PI ) * atan( 2 * ( score[i] - cc ) / ( cc1 - cc2 ) ) + 50;
        //cout << score[i] << endl;
    }
}

double *RNA::compute_score( double ** humanprotein , double *ten , double *w , double *fscore_tmp  , int len )
{
    fscore_tmp = new double[len];
    int count = 0 ;
    double score_tmp;
    double tmp[100];
    for( int i = 0 ; i < len ; ++i )
    {
        count = 0 ;
        for( int j = 0 ; j < 10 ; ++j )
        {
            for( int k = 0 ; k < 10 ; ++k )
            {
                tmp[count++] = humanprotein[i][j] * ten[k];
            }
        }

        score_tmp = 0 ;
        for( int j = 0 ; j < 100 ; ++j )
        {
            score_tmp += tmp[j] * w[j];
        }
        fscore_tmp[i] = score_tmp;
    }

    /*
    //output
    for( int i = 0 ; i < len ; ++i )
        cout << fscore_tmp[i] << endl;
    */
    return fscore_tmp;


}

void RNA::load_m(  string filename , double *m )
{
    ifstream in;
    char tmp;
    in.open( filename.c_str() , ifstream::in );
    for( int i = 0 ; i < 100 ; ++i )
    {
        in >> m[i];
        in >> tmp;
    }
    in.close();

    /*
    //output
    for( int i = 0 ; i < 100 ; ++i )
        cout << m [i] << endl;
    */
}

void RNA::load_w( string filename , double *w  )
{
    ifstream in;
    in.open( filename.c_str() , ifstream::in );
    for( int i = 0 ; i < 100 ; ++i )
        in >> w[i];
    in.close();

    /*
    //output
    for( int i = 0 ;i  < 100 ; ++i )
        cout << w[i] << endl;
    */
}

double **RNA::load_humanprotein( vector<double*> &other , double **humanprotein , int &len )
{
    len = other.size() ;

    humanprotein = new double*[len];
    for( int i = 0 ; i < len ; ++i )
        humanprotein[i] = new double[10];

    for( int i = 0 ; i < len ; ++i )
    {
        for( int j = 0 ; j < 10 ; ++j )
            humanprotein[i][j] = other[i][j];
    }


//    //output
//    cout << "len is " << len << endl;
//    for( int i = 0 ; i < len ; ++i )
//    {
//        for( int j = 0 ; j < 10 ; ++j )
//            cout << humanprotein[i][j] << " ";
//        cout << endl;
//    }
    return humanprotein;
}


void RNA::lq( vector<int> score , int len , double *ten )
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

    /*
    //output
    for( int i = 0 ; i < 10 ; ++i )
        cout << ten[i] << " ";
    cout << endl;
    */
}

void RNA::encode1()
{
    string cmd ;
    cmd = "echo " + RNA_string + " | RNAsubopt -p 6";
    //cout <<RNA_string <<endl;
	//cout <<len_RNA_string <<endl;
    char *a = new char[len_RNA_string];
score1.clear();
	//cout <<"memory a" <<endl;
    //score1 = new int[len_RNA_string];
	//cout <<"memory score1" <<endl;
    for( int i = 0 ; i < len_RNA_string ; ++i )
        score1.push_back( 0);
    //cout <<"assign score1";
    FILE *pp;
    if( ( pp = popen(cmd.c_str() , "r") ) == NULL )
    {
        cout << "popen error!\n";
    }
    for( int i = 0 ; i < 6 ; ++ i )
    {
        fscanf( pp , "%s" , a );
        //cout << a << endl;
        for( int j = 0 ; j < len_RNA_string ; ++j )
        {
            if( a[j] == '.')
                continue;
            score1[j] ++;
        }
    }
    pclose( pp );
//cout <<"subplot2";
    /*
    //output score1
    for( int i = 0 ; i < len_RNA_string ; ++i )
        cout << score1[i];
    cout << endl;
    */
}


//RNA_Hydro.cpp
void RNA::encode2()
{
    score2.clear();
    //score2 = new int[50000];
    for( int i = 0 ; i < len_RNA_string ; ++i )
    {
        if( RNA_string[i] == 'A')
        {
            score2.push_back(24);
            score2.push_back(17);
            score2.push_back(40);
            score2.push_back(10);
        }
        else if( RNA_string[i] == 'G' )
        {
            score2.push_back(21);
            score2.push_back(86);
            score2.push_back(17);
            score2.push_back(41);
            score2.push_back(29);
        }
        else if( RNA_string[i] == 'C' )
        {
            score2.push_back(49);
            score2.push_back(21);
            score2.push_back(26);

        }
        else
        {
            score2.push_back(24);
            score2.push_back(17);
            score2.push_back(22);
        }
    }
len_score2 = score2.size();
    /*
    //output
    for( int i = 0 ; i < len_score2 ; ++i )
    {
        cout << score2[i] << " ";
    }
    cout << endl << "len is : " << len_score2 << endl;
    */
}

//RNA_Van.cpp
void RNA::encode3()
{
    score3.clear();
    //score3 = new int[50000];
    for( int i = 0 ; i < len_RNA_string ; ++i )
    {
        if( RNA_string[i] == 'A')
        {
            score3.push_back(79);
            score3.push_back(98);
            score3.push_back(69);
            score3.push_back(40);
            score3.push_back(53);
            score3.push_back(37);
            score3.push_back(84);
            score3.push_back(62);
            score3.push_back(49);
            score3.push_back(28);
        }
        else if( RNA_string[i] == 'G' )
        {
            score3.push_back(26);
            score3.push_back(74);
            score3.push_back(24);
            score3.push_back(37);
            score3.push_back(22);
            score3.push_back(21);
            score3.push_back(19);
            score3.push_back(67);
            score3.push_back(48);
            score3.push_back(44);
            score3.push_back(21);
        }
        else if( RNA_string[i] == 'C' )
        {
            score3.push_back(14);
            score3.push_back(44);
            score3.push_back(98);
            score3.push_back(42);
            score3.push_back(30);
            score3.push_back(50);
            score3.push_back(39);
            score3.push_back(19);
        }
        else
        {
            score3.push_back(25);
            score3.push_back(42);
            score3.push_back(74);
            score3.push_back(53);
            score3.push_back(43);
            score3.push_back(67);
            score3.push_back(44);
            score3.push_back(24);
        }
    }
    len_score3 = score3.size();
    /*
    //output
    for( int i = 0 ; i < len_score3 ; ++i )
    {
        cout << score3[i] << " ";
    }
    cout << endl << "len is : " << len_score3 << endl;
    */
}

void RNA::output()
{
    for( int i = 0 ; i < len_humanprotein_1 ; ++i )
    {
        cout << protein_name[i] << "\t" << fscore[i] << "\r\n";
    }
}
