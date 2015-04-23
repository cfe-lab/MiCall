#include <string>

#ifdef __PYTHON__
    #include <Python.h>
#else
    #include "ruby.h"
#endif

using namespace std;

/*
I think this application should do a complete alignment.  Unfortunately aligning is way to slow in perl,
so I suspect the alignment, merging and possibly the gap widening should be done in c.  Another possibility
is to call the C functions from perl, which would simpify things quite a bit!  Unfortunately, I'm not entirely
confident in that perl can do this seamlessly(unlike nicer languages like ruby and python).
*/

void trim(string* seq);

static int nucMat[127][127]; // ASCII runs from 0 to 127
void init_pairscore(int matchscore, int mismatchPenalty)
{
	for (int i=0; i<127; i++)
	{
		for (int j=0; j<127; j++)
		{
			if (i==j)
			{
				nucMat[i][j]=matchscore;
			}
			else
			{
				nucMat[i][j]=-mismatchPenalty;
//				if ((char)i=='N' || (char)i=='n' || (char)j=='N' || (char)j=='n')
//				{
//					nucMat[i][j]=-mismatchPenalty;
//				}
			}
		}
	}

	// adjust naive assignments for case-insensitivity
	nucMat['a']['A']=nucMat['A']['a']=matchscore;
	nucMat['c']['C']=nucMat['C']['c']=matchscore;
	nucMat['g']['G']=nucMat['G']['g']=matchscore;
	nucMat['t']['T']=nucMat['T']['t']=nucMat['u']['U']=nucMat['U']['u']=matchscore;
	nucMat['t']['u']=nucMat['t']['U']=nucMat['T']['u']=nucMat['T']['U']=matchscore;
	nucMat['u']['t']=nucMat['t']['T']=nucMat['U']['t']=nucMat['U']['T']=matchscore;
	nucMat['N']['N']=nucMat['n']['N']=nucMat['N']['n']=0;


	//bi-mixtures
	nucMat['A']['R']=nucMat['R']['A']=matchscore;
	nucMat['G']['R']=nucMat['R']['G']=matchscore;

	nucMat['C']['Y']=nucMat['Y']['C']=matchscore;
	nucMat['T']['Y']=nucMat['Y']['T']=matchscore;

	nucMat['G']['K']=nucMat['K']['G']=matchscore;
	nucMat['T']['K']=nucMat['K']['T']=matchscore;

	nucMat['C']['M']=nucMat['M']['C']=matchscore;
	nucMat['A']['M']=nucMat['M']['A']=matchscore;

	nucMat['C']['S']=nucMat['S']['C']=matchscore;
	nucMat['G']['S']=nucMat['S']['G']=matchscore;

	nucMat['T']['W']=nucMat['W']['T']=matchscore;
	nucMat['A']['W']=nucMat['W']['A']=matchscore;

	//tri-mixtures
	nucMat['C']['B']=nucMat['B']['C']=matchscore;
	nucMat['G']['B']=nucMat['B']['G']=matchscore;
	nucMat['T']['B']=nucMat['B']['T']=matchscore;

	nucMat['A']['D']=nucMat['D']['A']=matchscore;
	nucMat['G']['D']=nucMat['D']['G']=matchscore;
	nucMat['T']['D']=nucMat['D']['T']=matchscore;

	nucMat['A']['H']=nucMat['H']['A']=matchscore;
	nucMat['C']['H']=nucMat['H']['C']=matchscore;
	nucMat['T']['H']=nucMat['H']['T']=matchscore;

	nucMat['A']['V']=nucMat['V']['A']=matchscore;
	nucMat['C']['V']=nucMat['V']['C']=matchscore;
	nucMat['G']['V']=nucMat['V']['G']=matchscore;

   //Wild cards
   nucMat['*']['A']=nucMat['*']['a']=nucMat['A']['*']=nucMat['a']['*']=matchscore;
	nucMat['*']['C']=nucMat['*']['c']=nucMat['C']['*']=nucMat['c']['*']=matchscore;
	nucMat['*']['T']=nucMat['*']['t']=nucMat['T']['*']=nucMat['t']['*']=matchscore;
	nucMat['*']['G']=nucMat['*']['g']=nucMat['G']['*']=nucMat['g']['*']=matchscore;

   nucMat['$']['$']=50;
//    nucMat['$']['A']=nucMat['$']['a']=nucMat['A']['$']=nucMat['a']['$']=0;
//	nucMat['$']['T']=nucMat['$']['t']=nucMat['T']['$']=nucMat['t']['$']=0;
//	nucMat['$']['G']=nucMat['$']['g']=nucMat['G']['$']=nucMat['g']['$']=0;

	//For those annoying duplicate phred values.
	nucMat['.']['A']=nucMat['.']['a']=nucMat['A']['.']=nucMat['a']['.']=-20;
	nucMat['.']['C']=nucMat['.']['c']=nucMat['C']['.']=nucMat['c']['.']=-20;
	nucMat['.']['T']=nucMat['.']['t']=nucMat['T']['.']=nucMat['t']['.']=-20;
	nucMat['.']['G']=nucMat['.']['g']=nucMat['G']['.']=nucMat['g']['.']=-20;

	nucMat['N']['A']=nucMat['N']['a']=nucMat['A']['N']=nucMat['a']['N']=-3;
	nucMat['N']['C']=nucMat['N']['c']=nucMat['C']['N']=nucMat['c']['N']=-3;
	nucMat['N']['T']=nucMat['N']['t']=nucMat['T']['N']=nucMat['t']['N']=-3;
	nucMat['N']['G']=nucMat['N']['g']=nucMat['G']['N']=nucMat['g']['N']=-3;
	
	//for easy alignment to a standard with gaps
	nucMat['X']['A']=nucMat['X']['a']=nucMat['A']['X']=nucMat['a']['X']=-6;
	nucMat['X']['C']=nucMat['X']['c']=nucMat['C']['X']=nucMat['c']['X']=-6;
	nucMat['X']['T']=nucMat['X']['t']=nucMat['T']['X']=nucMat['t']['X']=-6;
	nucMat['X']['G']=nucMat['X']['g']=nucMat['G']['X']=nucMat['g']['X']=-6;
	nucMat['X']['R']=nucMat['X']['r']=nucMat['R']['X']=nucMat['r']['X']=-6;
	nucMat['X']['Y']=nucMat['X']['y']=nucMat['Y']['X']=nucMat['y']['X']=-6;
	nucMat['X']['K']=nucMat['X']['k']=nucMat['K']['X']=nucMat['k']['X']=-6;
	nucMat['X']['M']=nucMat['X']['m']=nucMat['M']['X']=nucMat['m']['X']=-6;
	nucMat['X']['S']=nucMat['X']['s']=nucMat['S']['X']=nucMat['s']['X']=-6;
	nucMat['X']['W']=nucMat['X']['w']=nucMat['W']['X']=nucMat['w']['X']=-6;
	nucMat['X']['B']=nucMat['X']['b']=nucMat['B']['X']=nucMat['b']['X']=-6;
	nucMat['X']['D']=nucMat['X']['d']=nucMat['D']['X']=nucMat['d']['X']=-6;
	nucMat['X']['H']=nucMat['X']['h']=nucMat['H']['X']=nucMat['h']['X']=-6;
	nucMat['X']['V']=nucMat['X']['v']=nucMat['V']['X']=nucMat['v']['X']=-6;
	nucMat['X']['-']=nucMat['X']['-']=3;
}


void init_pairscore_aa(int matchscore, int mismatchPenalty)
{
	for (int i=0; i<127; i++)
	{
		for (int j=0; j<127; j++)
		{
			if(i==j)
			{
				nucMat[i][j]=matchscore;
			}
			else
			{
				nucMat[i][j]=-mismatchPenalty;
               if((char)i=='X' || (char)j=='X')
               {
                   nucMat[i][j]=-4;    
               }
			}
		}
	}

	nucMat['Z']['Z']=nucMat['z']['Z']=nucMat['Z']['z']=0;
   nucMat['X']['-']=nucMat['-']['X']=matchscore;
}


/*
    Empirical score matrix based on 25% divergent HIV sequences
    See Nickle, David C., et al. "HIV-specific probabilistic models of protein evolution."
        PLoS One 2.6 (2007): e503.
*/
static int empirical_hiv25[24][24] = {\
{7,-7,-7,-4,-10,-11,-4,-3,-10,-6,-9,-9,-7,-13,-3,-2,1,-16,-15,0,-5,-5,-3,-17},\
{-7,7,-5,-11,-8,-2,-7,-2,0,-6,-6,2,-3,-12,-4,-2,-2,-5,-9,-10,-7,-3,-3,-17},\
{-7,-5,8,2,-9,-6,-6,-7,0,-6,-12,0,-10,-12,-9,1,0,-17,-3,-10,6,-6,-3,-17},\
{-4,-11,2,8,-14,-10,0,-2,-3,-11,-15,-7,-13,-15,-13,-5,-6,-16,-6,-5,7,0,-3,-17},\
{-10,-8,-9,-14,11,-16,-15,-5,-7,-11,-9,-13,-14,0,-12,-1,-6,-2,0,-8,-10,-16,-5,-17},\
{-11,-2,-6,-10,-16,8,-2,-10,0,-12,-4,0,-8,-12,-1,-9,-8,-14,-9,-13,-7,6,-4,-17},\
{-4,-7,-6,0,-15,-2,7,-1,-9,-12,-15,-1,-10,-17,-13,-11,-8,-15,-12,-5,0,6,-4,-17},\
{-3,-2,-7,-2,-5,-10,-1,7,-10,-11,-14,-6,-12,-9,-11,-1,-7,-5,-14,-5,-4,-3,-4,-17},\
{-10,0,0,-3,-7,0,-9,-10,10,-10,-4,-5,-10,-6,-3,-6,-6,-11,2,-14,-1,-2,-3,-17},\
{-6,-6,-6,-11,-11,-12,-12,-11,-10,7,0,-7,0,-2,-10,-4,0,-14,-9,2,-7,-12,-2,-17},\
{-9,-6,-12,-15,-9,-4,-15,-14,-4,0,6,-10,0,0,-3,-5,-8,-6,-8,-4,-13,-6,-4,-17},\
{-9,2,0,-7,-13,0,-1,-6,-5,-7,-10,7,-4,-14,-9,-5,-1,-12,-13,-9,-1,-1,-2,-17},\
{-7,-3,-10,-13,-14,-8,-10,-12,-10,0,0,-4,10,-7,-11,-9,-1,-11,-15,0,-11,-9,-3,-17},\
{-13,-12,-12,-15,0,-12,-17,-9,-6,-2,0,-14,-7,10,-11,-5,-10,-5,1,-5,-13,-14,-3,-17},\
{-3,-4,-9,-13,-12,-1,-13,-11,-3,-10,-3,-9,-11,-11,8,-1,-3,-13,-11,-12,-10,-3,-5,-17},\
{-2,-2,1,-5,-1,-9,-11,-1,-6,-4,-5,-5,-9,-5,-1,8,0,-12,-6,-9,0,-10,-3,-17},\
{1,-2,0,-6,-6,-8,-8,-7,-6,0,-8,-1,-1,-10,-3,0,7,-16,-10,-4,-2,-8,-2,-17},\
{-16,-5,-17,-16,-2,-14,-15,-5,-11,-14,-6,-12,-11,-5,-13,-12,-16,10,-4,-16,-16,-14,-8,-17},\
{-15,-9,-3,-6,0,-9,-12,-14,2,-9,-8,-13,-15,1,-11,-6,-10,-4,10,-12,-4,-10,-4,-17},\
{0,-10,-10,-5,-8,-13,-5,-5,-14,2,-4,-9,0,-5,-12,-9,-4,-16,-12,7,-7,-7,-3,-17},\
{-5,-7,6,7,-10,-7,0,-4,-1,-7,-13,-1,-11,-13,-10,0,-2,-16,-4,-7,7,-2,-4,-17},\
{-5,-3,-6,0,-16,6,6,-3,-2,-12,-6,-1,-9,-14,-3,-10,-8,-14,-10,-7,-2,6,-4,-17},\
{-3,-3,-3,-3,-5,-4,-4,-4,-3,-2,-4,-2,-3,-3,-5,-3,-2,-8,-4,-3,-4,-4,-3,-17},\
{-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,1}};

void init_pairscore_hiv25(void) {
    // ASCII codes for protein alphabet ARNDCQEGHILKMFPSTWYVBZ?*
    int aa_to_ascii[24] = { 65, 82, 78, 68, 67, 81, 69, 71, 72, 73, 76, 75, 77, 70, 80, 83, 84,
    87, 89, 86, 66, 90, 63, 42 };
    int i2, j2;

    // reset score matrix to be safe
    for (int i=0; i<127; i++) {
        for (int j=0; j<127; j++) {
            nucMat[i][j] = 0;
        }
    }

    // map HIV 25% empirical matrix to score matrix
    for (int i=0; i<24; i++) {
        i2 = aa_to_ascii[i];
        for (int j=0; j<24; j++) {
            j2 = aa_to_ascii[j];
            // also map to lowercase
            nucMat[i2+32][j2+32] = nucMat[i2+32][j2] = nucMat[i2][j2+32] = nucMat[i2][j2] = empirical_hiv25[i][j];
        }
    }
}


extern int pairscore(char a, char b)
{
	return nucMat[a][b];
}

void reverse(string* seq)
{
	string tmp = "";
	for(int i = seq->size() - 1; i >= 0; --i)
	{
		tmp += (*seq)[i];
	}
	*seq = tmp;
}


//Error must be somewhere in here.  Ug...
int align(string* seqa, string* seqb, string* newseqa, string* newseqb,
          int gip, int gep, int use_terminal_gap_penalty)
{
    /*
     Pairwise alignment with affine gap penalty.
     see Gotoh, Osamu. "Optimal alignment between groups of sequences and its application
     to multiple sequence alignment." Computer applications in the biosciences: CABIOS 9.3
     (1993): 361-370.

     Gap open and extension penalties [gip] and [gep] are assumed to take positive values.
    */

	int M = seqa->size(); // first group of pre-aligned sequences
	int N = seqb->size(); // second group

	int i, j;

    // not all elements of D, P, and Q need to be stored - vectors are adequate
	int *SS=new int[N+1];  // D(i, .)
	int *oldSS=new int[N+1];  // D(i-1, .)
	int *PP = new int[N+1];  // P(i, .)

    // Gotoh traceback matrices
    int **piSS = new int*[M+1];
	int **pjSS = new int*[M+1];

	int u = -gip; // affine gap initiation penalty
	int v = -gep; // affine gap extension penalty

	int w1 = u + v; // gap weight w_k = v * k + u for k = 1
	int t = u;
	int s, q;

	// initialize vectors
	for (j=0; j<N+1; j++)
	{
		SS[j]=0;
		oldSS[j]=0;
		PP[j]=0;
	}

    // initialize traceback matrices
	piSS[0] = new int[N+1];
	pjSS[0] = new int[N+1];
	piSS[1] = new int[N+1];
	pjSS[1] = new int[N+1];
	piSS[1][0] = 0;
	pjSS[1][0] = 0;
	piSS[0][1] = 0;
	pjSS[0][1] = 0;

	int maxiS = -100000;
	int maxjS = -100000;
	int maxij, maxji;

	for (i=1; i < M+1; i++)
	{
		t += v; // update gap extension
		s = t;
		SS[0]=0;
		q = t + u;

        // add new rows
		if (i>1)
		{
			piSS[i] = new int[N+1];
			pjSS[i] = new int[N+1];
		}

		for (j = 1; j < N + 1; j++)
		{
		    // recursive calculation of Q
			if (q >= s + u )
				q += v; // extension
			else
				q = s + u + v; // open

            // recursive calculation of P
			if ((oldSS[j] + w1) > (PP[j] + v))
				PP[j] = oldSS[j] + w1;
			else
				PP[j] += v;

			int tmp_pp = PP[j];

			// D(i-1, j-1) + d(a_i, b_j)
			int pscore = oldSS[j - 1] + pairscore((*seqa)[i - 1], (*seqb)[j - 1]);

            //no idea if this will work, but its supposed to be a stop codon aligner
            //the bonus is  assigned on the last codon, if the codons between don't make a big difference it'll be wrong.  Hrm.

            if((*seqa)[i-3] == '$' && (*seqa)[i-2] == '$' && (*seqa)[i-1] == '$'  &&
               (((*seqb)[j-3] == 'T' && (*seqb)[j-2] == 'A' && (*seqb)[j-1] == 'G') ||
               ((*seqb)[j-3] == 'T' && (*seqb)[j-2] == 'A' && (*seqb)[j-1] == 'A') ||
               ((*seqb)[j-3] == 'T' && (*seqb)[j-2] == 'G' && (*seqb)[j-1] == 'A') ))
            {
               pscore += 6;
            }
            if((*seqa)[i-2] == '$' && (*seqa)[i-1] == '$' && (*seqa)[i-0] == '$'  &&
               (((*seqb)[j-2] == 'T' && (*seqb)[j-1] == 'A' && (*seqb)[j-0] == 'G') ||
               ((*seqb)[j-2] == 'T' && (*seqb)[j-1] == 'A' && (*seqb)[j-0] == 'A') ||
               ((*seqb)[j-2] == 'T' && (*seqb)[j-1] == 'G' && (*seqb)[j-0] == 'A') ))
            {
               pscore += 6;
            }
            if((*seqa)[i-1] == '$' && (*seqa)[i-0] == '$' && (*seqa)[i+1] == '$'  &&
               (((*seqb)[j-1] == 'T' && (*seqb)[j-0] == 'A' && (*seqb)[j+1] == 'G') ||
               ((*seqb)[j-1] == 'T' && (*seqb)[j-0] == 'A' && (*seqb)[j+1] == 'A') ||
               ((*seqb)[j-1] == 'T' && (*seqb)[j-0] == 'G' && (*seqb)[j+1] == 'A') ))
            {
               pscore += 6;
            }

            /*
             D(i,j) = Min { D(i-1, j-1) + d(a_i, b_j), P(i,j), Q(i,j) }
			 where P(i,j) = Min { D(i-k, j) + w_k } for k = 1, .., i
			   and Q(i,j) = Min { D(i, j-k) + w_k } for k = 1, ..., j

			 i.e., three options are:
			  1. match/mismatch,
			  2. gap open/extension in sequence (a),
			  3. gap open/extension in sequence (b)

			 pscore = D(i-1, j-1) + d(a_i, b_j)
			 tmp_pp = P(i,j)
			 q = Q(i,j)
            */

			//maybe just >?
			if (tmp_pp >= pscore)
			{
				if (tmp_pp > q)
				{
				    // gap open / extension in (a)
					s = tmp_pp;
					piSS[i][j] = i - 1;
					pjSS[i][j] = j;
				}
				else // q > tmp_pp > pscore
				{
				    // gap open / extension in (b)
					s = q;
					piSS[i][j] = i;
					pjSS[i][j] = j - 1;
				}
			}
			else // pscore > tmp_pp)
			{
				if (pscore > q)
				{
				    // match / mismatch
					s = pscore;
					piSS[i][j] = i - 1;
					pjSS[i][j] = j - 1;
				}
				else // q > pscore > tmp_pp
				{
				    // gap open / extension in (b)
					s = q;
					piSS[i][j] = i;
					pjSS[i][j] = j - 1;
				}
			}

			SS[j] = s;

			if (i == M && SS[j] >= maxiS)
			{
				maxiS = SS[j];
				maxij = j;
			}
		}

		if (SS[N] >= maxjS)
		{
			maxjS = SS[N];
			maxji = i;
		}

		for (j = 0; j < N + 1; j++)
		{
			oldSS[j] = SS[j];
		}
	}

	if (maxij>N)
		maxij=N;
	if (maxji>M)
		maxji=M;
	if (maxij<0)
		maxij=0;
	if (maxji<0)
		maxji=0;

	//add starting -'s
	int alignment_score;
	if (maxiS > maxjS)
	{
	    alignment_score = maxiS;
		i = M;
		j = maxij;
		for (int kk = N; kk > maxij; kk--)
		{
			*newseqb += (*seqb)[kk - 1];
			*newseqa += '-';
		}
	}
   else
	{
	    alignment_score = maxjS;
		i = maxji;
		j = N;
		for (int kk = M; kk > maxji; kk--)
		{
			*newseqa += (*seqa)[kk - 1];
			*newseqb += '-';
		}
	}

	bool decI = false;
	bool decJ = false;
	//inserting -'s in the middle!
	while(i >= 1 && j >= 1)
	{
		decI=false;
		decJ=false;
		if (piSS[i][j] < i)
		{
			*newseqa += (*seqa)[i - 1];
			decI = true;
		}
		else
		{
			*newseqa += '-';
		}

		if (pjSS[i][j] < j)
		{
			*newseqb += (*seqb)[j - 1];
			decJ=true;
		}
		else
		{
			*newseqb += '-';
		}

		if (decI)
		{
			i--;
		}
		if (decJ)
		{
			j--;
		}
	}

	//add extra trailing -'s
	//forgive terminal gap penalties if user specifies this option
	if (i < j)
	{
		for (int jj = j; jj >= 1; jj--)
		{
			*newseqb += (*seqb)[jj - 1];
			*newseqa += '-';
			if (use_terminal_gap_penalty==0) alignment_score += gep;
		}
		if (use_terminal_gap_penalty==0) alignment_score += gip;
	}
	else if(i > j)
	{
		for (int ii = i; ii >= 1; ii--)
		{
			*newseqa += (*seqa)[ii - 1];
			*newseqb += '-';
			if (use_terminal_gap_penalty==0) alignment_score += gep;
		}
		if (use_terminal_gap_penalty==0) alignment_score += gip;
	}

	reverse(newseqa);
	reverse(newseqb);

	for (i = 0; i < M + 1; i++)
	{
		delete []piSS[i];
		delete []pjSS[i];
	}

	delete []SS;
	delete []oldSS;
	delete []piSS;
	delete []pjSS;
	delete []PP;
	return alignment_score;
}

void degap(string* seq)
{
    /*
    Remove pre-existing gap characters from sequences prior to alignment.
    */
	unsigned int pos = 0;
	while(pos != -1)
	{
		pos = seq->find('-', 0);
		if(pos != -1)
		{
			seq->erase(pos, 1);
		}
	}
}

void trim(string* seq)
{
    /*
    Remove trailing whitespace from sequences.
    */
	while((*seq)[0] == ' ' || (*seq)[0] == '\t' || (*seq)[0] == '\n' || (*seq)[0] == '\r')
	{
		seq->erase(0, 1);
	}

	while((*seq)[seq->size() - 1] == ' ' || (*seq)[seq->size() - 1] == '\t' || (*seq)[seq->size() - 1] == '\n'  || (*seq)[seq->size() - 1] == '\r')
	{
		seq->erase(seq->size() - 1, 1);
	}
}


void widen_gaps(string* seq)
{
	int size = seq->size();
	for(int i = 0; i < size; i++)
	{
		if((*seq)[i] == '-')
		{ //start searching for gaps to cluster

			//backwards, seqa
			unsigned int j = i - 1;
			int letter = (*seq)[j];
			j--;
			while(j >= 0)
			{
				if((*seq)[j] == '-')
				{
					//woo, swap this with i - 1
					(*seq)[j] = letter;
					(*seq)[i - 1] = '-';
					break;
				}
				else if((*seq)[j] == letter)
				{
					//nothing really
				}
				else if((*seq)[j] != letter)
				{
					break;
				}
				j--;
			}


			//forward, seqa
			j = i + 1;
			letter = (*seq)[j];
			j++;
			while(j < seq->size())
			{
				if((*seq)[j] == '-')
				{
					//woo, swap this with i + 1
					(*seq)[j] = letter;
					(*seq)[i + 1] = '-';
					break;
				}
				else if((*seq)[j] == letter)
				{
					//nothing really
				}
				else if((*seq)[j] != letter)
				{
					break;
				}
				j++;
			}
		}
	}
}

#ifdef __PYTHON__
    /* Python wrapper functions */
    static PyObject * align_it(PyObject * self, PyObject * args)
    {
        const char * standard;
        const char * seq;
        int gap_init_penalty;
        int gap_extend_penalty;
        int use_terminal_gap_penalty;
        int score;

        if (!PyArg_ParseTuple(args, "ssiii", &standard, &seq, &gap_init_penalty, &gap_extend_penalty, &use_terminal_gap_penalty)) {
            return NULL;
        }

        init_pairscore(5, 4); // match, mismatch scores +5, -4 respectively (HyPhy defaults)

        string* seqa = new string(standard);
        string* seqb = new string(seq);
        trim(seqa);
        trim(seqb);
        //degap(seqa);
        //degap(seqb);
        string* newseqa = new string();
        string* newseqb = new string();

        score = align(seqa, seqb, newseqa, newseqb, gap_init_penalty, gap_extend_penalty, use_terminal_gap_penalty);

        PyObject * retval = Py_BuildValue("ssi", newseqa->c_str(), newseqb->c_str(), score);

        delete seqa;
        delete seqb;
        delete newseqa;
        delete newseqb;

        return retval;
    }

    static PyObject * align_it_aa(PyObject * self, PyObject * args)
    {
        const char * standard;
        const char * seq;
        int gap_init_penalty;
        int gap_extend_penalty;
        int use_terminal_gap_penalty;
        int score;

        if (!PyArg_ParseTuple(args, "ssiii", &standard, &seq, &gap_init_penalty, &gap_extend_penalty, &use_terminal_gap_penalty)) {
            return NULL;
        }

        init_pairscore_hiv25();

        string* seqa = new string(standard);  // reference
        string* seqb = new string(seq);  // query
        trim(seqa);
        trim(seqb);
        //degap(seqa);  // HyPhy behaviour is to not remove gaps
        //degap(seqb);
        string* newseqa = new string();
        string* newseqb = new string();

        score = align(seqa, seqb, newseqa, newseqb, gap_init_penalty, gap_extend_penalty, use_terminal_gap_penalty);

        PyObject * retval = Py_BuildValue("ssi", newseqa->c_str(), newseqb->c_str(), score);
        delete seqa;
        delete seqb;
        delete newseqa;
        delete newseqb;

        return retval;
    }

    static PyMethodDef AlignmentMethods [] =
    {
        {"align_it", align_it, METH_VARARGS, "Pairwise alignment of nucleotide sequences."},
        {"align_it_aa", align_it_aa, METH_VARARGS, "Pairwise alignment of protein sequences using empirical HIV 25% score matrix."},
        {NULL, NULL, 0, NULL}
    };

    PyMODINIT_FUNC initalignment (void) {
        (void) Py_InitModule("alignment", AlignmentMethods);
    }

#else
    /* Ruby wrapper functions */
    extern "C" VALUE align_it(VALUE self, VALUE standard, VALUE seq, VALUE gap_init_penalty, VALUE gap_extend_penalty)
    {
       init_pairscore(1, 1);

       string* seqa = new string(RSTRING_PTR(standard));
       string* seqb = new string(RSTRING_PTR(seq));
       trim(seqa);
        trim(seqb);
       degap(seqa);
        degap(seqb);
       string* newseqa = new string();
       string* newseqb = new string();
       align(seqa, seqb, newseqa, newseqb, NUM2INT(gap_init_penalty), NUM2INT(gap_extend_penalty), 0);

       VALUE ret = rb_ary_new3(2, rb_str_new2(newseqa->c_str()),rb_str_new2(newseqb->c_str()));

       delete seqa;
       delete seqb;
       delete newseqa;
       delete newseqb;

       return ret;
    }

    extern "C" VALUE align_it_aa(VALUE self, VALUE standard, VALUE seq, VALUE gap_init_penalty, VALUE gap_extend_penalty)
    {
       init_pairscore_aa(4, -2);

       string* seqa = new string(RSTRING_PTR(standard));
       string* seqb = new string(RSTRING_PTR(seq));
       trim(seqa);
        trim(seqb);
       degap(seqa);
        degap(seqb);
       string* newseqa = new string();
       string* newseqb = new string();
       align(seqa, seqb, newseqa, newseqb, NUM2INT(gap_init_penalty), NUM2INT(gap_extend_penalty), 0);

       VALUE ret = rb_ary_new3(2, rb_str_new2(newseqa->c_str()),rb_str_new2(newseqb->c_str()));

       delete seqa;
       delete seqb;
       delete newseqa;
       delete newseqb;

       return ret;
    }


    extern "C" void Init_alignment()
    {
       rb_define_global_function("align_it", (VALUE(*)(...))align_it, 4);
       rb_define_global_function("align_it_aa", (VALUE(*)(...))align_it_aa, 4);
    }

#endif
