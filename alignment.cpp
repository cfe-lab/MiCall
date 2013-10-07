#include <string>
#include "ruby.h"
using namespace std;

/*
I think this application should do a complete alignment.  Unfortunately aligning is way to slow in perl,
so I suspect the alignment, merging and possibly the gap widening should be done in c.  Another possibility
is to call the C functions from perl, which would simpify things quite a bit!  Unfortunately, I'm not entirely
confident in that perl can do this seamlessly(unlike nicer languages like ruby and python).
*/

void trim(string* seq);

static int nucMat[127][127];
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
int align(string* seqa, string* seqb, string* newseqa, string* newseqb, int gip, int gep)
{
	int M = seqa->size();
	int N = seqb->size();

	int i, j;

	int *SS=new int[N+1];
	int *oldSS=new int[N+1];
	int **piSS = new int*[M+1];
	int **pjSS = new int*[M+1];
	int *PP = new int[N+1];


	int u = -gip; // affine gap initiation penalty
	int v = -gep; // affine gap extension penalty

	int w1 = u + v;
	int t = u;
	int s, q;
	for (j=0; j<N+1; j++)
	{
		SS[j]=0;
		oldSS[j]=0;
		PP[j]=0;
	}

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
		t += v;
		s = t;
		SS[0]=0;
		q = t + u;

		if (i>1)
		{
			piSS[i] = new int[N+1];
			pjSS[i] = new int[N+1];
		}

		for (j = 1; j < N + 1; j++)
		{
			if (q >= s + u )
				q += v;
			else
				q = s + u + v;

			if ((oldSS[j] + w1) > (PP[j] + v))
				PP[j] = oldSS[j] + w1;
			else
				PP[j] += v;

			int tmp_pp = PP[j];
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


			//maybe just >?
			if (tmp_pp >= pscore)
			{
				if (tmp_pp > q)
				{
					s = tmp_pp;
					piSS[i][j] = i - 1;
					pjSS[i][j] = j;
				}
				else
				{
					s = q;
					piSS[i][j] = i;
					pjSS[i][j] = j - 1;
				}
			}
			else
			{
				if (pscore > q)
				{
					s = pscore;
					piSS[i][j] = i - 1;
					pjSS[i][j] = j - 1;
				}
				else
				{
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
	if (maxiS > maxjS)
	{
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
	if (i < j)
	{
		for (int jj = j; jj >= 1; jj--)
		{
			*newseqb += (*seqb)[jj - 1];
			*newseqa += '-';
		}
	}
	else if(i > j)
	{
		for (int ii = i; ii >= 1; ii--)
		{
			*newseqa += (*seqa)[ii - 1];
			*newseqb += '-';
		}
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
	return 0;
}

void degap(string* seq)
{
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

extern "C" VALUE align_it(VALUE self, VALUE standard, VALUE seq, VALUE gap_init_penalty, VALUE gap_extend_penalty)
{
   init_pairscore(1, 1);

   string* seqa = new string(RSTRING(standard)->ptr);
   string* seqb = new string(RSTRING(seq)->ptr);
   trim(seqa);
	trim(seqb);
   degap(seqa);
	degap(seqb);
   string* newseqa = new string();
   string* newseqb = new string();
   align(seqa, seqb, newseqa, newseqb, NUM2INT(gap_init_penalty), NUM2INT(gap_extend_penalty));

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

   string* seqa = new string(RSTRING(standard)->ptr);
   string* seqb = new string(RSTRING(seq)->ptr);
   trim(seqa);
	trim(seqb);
   degap(seqa);
	degap(seqb);
   string* newseqa = new string();
   string* newseqb = new string();
   align(seqa, seqb, newseqa, newseqb, NUM2INT(gap_init_penalty), NUM2INT(gap_extend_penalty));

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
