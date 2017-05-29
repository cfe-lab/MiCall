"""
This class loads an ASI2 algorithm XML file and uses the loaded algorithm to
interpret amino acid strings into resistance scores.

#loading algorithms:
algorithm_hivdb = AsiAlgorithm("HIVDB_6.1.1F.xml")
algorithm_anrs = AsiAlgorithm("ANRS_max.xml")
algorithm_rega = AsiAlgorithm("RegaInst_max.xml")

#interpretting results:
res = algorithm_hivdb.interpret(aa_seq, 'RT')

The interpret call expects an amino acid sequence and a region(PR/RT/IN).  Each
region must be processed seperately.  The Amino acid string must be list of lists, like:

[['P'], ['I'], ['S'], ['P','K'], ['I'], ['E'], ['T'], ['V','T'], ['P'], ['V'], ['K'], ['L'], ['K'], ['P'], ...]

All amino acids MUST be uppercase.  Deletions are represented by a lowercase ['d'] and
insertions must be represented by an lowercase ['i'].

The algorithm returns a result class, which provides the following:

result.alg_name     #name of the algorithm
result.alg_version  #version of the algorithm
result.mutation_comments #list of comments
result.drugs  #a list of drugs

drug.name #drug name
drug.code  #drug code (like 3TC)
drug.drug_class #drug class (PI, NRTI, NNRTI, INT)
drug.score  #score the algorithm assigned
drug.level   #resistance level the algorithm assigned
drug.comments  #list of comments associated with this drug
"""

import xml.dom.minidom as minidom  # Bah Weep Granah Weep Mini Dom?
import re

# utility code ---------------------------------------------------------------
# random useful junk that should probably be in a util file.  Mostly translated
# over from my ruby code.

# nuc to amino hash
AA_HASH = {
    'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
    'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
    'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
    'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',

    'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
    'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
    'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
    'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',

    'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
    'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
    'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
    'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',

    'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
    'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
    'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
    'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G',
    '---': '-',
}

# Hash of Mixture nuc to list of un-mixed nucs
AMBIG = {
    'A': ['A'], 'G': ['G'], 'T': ['T'], 'C': ['C'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'K': ['G', 'T'],
    'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'S': ['C', 'G'],
    'W': ['A', 'T'], 'N': ['A', 'C', 'G', 'T'], 'X': ['X'], '-': ['-']
}


# Generates all possible non-ambigious nucleotides from an ambiguous 3 nucleotide sequence
def generate(nuc):
    try:
        posa = AMBIG[nuc[0]]
        posb = AMBIG[nuc[1]]
        posc = AMBIG[nuc[2]]

        if nuc == 'XXX':
            return ['X']

        nuclist = []
        for a in posa:
            for b in posb:
                for c in posc:
                    nuclist.append(a + b + c)

        return nuclist
    except Exception as ex:
        print
        ex
        print
        nuc
        return None


# Turns nucleotide sequence into an amino acid array
def translate_complete_to_array(sequence):
    protseq = []
    for i in range(0, len(sequence) // 3):
        alist = generate(sequence[i * 3: i * 3 + 3])
        newlist = []
        for entry in alist:
            if entry == 'X':
                newlist.append('X')
            else:
                newlist.append(raw_translator(entry))

        alist = set(newlist)  # dunno if this works?

        if len(alist) > 1:
            protseq.append(list(alist))
        elif len(alist) == 1:
            protseq.append([list(alist)[0]])
        else:
            protseq.append(['*'])
    return protseq


def raw_translator(str, keepdashs=False):
    aa = ''
    for i in range(0, (len(str) // 3)):
        if keepdashs and str[i * 3: i * 3 + 3] == '---':
            aa += 'd'  # lowercase d for deletion.
            continue

        x = AA_HASH[str[i * 3: i * 3 + 3].lower()]
        if x == None:
            x = 'X'
        aa += x
    return aa


# BNF Support code: ----------------------------------------------------------
class BNFVal:
    def __init__(self, cond, truth=False, score=0):
        self.truth = truth
        self.cond = cond
        self.score = score
        self.logic = None  # will be 'AND' or 'OR' if its from a condition2


# Now the actual code:----------------------------------------------------------
class AsiAlgorithm:
    def __init__(self, filename):
        # Extra hacks to conform to sierra's extra hacks
        self.pr_std = 'PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF'
        self.rt_std = 'PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIKVRQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWET'
        self.int_std = 'FLDGIDKAQDEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTGATVRAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRNPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED'

        # Algorithm info
        self.alg_version = ''
        self.alg_name = ''

        # definitions
        self.gene_def = []  # [['PR', ['PI']], ['RT',['NNRTI','NRTI']], ...]
        self.level_def = []  # [['1', 'Susceptible', 'S'], ...]
        self.drug_class = []  # [ ['PI', ['FPV/r', 'IDV/r', ...] ], ...]
        self.global_range = []  # [ ['-INF', '9', '1'] , ...]  #first two are the range, the third one is the res level
        self.comment_def = []  # something...

        self.drugs = []  # sub objects?  Hashes?
        self.mutation_comments = []  # maybe skip for now?  We don't really use this atm.
        self.mutations = []  # only for hivdb.  This is technically a hack that isn't part of the alg

        dom = minidom.parse(filename)

        # algorithm info
        for node in dom.getElementsByTagName('ALGNAME'):
            self.alg_name = node.childNodes[0].nodeValue
        for node in dom.getElementsByTagName('ALGVERSION'):
            self.alg_version = node.childNodes[0].nodeValue

        # definitions
        defs = dom.getElementsByTagName('DEFINITIONS')[0]
        for node in defs.getElementsByTagName('GENE_DEFINITION'):
            a = node.getElementsByTagName('NAME')[0].childNodes[0].nodeValue
            b = node.getElementsByTagName('DRUGCLASSLIST')[0].childNodes[0].nodeValue.split(',')
            b = map((lambda e: e.strip()), b)
            self.gene_def.append([a.strip(), b])

        for node in defs.getElementsByTagName('LEVEL_DEFINITION'):
            a = node.getElementsByTagName('ORDER')[0].childNodes[0].nodeValue
            b = node.getElementsByTagName('ORIGINAL')[0].childNodes[0].nodeValue
            c = node.getElementsByTagName('SIR')[0].childNodes[0].nodeValue
            self.level_def.append([a.strip(), b.strip(), c.strip()])

        for node in defs.getElementsByTagName('DRUGCLASS'):
            a = node.getElementsByTagName('NAME')[0].childNodes[0].nodeValue
            b = node.getElementsByTagName('DRUGLIST')[0].childNodes[0].nodeValue.split(',')
            b = map((lambda e: e.strip()), b)
            self.drug_class.append([a.strip(), b])

        if defs.getElementsByTagName('GLOBALRANGE') != []:
            self.global_range = defs.getElementsByTagName('GLOBALRANGE')[0].childNodes[0].nodeValue
            self.global_range = self.global_range.strip(")( \n").split(',')
            self.global_range = map(
                lambda a: (list(re.match('\s*(\S+)\s*TO\s*(\S+)\s*=>\s*(\S+)\s*', a.strip("\n \t")).groups())),
                self.global_range)
        # print repr(self.global_range)
        if defs.getElementsByTagName('COMMENT_DEFINITIONS') != []:
            comment_defs = defs.getElementsByTagName('COMMENT_DEFINITIONS')[0]
            for node in comment_defs.getElementsByTagName('COMMENT_STRING'):
                a = node.attributes.item(0).value
                b = node.getElementsByTagName('TEXT')[0].childNodes[0].nodeValue
                c = node.getElementsByTagName('SORT_TAG')[0].childNodes[0].nodeValue

                # mutations look like [REGION, POS, AAS, TYPE]
                # tmp = re.match('^([^_]+)_POS(\d+)([a-z]+)_(.+)$', a, flags=re.I)
                # if tmp:
                #  self.mutations.append([tmp.group(1), tmp.group(2), tmp.group(3), tmp.group(4)])

                self.comment_def.append([a, b, c])

        # done with the definitions, on to more interesting pastures...
        for drug_node in dom.getElementsByTagName('DRUG'):
            name = drug_node.getElementsByTagName('NAME')[0].childNodes[0].nodeValue
            fullname = ''
            if drug_node.getElementsByTagName('FULLNAME'):
                fullname = drug_node.getElementsByTagName('FULLNAME')[0].childNodes[0].nodeValue
            rules = []
            for rule_node in drug_node.getElementsByTagName('RULE'):
                # Yeah, rules!
                condition = rule_node.getElementsByTagName('CONDITION')[0].childNodes[0].nodeValue
                condition = re.sub('\s+', ' ', condition)
                actions = []
                # Need to turn condition into science?

                for action_node in rule_node.getElementsByTagName('ACTIONS'):
                    if action_node.getElementsByTagName('LEVEL') != []:
                        level = action_node.getElementsByTagName('LEVEL')[0].childNodes[0].nodeValue  # easyish
                        actions.append(['level', int(level)])

                    # Double check
                    if action_node.getElementsByTagName('COMMENT') != []:
                        comment = action_node.getElementsByTagName('COMMENT')[0]
                        actions.append(['comment', comment.attributes.item(0).value])

                    if action_node.getElementsByTagName('SCORERANGE') != [] and action_node.getElementsByTagName(
                            'USE_GLOBALRANGE') != []:
                        actions.append(['scorerange', 'useglobalrange'])
                    elif action_node.getElementsByTagName('SCORERANGE') != []:
                        srange = action_node.getElementsByTagName('SCORERANGE')[0].childNodes[0].nodeValue
                        srange = srange.strip(")( \n").split(',')
                        srange = map(lambda a: (
                        list(re.match('\s*(\S+)\s*TO\s*(\S+)\s*=>\s*(\S+)\s*', a.strip("\n \t")).groups())), srange)
                        actions.append(['scorerange', srange])

                rules.append([condition, actions])  # hrmmm
            self.drugs.append([name, fullname, rules])

        # and now we do the comments!
        if dom.getElementsByTagName('MUTATION_COMMENTS') != []:
            for gene_node in dom.getElementsByTagName('MUTATION_COMMENTS')[0].getElementsByTagName('GENE'):
                gene_name = gene_node.getElementsByTagName('NAME')[0].childNodes[0].nodeValue
                # rules and such again
                rules = []
                for rule_node in gene_node.getElementsByTagName('RULE'):
                    # Yeah, rules!
                    condition = rule_node.getElementsByTagName('CONDITION')[0].childNodes[0].nodeValue
                    condition = re.sub('\s+', ' ', condition)
                    actions = []
                    # Need to turn condition into science?

                    for action_node in rule_node.getElementsByTagName('ACTIONS'):
                        if action_node.getElementsByTagName('LEVEL') != []:
                            level = action_node.getElementsByTagName('LEVEL')[0].childNodes[0].nodeValue  # easyish
                            actions.append(['level', int(level)])

                        if action_node.getElementsByTagName('COMMENT') != []:
                            comment = action_node.getElementsByTagName('COMMENT')[0]
                            actions.append(['comment', comment.attributes.item(0).value])

                        if action_node.getElementsByTagName('SCORERANGE') != [] and action_node.getElementsByTagName(
                                'USE_GLOBALRANGE') != []:
                            actions.append(['scorerange', 'useglobalrange'])
                        elif action_node.getElementsByTagName('SCORERANGE') != []:
                            srange = action_node.getElementsByTagName('SCORERANGE')[0].childNodes[0].nodeValue
                            srange = srange.strip(")( \n").split(',')
                            srange = map(lambda a: (
                            list(re.match('\s*(\S+)\s*TO\s*(\S+)\s*:\s*(\S+)\s*', a.strip("\n \t")).groups())), srange)
                            actions.append(['scorerange', srange])
                    rules.append([condition, actions])  # hrmm
                self.mutation_comments.append([gene_name, rules])

    # This is going to be harder than I thought.  Darn you BNF!
    # starts the crazy BNF Parsing
    def interp_condition(self, cond, aaseq):
        bnf = self.bnf_statement(cond + '|', aaseq)
        if bnf.cond == False:
            print
            "---Could not parse algorithm condition:  " + cond
            exit()
            return None
        elif bnf.truth:
            # print "RESULT IS TRUE, SCORE: " + str(bnf.score)
            return [True, bnf.score]
        else:
            # print "RESULT IS FALSE, SCORE: " + str(bnf.score)
            return [False, bnf.score]

    # Hmm, what is this?  Oh not much, just a Backus-Naur Form parser.
    # booleancondition | scorecondition
    def bnf_statement(self, cond, aaseq):
        for func in [self.bnf_booleancondition, self.bnf_scorecondition]:
            bnf = func(cond, aaseq)
            if bnf.cond:
                return bnf
        return BNFVal(False)

    # condition condition2*;
    def bnf_booleancondition(self, cond, aaseq):
        bnflist = []
        bnf = self.bnf_condition(cond, aaseq)
        if bnf.cond:
            bnflist.append(bnf)
            while True:
                bnf = self.bnf_condition2(bnf.cond, aaseq)
                if bnf.cond:
                    bnflist.append(bnf)
                else:
                    break
            # Use the logic
            left_truth = None
            for bnf in bnflist:
                # print "boolean_truth " + str(bnf.truth)
                if left_truth == None:
                    left_truth = bnf.truth
                else:  # Not quite proper logic order, but I don't think anybody is randomly mixing OR's and AND's so it should be okay
                    if bnf.logic == 'AND' and left_truth and bnf.truth:
                        left_truth = True
                    elif bnf.logic == 'AND':
                        left_truth = False
                    elif bnf.logic == 'OR' and (left_truth or bnf.truth):
                        left_truth = True
                    elif bnf.logic == 'OR':
                        left_truth = False
            return BNFVal(bnflist[-1].cond, left_truth)
        else:
            return BNFVal(False)

    # l_par booleancondition r_par | residue | excludestatement | selectstatement
    def bnf_condition(self, cond, aaseq):
        for func in [self.bnf_booleancondition, self.bnf_residue, self.bnf_excludestatement, self.bnf_selectstatement]:
            if func == self.bnf_booleancondition:
                lpi = -1
                rpi = -1
                cnt = 0
                # immediaate break if it doesn't match this regexp:
                if not re.match('^\s*\(', cond):
                    continue

                # Search for parens
                for i in range(0, len(cond)):
                    if cond[i] == '(':
                        if lpi == -1:
                            lpi = i
                        cnt += 1
                    elif cond[i] == ')':
                        cnt -= 1
                        if cnt == 0:
                            rpi = i
                            break
                            # elif cond[i] == ' ':
                    # continue
                    else:
                        continue
                        # break #no good

                if lpi == -1 or rpi == -1:
                    continue

                bnf = func(cond[lpi + 1: rpi] + ' |', aaseq)
                if bnf.cond:
                    return bnf
            else:
                bnf = func(cond, aaseq)
                if bnf.cond:
                    return bnf
        return BNFVal(False)

    # logicsymbol condition;
    def bnf_condition2(self, cond, aaseq):
        bnf_logic = self.bnf_logicsymbol(cond, aaseq)
        if bnf_logic.cond:
            bnf = self.bnf_condition(bnf_logic.cond, aaseq)
            if bnf.cond:
                bnf.logic = bnf_logic.logic
                return bnf
        return BNFVal(False)

    # and | or
    def bnf_logicsymbol(self, cond, aaseq):
        if re.search('^\s*AND\s*', cond, flags=re.I):
            bnf = BNFVal(re.sub('^\s*AND\s*', '', cond))
            bnf.logic = 'AND'
            return bnf
        elif re.search('^\s*OR\s*', cond, flags=re.I):
            bnf = BNFVal(re.sub('^\s*OR\s*', '', cond))
            bnf.logic = 'OR'  # I think this works????
            return bnf
        return BNFVal(False)

    # [originalaminoacid]:amino_acid? integer [mutatedaminoacid]:amino_acid+ |
    # not [originalaminoacid]:amino_acid? Integer [mutatedaminoacid]:amino_acid+ |
    # [originalaminoacid]:amino_acid? integer l_par not [mutatedaminoacid]:amino_acid+ r_par
    def bnf_residue(self, cond, aaseq):  # this looks hard yo
        # I think we'll have to go the regexp route here.  Haha.
        mo_a = re.search('^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)', cond)
        mo_b = re.search('^\s*NOT\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)', cond)
        mo_c = re.search('^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*\(\s*NOT\s*([ARNDCEQGHILKMFPSTWYVid]+)\s*\)', cond)
        truth = False
        if mo_a:
            loc = int(mo_a.group(1))
            aas = mo_a.group(2)
            if (len(aaseq) > loc):
                for aa in aaseq[loc - 1]:
                    if aa in aas:
                        truth = True
            # print '--test' + re.sub('^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)', '', cond)
            bnf = BNFVal(re.sub('^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)', '', cond), truth)
            return bnf
        elif mo_b:
            loc = int(mo_b.group(1))
            aas = mo_b.group(2)
            truth = True
            if (len(aaseq) > loc):
                for aa in aaseq[loc - 1]:
                    if aa in aas:
                        truth = False
            else:
                truth = False  # ????UNKNOWN
            bnf = BNFVal(re.sub('^\s*NOT\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)', '', cond),
                         truth)
            return bnf
        elif mo_c:
            loc = int(mo_c.group(1))
            aas = mo_c.group(2)
            truth = True
            if (len(aaseq) > loc):
                for aa in aaseq[loc - 1]:
                    if aa in aas:
                        truth = False  # I think this makes sense
                    elif aa != '*':  # TODO, I'm not 100% sure this is correct.  (Seems to be right)
                        truth = True  # TODO
                        break  # TODO
            else:
                truth = False  # ????UNKNOWN
            bnf = BNFVal(
                re.sub('^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*\(\s*NOT\s*([ARNDCEQGHILKMFPSTWYVid]+)\s*\)', '', cond),
                truth)
            return bnf
        return BNFVal(False)

    # exclude residue
    def bnf_excludestatement(self, cond, aaseq):
        if re.search('^\s*EXCLUDE\s*', cond, flags=re.I):
            bnf = self.bnf_residue(re.sub('^\s*EXCLUDE\s*', '', cond), aaseq)
            if bnf.cond:
                bnf.truth = not bnf.truth
                return bnf
        return BNFVal(False)

    # select selectstatement2
    def bnf_selectstatement(self, cond, aaseq):
        if re.search('^\s*SELECT\s*', cond, flags=re.I):
            bnf = self.bnf_selectstatement2(re.sub('^\s*SELECT\s*', '', cond), aaseq)
            if bnf.cond:
                return bnf
        return BNFVal(False)

    # exactly integer from l_par selectlist r_par |
    # atleast integer from l_par selectlist r_par |
    # notmorethan integer from l_par selectlist r_par |
    # atleast [atleastnumber]:integer logicsymbol notmorethan [notmorethannumber]:integer from l_par selectlist r_par
    def bnf_selectstatement2(self, cond, aaseq):
        lparen = -1
        rparen = -1
        cnt = 0
        for i in range(0, len(cond)):
            if cond[i] == '(' and cnt == 0:
                lparen = i
                cnt += 1
            elif cond[i] == '(':
                cnt += 1
            elif cond[i] == ')' and cnt == 1:
                rparen = i
                cnt -= 1
                break
            elif cond[i] == ')':
                cnt -= 1

        if (lparen == -1 or rparen == -1):
            return BNFVal(False)

        # get the list items
        bnflist = self.bnf_selectlist(cond[lparen + 1:rparen] + ' |', aaseq)

        mo_a = re.search(
            '^\s*(EXACTLY\s*(\d+)|ATLEAST\s*(\d+)\s*(AND|OR)\s*NOTMORETHAN\s*(\d+)|ATLEAST\s*(\d+)|NOTMORETHAN\s*(\d+))\s*from\s*\(\s*(.+)\s*\)',
            cond, flags=re.I)
        atleastn = -1
        atmostn = -1
        exactlyn = -1
        cnt = 0
        logic = None
        type = mo_a.group(1)
        for bnf in bnflist:
            if bnf.cond and bnf.truth:
                cnt += 1
        if re.search('^\s*ATLEAST\s*(\d+)\s*(AND|OR)\s*NOTMORETHAN', type, flags=re.I):
            atleastn = int(mo_a.group(3))
            logic = mo_a.group(4)
            atmostn = int(mo_a.group(5))
            if cnt >= atleastn and cnt <= atmostn:
                return BNFVal(cond[rparen + 1:], truth=True)
            else:
                return BNFVal(cond[rparen + 1:], truth=False)
        elif re.search('^\s*EXACTLY', type, flags=re.I):
            # print 'exactly'
            exactlyn = int(mo_a.group(2))
            if cnt == exactlyn:
                return BNFVal(cond[rparen + 1:], truth=True)
            else:
                return BNFVal(cond[rparen + 1:], truth=False)
        elif re.search('^\s*ATLEAST', type, flags=re.I):
            # print 'atleastn'
            atleastn = int(mo_a.group(6))
            if cnt >= atleastn:
                return BNFVal(cond[rparen + 1:], truth=True)
            else:
                return BNFVal(cond[rparen + 1:], truth=False)
        elif re.search('^\s*NOTMORETHAN', type, flags=re.I):
            # print 'notmorethan'
            atmostn = int(mo_a.group(7))
            if cnt <= atmostn:
                return BNFVal(cond[rparen + 1:], truth=True)
            else:
                return BNFVal(cond[rparen + 1:], truth=False)

        return BNFVal(False)

    # residue listitems*
    def bnf_selectlist(self, cond, aaseq):
        bnflist = []
        bnf = self.bnf_residue(cond, aaseq)
        if bnf.cond:
            newcond = bnf.cond
            bnflist.append(bnf)
            while True:
                tmp = re.search('^\s*,\s*', newcond)
                newcond = re.sub('^\s*,\s*', '', newcond)
                bnf = self.bnf_residue(newcond, aaseq)
                if bnf.cond:
                    bnflist.append(bnf)
                    newcond = bnf.cond
                else:
                    break
            return bnflist

        return [BNFVal(False)]

    # score from l_par scorelist r_par
    def bnf_scorecondition(self, cond, aaseq):
        tmp = re.search('^\s*score\s*from\s*\((.+)\)\s*|', cond, flags=re.I)
        if tmp and tmp.group(1):
            score = 0.0
            bnf_list = self.bnf_scorelist(tmp.group(1) + '|', aaseq)
            if bnf_list[0].cond:
                newcond = cond
                for bnf in bnf_list:
                    if bnf.cond:
                        newcond = bnf.cond
                        score += bnf.score
                    else:
                        break
            else:
                return BNFVal(False)
            # I guess evalate the truth values?
            return BNFVal(newcond, score=score)
            # return bnf
        return BNFVal(False)

    # Should actually have a comma before each scoreitem*
    # scoreitem scoreitems*
    def bnf_scorelist(self, cond, aaseq):
        bnflist = []
        bnf = self.bnf_scoreitem(cond, aaseq)
        if bnf.cond:
            bnflist.append(bnf)
            newcond = bnf.cond
            while (newcond):
                # check for comma?
                # print "test:  " + newcond
                tmp = re.search('^\s*,\s*', newcond)
                if tmp:
                    newcond = re.sub('^\s*,\s*', '', newcond)
                    bnf = self.bnf_scoreitem(newcond, aaseq)
                    if bnf.cond == False:
                        break
                    newcond = bnf.cond
                    bnflist.append(bnf)
                else:
                    break
            return bnflist
        return [BNFVal(False)]

    # booleancondition mapper min? number |
    # max l_par scorelist r_par
    def bnf_scoreitem(self, cond, aaseq):
        # mo_a has 4 groups, the booleanconditon, an optional pointless MIN, and the score, and then the rest of the string
        # I think we need to match a dash for negative numbers yo
        mo_a = re.search('^\s*([^=>]+)\s*=>\s*(min)?\s*(-?\d+\.?\d*)\s*(.+)$', cond, flags=re.I)

        # Trickier than we think, as the subexpressions could have parens.  Need to do the counting game.  Sadly.
        mo_b = re.search('^\s*MAX\s*\(', cond, flags=re.I)
        if mo_b:
            # print "to mob, or not to mob?"
            lparen = -1
            rparen = -1
            cnt = 0
            for i in range(0, len(cond)):
                if cond[i] == '(' and cnt == 0:
                    lparen = i
                    cnt += 1
                elif cond[i] == '(':
                    cnt += 1
                elif cond[i] == ')' and cnt == 1:
                    rparen = i
                    cnt -= 1
                    break
                elif cond[i] == ')':
                    cnt -= 1

            if (lparen == -1 or rparen == -1):
                return BNFVal(False)
            newcond = cond[lparen + 1: rparen]
            bnflist = self.bnf_scorelist(newcond, aaseq)
            score = -999  # close enough to infinity.

            if bnflist[0].cond:
                for bnf in bnflist:
                    if bnf.cond:
                        newcond = bnf.cond
                        if bnf.score > score and bnf.score != 0.0:
                            score = bnf.score
                if (score == -999):
                    score = 0.0
                return BNFVal(cond[rparen + 1:], score=score)
        elif mo_a:
            # print "to mo a, or not to mo a? "
            # print mo_a.group(1)
            # print mo_a.group(3)
            bnf = self.bnf_booleancondition(mo_a.group(1), aaseq)
            bnf_score = float(mo_a.group(3))
            if bnf.cond:
                if bnf.truth:
                    bnf.score = bnf_score
                bnf.cond = mo_a.group(4)
                return bnf
        return BNFVal(False)

    # handles a couple undocumented comment filtering things.
    def comment_filter(self, comment, aaseq, region=None):
        # listMutsIn
        tmp = re.match('^.*\$listMutsIn\{([^\}]+)\}.*$', comment)
        if tmp:
            tmporig = tmp.group(1)
            tmporig = tmporig.replace('(', '\(').replace(')', '\)')
            muts = tmp.group(1).split(',')
            final = []
            for mut in muts:
                # If it matches \d+\(NOT \w+\), then we got to do something fancy
                tmpa = re.match('[a-z]?(\d+)\(NOT\s+([a-z]+)\)', mut, flags=re.I)
                tmpb = re.match('[a-z]?(\d+)([a-z]+)', mut, flags=re.I)
                match = ''
                loc = ''
                if tmpa:
                    loc = tmpa.group(1)
                    match = 'ARNDCEQGHILKMFPSTWYVid'
                    for ch in tmpa.group(2):
                        match = match.replace(ch, '')
                elif tmpb:
                    loc = tmpb.group(1)
                    match = tmpb.group(2)

                aas = aaseq[int(loc) - 1]
                subs = ''
                for aa in aas:
                    if aa in match:
                        subs += aa
                if subs != '':
                    subs = list(subs)
                    subs.sort()
                    subs = ''.join(subs)
                    if region == 'PR':
                        final.append(self.pr_std[int(loc) - 1] + loc + subs)
                    elif region == 'RT':
                        final.append(self.rt_std[int(loc) - 1] + loc + subs)
                    elif region == 'IN':
                        final.append(self.int_std[int(loc) - 1] + loc + subs)
                    else:
                        final.append(loc + subs)

            comment = re.sub('\$listMutsIn\{' + tmporig + '\}', ', '.join(final), comment)
            comment = re.sub(' \(\)', '', comment)  # get rid of empty brackets.

        # numberOfMutsIn
        tmp = re.match('^.+\$numberOfMutsIn\{([^\}]+)\}.+$', comment)
        if tmp:
            tmpmatch = tmp.group(1)
            muts = tmp.group(1).split(',')
            cnt = 0
            for mut in muts:
                tmp = re.match('^(\d+)([a-z]+)$', mut, flags=re.I)
                aas = aaseq[int(tmp.group(1)) - 1]
                for aa in aas:
                    if aa in tmp.group(2):
                        cnt += 1
                        # break

            comment = re.sub('\$numberOfMutsIn\{' + tmpmatch + '\}', str(cnt), comment)

        comment = re.sub('  ', ' ', comment)  # Make spacing more like sierra
        unicod = u'\uf0b1'
        comment = re.sub(unicod, '+/-', comment)  # Fixing crazy unicode characters

        return comment

    # Most important method.
    def interpret(self, aaseq, region):
        result = AsiResult()
        result.alg_name = self.alg_name
        result.alg_version = self.alg_version
        genes = list(filter((lambda e: e[0] == region), self.gene_def))[0][1]
        drs = filter((lambda e: e[0] in genes), self.drug_class)
        # print drs
        for drcls in drs:
            cls = drcls[0]
            for drname in drcls[1]:
                # Request the actual drug thingy
                #        if drname != 'ABC':
                #          continue
                drug = list(filter((lambda e: e[0] == drname), self.drugs))[0]
                drug_result = AsiDrugResult()
                drug_result.code = drug[0]
                drug_result.name = drug[1]
                drug_result.drug_class = cls

                for rule in drug[2]:
                    cond = rule[0]
                    actions = rule[1]
                    interp = self.interp_condition(cond, aaseq)

                    if interp:
                        score = interp[1]
                        truth = interp[0]
                        if truth == False and score == 0.0:
                            continue
                        for act in actions:
                            # print act
                            if act[0] == 'level':
                                if int(act[1]) > drug_result.level:
                                    drug_result.level = int(act[1])
                            elif act[0] == 'comment':
                                comm = act[1]
                                comm = filter(lambda e: (e[0] == act[1]), self.comment_def)[0]
                                comment = comm[1]
                                while (re.search('\$numberOfMutsIn{', comment) or re.search('\$listMutsIn{', comment)):
                                    comment = self.comment_filter(comment, aaseq, region)
                                drug_result.comments.append(comment)
                            elif act[0] == 'scorerange':
                                drug_result.score = score
                                scorerange = act[1]
                                if scorerange == 'useglobalrange':
                                    scorerange = self.global_range
                                else:
                                    pass
                                    # I guess this was already done!
                                    # print scorerange
                                    # scorerange = scorerange.strip(")( \n").split(',')
                                    # scorerange = map(lambda a: ( list(re.match('\s*(\S+)\s*TO\s*(\S+)\s*=>\s*(\S+)\s*', a.strip("\n \t")).groups()) ), scorerange)
                                # use score range to determine level

                                for rng in scorerange:
                                    if rng[0] == '-INF':
                                        rng[0] = -99999  # that is close enough to negative infinity.
                                    else:
                                        try:
                                            rng[0] = float(rng[0])
                                        except:
                                            rng[0] = float(rng[0])

                                    if rng[1] == 'INF':
                                        rng[1] = 99999  # that is close enough to infinity.
                                    else:
                                        try:
                                            rng[1] = float(rng[1])
                                        except:
                                            rng[1] = float(rng[1])

                                    if drug_result.score >= rng[0] and drug_result.score <= rng[1]:
                                        #                    print drug_result.score
                                        #                    print drug_result.level
                                        #                    print str(rng)
                                        if int(rng[2]) > drug_result.level:
                                            drug_result.level = int(rng[2])
                                        break
                                if (drug_result.level == None):
                                    raise "drug score range level parsing error"
                    elif interp == None:
                        print
                        "ERROR in condition: " + cond
                result.drugs.append(drug_result)

        result.drugs.sort(key=lambda e: e.code)
        result.drugs.sort(key=lambda e: e.drug_class, reverse=True)
        # comments
        for gene in self.mutation_comments:
            if gene[0] != region:
                continue
            for mut in gene[1]:
                cond = mut[0]
                actions = mut[1]

                interp = self.interp_condition(cond, aaseq)
                if interp and interp[0]:
                    for act in actions:
                        comm = list(filter(lambda e: (e[0] == act[1]), self.comment_def))[0]
                        comment = self.comment_filter(comm[1], aaseq, region)
                        result.mutation_comments.append(comment)
                elif interp == None:
                    print
                    "ERROR in condition: " + cond

                    # mutations
                    # ~ for mut in self.mutations:
                    # ~ if mut[0] == region:
                    # ~ aas = aaseq[int(mut[1]) - 1]
                    # ~ subs = ''
                    # ~ for aa in aas:
                    # ~ if aa in mut[2]:
                    # ~ subs += aa
                    # ~ if subs != '':
                    # ~ subs = list(subs)
                    # ~ subs.sort()
                    # ~ subs = ''.join(subs)
                    # ~ moot = ''
                    # ~ if region == 'PR':
                    # ~ moot = self.pr_std[int(mut[1]) - 1] + mut[1] + subs
                    # ~ elif region == 'RT':
                    # ~ moot = self.rt_std[int(mut[1]) - 1] + mut[1] + subs
                    # ~ elif region == 'IN':
                    # ~ moot = self.int_std[int(mut[1]) - 1] + mut[1] + subs
                    # ~ else:
                    # ~ moot = loc + subs

                    # ~ if mut[3] == 'Other':
                    # ~ result.mutations_other.append(moot)
                    # ~ else:
                    # ~ result.mutations.append(moot)

        return result


class AsiDrugResult:
    def __init__(self):
        self.name = ''
        self.code = ''
        self.drug_class = ''
        self.score = 0.0
        self.level = 1
        self.comments = []  # Do we need this?


class AsiResult:
    # tempting to base off the sierra result...
    def __init__(self):
        self.alg_name = ''
        self.alg_version = ''

        self.mutation_comments = []
        self.mutations = []  # mut hash yo!  Don't forget yer proteins!  (Only used for HIVDB, and its a bit of a hack)
        self.mutations_other = []  # mut hash yo!  Don't forget yer proteins!  (Only used for HIVDB, and its a bit of a hack)
        self.drugs = []  # drug hash yo!  Don't forget your scores!


"""
#Testing
algorithm_hivdb = AsiAlgorithm("HIVDB_6.1.1F.xml")
#algorithm_anrs = AsiAlgorithm("ANRS_max.xml")
#algorithm_rega = AsiAlgorithm("RegaInst_max.xml")
#print algorithm_hivdb.drugs[0][2]
#print algorithm_hivdb.mutation_comments[0]


int_seq = 'TTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCATGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCGGTAAAAACAATACATACAGACAATGGCAGCAATTACACCAGTGCTACGGTTAAGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAATAGAATCTATGAATAAAGAATTAAAGAAAATTATAAGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTRACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGAT'
pr_seq = 'CCTCAAATCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGRGGGCAGCTAAMGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGACATGGARTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGRTACCCATAGAAATTTGTGGACAYAAAACTATAGGTWCAGTATTAATAGGACCTACACCWGTTAACATAATTGGAAGAAATCTGATGAYTCAGCTTGGTTGCACTTTAAATTTT'
rt_seq = 'CCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAGGTYAARCAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGAAGGRAAGATTTCAAAAATTGGACCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAARAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCYGCAGGGTTAAAAAAGAAMAAGTCAGTAACAGTACTRGATGTGGGTGATGCATATTTTTCAGTTCCCTTATATGAAGACTTCAGGAAGTATACTGCATTCACCATACCTAGYACAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTGCCACAAGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGATAAAAATCTTAGAGCCTTTCAGAAAACAAAATCCAGARATAGTCATCTATCAATACGTGGATGATTTGTATGTAGSATCTGACTTAGAAATAGGGCAGCATAGAACAAAGATAGAGGAACTGAGAGCACATCTRTTRAAGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAGCCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACR'
aa_seq = translate_complete_to_array(rt_seq)
print aa_seq[0:15]

#hrm.

res = algorithm_hivdb.interpret(aa_seq, 'RT')
print res
for drug in res.drugs:
  print drug.code + ":  SCORE: " + str(drug.score) + " , LEVEL: " + str(drug.level)
  for com in drug.comments:
    print com
print "Mutation Comments:"
for com in res.mutation_comments:
  print com


"""
