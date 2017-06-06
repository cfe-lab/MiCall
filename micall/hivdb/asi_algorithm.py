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

from collections import defaultdict
import re
import xml.dom.minidom as minidom

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
    'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
}

# Hash of Mixture nuc to list of un-mixed nucs
AMBIG = {
    'A': ['A'], 'G': ['G'], 'T': ['T'], 'C': ['C'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'K': ['G', 'T'],
    'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'S': ['C', 'G'],
    'W': ['A', 'T'], 'N': ['A', 'C', 'G', 'T']
}


# Generates all possible non-ambigious nucleotides from an ambiguous 3 nucleotide sequence
def generate(nuc):
    posa = AMBIG[nuc[0]]
    posb = AMBIG[nuc[1]]
    posc = AMBIG[nuc[2]]

    nuclist = []
    for a in posa:
        for b in posb:
            for c in posc:
                nuclist.append(a + b + c)

    return nuclist


# Turns nucleotide sequence into an amino acid array
def translate_complete_to_array(sequence):
    protseq = []
    for i in range(0, len(sequence) // 3):
        expanded_codon = generate(sequence[i * 3: i * 3 + 3])
        aminos = {raw_translator(codon) for codon in expanded_codon}

        protseq.append(list(aminos))
    return protseq


def raw_translator(nucs):
    aa = ''
    for i in range(0, (len(nucs) // 3)):
        x = AA_HASH[nucs[i * 3: i * 3 + 3].lower()]
        aa += x
    return aa


# BNF Support code: ----------------------------------------------------------
class BNFVal:
    def __init__(self, cond, truth=False, score=0, mutations=None):
        self.truth = truth
        self.cond = cond
        self.score = score
        self.logic = None  # will be 'AND' or 'OR' if its from a condition2
        if mutations is None:
            self.mutations = set()  # mutations that triggered this, like 'M41L'
        else:
            self.mutations = mutations

    def __repr__(self):
        if not self.cond or len(self.cond) <= 21:
            cond = self.cond
        else:
            cond = self.cond[:9] + '...' + self.cond[-9:]
        return 'BNFVal({!r}, {!r}, {!r}, {!r})'.format(cond,
                                                       self.truth,
                                                       self.score,
                                                       self.mutations)


# Now the actual code:----------------------------------------------------------
class AsiAlgorithm:
    def __init__(self, file):
        """ Load ASI rules from a file or file object. """

        # Extra hacks to conform to sierra's extra hacks
        self.stds = {
            'PR': 'PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQY'
                  'DQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF',
            'RT': 'PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTP'
                  'VFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSV'
                  'PLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPD'
                  'IVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHP'
                  'DKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIKVRQLCKLLRGTKALTEVIPL'
                  'TEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKY'
                  'ARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWET',
            'IN': 'FLDGIDKAQDEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPG'
                  'IWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNG'
                  'SNFTGATVRAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQ'
                  'MAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRNPLWK'
                  'GPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED'}

        # Algorithm info
        self.alg_version = ''
        self.alg_name = ''

        # definitions
        self.gene_def = {}  # {code: [drug_class_code]}
        self.level_def = {}  # {'1': 'Susceptible'}
        self.drug_class = {}  # {code: [drug_code]}
        self.global_range = []  # [ ['-INF', '9', '1'] , ...]  #first two are the range, the third one is the res level
        self.comment_def = {}  # {code: comment_text}

        self.drugs = {}  # {code: (name, [condition, [(action_type, action_value)]])}
        self.mutation_comments = []  # maybe skip for now?  We don't really use this atm.
        self.mutations = []  # only for hivdb.  This is technically a hack that isn't part of the alg

        dom = minidom.parse(file)

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
            self.gene_def[a.strip()] = [e.strip() for e in b]

        for node in defs.getElementsByTagName('LEVEL_DEFINITION'):
            a = node.getElementsByTagName('ORDER')[0].childNodes[0].nodeValue
            b = node.getElementsByTagName('ORIGINAL')[0].childNodes[0].nodeValue
            self.level_def[a.strip()] = b.strip()

        for node in defs.getElementsByTagName('DRUGCLASS'):
            a = node.getElementsByTagName('NAME')[0].childNodes[0].nodeValue
            b = node.getElementsByTagName('DRUGLIST')[0].childNodes[0].nodeValue.split(',')
            self.drug_class[a.strip()] = [e.strip() for e in b]

        if defs.getElementsByTagName('GLOBALRANGE'):
            global_range_text = defs.getElementsByTagName('GLOBALRANGE')[0].childNodes[0].nodeValue
            global_range_items = global_range_text.strip(")( \n").split(',')
            self.global_range = list(map(
                lambda item: (list(re.match('\s*(\S+)\s*TO\s*(\S+)\s*=>\s*(\S+)\s*',
                                            item.strip("\n \t")).groups())),
                global_range_items))
        if defs.getElementsByTagName('COMMENT_DEFINITIONS'):
            comment_defs = defs.getElementsByTagName('COMMENT_DEFINITIONS')[0]
            for node in comment_defs.getElementsByTagName('COMMENT_STRING'):
                comment_id = node.attributes.item(0).value
                b = node.getElementsByTagName('TEXT')[0].childNodes[0].nodeValue
                c = node.getElementsByTagName('SORT_TAG')[0].childNodes[0].nodeValue

                # mutations look like [REGION, POS, AAS, TYPE]
                # tmp = re.match('^([^_]+)_POS(\d+)([a-z]+)_(.+)$', a, flags=re.I)
                # if tmp:
                #  self.mutations.append([tmp.group(1), tmp.group(2), tmp.group(3), tmp.group(4)])

                self.comment_def[comment_id] = (b, c)

        # done with the definitions, on to more interesting pastures...
        for drug_node in dom.getElementsByTagName('DRUG'):
            name = drug_node.getElementsByTagName('NAME')[0].childNodes[0].nodeValue
            fullname = ''
            if drug_node.getElementsByTagName('FULLNAME'):
                fullname = drug_node.getElementsByTagName('FULLNAME')[0].childNodes[0].nodeValue
            rules = []
            for rule_node in drug_node.getElementsByTagName('RULE'):
                condition = rule_node.getElementsByTagName('CONDITION')[0].childNodes[0].nodeValue
                condition = re.sub('\s+', ' ', condition)
                actions = []

                for action_node in rule_node.getElementsByTagName('ACTIONS'):
                    if action_node.getElementsByTagName('LEVEL'):
                        level = action_node.getElementsByTagName('LEVEL')[0].childNodes[0].nodeValue  # easyish
                        actions.append(('level', level))

                    if action_node.getElementsByTagName('COMMENT'):
                        comment = action_node.getElementsByTagName('COMMENT')[0]
                        actions.append(('comment', comment.attributes.item(0).value))

                    if (action_node.getElementsByTagName('SCORERANGE') and
                            action_node.getElementsByTagName('USE_GLOBALRANGE')):
                        actions.append(('scorerange', 'useglobalrange'))
                    elif action_node.getElementsByTagName('SCORERANGE'):
                        srange = action_node.getElementsByTagName('SCORERANGE')[0].childNodes[0].nodeValue
                        srange = srange.strip(")( \n").split(',')
                        srange = [re.match('\s*(\S+)\s*TO\s*(\S+)\s*=>\s*(\S+)\s*',
                                           item.strip("\n \t")).groups()
                                  for item in srange]
                        actions.append(('scorerange', srange))

                rules.append((condition, actions))  # hrmmm
            self.drugs[name] = (fullname, rules)

        # and now we do the comments!
        if dom.getElementsByTagName('MUTATION_COMMENTS'):
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
                        if action_node.getElementsByTagName('COMMENT'):
                            comment = action_node.getElementsByTagName('COMMENT')[0]
                            actions.append(['comment', comment.attributes.item(0).value])
                    rules.append([condition, actions])  # hrmm
                self.mutation_comments.append([gene_name, rules])

    # This is going to be harder than I thought.  Darn you BNF!
    # starts the crazy BNF Parsing
    def interp_condition(self, cond, aaseq):
        bnf = self.bnf_statement(cond + '|', aaseq)
        assert bnf.cond
        if bnf.truth:
            return [True, bnf.score, bnf.mutations]
        return [False, bnf.score, bnf.mutations]

    # Hmm, what is this?  Oh not much, just a Backus-Naur Form parser.
    # booleancondition | scorecondition
    def bnf_statement(self, cond, aaseq):
        for func in [self.bnf_booleancondition, self.bnf_scorecondition]:
            bnf = func(cond, aaseq)
            if bnf.cond:
                return bnf

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
            mutations = None
            for bnf in bnflist:
                # print "boolean_truth " + str(bnf.truth)
                if left_truth is None:
                    left_truth = bnf.truth
                    mutations = bnf.mutations
                else:
                    # Not quite proper logic order, but I don't think anybody
                    # is randomly mixing OR's and AND's so it should be okay
                    if bnf.logic == 'AND':
                        left_truth = left_truth and bnf.truth
                        if left_truth:
                            mutations |= bnf.mutations
                        else:
                            mutations.clear()
                    elif bnf.logic == 'OR':
                        left_truth = left_truth or bnf.truth
                        mutations |= bnf.mutations
            return BNFVal(bnflist[-1].cond, left_truth, mutations=mutations)
        else:
            return BNFVal(False)

    # l_par booleancondition r_par | residue | excludestatement | selectstatement
    def bnf_condition(self, cond, aaseq):
        for func in [self.bnf_booleancondition, self.bnf_residue, self.bnf_excludestatement, self.bnf_selectstatement]:
            if func == self.bnf_booleancondition:
                lpi = -1
                rpi = -1
                cnt = 0
                # immediate break if it doesn't match this regexp:
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

                assert lpi >= 0 and rpi >= 0, (lpi, rpi)

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
        bnf_logic = self.bnf_logicsymbol(cond)
        if bnf_logic.cond:
            bnf = self.bnf_condition(bnf_logic.cond, aaseq)
            if bnf.cond:
                bnf.logic = bnf_logic.logic
                return bnf
        return BNFVal(False)

    # and | or
    @staticmethod
    def bnf_logicsymbol(cond):
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
    @staticmethod
    def bnf_residue(cond, aaseq):  # this looks hard yo
        # I think we'll have to go the regexp route here.  Haha.
        mo_a = re.search(
            '^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)(.*)',
            cond)
        mo_b = re.search(
            '^\s*(?:NOT|EXCLUDE)\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYVid]+)(.*)',
            cond)
        mo_c = re.search(
            '^\s*[ARNDCEQGHILKMFPSTWYVid]?\s*(\d+)\s*\(\s*NOT\s*([ARNDCEQGHILKMFPSTWYVid]+)\s*\)(.*)',
            cond)
        truth = False
        mutations = None
        if mo_a:
            loc = int(mo_a.group(1))
            aas = mo_a.group(2)
            if len(aaseq) >= loc:
                mutations = {str(loc) + aa
                             for aa in aaseq[loc - 1]
                             if aa in aas}
                truth = bool(mutations)
            bnf = BNFVal(mo_a.group(3), truth, mutations=mutations)
            return bnf
        elif mo_b or mo_c:
            match = mo_b or mo_c
            loc = int(match.group(1))
            aas = match.group(2)
            if len(aaseq) >= loc:
                truth = not any(aa in aas for aa in aaseq[loc - 1])
            else:
                truth = False  # ????UNKNOWN
            bnf = BNFVal(match.group(3), truth)
            return bnf
        return BNFVal(False)

    # exclude residue
    def bnf_excludestatement(self, cond, aaseq):
        return self.bnf_residue(cond, aaseq)

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
        """ Parse the details of a select statement.

        :param str cond: the details of the select statement to parse.
        :param aaseq: the list of amino acid lists for each position
        """
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

        if lparen == -1 or rparen == -1:
            return BNFVal(False)

        # get the list items
        bnflist = self.bnf_selectlist(cond[lparen + 1:rparen] + ' |', aaseq)

        mo_a = re.search(
            '^\s*(EXACTLY\s*(\d+)|ATLEAST\s*(\d+)\s*(AND|OR)\s*'
            'NOTMORETHAN\s*(\d+)|ATLEAST\s*(\d+)|NOTMORETHAN\s*(\d+))\s*'
            'from\s*\(\s*(.+)\s*\)',
            cond, flags=re.I)
        cnt = 0
        select_type = mo_a.group(1)
        result = BNFVal(cond[rparen + 1:])
        for bnf in bnflist:
            if bnf.cond and bnf.truth:
                cnt += 1
                result.mutations |= bnf.mutations

        if re.search('^\s*ATLEAST\s*(\d+)\s*(AND|OR)\s*NOTMORETHAN',
                     select_type,
                     flags=re.I):
            atleastn = int(mo_a.group(3))
            logic = mo_a.group(4)
            atmostn = int(mo_a.group(5))
            if logic.upper() == 'AND':
                result.truth = atleastn <= cnt <= atmostn
            else:
                result.truth = atleastn <= cnt or cnt <= atmostn
        elif re.search('^\s*EXACTLY', select_type, flags=re.I):
            # print 'exactly'
            exactlyn = int(mo_a.group(2))
            result.truth = cnt == exactlyn
        elif re.search('^\s*ATLEAST', select_type, flags=re.I):
            # print 'atleastn'
            atleastn = int(mo_a.group(6))
            result.truth = cnt >= atleastn
        elif re.search('^\s*NOTMORETHAN', select_type, flags=re.I):
            # print 'notmorethan'
            atmostn = int(mo_a.group(7))
            result.truth = cnt <= atmostn

        return result

    # residue listitems*
    def bnf_selectlist(self, cond, aaseq):
        bnflist = []
        bnf = self.bnf_residue(cond, aaseq)
        if bnf.cond:
            newcond = bnf.cond
            bnflist.append(bnf)
            while True:
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
            mutations = set()
            if bnf_list[0].cond:
                newcond = cond
                for bnf in bnf_list:
                    if bnf.cond:
                        newcond = bnf.cond
                        score += bnf.score
                        mutations |= bnf.mutations
                    else:
                        break
            else:
                return BNFVal(False)
            # I guess evalate the truth values?
            return BNFVal(newcond, score=score, mutations=mutations)
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
            while newcond:
                # check for comma?
                # print "test:  " + newcond
                tmp = re.search('^\s*,\s*', newcond)
                if tmp:
                    newcond = re.sub('^\s*,\s*', '', newcond)
                    bnf = self.bnf_scoreitem(newcond, aaseq)
                    if not bnf.cond:
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
        # mo_a has 4 groups, the booleanconditon, an optional pointless MIN,
        # and the score, and then the rest of the string
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

            if lparen == -1 or rparen == -1:
                return BNFVal(False)
            newcond = cond[lparen + 1: rparen]
            bnflist = self.bnf_scorelist(newcond, aaseq)
            score = -999  # close enough to infinity.

            if bnflist[0].cond:
                mutations = set()
                for bnf in bnflist:
                    if bnf.cond:
                        mutations |= bnf.mutations
                        if bnf.score > score and bnf.score != 0.0:
                            score = bnf.score
                if score == -999:
                    score = 0.0
                return BNFVal(cond[rparen + 1:],
                              score=score,
                              mutations=mutations)
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
                    std = self.stds.get(region)
                    if std is not None:
                        final.append(std[int(loc) - 1] + loc + subs)
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
        raw_mutations = defaultdict(set)
        result.alg_name = self.alg_name
        result.alg_version = self.alg_version
        drug_classes = self.gene_def[region]
        default_level = 1
        default_level_name = self.level_def['1']
        for drug_class in drug_classes:
            for drug_code in self.drug_class[drug_class]:
                drug_name, drug_rules = self.drugs[drug_code]
                drug_result = AsiDrugResult()
                drug_result.code = drug_code
                drug_result.name = drug_name
                drug_result.level = default_level
                drug_result.level_name = default_level_name
                drug_result.drug_class = drug_class

                for rule in drug_rules:
                    cond = rule[0]
                    actions = rule[1]
                    interp = self.interp_condition(cond, aaseq)

                    score = interp[1]
                    truth = interp[0]
                    raw_mutations[drug_class] |= interp[2]
                    if not truth and score == 0.0:
                        continue
                    for act in actions:
                        if act[0] == 'level':
                            if int(act[1]) > drug_result.level:
                                drug_result.level = int(act[1])
                                drug_result.level_name = self.level_def[act[1]]
                        elif act[0] == 'comment':
                            comment, _ = self.comment_def[act[1]]
                            while (re.search('\$numberOfMutsIn{', comment) or
                                   re.search('\$listMutsIn{', comment)):
                                comment = self.comment_filter(comment, aaseq, region)
                            drug_result.comments.append(comment)
                        elif act[0] == 'scorerange':
                            drug_result.score = score
                            scorerange = act[1]
                            if scorerange == 'useglobalrange':
                                scorerange = self.global_range

                            # use score range to determine level
                            for low_score, high_score, level in scorerange:
                                if low_score == '-INF':
                                    low_score = -99999  # that is close enough to negative infinity.
                                else:
                                    low_score = float(low_score)

                                if high_score == 'INF':
                                    high_score = 99999  # that is close enough to infinity.
                                else:
                                    high_score = float(high_score)

                                if low_score <= drug_result.score <= high_score:
                                    if int(level) > drug_result.level:
                                        drug_result.level = int(level)
                                        drug_result.level_name = self.level_def[level]
                                    break
                result.drugs.append(drug_result)

        for cls, cls_mutations in raw_mutations.items():
            mutation_parts = sorted((int(mutation[:-1]), mutation[-1])
                                    for mutation in cls_mutations)
            std = self.stds[region]
            result.mutations[cls] = ['{}{}{}'.format(std[pos-1], pos, amino)
                                     for pos, amino in mutation_parts]
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
                if interp[0]:
                    for act in actions:
                        comment_template, _ = self.comment_def[act[1]]
                        comment = self.comment_filter(comment_template, aaseq, region)
                        result.mutation_comments.append(comment)

        return result


class AsiDrugResult:
    def __init__(self):
        self.name = ''
        self.code = ''
        self.drug_class = ''
        self.score = 0.0
        self.level = self.level_name = None
        self.comments = []  # Do we need this?


class AsiResult:
    # tempting to base off the sierra result...
    def __init__(self):
        self.alg_name = ''
        self.alg_version = ''

        self.mutation_comments = []
        self.mutations = {}
        self.drugs = []  # drug hash yo!  Don't forget your scores!
