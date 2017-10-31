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
from pyvdrm.asi2 import ASI2
from pyvdrm.vcf import MutationSet

def call_mutations(aaseq):
    '''shim to shoehorn new DRM evaluator -- replace'''
    muset = set([])
    for pos, mus in zip(range(1, len(aaseq) + 1), aaseq):
        mus = ''.join(mus)
        if mus == '':
            continue
        muset |= MutationSet.from_string("{}{}".format(pos, mus)).mutations
    return muset


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

    def get_gene_positions(self, gene):
        positions = set([])
        for drug_class in self.gene_def[gene]:
            for drug_code in self.drug_class[drug_class]:
                drug_config = self.drugs[drug_code]
                rules = drug_config[1]
                for condition, _ in rules:
                    for match in re.finditer(r'(\d+)[A-Zid]', condition):
                        positions.add(int(match.group(1)))
        return positions

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

        mutations = call_mutations(aaseq)

        for drug_class in drug_classes:
            for drug_code in self.drug_class[drug_class]:
                drug_name, drug_rules = self.drugs[drug_code]
                drug_result = AsiDrugResult()
                drug_result.code = drug_code
                drug_result.name = drug_name
                drug_result.level = default_level
                drug_result.level_name = default_level_name
                drug_result.drug_class = drug_class

                for condition, actions in drug_rules:

                    rule = ASI2(condition)
                    rule_result = rule.dtree(mutations)

                    if rule_result is None:
                        continue

                    score = float(rule_result.score)
                    m = rule_result.residues
                    raw_mutations[drug_class] |= m

                    for action, comment in actions:
                        if action == 'level':
                            if int(comment) > drug_result.level:
                                drug_result.level = int(comment)
                                drug_result.level_name = self.level_def[comment]
                        elif action == 'comment':
                            comment, _ = self.comment_def[comment]
                            while (re.search('\$numberOfMutsIn{', comment) or
                                   re.search('\$listMutsIn{', comment)):
                                comment = self.comment_filter(comment, aaseq, region)
                            drug_result.comments.append(comment)
                        elif action == 'scorerange':
                            drug_result.score = score
                            scorerange = comment
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
            for m in cls_mutations:
                m.setWt(self.stds[region])
            result.mutations[cls] = [str(m) for m in cls_mutations]

        result.drugs.sort(key=lambda e: e.code)
        result.drugs.sort(key=lambda e: e.drug_class, reverse=True)
        # comments
        for target_region, results in self.mutation_comments:
            if target_region != region:
                continue
            for cond, actions in results:

                # Why are we evaluating the mutations twice???
                # this code should be considered broken
                rule = ASI2(cond)
                scoring = rule(mutations)

                if scoring:
                    for _, act in actions:
                        comment_template, _ = self.comment_def[act]
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
