=begin
Iterate through a CSV (produced by remap.py) where each row contains
 the contents of a SAM file.
Sequences have already been restricted to HIV env V3 if the reads were
 mapped to the 'V3LOOP' region.
Do not attempt to clean sequence based on CIGAR string.
Apply Conan's implementation of the geno2pheno (g2p) coreceptor
 tropism prediction algorithm to each entry.  Output tuples of header,
 g2p score, FPR, and g2p-aligned protein sequence.

Call:
ruby sam_to_g2p.rb <INPUT SAM csv> <OUTPUT scored csv>

Dependencies:
BioRuby
pssm_lib.rb
g2p.matrix
CSV [Ruby module]
=end

require 'csv'
require './pssm_lib'
require 'bio'

QMIN = 20
QCUT = 10
QDELTA = 5

def apply_cigar (cigar, seq, qual)
    newseq = ''
    newqual = ''
    tokens = cigar.scan(/[0-9]+[MIDNSHPX=]/)
    if tokens.length == 0
        return Nil, Nil, Nil
    end
    shift = 0
    if tokens[0][-1,1] == 'S'
        shift = tokens[0][0..-2].to_i
    end
    left = 0
    tokens.each { |token|
        length = token[0..-2].to_i
        right = left + length - 1
        operation = token[-1,1]
        if operation == 'M'
            newseq += seq[left..right]
            newqual += qual[left..right]
            left += length
        elsif operation == 'D'
            newseq += '-'*length
            newqual += ' '*length
        elsif operation == 'I'
            if length % 3 > 0
                left += length
                next
            end
            its_quals = qual[left..right].unpack('C*') # array of ASCII codes
            if its_quals.map{|asc| asc-33 >= QMIN}.all?
                newseq += seq[left..right]
                newqual += qual[left..right]
            end
            left += length
            next
        elsif operation == 'S'
            left += length
            next
        else
            puts 'Unable to handle operation #{operation}'
        end
    }
    return shift, newseq, newqual
end

def merge_pairs (seq1, seq2, qual1, qual2)
    mseq = ''
    if seq2.length < seq1.length # make sure seq2 is longer
      seq1, seq2 = seq2, seq1 
      qual1, qual2 = qual2, qual1
    end
    i = 0
    qual1_ints = qual1.unpack('C*')
    qual2_ints = qual2.unpack('C*')
    qual1_ints.map! {|x| x - 33}
    qual2_ints.map! {|x| x - 33}
    
    seq2.each_byte { |b2|
        c2 = b2.chr
        q2 = qual2_ints[i]
        if i < seq1.length
            c1 = seq1[i,1]
            q1 = qual1_ints[i]
            if c1 == '-' and c2 == '-'
                mseq += '-'
                next
            end
            if c1 == c2
                if q1 > QCUT or q2 > QCUT
                    mseq += c1
                else
                    mseq += 'N'
                end
            else
                if (q2 - q1).abs >= QDELTA
                    if q1 > [q2, QCUT].max
                        mseq += c1
                    elsif q2 > [q1, QCUT].max
                        mseq += c2
                    else
                        mseq += 'N'
                    end
                else
                    mseq += 'N'
                end
            end
        else
            # past the end of read 1
            if c2 == '-' and q2 == 0
                mseq += 'n' # interval between reads
            else
                if q2 > QCUT
                    mseq += c2
                else
                    mseq += 'N'
                end
            end
        end
        i += 1 # increment at end of loop
    }
    return mseq
end




# convert file contents into hash
pairs = Hash.new  # cache read for pairing
merged = Hash.new # tabulate merged sequence variants
sample_name = ''

CSV.foreach(ARGV[0]) do |row|
    sample_name, qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = row
    if rname != 'V3LOOP'
        # uninteresting region or the header row
        next
    end
    shift, seq1, qual1 = apply_cigar(cigar, seq, qual)
    pos1 = pos.to_i-1
    seq2 = '-'*pos1 + seq1 # pad the sequence
    qual2 = '!'*pos1 + qual1
    
    if pairs.has_key?(qname) # merge reads
        seq1 = pairs[qname]['seq'] # retrieve from cache
        qual1 = pairs[qname]['qual']
        
        mseq = merge_pairs(seq1, seq2, qual1, qual2)
        
        seqlen = (mseq.delete '-').length
        if (mseq.count 'N') > (0.5*seqlen)
            # if more than 50% of sequence is garbage
            next
        end
        
        if merged.has_key?(mseq)
            merged[mseq] += 1
        else
            merged[mseq] = 1
        end
    else
        pairs[qname] = {'seq'=>seq2, 'qual'=>qual2} # cache this read
    end
end  

sorted = merged.sort_by {|k,v| v}.reverse  # sort by count in descending order

# apply g2p algorithm to merged reads
f = File.open(ARGV[1], mode='w') # write-only
f.write("sample,rank,count,g2p,fpr,aligned,error\n")  # CSV header

rank = 0
sorted.each do |s, count|
    rank += 1
    prefix = "#{sample_name},#{rank},#{count}"
    seq = s.delete '-'
    if s.length % 3 != 0
        f.write("#{prefix},,,,notdiv3\n")
        next
    end
    
    if seq.length == 0
      f.write("#{prefix},,,,zerolength\n")
      next
    end
    
    dna = Bio::Sequence.auto(seq)
    prot = dna.translate

    # sanity check 1 - bounded by cysteines
    if !prot.match(/^C/) || !prot.match(/C$/)
        f.write("#{prefix},,,#{prot},cysteines\n")
        next
    end

    # sanity check 2 - no ambiguous codons
    if prot.count('X') > 0
        f.write("#{prefix},,,#{prot},ambiguous\n")
        next
    end

    # sanity check 3 - no stop codons
    if prot.count('*') > 0
        f.write("#{prefix},,,#{prot},stop codons\n")
        next
    end

    # sanity check 4 - V3 length in range 32-40 inclusive
    if prot.length < 32 || prot.length > 40
        f.write("#{prefix},,,#{prot},length\n")
        next
    end

    score, aligned = run_g2p(seq, $std_v3_g2p, load_matrix('g2p'))
    aligned2 = aligned.flatten * ""
    fpr = g2p_to_fpr(score)
    f.write("#{prefix},#{score},#{fpr},#{aligned2},\n")
end

f.close()

