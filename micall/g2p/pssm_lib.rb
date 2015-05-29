=begin
pssm_lib.rb
Written by Conan
=end

require 'tempfile'
require '../alignment/alignment' #recall alignment alg
#require 'net/http'
#require 'uri'
#require 'ckwlib/ckw'
#require 'ckwlib/cfe_scores.rb'

$std_v3_pssm = 'CTRPNNNTRKGIHIGPGRAFYATGEIIGDIRQAHC'
$std_v3_g2p = 'CTRPNXNNTXXRKSIRIXXXGPGQXXXAFYATXXXXGDIIGDIXXRQAHC'
$std_v3_nuc = 'TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT'

#Returns FPR and score of an AA string
def run_g2p_on_aa(aa_seq)
  score = run_g2p([aa_seq.split('').map{|a| [a]}], $std_v3_g2p, load_matrix('g2p'))  
  return score[0], g2p_to_fpr(score[0])
end

#Returns a score and an alignment
def g2p(aa_seq, matrix)
    rho = 1.33153
    probA = -2.31191
    probB = 0.244784

    sum = 0.to_f
    0.upto(49) do |i|
        w_arr = []
        0.upto(aa_seq[i].size - 1) do |j|
            w_arr << matrix[i][aa_seq[i][j]]
        end
        s = 0
        w_arr.each {|v| s += v }
        s = s / w_arr.size #think that makes sense
        sum += s #It gets inverted funny at some point, so min actually will result in more x4
    end

    dv = rho - sum
    fapb = (dv * probA) + probB
    score = (1 / (1 + Math.exp(fapb-0.5)))
    return [rho,probA,probB,sum,dv,fapb,score]
end

def g2p_worstcase(aa_seq, matrix)
    rho = 1.33153
    probA = -2.31191
    probB = 0.244784

    sum = 0.to_f
    0.upto(49) do |i|
        w_arr = []
        0.upto(aa_seq[i].size - 1) do |j|
            w_arr << matrix[i][aa_seq[i][j]]
        end
        sum += w_arr.min #It gets inverted funny at some point, so min actually will result in more x4
    end

    dv = rho - sum
    fapb = (dv * probA) + probB
    score = (1 / (1 + Math.exp(fapb-0.5)))
    return [rho,probA,probB,sum,dv,fapb,score]
end

#Returns a score.  
def pssm(aa_seq, matrix=nil)
    min_seq_score = 0
    max_seq_score = 0
    0.upto(aa_seq.size - 1) do |i|
        aa_min = 0
        aa_max = 0
        if(aa_seq[i].size > 1)
            0.upto(aa_seq[i].size - 1) do |j|
                aa_score = matrix[i][aa_seq[i][j]]
                if(j == 0)
                    aa_min = aa_score 
                    aa_max = aa_score
                elsif(aa_score < aa_min)
                    aa_min = aa_score
                elsif(aa_score > aa_max)
                    aa_max = aa_score
                end
            end
            min_seq_score += aa_min
            max_seq_score += aa_max
        else
            aa_score = matrix[i][aa_seq[i][0]]
            if(aa_score)
                min_seq_score += aa_score
                max_seq_score += aa_score
            end
        end
      end
    return [min_seq_score, max_seq_score]
  end
  
  def pssm_avg(aa_seq, matrix=nil)
    final_score = 0
    0.upto(aa_seq.size - 1) do |i|
      aa_score = 0
      aa_cnt = 0
      if(aa_seq[i].size > 1)
        0.upto(aa_seq[i].size - 1) do |j|
          aa_score += matrix[i][aa_seq[i][j]]
          aa_cnt +=1
        end
      else
        aa_score = matrix[i][aa_seq[i][0]]
        aa_cnt = 1
      end
      final_score += (aa_score / aa_cnt)
    end
    return final_score
  end
  
def pssm_nuc(nuc_seq, matrix=nil)
  final_score = 0
  nuc_seq = nuc_seq.split('').map{|a| AMBIG[a]}
    
  0.upto(nuc_seq.size - 1) do |i|
    nuc_score = 0.0
    nuc_cnt = 0.0
    if(nuc_seq[i].size > 1) #expanded?
      0.upto(nuc_seq[i].size - 1) do |j|
        nuc_score += matrix[i][nuc_seq[i][j]]
        nuc_cnt += 1
      end
    else
      nuc_score = matrix[i][nuc_seq[i][0]]
      nuc_cnt = 1
    end
    final_score += (nuc_score / nuc_cnt)
  end
  return final_score 
end

#If there is an insertion (or deletion?) or a mixture, then modify.
#indel: -0.06 binary, exists are doesn't
#mix: 0.57 binary, exists or doesn't.  Also includes deletions for some reason...
=begin
def cfe_score(aa_seq, has_indel = false)
    min_seq_score = 0
    max_seq_score = 0
    has_mix = false

    0.upto(aa_seq.size - 1) do |i|
        aa_min = 0
        aa_max = 0
        has_mix = true if(aa_seq[i] == ['-'])
        if(aa_seq[i].size > 1)
                has_mix = true
                0.upto(aa_seq[i].size - 1) do |j|
                #aa_score = matrix[i][aa_seq[i][j]]
                aa_score = $cfe_score_hash[(i + 1).to_s + aa_seq[i][j]]
#                puts (i + 1).to_s + aa_seq[i][j] + "\t" + $cfe_score_hash[(i + 1).to_s + aa_seq[i][j]].to_s
                if(j == 0)
                    aa_min = aa_score 
                    aa_max = aa_score
                elsif(aa_score < aa_min)
                    aa_min = aa_score
                elsif(aa_score > aa_max)
                    aa_max = aa_score
                end
            end
            min_seq_score += aa_min
            max_seq_score += aa_max
        else
            #aa_score = matrix[i][aa_seq[i][0]]
            aa_score = $cfe_score_hash[(i + 1).to_s + aa_seq[i][0]]
#            puts (i + 1).to_s + aa_seq[i][0] + "\t" + $cfe_score_hash[(i + 1).to_s + aa_seq[i][0]].to_s
            if(aa_score)
                min_seq_score += aa_score
                max_seq_score += aa_score
            end
        end
    end
    min_seq_score += 0.57 if(has_mix)
    max_seq_score += 0.57 if(has_mix)
    min_seq_score += -0.06 if(has_indel)
    max_seq_score += -0.06 if(has_indel)
    return [min_seq_score, max_seq_score]
end
=end
#Returns an interpretation
#Matrix includes std
def interpret_g2p_score(score)
    return (score > 0.156) ? 'X4' : 'R5'
end

def interpret_g2p_score_5perc(score)
    return (score > 0.25868) ? 'X4' : 'R5'
end

def interpret_sinsi_score(score)
    return (score > -5.40) ? 'X4' : 'R5'
end

=begin
def interpret_cfe_score(score)
#    return (score > -5.40) ? 'X4' : 'R5'
    return "No idea"
end
=end

#Returns an interpretation
#Matrix includes std
def interpret_pssm_score(score)
    return (score > -6.96) ? 'X4' : 'R5'
end


#Returns all possible equivilent alignments to an existing alignment
#WARNING, this may go slowly if lots of gaps.
#1) finds all possible alignments of the first gap
#2) for all those, change the gap into a period or something
#3) Recursively call malign on each of THESE alignments to start on the next gap.
#4) Return everything
# note, may actually be easier if you align all gaps to the left first.
# Hrm, I actually suspect that we should probably do this on the nucleotide sequence first...
# Big gaps should be aligned as a whole.
def malign(std, orig, first=true)
#    puts std
#    puts orig.join
    seqs = []
    bestscore = score_seq(std, orig)
    dex = orig.index(['-'])
    
    return [Array.new(orig)] if(dex == nil)
    
    new = Array.new(orig)
    new[dex] = ['.']
    while(dex >= 1)
        new[dex] = new[dex - 1]
        new[dex - 1] = ['.']
        dex-=1
        if(score_seq(std, new) < bestscore)
            new[dex] = new[dex + 1]
            new[dex + 1] = ['.']
            dex+=1
            break
        end
    end
    seqs += malign(std, Array.new(new), false)
    
    #Ok, we are at the the first possible gap location, lets go forward and analyze each.
    while(dex < std.size - 1)
        new[dex] = new[dex + 1]
        new[dex + 1] = ['.']
        dex += 1
        if(score_seq(std, new) < bestscore)
            break
        end
        seqs += malign(std, Array.new(new),false) #Can go to recursive loop here?
    end
#    seqs.each do |s|
#        puts s.map{|a| a[0]}.join('')
#    end
#    puts '--'
    if(first)
        seqs.each do |s|
            s.each do |aa|
                aa[0] = '-' if(aa == ['.'])
            end
        end
    end
    return seqs.uniq
#    return seqs.map{|s| s.gsub('.', '-')}
end

#Only works proper like on aminos
#Should handle mixtures  
def score_seq(std, seq)
#    begin
    score = 0
    std = std.split('')
    0.upto(std.size - 1) do |i|
        #if(alignment[0][i,1] == alignment[1][i,1])
        if(seq[i].include?(std[i]) )
            score += 1
        end
    end
    return score
#    rescue
#        puts std.inspect
        #puts seq.inspect
        #exit
    #end
end


#Aligns the amino acids to a standard
#Method can be :recall, :mynap, :muscle (In reality, only :recall at the moment)
#Expects dashes, not X's
def align_aminos(std, seq, method=:recall, removeinserts=true, qachecks=true)  #Hmmm?
    if(qachecks)
      return -1 if(seq == nil)
      return -1 if(seq.size % 3 != 0 or seq.size < 99)
      return -1 if(seq =~ /^----/ or seq =~ /----$/)
    end
    indels = false
    #Hmmm!  
    
    aa_seq = nil
    if(seq.class == Array)
      aa_seq = seq
    else
      aa_seq = Pssm.translate_complete_to_array(seq)
    end
    
    aa_seq.each do |a|
        0.upto(a.size - 1) { |i| a[i] = '-' if(a[i] == 'X') }
        a.delete_if {|aa| aa == '*'} if(a.size > 1)
    end
    aa_seq.delete_if{|s| s == ['-']}
    aa_seq_s = aa_seq.map {|a| a.sort[0]}.join('')
    elem = nil

    if(qachecks)
      return -1 if(aa_seq.any? {|a| a.include?('*') }) #Has stop codons
    end
    
    if(method == :none)
        elem = method_none(std, aa_seq_s)
    elsif(method == :recall1)
        std.gsub!('-', 'X')
        elem = method_recall(std, aa_seq_s,1)
        elem[0].gsub!('X', '-')
        std.gsub!('X', '-')
    elsif(method == :recall)
        std.gsub!('-', 'X')
        elem = method_recall(std, aa_seq_s)
        elem[0].gsub!('X', '-')
        std.gsub!('X', '-')
    elsif(method == :recall6)
        std.gsub!('-', 'X')
        elem = method_recall(std, aa_seq_s, 6)
        elem[0].gsub!('X', '-')
        std.gsub!('X', '-')
    elsif(method == :muscle)
        elem = method_muscle(std, aa_seq_s)
        if(elem[0] != std)
            0.upto(std.size - 1) do |i|
                if(std[i, 1] != elem[0][i,1] and std[i,1] == '-')
                    elem[0].insert(i, '-')
                    elem[1].insert(i, '-') 
                    break
                end
            end
        end
    elsif(method == :nucleotide) #I dunno.....
        elem = method_nucleotide(std, aa_seq_s)
    end

   
    aa_seq_s = elem[1]
    aa_seq.delete_if {|a| a == ['-'] }

    0.upto(aa_seq_s.size - 1) do |i|
        aa_seq.insert(i, ['-']) if(aa_seq_s[i, 1] == '-')
    end
      
    if(removeinserts and elem[0].include?('-')) #cut insertions
        #wait, does this interfere with the malign?
        (elem[0].size - 1).downto(0) do |i|
            if(elem[0][i,1] == '-')
                elem[0][i,1] = ''
                elem[1][i,1] = ''
                aa_seq.delete_at(i)
                indels = true
            end
        end
    else
        #puts elem[1]
        #puts elem[0]
        #puts std
        return -2 if(elem[0] != std)
    end
    
    return aa_seq, indels
  end
  
  def align_nucs(std, seq, method=:recall, removeinserts=true)
    elem = align_it(std, seq, 3, 1)
    0.upto(elem[0].size - 1) do |i|
      if(elem[0][i,1] == '-')
        elem[0][i] ='X'
        elem[1][i] ='X'
      end
    end
    
#    puts "--------------------------------"
#    puts elem[0]
#    puts elem[1]
#    puts
    elem[0].gsub!('X','')
    elem[1].gsub!('X','')
#    puts elem[0]
#    puts elem[1]
    return elem
  end

def method_none(std, seq)
    return [std, seq]
end

#currently only works in windows, and also we aren't really using it.
#Also, I don't think it actually works.  I think its writing to the wrong spot(or not at all)
def method_muscle(std, seq)
    file = Tempfile.new('muscle_align')
    filename = file.path
    file.puts ">std\n#{std}\n>seq\n#{seq}"
    file.close 
    system("./muscle.exe -quiet -stable -in #{filename} -gapopen -3 -gapextend -1 -out #{filename}")

    newstd = ''
    newseq = ''
    state = 0
    File.open(filename) do |file|
        file.each_line do |line|
            if(line =~ /^>std/)
                state = 0
            elsif(line =~ /^>seq/)
                state = 1
            elsif(state == 0)
                newstd += line.strip
            elsif(state == 1)
                newseq += line.strip
            end
        end
    end
#    puts newstd
#    puts newseq
    return [newstd,newseq] #Need to set this up...
end

def method_recall(std, seq, gap_init=3,gap_ext=1)
    return align_it_aa(std, seq, gap_init, gap_ext)
end

#sinsi, g2p, x4r5
def load_matrix(label)
  #__FILE__
    filename = $:.map do |path|
      File.exists?("#{__FILE__.gsub('pssm_lib.rb', '')}/#{label}.matrix") ? "#{__FILE__.gsub('pssm_lib.rb', '')}/#{label}.matrix" : nil
    end
    filename = filename.delete_if{|a| a==nil}.sort[0]
    if(filename == nil)
      raise "File not found in path: #{label}.matrix"
    end
    matrix = nil
    File.open(filename) do |file|
        file.each_line do |line|
            next if(line =~ /^\#/ or line =~ /^\s*$/)
            vals = line.strip.split(/\t/)
            if(matrix == nil)
                matrix = Array.new(vals.size - 1) {Hash.new}
                matrix.each do |h|
                    h.default = 0.0
                end
            end
            matrix.each_with_index do |h, i|
                h[vals[0]] = vals[i + 1].to_f
            end
        end
    end
    
    return matrix
end

#This should do it
def run_pssm(seqs, std, matrix)
    aa_aligned = ''
    is_array = seqs.is_a?(Array)
    scores = []
    seqs = [seqs] if(!is_array)
    seqs.each do |seq|
        #puts std
        #puts seq
        aa, indels = align_aminos(std, seq, :recall)
        dels = 0
        score = nil
        dels = aa.find_all{|a| a == ['-']}.size if(aa != -1 and aa != -2)
        if(dels < 5 and aa != -1 and aa != -2)
            subseqs = malign(std, aa)
            score = subseqs.map {|s| [pssm(s, matrix)[1], s]}.max{|a,b| a[0] <=> b[0]}
            aa_aligned = score[1]
            score = score[0]
        else
            score = pssm(aa, matrix)[1] if(aa != -1 and aa != -2)
            aa_aligned = aa
        end
        scores.push(score)
    end
    if(!is_array)
        return scores[0], aa_aligned
    else
        return scores
    end
end
  
  
  

def run_g2p(seqs, std, matrix)
    aa_aligned = ''
    is_array = seqs.is_a?(Array)
    scores = []
    seqs = [seqs] if(!is_array)
    seqs.each do |seq|
      aa, indels = nil, nil
      if(seq.class == Array)
        aa, indels = align_aminos(std, seq, :recall, false, false)
        aa, indels = align_aminos(std, seq, :recall6, false, false) if(aa == -2 or aa == -1)
      else
        aa, indels = align_aminos(std, seq, :recall, false)
        aa, indels = align_aminos(std, seq, :recall6, false) if(aa == -2 or aa == -1)
      end
        #puts std
        #puts aa.join
        score = nil
        score = g2p(aa, matrix)[6] if(aa != -1 and aa != -2)
        scores.push(score)
        if(aa == -1)
          aa, indels = align_aminos(std, seq, :recall6, false, false) 
        end
        aa_aligned = aa
    end
    if(!is_array)
        return scores[0], aa_aligned
    else
        return scores
    end
  end
  
  def run_valid(seqs, std, matrix)
    aa_aligned = ''
    is_array = seqs.is_a?(Array)
    scores = []
    seqs = [seqs] if(!is_array)
    seqs.each do |seq|
        aa, indels = align_aminos(std,seq, :recall1, false)
        aa, indels = align_aminos(std,seq, :recall, false) if(aa == -2 or aa == -1)
        aa, indels = align_aminos(std,seq, :recall6, false) if(aa == -2 or aa == -1)
        score = nil
        #puts std
        #puts aa.join
        score = pssm_avg(aa, matrix) if(aa != -1 and aa != -2)
        #puts score
        scores.push(score)
        aa_aligned = aa
    end
    if(!is_array)
        return scores[0], aa_aligned
    else
        return scores
    end
  end
  
  #aaseq must be a string, no mixtures
  def score_valid_aa(aaseq, matrix)
    return pssm_avg(aaseq.split('').map{|a| [a]}, matrix)
  end
  
  def run_valid_nuc(seqs, std, matrix)
    nuc_aligned = ''
    is_array = seqs.is_a?(Array)
    scores = []
    seqs = [seqs] if(!is_array)
    seqs.each do |seq|
        begin
          elem = align_nucs(std,seq)
        #aa, indels = align_aminos(std,seq, :recall, false) if(aa == -2 or aa == -1)
        #aa, indels = align_aminos(std,seq, :recall6, false) if(aa == -2 or aa == -1)
          score = nil
          score = pssm_nuc(elem[1], matrix)
          scores.push(score)
          nuc_aligned = elem[1]
        rescue
          return [-1]
        end
    end
    if(!is_array)
        return scores[0], nuc_aligned
    else
        return scores
    end
  end


$g2p_fpr_data = []

def g2p_to_fpr(g2p)
    return nil if(g2p == nil)
    if($g2p_fpr_data == [])
	$: << '.'
        filename = $:.map do |path|
            File.exists?(path + "/g2p_fpr.txt") ? path + "/g2p_fpr.txt" : nil
        end

        raise "missing fpr file" if(filename == nil or filename.size == 0)
        filename = filename.delete_if{|a| a==nil}.sort[0]

        File.open(filename) do |file|
            file.each_line do |line|
                row = line.strip.split(',')
                row[0] = row[0].to_f
                row[1] = row[1].to_f
                $g2p_fpr_data << row
            end
        end
    end
    closest = nil
    $g2p_fpr_data.each do |d|
        closest = d if(closest == nil or (g2p - d[0]).abs < (g2p - closest[0]).abs)
    end
    return closest[1]
end

module Pssm

#returns an array of arrays,
#translate_complete_to_array('AAARRR') => [["K"], ["K", "R", "E", "G"]]
def Pssm.translate_complete_to_array(sequence)
	protseq = [];
	0.upto((sequence.size() / 3) - 1) do |i|
		list = Pssm.generate(sequence[i * 3, 3])
		list.map! { |entry| 
				if(entry == 'X')
					'X'
				else
					raw_translator(entry)
				end
		}.uniq!
		
		if(list.size > 1)
			tmp = []
			list.each do |entry|
				tmp.push(entry)
			end
			
			protseq.push(tmp)
		elsif(list.size == 1)
			protseq.push([list[0]])
		else
			protseq.push(['*'])
		end
	end

	return protseq;
end

def Pssm.raw_translator(str, keepdashs = false)
	aa = ''
	0.upto(((str.length - 1) / 3)) do |i|
        if(keepdashs and str[i * 3, 3] == '---')
            aa += '-'
            next
        end
		x = Pssm::AA_HASH[str[i * 3, 3].downcase]
		x = 'X' if(x == nil)
		aa += x
	end
	return aa;
end

def Pssm.generate(nuc)
	posa = Pssm::AMBIG[nuc[0,1]]
	posb = Pssm::AMBIG[nuc[1,1]]
	posc = Pssm::AMBIG[nuc[2,1]]
	
	if(nuc == 'XXX')
		return ['X']
	end
	nuclist = []
	posa.each do |a|
		posb.each do |b|
			posc.each do |c|
				nuclist.push(a + b + c)
			end
		end
	end
	return nuclist
rescue StandardError => error
	puts error
	puts nuc
	return nil
end



AA_HASH = {
      'ttt' => 'F', 'tct' => 'S', 'tat' => 'Y', 'tgt' => 'C',
      'ttc' => 'F', 'tcc' => 'S', 'tac' => 'Y', 'tgc' => 'C',
      'tta' => 'L', 'tca' => 'S', 'taa' => '*', 'tga' => '*',
      'ttg' => 'L', 'tcg' => 'S', 'tag' => '*', 'tgg' => 'W',

      'ctt' => 'L', 'cct' => 'P', 'cat' => 'H', 'cgt' => 'R',
      'ctc' => 'L', 'ccc' => 'P', 'cac' => 'H', 'cgc' => 'R',
      'cta' => 'L', 'cca' => 'P', 'caa' => 'Q', 'cga' => 'R',
      'ctg' => 'L', 'ccg' => 'P', 'cag' => 'Q', 'cgg' => 'R',

      'att' => 'I', 'act' => 'T', 'aat' => 'N', 'agt' => 'S',
      'atc' => 'I', 'acc' => 'T', 'aac' => 'N', 'agc' => 'S',
      'ata' => 'I', 'aca' => 'T', 'aaa' => 'K', 'aga' => 'R',
      'atg' => 'M', 'acg' => 'T', 'aag' => 'K', 'agg' => 'R',

      'gtt' => 'V', 'gct' => 'A', 'gat' => 'D', 'ggt' => 'G',
      'gtc' => 'V', 'gcc' => 'A', 'gac' => 'D', 'ggc' => 'G',
      'gta' => 'V', 'gca' => 'A', 'gaa' => 'E', 'gga' => 'G',
      'gtg' => 'V', 'gcg' => 'A', 'gag' => 'E', 'ggg' => 'G',
}

AMBIG = { 'A'=> ['A'], 'G'=> ['G'], 'T'=> ['T'], 'C'=> ['C'],
'R' => ['A', 'G'], 'Y' => ['C', 'T'], 'K' => ['G', 'T'],
'M' => ['A', 'C'], 'B' => ['C', 'G', 'T'], 'D' => ['A', 'G', 'T'],
'H' => ['A', 'C', 'T'], 'V' => ['A', 'C', 'G'], 'S' => ['C', 'G'],
'W' => ['A', 'T'], 'N' => ['A', 'C', 'G', 'T'], 'X' => ['X'], '-' => ['-'] }

end
