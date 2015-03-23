#! python3

# Scoring assumptions: The scoring gives equal weight to each of the three
# scoring criteria. The higher the score the closer to ideal the sequence
# is. 


# check_homopolymers: uses a 5 nucleotide sliding window to check for
# homopolymers if any are found in the sequence True is returned and
# the sequence is disregarded.
def check_homopolymers(sequence):
    homopolymers = False
    four_bp_window_list = [sequence[i:i+5] for i in range(len(sequence)-4)]
    for four_bp_window in four_bp_window_list:
        bp_list = []
        for chr in four_bp_window:
            bp_list.append(chr)
            
        if((all(x == bp_list[0] for x in bp_list)) is True):
            homopolymers = True
    return homopolymers

# check_for_ATF_motif: Checks if the motif is in the sequence at all
# if found True is returned and the sequence is disregarded.
def check_for_ATG_motif(sequence):
    if('ATG' in sequence):
        return True
    else:
        return False

# get_GC_content: returns the percent GC of the given sequence.
def get_GC_content(sequence):
    GC_sum = float((sequence.count('G')+sequence.count('C')))
    GC_content = GC_sum/float(len(sequence))*100
    return GC_content

# check_GC_content: Checks that the GC content of the sequence is in
# the range from 20% to 80%.
def check_GC_content(sequence):
    GC_content = get_GC_content(sequence)
    if(GC_content >= 20 and GC_content <= 80):
        return True
    else:
        return False

# check_twentyfirst_nt: checks if the 21st nucleotide is a G or C.
# Uses 0 base numbering for the actual comparison. If G or C is
# present at this position a score of 100 is returned.
def check_twentyfirst_nt(sequence):
    if(sequence[20] == 'G' or sequence[20] == 'C'):
        return 100
    else:
        return 0
    
# GC_score_entire_sequence: Determines how close a sequence is to 55%
# GC content. The max score for a perfect 55% GC content would be 100.
def GC_score_entire_sequence(sequence):
    return 100-abs(55-get_GC_content(sequence))

# GC_score_in_range: score based on the nucleotides in positions 15-20.
# Using 0 based numbering for analysis.
def GC_score_in_range(sequence):
    return (get_GC_content(sequence[14:19]))

#Input is assumed to be a valid DNA sequence containing only ATCG
def main():
##    sequence = ''
##    sequence = raw_input('DNA sequence: ') #python2 input
    sequence = input('DNA sequence: ') #python3 input
    sequence = sequence.upper()
    sequences_ending_in_GG = []
    #sliding window of width 23 used to find potential subsequences.
    twentythree_nt_windows = [sequence[i:i+23] for i in range(len(sequence)-22)]
    #determination of if the subsequence ends in GG
    for i in twentythree_nt_windows:
        if(i[-2] == 'G' and i[-1] == 'G'):
            sequences_ending_in_GG.append(i)
    print ('23-mer: %-17s score: %-5s'%('', ''))
    #filering: checks the 3 filtering criteria before scoring.
    for twentythree_nt_sequence in sequences_ending_in_GG:
        homopolymers = check_homopolymers(twentythree_nt_sequence)
        ATG_motif = check_for_ATG_motif(twentythree_nt_sequence)
        GC_content = check_GC_content(twentythree_nt_sequence)
        if(homopolymers is False and ATG_motif is False and GC_content is True):
            #scoring: sums up the score from each of the three criteria
            score = 0
            score = check_twentyfirst_nt(twentythree_nt_sequence)
            score = score + GC_score_entire_sequence(twentythree_nt_sequence)
            score = score + GC_score_in_range(twentythree_nt_sequence)
            print ('%-25s %-5s'%(twentythree_nt_sequence, str(round(score, 2))))

main()
