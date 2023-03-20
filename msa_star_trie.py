import copy, math
import random
import time
from aho_corasick import AhoCorasick

SEGMENT_LENGTH = 3
S_MATCH = 3
S_MISMATCH = -1
S_GAP = -2
S_GAP_GAP = 0
gap_indices = []
seq_num = 0

def global_align(x, y, s_match, s_mismatch, s_gap):
    global gap_indices
    A = []
    for i in range(len(y) + 1):
        A.append([0] * (len(x) + 1))
    for i in range(len(y) + 1):
        A[i][0] = s_gap * i
    for i in range(len(x) + 1):
        A[0][i] = s_gap * i
    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            A[i][j] = max(
                A[i][j - 1] + s_gap,
                A[i - 1][j] + s_gap,
                A[i - 1][j - 1] + (s_match if (y[i - 1] == x[j - 1] and y[i - 1] != '-') else 0) + (
                    s_mismatch if (y[i - 1] != x[j - 1] and y[i - 1] != '-' and x[j - 1] != '-') else 0) + (
                    s_gap if (y[i - 1] == '-' or x[j - 1] == '-') else 0)
            )
    align_X = ""
    align_Y = ""
    i = len(x)
    j = len(y)

    gap_indices = []
    while i > 0 or j > 0:
        current_score = A[j][i]
        if i > 0 and j > 0 and (
                ((x[i - 1] == y[j - 1] and y[j - 1] != '-') and current_score == A[j - 1][i - 1] + s_match) or
                ((y[j - 1] != x[i - 1] and y[j - 1] != '-' and x[i - 1] != '-') and current_score == A[j - 1][i - 1] + s_mismatch) or
                ((y[j - 1] == '-' or x[i - 1] == '-') and current_score == A[j - 1][i - 1] + s_gap)
        ):
            align_X = x[i - 1] + align_X
            align_Y = y[j - 1] + align_Y
            i = i - 1
            j = j - 1
        elif i > 0 and (current_score == A[j][i - 1] + s_gap):
            align_X = x[i - 1] + align_X
            if len(gap_indices) != 0:
                gap_indices = [x+1 for x in gap_indices]
                gap_indices.append(j)
            else:
                gap_indices.append(j)
            align_Y = "-" + align_Y
            i = i - 1
        else:
            align_X = "-" + align_X
            align_Y = y[j - 1] + align_Y
            j = j - 1
    return (align_X, align_Y, A[len(y)][len(x)])

if __name__=="__main__":
    
    seq_num = int(input())
    seqs = []

    for _ in range(seq_num):
        seqs.append(input())

    # time to align the sequences
    start_time = time.time()

    # make all sequences lower case
    seqs = [x.lower() for x in seqs]

    # first break each sequence into k-mers of length k

    seqs_kmers = []
    for i in range(seq_num):
        seqs_kmers.append([seqs[i][j:j+SEGMENT_LENGTH] for j in range(0, int(len(seqs[i])/SEGMENT_LENGTH)*SEGMENT_LENGTH, SEGMENT_LENGTH)])

    # now construct trie trees for each sequence
    trie_trees = []
    for i in range(seq_num):
        trie_trees.append(AhoCorasick(seqs_kmers[i]))

    # calculate the score of each sequence
    seqs_score = []
    for i in range(seq_num):
        score = 0;
        for j in range(seq_num):
            if i != j:
                result = trie_trees[j].search_words(seqs[i])
                for word in result:
                    score += len(result[word])
        seqs_score.append(score)

    # get the centre sequence
    centre_seq_index = seqs_score.index(max(seqs_score))
    centre_seq = seqs[centre_seq_index]
    # print("Centre Sequence: " + centre_seq)

    # compute positions of matches in centre sequence
    aho_result = []
    for i in range(seq_num):
        if(i!=centre_seq_index):
            ahocorasick_result = trie_trees[i].search_words(centre_seq)
            aho_result.append(ahocorasick_result)
        else:
            aho_result.append(None)

    alignments_with_center = {}
    for i in range(seq_num):
        if i != centre_seq_index:

            for j in range(len(seqs_kmers[i])):
                seg = seqs_kmers[i][j]
                min_index = 100000000
                for index in aho_result[i][seg]:
                    if index < min_index:
                        min_index = index
                aho_result[i][seg] = [min_index]

            align_X = ""
            # align_Y will be the center sequence
            align_Y = "" 

            x_index = 0
            y_index = 0

            for j in range(len(seqs_kmers[i])):
                seg = seqs_kmers[i][j]
                index = aho_result[i][seg][0]
                if index == 100000000:
                    continue
                else:
                    if(index < y_index):
                        continue
                    unmatched_X = seqs[i][x_index:SEGMENT_LENGTH*j]
                    unmatched_Y = centre_seq[y_index:index]
                    matched_X, matched_Y, s = global_align(unmatched_X, unmatched_Y, S_MATCH, S_MISMATCH, S_GAP)
                    align_X += matched_X+seg
                    align_Y += matched_Y+seg
                    x_index = SEGMENT_LENGTH*j + SEGMENT_LENGTH
                    y_index = index + SEGMENT_LENGTH

            unmatched_X = seqs[i][x_index:]
            unmatched_Y = centre_seq[y_index:]
            matched_X, matched_Y, s = global_align(unmatched_X, unmatched_Y, S_MATCH, S_MISMATCH, S_GAP)
            align_X += matched_X
            align_Y += matched_Y

            alignments_with_center[i] = [align_X, align_Y, s]

    # need to allign all centre sequences now 
    curr_pos=0;
    chars = len(centre_seq)
    while(chars>0):
        gap_flag = False
        for i in range(seq_num):
            if i != centre_seq_index:
                if alignments_with_center[i][1][curr_pos] == '-':
                    gap_flag = True
                    break

        if gap_flag == True:
            for i in range(seq_num):
                if i != centre_seq_index:
                    if(alignments_with_center[i][1][curr_pos] != '-'):
                        alignments_with_center[i][0] = alignments_with_center[i][0][:curr_pos] + '-' + alignments_with_center[i][0][curr_pos:]
                        alignments_with_center[i][1] = alignments_with_center[i][1][:curr_pos] + '-' + alignments_with_center[i][1][curr_pos:]
        else:
            chars = chars - 1
        curr_pos += 1
    
    end_time = time.time()
    print("")
    print("Time taken to align sequences: " + str(end_time-start_time))

    print("")
    print("The aligned sequences are: ")
    print("")
    
    for i in range(seq_num):
        if i!=centre_seq_index:
            print(alignments_with_center[i][0].upper())
        else:
            index = centre_seq_index
            while(index == centre_seq_index):
                index = random.randint(0, seq_num-1)
            print(alignments_with_center[index][1].upper())


