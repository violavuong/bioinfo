# implementation of the vonHeijne algorithm
# Cross/validation of the training set. every protein in the training set is labeled with a class, ie. 0, 1, 2, 3, 4 as we perform a 5-cross validation. So for run 0, the training set are 1, 2, 3, 4 while the testing one is 0 and so on. Each training set has its own threshold which is used to test the testing set from which a confusion matrix is build, ie. each run has its own threshold and confusion matrix. The avg threshold of all threshold obtained will be use to evaluate the benchmark set later on. Remember: each set of trainin sets should have its own position-specific weight matrix; also the threshold should not be optimized on the testing data, as you cannot obtain then an avg performance of your whole data.

def get_data(f):
    '''function that takes as input a file, reads it and returns the corresponding pandas df'''
    with open(f, 'r') as seq_df:
        return pd.read_csv(f, sep='\t', header=0)

def get_sets_by_class(seq_df, num_fold):
    '''function that takes as input the seqs df, divided it into testing and training sets according to the fold number'''
    test_set = seq_df[seq_df['Cross-validation fold'] == num_fold][['Class', 'Sequence (first 50 N-terminal residues)', 'SP cleavage-site annotation']]
    train_set = seq_df[seq_df['Cross-validation fold'] != num_fold][['Class', 'Sequence (first 50 N-terminal residues)', 'SP cleavage-site annotation']]
    return (test_set, train_set)

def get_SP_cleaved_sites(train_set):
    '''function that takes as input the training set and returns a list containing the SP proteins cleaved sites'''
    sp_len = list(train_set.query('Class=="SP"')['SP cleavage-site annotation'].str.count('S')) #creating a list containing the SP length of the seqs read
    sp_seq = list(train_set.query('Class=="SP"')['Sequence (first 50 N-terminal residues)']) #creating another list containing all SP seqs
    sp_sites = []
    for length, seq in zip(sp_len, sp_seq): #filling the above list with the SP seqs cut according to the cleavage site
        sp_sites.append(seq[length - 13: length + 2])
    return sp_sites

def create_pswm(sp_sites, aa_bg_comp):
    '''function that takes as input a list containing all SP proteins cleaved sites and returns the corresponding PSWM matrix'''
    bg_keys_ls = []
    for key, value in aa_bg_comp.items(): #extracting the key of the dict, these will become the columns
        bg_keys_ls.append(key)
    matrix = np.ones([len(sp_sites[0]), len(bg_keys_ls)])
    for seq in sp_sites: #filling the matrix
        for row, aa in enumerate(seq):
            matrix[row][bg_keys_ls.index(aa)] += 1
    matrix = matrix / (len(bg_keys_ls) + len(sp_sites)) #building the PSPM
    for res in range(0, len(sp_sites[0])): #dividing by the bg freq
        for key, value in aa_bg_comp.items():
            matrix[res][bg_keys_ls.index(key)] /= value
    matrix = np.log2(matrix) #computing the log to obtain the PSWM
    return matrix

def get_score(train_set, matrix, aa_order):
    '''function that takes every seq in the set, slice it into 15-res long subseqs, computes their score, find the highest one and assign this max score to the entire original seq. returns a dict which keys are the seq and values the corresponding score'''
    window_len = 15
    seqs_score = {}
    for seq in list(train_set['Sequence (first 50 N-terminal residues)']):
        window_seq = '' #empty string to store the subseq to be analyzed
        scores = []
        for i in range(36): #slicing the seq: 50 - 15 + 1 
            window_seq = seq[i:i + window_len]
            score = 0
            for j in range(len(window_seq)):
                aa_index = aa_order.find(window_seq[j]) #finding the pos of the res in the aa order
                score += matrix[j][aa_index]
            scores.append(score)
        seqs_score[seq] = max(scores)
    return seqs_score

def create_binary_vect(df):
    '''returns a binary vect of the input dataframe, ie SP=1, NO_SP=0'''
    return list(df['Class'].replace(to_replace = ['SP', 'NO_SP'], value = [1, 0]))

def get_thr(binary_vect, train_set_scores):
    '''function that takes as input the list of binary classification and the dict containing, for each seq in the training set, the best window of 15 res with the highest score and finds the best threshold according to the following optimized metrics (F1 score, precision, recall)'''
    precision, recall, thresholds = precision_recall_curve(binary_vect, list(train_set_scores.values())) #precision: contains precision scores at varying threshold; recall: contains recall scores at varying thresholds; thresholds: the thresholds values
    fscore = (2 * precision * recall) / (precision + recall) #compute f-scores at varying thresholds
    index = np.argmax(fscore) #get the index of the maximum value of the f-score
    opt_thr = thresholds[index] #retrieve the threshold value corresponding to the max f-score computed above
    return opt_thr

def get_classification(bench_df, true, pred):
    '''function that classify each entry in dataset as either TN, FP, FN, TP. returns the dataset with the column specifying this'''
    class_order = []
    for true_value, pred_value in zip(true, pred):
        if true_value == pred_value == 1:
            class_order.append('TP')
        elif true_value == pred_value == 0:
            class_order.append('TN')
        elif true_value == 0 and pred_value == 1:
            class_order.append('FP')
        else:
            class_order.append('FN')
    bench_df = bench_df.assign(Classification=class_order)
    bench_df.to_csv('benchmark_set_vh_class.tsv', index=False, sep='\t')
    return ('The VH classification TSV file of the benchmark entries has been produced.')

def get_metrics(true, pred):
    '''returns precision, recall, accuracy and MCC'''
    accuracy = accuracy_score(true, pred)
    f_score = f1_score(true, pred)
    mcc = matthews_corrcoef(true, pred)
    precision = precision_score(true, pred)
    recall = recall_score(true, pred)
    return accuracy, f_score, mcc, precision, recall

def get_se(values_ls, num_class):
    '''returns the standard error'''
    return np.std(values_ls) / sqrt(num_class)
  
if __name__ == '__main__':
    import numpy as np
    import pandas as pd
    import sys
    from math import sqrt
    from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, matthews_corrcoef, precision_recall_curve, precision_score, recall_score
    from statistics import mean
    
    train_df = get_data(sys.argv[1])
    bench_df = get_data(sys.argv[2])

    aa_bg_comp = {'A': 0.0825, 'Q': 0.0393, 'L': 0.0965, 'S': 0.0664, 'R': 0.0553, 'E': 0.0672, 'K': 0.058, 'T': 0.0535, 'N': 0.0406, 'G': 0.0707, 'M': 0.0241, 'W': 0.011, 'D': 0.0546, 'H': 0.0227, 'F': 0.0386, 'Y': 0.0292, 'C': 0.0138, 'I': 0.0591, 'P': 0.0474, 'V': 0.0686}
    aa_order = 'AQLSREKTNGMWDHFYCIPV'
    num_class = 5
    accuracies, f_scores, mccs, precisions, recalls, thrs = [], [], [], [], [], []

    for num_fold in range(num_class):
        test_set, train_set = get_sets_by_class(train_df, num_fold)
        matrix = create_pswm(get_SP_cleaved_sites(train_set), aa_bg_comp)
        test_set_scores, train_set_scores = get_score(test_set, matrix, aa_order), get_score(train_set, matrix, aa_order)
        test_set_vect, train_set_vect = create_binary_vect(test_set), create_binary_vect(train_set)
        opt_thr = get_thr(create_binary_vect(train_set), train_set_scores)
        test_set_pred = [int(score >= opt_thr) for score in list(test_set_scores.values())]
        #computing the necessary metrics
        thrs.append(opt_thr)
        accuracy, f_score, mcc, precision, recall = get_metrics(test_set_vect, test_set_pred)
        accuracies.append(accuracy), f_scores.append(f_score), mccs.append(mcc), precisions.append(precision), recalls.append(recall)
        accuracy_se, fscore_se, mcc_se, precision_se, recall_se = get_se(accuracies, num_class), get_se(f_scores, num_class), get_se(mccs, num_class), get_se(precisions, num_class), get_se(recalls, num_class)
        print('Optimal threshold for testing set', num_fold,':', opt_thr)
        print('Confusion matrix of testing set', num_fold)
        print(confusion_matrix(test_set_vect, test_set_pred))
        tn, fp, fn, tp = confusion_matrix(test_set_vect, test_set_pred).ravel()
        print('TN:', tn, 'FP:', fp, 'FN:', fn, 'TP:', tp)
        print()
        
    print('Thresholds:', thrs)
    print('Average threshold:', mean(thrs), '| Approximation: %.2f' %mean(thrs))
    print('Average accuracy:', mean(accuracies), '| Approximation: %.2f' %mean(accuracies))
    print('Average F-score:', mean(f_scores), '| Approximation: %.2f' %mean(f_scores))
    print('Average MCC:', mean(mccs), '| Approximation: %.2f' %mean(mccs))
    print('Average precision:', mean(precisions), '| Approximation: %.2f' %mean(precisions))
    print('Average recalls:', mean(recalls), '| Approximation: %.2f' %mean(recalls))
    print()
    print('Accuracy standard error:', accuracy_se, '| Approximation: %.2f' %accuracy_se)
    print('F1 score standard error:', fscore_se, '| Approximation: %.2f' %fscore_se)
    print('MCC standard error:', mcc_se, '| Approximation: %.2f' %mcc_se)
    print('Precision standard error:', precision_se, '| Approximation: %.2f' %precision_se)
    print('Recall standard error:', recall_se, '| Approximation: %.2f' %recall_se)
    print()

    bench_matrix = create_pswm(get_SP_cleaved_sites(train_df), aa_bg_comp)
    bench_df_scores = get_score(bench_df, bench_matrix, aa_order)
    bench_df_vect = create_binary_vect(bench_df)
    bench_df_pred = [int(score >= opt_thr) for score in list(bench_df_scores.values())]
    
    tn, fp, fn, tp = confusion_matrix(bench_df_vect, bench_df_pred).ravel()
    accuracy, f_score, mcc, precision, recall = get_metrics(bench_df_vect, bench_df_pred)
    print('Benchmark accuracy:', accuracy, '| Approximation: %.2f' %accuracy)
    print('Benchmark F-score:', f_score, '| Approximation: %.2f' %f_score)
    print('Benchmark MCC:', mcc, '| Approximation: %.2f' %mcc)
    print('Benchmark precision:', precision, '| Approximation: %.2f' %precision)
    print('Benchmark recall:', recall, '| Approximation: %.2f' %recall)
    print('Benchmark Set confusion matrix:')
    print(confusion_matrix(bench_df_vect, bench_df_pred))
    print('TN:', tn, 'FP:', fp, 'FN:', fn, 'TP:', tp)
    print(get_classification(bench_df, bench_df_vect, bench_df_pred))
    
    all_train_df_scores = get_score(train_df, matrix, aa_order)
    all_train_df_vect = create_binary_vect(train_df)
    all_train_df_pred = [int(score >= opt_thr) for score in list(all_train_df_scores.values())]



