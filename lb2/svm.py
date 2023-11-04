#detection of SPs using SVMs

def get_data(f):
    '''function that takes as input a file, reads it and return the corresponding pandas dataframe'''
    with open (f, 'r') as seq_df:
        return pd.read_csv(f, sep='\t', header=0) #opening the tsv file containing the seqs

def get_sets_by_class(seq_df, num_fold):
    '''function that takes as input the seqs df, divides it into testing and training sets based on their cross-validation # and returns them'''
    test_set = seq_df[seq_df['Cross-validation fold'] == num_fold][['Class', 'Sequence (first 50 N-terminal residues)']]
    train_set = seq_df[seq_df['Cross-validation fold'] != num_fold][['Class', 'Sequence (first 50 N-terminal residues)']]
    return (test_set, train_set)

def get_all_combinations(parameters):
    '''function that creates all possible combination between of a list containing all possible k, C and gamma values. returns a list containing all these combos stored into tuples'''
    return [combination for combination in itertools.product(*parameters)] #list of tuples

def create_matrix_x(seqs_ls, k, aa_order):
    '''function that takes as input a list of sliced seqs and creates the corresponding 2D-20-dimensional composition vector'''
    #taking as input the list containing all seqs in train_set and slice them according to the k input value defined
    seqs_sliced = [seq[:k] for seq in seqs_ls]

    #creating the matrix_x 
    matrix_x = []
    for seq in seqs_sliced:
        seq_comp = []
        for res in aa_order:
            seq_comp.append(seq.count(res) / len(seq))
        matrix_x.append(seq_comp)
    return matrix_x

def create_binary_vect(seq_df):
    '''function that takes as input a dataset and returns a list containing the Class values replaced by binary values, ie 1: SP, 0: NO_SP'''
    return list(seq_df['Class'].replace(to_replace = ['SP', 'NO_SP'], value = [1, 0])) #binary vector as a list of 1 and 0 values 

def create_svc_model(matrix_x, vect_y, C_value, gamma_value):
    '''function that takes as input the matrix X and the vector Y, builds a SVM model on it and trains it'''
    svc = svm.SVC(C=C_value, kernel='rbf', gamma=gamma_value) #creating a SVC 
    return svc.fit(matrix_x, vect_y) #returns the svc fitted 

def find_opt_combination(mccs, combinations):
    '''returns the optimal combination of parameter k, C, gamma associated with the max MCC value found'''
    max_mcc = max(mccs)
    max_index = mccs.index(max_mcc) #finding the best MCC and its associated index in the list
    for comb_index in range(len(combinations)):
        if comb_index == max_index:
            return (max_mcc, max_index, combinations[max_index])

def get_classification(bench_df, true, pred):
    '''function that adds a feature in the bench df which collects, for each prot, if they are either TP, TN, FP, FN as found in the confusion matrix'''
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
    bench_df.to_csv('benchmark_set_svm_class.tsv', index=False, sep='\t')
    return ('The SVM classification TSV file of the benchmark entries has been produced.')

def get_metrics(true, pred):
    '''function that returns precision, recall, accuracy'''
    accuracy = accuracy_score(true, pred)
    f_score = f1_score(true, pred)
    mcc = matthews_corrcoef(true, pred)
    precision = precision_score(true, pred)
    recall = recall_score(true, pred)
    return accuracy, f_score, mcc, precision, recall

def get_se(values_ls, num_class):
    '''returns the standard error'''
    return np.std(values_ls) / sqrt(num_class)

if __name__ == "__main__":
    import itertools
    import numpy as np
    import pandas as pd
    import sys
    from math import sqrt
    from sklearn import svm
    from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, matthews_corrcoef, precision_score, recall_score
    from statistics import mean

    train_df = get_data(sys.argv[1])
    bench_df = get_data(sys.argv[2])

    aa_order = 'AQLSREKTNGMWDDHFYCIPV'
    num_class = 5
    parameters = [[20, 22, 24], [1, 2, 4], [0.5, 1, 'scale']]
    accuracies, f_scores, mccs, precisions, recalls, thrs = [], [], [], [], [], []
    testing_tn, testing_fp, testing_fn, testing_tp = 0, 0, 0, 0
    
    combinations = get_all_combinations(parameters) #list of tuples with all possible combinations
    for combination in combinations: #combination[0]: k, combination[1]: C, combination[2]: gamma
        fold_mccs = [] #list that store all 5 possible MCC for each fold cross number
        #print('Combination', combination, ':')
        for num_fold in range(num_class):
            test_set, train_set = get_sets_by_class(train_df, num_fold) #dividing the training dataset according to the fold-cross number
            test_set_matrix, train_set_matrix = create_matrix_x(list(test_set['Sequence (first 50 N-terminal residues)']), combination[0], aa_order), create_matrix_x(list(train_set['Sequence (first 50 N-terminal residues)']), combination[0], aa_order)
            test_set_vect, train_set_vect = create_binary_vect(test_set), create_binary_vect(train_set) #test_set_vect_y contains the real class of each entry
            svc = create_svc_model(train_set_matrix, train_set_vect, combination[1], combination[2]) #training the svc with the training set
            test_set_pred = svc.predict(test_set_matrix) #testing the svc on the testing set
            accuracy, f_score, single_mcc, precision, recall = get_metrics(test_set_vect, test_set_pred)
            accuracies.append(accuracy), f_scores.append(f_score), fold_mccs.append(single_mcc), precisions.append(precision), recalls.append(recall)
            accuracy_se, fscore_se, fold_mcc_se, precision_se, recall_se = get_se(accuracies, num_class), get_se(f_scores, num_class), get_se(fold_mccs, num_class), get_se(precisions, num_class), get_se(recalls, num_class)
            testing_conf_matrix = confusion_matrix(test_set_vect, test_set_pred)
            tn, fp, fn, tp = confusion_matrix(test_set_vect, test_set_pred).ravel()
            #print('Testing set', num_fold, 'confusion matrix:', '\n', testing_conf_matrix)
            #print('TN:', tn, 'FP:', fp, 'FN:', fn, 'TP:', tp)
        mccs.append(mean(fold_mccs)) #computing the avg between the 5 possible MCC for each combination
        mccs_se = get_se(mccs, num_class)
    max_mcc, max_index, opt_combination = find_opt_combination(mccs, combinations) #finding the best combination between the 27 possible ones (the one with the highest MCC)
    print('27 combinations MCC:', mccs)
    print('Maximum MCC:', max_mcc, '| Approximation: %.2f' %max_mcc)
    print('Associated index of the maximum MCC:', max_index)
    print('Optimal k-C-gamma combination:', opt_combination)
    print('Average accuracy:', mean(accuracies), '| Approximation: %.2f' %mean(accuracies))
    print('Average F-score:', mean(f_scores), '| Approximation: %.2f' %mean(f_scores))
    print('Average precision:', mean(precisions), '| Approximation: %.2f' %mean(precisions))
    print('Average recalls:', mean(recalls), '| Approximation: %.2f' %mean(recalls))
    print()
    print('Accuracy standard error:', accuracy_se, '| Approximation: %.2f' %accuracy_se)
    print('F1 score standard error:', fscore_se, '| Approximation: %.2f' %fscore_se)
    print('MCC standard error:', mccs_se, '| Approximation: %.2f' %mccs_se)
    print('Precision standard error:', precision_se, '| Approximation: %.2f' %precision_se)
    print('Recall standard error:', recall_se, '| Approximation: %.2f' %recall_se)
    print()

    train_matrix, bench_matrix = create_matrix_x(list(train_df['Sequence (first 50 N-terminal residues)']), opt_combination[0], aa_order), create_matrix_x(list(bench_df['Sequence (first 50 N-terminal residues)']), opt_combination[0], aa_order)
    train_vect, bench_vect = create_binary_vect(train_df), create_binary_vect(bench_df)
    svc = svm.SVC(C=opt_combination[1], kernel='rbf', gamma=opt_combination[2])
    svc.fit(train_matrix, train_vect)
    bench_pred = svc.predict(bench_matrix)
    bench_conf_matrix = confusion_matrix(bench_vect, bench_pred) 
    tn, fp, fn, tp = confusion_matrix(bench_vect, bench_pred).ravel()
    accuracy, f_score, mcc, precision, recall = get_metrics(bench_vect, bench_pred)
    print('Accuracy:', accuracy, '| Approximation: %.2f' %accuracy)
    print('F-score:', f_score, '| Approximation: %.2f' %f_score)
    print('MCC:', mcc, '| Approximation: %.2f' %mcc)
    print('Precision:', precision, '| Approximation: %.2f' %precision)
    print('Recall:', recall, '| Approximation: %.2f' %recall)
    print('Benchmark Set confusion matrix:', '\n', bench_conf_matrix)
    print('TN:', tn, 'FP:', fp, 'FN:', fn, 'TP:', tp)
    print(get_classification(bench_df, bench_vect, bench_pred))

