# python program that evaluates and discuss the FP and FN found in both methods (von Heijne algo, SVM)
# gonna conduct the analysis only for the benchmarking dataset

def get_fpr(f1, f2):
    '''returns the fpr rate'''
    with open(f1, 'r') as txt_fp:
        fps = len(txt_fp.readlines())
    with open(f2, 'r') as txt_tn:
        tns = len(txt_tn.readlines())
    return (fps / (fps + tns))

def get_fp_tm_ratio(f1, f2, ecos):
    '''function that returns the FP-TM ratio'''
    fp_tm, tm = 0, 0
    fh = open(f1, 'r')
    fn = open(f2, 'r')

    for line in fh:
        line = line.rstrip()
        for eco in ecos:
            if eco in line and re.match('[1-9]\.|[1234][0-9]|50', line[23:25]):
                fp_tm += 1
                break

    for line in fn:
        line = line.rstrip()
        for eco in ecos: 
            if eco in line and re.match('[1-9]\.|[1234][0-9]|50', line[23:25]):
                tm += 1
                break

    return (fp_tm, tm, fp_tm / tm)

def get_ts(f1, f2, ecos):
    '''function that finds and returns the # of mitho, chloro and perox seqs'''
    fp_chloro, fp_mitho, fp_perox, fp_ts_tot = 0, 0, 0, 0
    chloro, mitho, perox, ts_tot = 0, 0, 0, 0
    fh = open(f1, 'r')
    fn = open(f2, 'r')

    for line in fh:
        line = line.rstrip()
        fp_ts_tot += 1
        for eco in ecos:
            if 'Chloroplast' in line and eco in line:
                fp_chloro += 1
            elif 'Mitochondrion' in line and eco in line: 
                fp_mitho += 1
            elif 'Peroxisome' in line and eco in line:
                fp_perox += 1
    print('FP Chloro | FP Mitho | FP Perox | FP Tot')
    print(fp_chloro, '|', fp_mitho, '|', fp_perox, '|', fp_ts_tot)
    
    for line in fn:
        ts_tot += 1
        for eco in ecos:
            if 'Chloroplast' in line and eco in line:
                chloro += 1
            elif 'Mitochondrion' in line and eco in line:
                mitho += 1
            elif 'Peroxisome' in line and eco in line:
                perox += 1
    print('TN Chloro | TN Mitho | TN Perox | TN Tot')
    print(chloro, '|', mitho, '|', perox, '|', ts_tot)
    try:
        fpr_c = fp_chloro / (fp_chloro + chloro)
        fpr_m = fp_mitho / (fp_mitho + mitho)
        fpr_p = fp_perox / (fp_perox + perox)
        fpr_ts = fp_ts_tot / (fp_ts_tot + ts_tot)
    except ZeroDivisionError:
        fpr_c = 0
        fpr_m = 0
        fpr_p = 0
        fpr_ts = 0
    return (fpr_c, fpr_m, fpr_p, fpr_ts)

if __name__ == '__main__':
    import re
    import sys

    ecos = ['ECO:0000269','ECO:0000303','ECO:0000305','ECO:0000250','ECO:0000255','ECO:0000312','ECO:0007744']
    
    #vonHeijne 
    vh_fp = sys.argv[1]
    vh_tn = sys.argv[2]
    vh_fp_tm = sys.argv[3]
    vh_tn_tm = sys.argv[4]
    vh_fp_ts = sys.argv[5]
    vh_tn_ts = sys.argv[6]
    vh_fp_tm, vh_tm, vh_ratio = get_fp_tm_ratio(vh_fp_tm, vh_tn_tm, ecos)
    vh_fpr_c, vh_fpr_m, vh_fpr_p, vh_fpr_ts = get_ts(vh_fp_ts, vh_tn_ts, ecos)
    print('vh FPR:', get_fpr(vh_fp, vh_tn), '| Approximation %.2f' %get_fpr(vh_fp, vh_tn))
    print('vh FP TM:', vh_fp_tm)
    print('vh tot TM:', vh_tm)
    print('vh FP_TM-TM ratio: ', vh_ratio, '| Approximation %.2f' %vh_ratio)
    print('vh chloroplasts FPR:', vh_fpr_c, '| Approximation %.2f' %vh_fpr_c)
    print('vh mithocondrion FPR:', vh_fpr_m, '| Approximation %.2f' %vh_fpr_m)
    print('vh peroxisomes FPR:', vh_fpr_p, '| Approximation %.2f' %vh_fpr_p)
    print('vh ts FPR:', vh_fpr_ts, '| Approximation %.2f' %vh_fpr_ts)
    print()

    #SVM
    svm_fp = sys.argv[7]
    svm_tn = sys.argv[8]
    svm_fp_tm = sys.argv[9]
    svm_tn_tm = sys.argv[10]
    svm_fp_ts = sys.argv[11]
    svm_tn_ts = sys.argv[12]
    svm_fp_tm, svm_tm, svm_ratio = get_fp_tm_ratio(svm_fp_tm, svm_tn_tm, ecos)
    svm_fpr_c, svm_fpr_m, svm_fpr_p, svm_fpr_ts = get_ts(svm_fp_ts, svm_tn_ts, ecos)
    print('svm FPR:', get_fpr(svm_fp, svm_tn), '| Approximation %.2f' %get_fpr(svm_fp, svm_tn))
    print('svm FP TM:', svm_fp_tm)
    print('svm tot TM:', svm_tm)
    print('svm FP_TM-TM ratio: ', svm_ratio, '| Approximation %.2f' %svm_ratio)
    print('svm chloroplasts FPR:', svm_fpr_c, '| Approximation %.2f' %svm_fpr_c)
    print('svm mithocondrion FPR:', svm_fpr_m, '| Approximation %.2f' %svm_fpr_m)
    print('svm peroxisomes FPR:', svm_fpr_p, '| Approximation %.2f' %svm_fpr_p)
    print('svm ts FPR:', svm_fpr_ts, '| Approximation %.2f' %svm_fpr_ts)
