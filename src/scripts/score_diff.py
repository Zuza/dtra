def create_cftr():
    cftr = [[0 for x in xrange(4)] for x in xrange(4)] 
    cftr[0][0] = 67
    cftr[0][1] = -96
    cftr[0][2] = -20
    cftr[0][3] = -117
    cftr[1][0] = -96
    cftr[1][1] = 100
    cftr[1][2] = -79
    cftr[1][3] = -20
    cftr[2][0] = -20
    cftr[2][1] = -79
    cftr[2][2] = 100
    cftr[2][3] = -96
    cftr[3][0] = -117
    cftr[3][1] = -20
    cftr[3][2] = -96
    cftr[3][3] = 67
    return cftr

def create_hoxd():
    hoxd = [[0 for x in xrange(4)] for x in xrange(4)] 
    hoxd[0][0] = 91
    hoxd[0][1] = -114
    hoxd[0][2] = -31
    hoxd[0][3] = -123
    hoxd[1][0] = -114
    hoxd[1][1] = 100
    hoxd[1][2] = -125
    hoxd[1][3] = -31
    hoxd[2][0] = -31
    hoxd[2][1] = -125
    hoxd[2][2] = 100
    hoxd[2][3] = -114
    hoxd[3][0] = -123
    hoxd[3][1] = -31
    hoxd[3][2] = -114
    hoxd[3][3] = 91
    return hoxd

def create_hum16pter():
    hum16pter = [[0 for x in xrange(4)] for x in xrange(4)] 
    hum16pter[0][0] = 100
    hum16pter[0][1] = -123
    hum16pter[0][2] = -28
    hum16pter[0][3] = -109
    hum16pter[1][0] = -123
    hum16pter[1][1] = 91
    hum16pter[1][2] = -140
    hum16pter[1][3] = -28
    hum16pter[2][0] = -28
    hum16pter[2][1] = -140
    hum16pter[2][2] = 91
    hum16pter[2][3] = -123
    hum16pter[3][0] = -109
    hum16pter[3][1] = -28
    hum16pter[3][2] = -123
    hum16pter[3][3] = 100
    return hum16pter

def get_score(similarity_matrix, a, b):
    actg = "ACGT"
    return similarity_matrix[actg.find(a)][actg.find(b)]

def ideal_score(similarity_matrix, prob):
    return prob["A"] * get_score(similarity_matrix, "A", "A") +\
        prob["C"] * get_score(similarity_matrix, "C", "C")+\
        prob["G"] * get_score(similarity_matrix, "G", "G")+\
        prob["T"] * get_score(similarity_matrix, "T", "T")

def error_score(similarity_matrix, prob, error_percentage):
    s_same = 0 
    for i in ["A", "C", "G", "T"]:
        s_same += prob[i] * get_score(similarity_matrix, i, i)

    s_diff = 0
    prob_norm = 0
    for i in ["A", "C", "G", "T"]:
        for j in ["A", "C", "G", "T"]:
            if i != j:
                s_diff += prob[i]*prob[j]*\
                    get_score(similarity_matrix, i,j)
                prob_norm += prob[i]*prob[j]
    s_diff /= prob_norm

    return (1-error_percentage)*s_same + error_percentage*s_diff

def main():
    prob = {}
    prob["A"] = 0.3
    prob["C"] = 0.2
    prob["G"] = 0.2
    prob["T"] = 0.3
    error_percentage = 0.20

    print "error fraction %0.2f" % (error_percentage)

    print "CFTR"
    cftr = create_cftr()
    # najveci score ako je greska besplatna
    cftr_ideal_score = ideal_score(cftr, prob) * (1-error_percentage) 
    cftr_error_score = error_score(cftr, prob, error_percentage)
    cftr_diff = ((float)(cftr_ideal_score - cftr_error_score)) / \
        cftr_ideal_score
    print "ideal=%0.2f" % (cftr_ideal_score)
    print "error=%0.2f" % (cftr_error_score)
    print "relative=%0.2f" % (cftr_diff)

    print "HOXD"
    hoxd = create_hoxd()
    hoxd_ideal_score = ideal_score(hoxd, prob) * (1-error_percentage)
    hoxd_error_score = error_score(hoxd, prob, error_percentage)
    hoxd_diff = ((float)(hoxd_ideal_score - hoxd_error_score)) / \
        hoxd_ideal_score
    print "ideal=%0.2f" % (hoxd_ideal_score)
    print "error=%0.2f" % (hoxd_error_score)
    print "relative=%0.2f" % (hoxd_diff)

    print "hum16pter"
    hum16pter = create_hum16pter()
    hum16pter_ideal_score = ideal_score(hum16pter, prob) * (1-error_percentage)
    hum16pter_error_score = error_score(hum16pter, prob, error_percentage)
    hum16pter_diff = ((float)(hum16pter_ideal_score -\
                                  hum16pter_error_score)) / \
                                  hum16pter_ideal_score
    print "ideal=%0.2f" % (hum16pter_ideal_score)
    print "error=%0.2f" % (hum16pter_error_score)
    print "relative=%0.2f" % (hum16pter_diff)


        
if __name__ == "__main__":
    main()
