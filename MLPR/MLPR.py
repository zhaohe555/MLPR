# coding: utf-8
import networkx as nx
import numpy as np
import os
import math
import pickle
import matplotlib.pyplot as plt
def result_analyse(G, protein_score_sorted, esspro_inppi, cc):
    if len(esspro_inppi) / 1000 <5:
        t = 100
    else:
        t = 1000
    q = t
    while q <= math.ceil(len(esspro_inppi)/100)*100:
        numberof_esspro = 0
        for i in protein_score_sorted[0: q]:
            if i[0] in esspro_inppi:
                numberof_esspro = numberof_esspro + 1
        #print('Top' + str(q) + 'have essential protein', numberof_esspro)
        q = q + t
    # Statistical measures
    TP = 0
    FN = 0
    for i in protein_score_sorted[0: int(len(G.nodes())*cc)]:
        if i[0] in esspro_inppi:
            TP = TP + 1
    FP = int(len(G.nodes())*cc) - TP
    for i in protein_score_sorted[int(len(G.nodes())*cc):]:
        if i[0] in esspro_inppi:
            FN = FN + 1
    TN = len(G.nodes()) - int(len(G.nodes())*cc) - FN
    SN = TP / (TP + FN)
    SP = TN / (TN + FP)
    PPV = TP / (TP + FP)
    NPV = TN / (TN + FN)
    Fmeasure = 2 * SN * PPV / (SN + PPV)
    ACC = (TP + TN) / (TP + TN + FP + FN)
    #print(TP, FP, TN, FN,'\n', "%.4f" % SN,'\n', "%.4f" % SP,'\n', "%.4f" % PPV,'\n', "%.4f" % NPV,'\n', "%.4f" % Fmeasure,'\n', "%.4f" % ACC)
    return ACC, TP

def proscore(G,r,ID_NUM):
    protein_score = {}
    r_pro = {}
    for each in G.nodes():
        r_pro.update({each: r[ID_NUM[each]]})
    for each in G.nodes():
        score = r_pro[each]
        protein_score.update({each: score})
    return protein_score

# load network
def load_network(netpath,netName):
    os.chdir(netpath)
    a = open(netName, "r")
    G = nx.Graph()
    for i in a:
        n = i.strip().split("\t")
        G.add_edge(n[0], n[1])
    a.close()
    return G

def load_subnetwork(edgesub):

    G = nx.Graph()
    for i in edgesub:

        G.add_edge(i[0], i[1])

    return G

# mapping
def mapping(G):
    ID_NUM = {}
    NUM_ID = {}
    numAllGene = 0
    for x in G.nodes():
        ID_NUM[x] = numAllGene
        NUM_ID[numAllGene] = x
        numAllGene = numAllGene + 1
    return ID_NUM, NUM_ID, numAllGene

def readEssentialProtein(G, esspro_path, esspro_filename):
    esspro = []
    esspro_inppi = []
    os.chdir(esspro_path)
    essentialprotein_file = open(esspro_filename)
    for eachline in essentialprotein_file:
        a = eachline.replace('\n', '')
        esspro.append(a)
    essentialprotein_file.close()
    for each in G.nodes():
        if each in esspro:
            esspro_inppi.append(each)
    esspro1 = list(set(G.nodes())&set(esspro))

    return esspro_inppi

# load homology information
def load_homo_info(homoPath, homoName):
    tupleHomoInfo =[]
    os.chdir(homoPath)
    f = open(homoName,'r')
    for line in f.readlines():
        list = line.strip().split('\t')
        tupleHomoInfo.append((list[0],list[1]))
    f.close()
    return(tupleHomoInfo)

def ini_homo(G1, G2, G3):
    homobias1 = {}
    homobias2 = {}
    homobias3 = {}
    for each in G1.nodes():
        homobias1[each] = 1
    for each in G2.nodes():
        homobias2[each] = 1
    for each in G3.nodes():
        homobias3[each] = 1
    return homobias1, homobias2, homobias3

def ini_homo_s(tupleHomoInfo, homobias1, homobias2, G1, G2):
    for each in tupleHomoInfo:
        a = each[0]
        b = each[1]
        if a in G1.nodes() and b in G2.nodes():
            homobias1[a] = homobias1[a] + 1
            homobias2[b] = homobias2[b] + 1

def multi_homo_s(tupleHomoInfo21, tupleHomoInfo31, tupleHomoInfo23, homobias1, homobias2, homobias3, G1, G2, G3):
    ini_homo_s(tupleHomoInfo21, homobias2, homobias1, G2, G1)
    ini_homo_s(tupleHomoInfo31, homobias3, homobias1, G3, G1)
    ini_homo_s(tupleHomoInfo23, homobias2, homobias3, G2, G3)
    maxhomo1 = max(homobias1.values())
    maxhomo2 = max(homobias2.values())
    maxhomo3 = max(homobias3.values())
    for each in homobias1.keys():
        homobias1[each] = homobias1[each] / maxhomo1
    for each in homobias2.keys():
        homobias2[each] = homobias2[each] / maxhomo2
    for each in homobias3.keys():
        homobias3[each] = homobias3[each] / maxhomo3
    nx.set_node_attributes(G1, homobias1, 'weight')
    nx.set_node_attributes(G2, homobias2, 'weight')
    nx.set_node_attributes(G3, homobias3, 'weight')
    biasR1 = []
    biasR2 = []
    biasR3 = []
    for each in homobias1.keys():
        biasR1.append(homobias1[each])
    for each in homobias2.keys():
        biasR2.append(homobias2[each])
    for each in homobias3.keys():
        biasR3.append(homobias3[each])

    return homobias1, biasR1, homobias2, biasR2, homobias3, biasR3

def homo_e(G,homobias):
    for (a, b) in G.edges():
        G[a][b]['homoweight'] = homobias[a]*homobias[b]
    # edge 'complexweight'
    complex_weights = [G[a][b]['homoweight'] for a, b in G.edges() if 'homoweight' in G[a][b]]

    max_homo_weight = max(complex_weights) if complex_weights else None
    for (a, b) in G.edges():
        G[a][b]['homoweight'] = G[a][b]['homoweight']/max_homo_weight

def GO_e(G1, GOname):
    pro_go = {}
    for eachpro in G1.nodes():
        pro_go.update({eachpro: []})
    with open(GOname, 'r') as f:
        for i in f:
            i = i.strip().split('\t')
            if i[0] in G1.nodes():
                pro_go[i[0]].append(i[1])
    for (a, b) in G1.edges():
        go_score = 0
        pro1_go = pro_go[a]
        pro2_go = pro_go[b]
        pro1_pro2_and = list(set(pro1_go) & set(pro2_go))
        if (len(pro1_go) * len(pro2_go)) != 0:
            go_score = (len(pro1_pro2_and) ** 2) / (len(pro1_go) * len(pro2_go))
        G1[a][b]['weight'] = go_score
    return pro_go

def subceller_select(G, subcellername, subceller_11, esspro_inppi):
    threshold_sub = len(esspro_inppi)/len(G.nodes())
    subceller_11_protein = {}
    for each in subceller_11.keys():
        subceller_11_protein.setdefault(each, [])
    subceller_11_essprotein = {}
    for each in subceller_11.keys():
        subceller_11_essprotein.setdefault(each, [])
    with open(subcellername, 'r') as f:
        for i in f:
            i = i.strip().split('\t')
            if i[0] in G.nodes() and i[1] in subceller_11_protein.keys():
                subceller_11_protein[i[1]].append(i[0])
                if i[0] in esspro_inppi:
                    subceller_11_essprotein[i[1]].append(i[0])
    subceller_selected = {}
    for each in subceller_11.keys():
        if len(subceller_11_essprotein[each])/len(subceller_11_protein[each]) >= threshold_sub:
            subceller_selected.setdefault(each, 0)
    return subceller_selected


def subceller_node(G1, subcellername, subceller_11):
    subceller_11_protein = {}
    for each in subceller_11.keys():
        subceller_11_protein.setdefault(each, [])

    pro_subceller_11 = {}
    for each in G1.nodes():
        pro_subceller_11.setdefault(each, [])

    with open(subcellername, 'r') as f:
        for i in f:
            i = i.strip().split('\t')
            if i[0] in G1.nodes() and i[1] in subceller_11_protein.keys():
                subceller_11_protein[i[1]].append(i[0])
                pro_subceller_11[i[0]].append(i[1])
    subceller_11_proteinscore = {key: len(value) for key, value in subceller_11_protein.items() if isinstance(value, list)}
    max_length = max(subceller_11_proteinscore.values())
    subceller_11 = {key: value / max_length for key, value in subceller_11_proteinscore.items()}
    score_subceller1 = {key: sum(subceller_11.get(item, 0) for item in value) for key, value in pro_subceller_11.items()}
    maxvalue = max(score_subceller1.values())
    score_subceller =  {key: value / maxvalue for key, value in score_subceller1.items()}
    score_sublist = []
    for each in score_subceller.keys():
        score_sublist.append(score_subceller[each])

    edge_sub = []
    for each in G1.edges():
        a = each[0]
        b = each[1]

        for sub in subceller_11.keys():
            if sub in pro_subceller_11[a] and sub in pro_subceller_11[b]:
                edge_sub.append([a, b])
                break
    return score_subceller, score_sublist, edge_sub

def complex_node(G, complexflename):
    proteincomplex = []
    with open(complexflename, 'r', encoding='utf-8') as file:
        for eachline in file:
            eachproteincomplex = eachline.split()
            proteincomplex.append(eachproteincomplex)
    protein_frequencycomplex = {}
    for each in G.nodes():
        protein_frequencycomplex.setdefault(each, 0)
    for eachcomplex in proteincomplex:
        for eachprotein in eachcomplex:
            if eachprotein in G.nodes():
                protein_frequencycomplex[eachprotein] = protein_frequencycomplex[eachprotein]+1
    maxvalue = max(protein_frequencycomplex.values())
    protein_frequencycomplex1 = {key: value / maxvalue for key, value in protein_frequencycomplex.items()}
    score_complexlist = []
    for each in protein_frequencycomplex1.keys():
        score_complexlist.append(protein_frequencycomplex1[each])

    for (a, b) in G.edges():
        G[a][b]['complexweight'] = protein_frequencycomplex1[a]*protein_frequencycomplex1[b]
    # edge 'complexweight'
    complex_weights = [G[a][b]['complexweight'] for a, b in G.edges() if 'complexweight' in G[a][b]]

    max_complex_weight = max(complex_weights) if complex_weights else None
    for (a, b) in G.edges():
        G[a][b]['complexweight'] = G[a][b]['complexweight']/max_complex_weight

    return protein_frequencycomplex1, score_complexlist

def eccMatrix(numAllGene1, G1, ID_NUM1, genescore):
    ecc = np.zeros(shape=(numAllGene1, numAllGene1))
    for (a, b) in G1.edges():
        '''
        i = ID_NUM1[a]
        j = ID_NUM1[b]
        geneexpressionscore = genescore[(a,b)]
        ecc[i][j] = G1[a][b]['weight'] * geneexpressionscore
        ecc[j][i] = G1[a][b]['weight'] * geneexpressionscore
        '''
        geneexpressionscore = genescore[(a,b)]
        Neighbori = list(G1.neighbors(a))
        Neighborj = list(G1.neighbors(b))
        coNeighbor = set(Neighbori) & set(Neighborj)
        i = ID_NUM1[a]
        j = ID_NUM1[b]
        ecc[i][j] = 1/numAllGene1
        ecc[j][i] = 1/numAllGene1
        if (len(Neighbori) > 1) & (len(Neighborj) > 1):
            coscore = 0
            for each in coNeighbor:
                coscore = coscore + G1.nodes[each]['weight']
            ecc[i][j] = ecc[i][j] + coscore / min(len(Neighbori)-1, len(Neighborj)-1)
            ecc[j][i] = ecc[i][j]
        ecc[i][j] = ecc[i][j]*(G1[a][b]['weight']+G1[a][b]['complexweight'])
        ecc[j][i] = ecc[j][i]*(G1[a][b]['weight']+G1[a][b]['complexweight'])
        # ecc[i][j] = ecc[i][j] * (G1[a][b]['weight'] + geneexpressionscore)
        # ecc[j][i] = ecc[j][i] * (G1[a][b]['weight'] + geneexpressionscore)
    return ecc

def eccMatrix1(numAllGene1, G1, ID_NUM1):
    ecc = np.zeros(shape=(numAllGene1, numAllGene1))
    for (a, b) in G1.edges():
        Neighbori = list(G1.neighbors(a))
        Neighborj = list(G1.neighbors(b))
        coNeighbor = set(Neighbori) & set(Neighborj)
        i = ID_NUM1[a]
        j = ID_NUM1[b]
        if (len(Neighbori) > 1) & (len(Neighborj) > 1):
            coscore = len(coNeighbor)
            ecc[i][j] = ecc[i][j] + coscore / min(len(Neighbori)-1, len(Neighborj)-1)
            ecc[j][i] = ecc[i][j]
        ecc[i][j] = ecc[i][j]*(G1[a][b]['weight']+G1[a][b]['complexweight'])
        ecc[j][i] = ecc[j][i]*(G1[a][b]['weight']+G1[a][b]['complexweight'])
    return ecc

def genescore_node(G, genescore):
    genescore_node1 = {}
    for each in G.nodes():
        genescore_node1.setdefault(each, 0)
    for (a, b) in G.edges():
        geneexpressionscore = genescore[(a, b)]
        genescore_node1[a] = genescore_node1[a]+geneexpressionscore
        genescore_node1[b] = genescore_node1[b]+geneexpressionscore
    maxvalue = max(genescore_node1.values())
    genescore_node = {key: value / maxvalue for key, value in genescore_node1.items()}
    score_genenodelist = []
    for each in genescore_node.keys():
        score_genenodelist.append(genescore_node[each])
    return genescore_node, score_genenodelist

def jackknifeCurve(cc, G, esspro_inppi, BETAMAX):
    # Jackknife curves
    variable = globals()
    N = int(len(G.nodes()) * cc)
    labels = ['λ=0.1', 'λ=0.2', 'λ=0.3', 'λ=0.4', 'λ=0.5', 'λ=0.6', 'λ=0.7', 'λ=0.8', 'λ=0.9']

    plt.figure(dpi=300)

    i = 0
    for j in range(1, 10):
        beta = j * 0.1
        beta = round(beta, 1)
        x = []
        y = []
        x.append(0)
        y.append(0)
        xnum = 0
        ynum = 0
        for each in variable['score2_sorted' + str(beta)][:N]:
            x.append(xnum + 1)
            if each[0] in esspro_inppi:
                ynum = ynum + 1
            y.append(ynum)
            xnum = xnum + 1
        if str(beta) != str(BETAMAX):
            plt.plot(x, y, label=labels[i])
        else:
            plt.plot(x, y, lw=6, color='black', label=labels[i])
        i = i + 1
    #plt.title('biogrid')
    plt.xlabel('Number of top ranked predicted essential protein')
    plt.ylabel('Number of essential protein')
    plt.grid()
    plt.legend(loc="upper left")
    plt.show()

def Standard_W(w_lr):
    ff = np.sum(w_lr, axis=0)
    w_lr = w_lr / ff
    #print(w_lr[:, 0], w_lr[:, 1])
    return w_lr


def Standard_W1(w_lr):
    ff = np.sum(w_lr, axis=0)
    for i in range(len(ff)):
        if ff[i] != 0:  # Normalize only columns with non-zero column sums
            w_lr[:, i] = w_lr[:, i] / ff[i]
    return w_lr

def initial_Wlayers1(G1, G2, tupleHomoInfo, numAllGene1, numAllGene2, ID_NUM1, ID_NUM2, pro_go1, pro_go2  ):
    W12 = np.zeros(shape=(numAllGene1, numAllGene2))
    W21 = np.zeros(shape=(numAllGene2, numAllGene1))

    for each in tupleHomoInfo:
        a = each[0]
        b = each[1]

        if a in G1.nodes() and b in G2.nodes():

            go_score = 0
            pro1_go = pro_go1[a]
            pro2_go = pro_go2[b]
            pro1_pro2_and = list(set(pro1_go) & set(pro2_go))
            if (len(pro1_go) * len(pro2_go)) != 0:
                go_score = len(pro1_pro2_and)

            wab = G1.nodes[a]['weight'] * G2.nodes[b]['weight']
            W12[ID_NUM1[a]][ID_NUM2[b]] = go_score * wab
            W21[ID_NUM2[b]][ID_NUM1[a]] = go_score * wab
    #print(np.max(W12),np.max(W21))
    W12 = W12/np.max(W12)
    W21 = W21/np.max(W21)


    return W12, W21

def rwr_homo(w_lr1, beta1, biasR1, score_sublist1, score_complexlist1,genescorelist):
    score_sublist1 = np.array(score_sublist1)
    genescorelist = np.array(genescorelist)
    biasR1 = np.array(biasR1)
    #biasR1 = biasR1*score_sublist1
    biasR1 = biasR1*score_sublist1*(score_complexlist1+genescorelist)
    #biasR1 = biasR1 * score_complexlist1#fruitfly
    biasR1 = biasR1/max(biasR1)
    #biasR1 = score_sublist1 * score_complexlist1
    #biasR1 = biasR1 * score_complexlist1
    #biasR1 = np.array(score_complexlist1)
    r1 = biasR1
    l = 0.0000001  #
    change1 = 100000
    t = 0  #
    while (change1 > l):
        r_new1 = (1 - beta1) * np.dot(w_lr1, r1) + beta1 * biasR1
        change1_list = r_new1 - r1
        change1 = max(map(abs, change1_list))
        r1 = r_new1
        t = t + 1
        #print('1The number of iterations is：%s' % str(t), r_new1, change1)
    r1 = list(r1)
    return r1

def output_rwr_homo(NetName1, G1, numAllGene1, ID_NUM1, esspro_inppi1, biasR1, cc1, ecc1, score_sublist1,score_complexlist1, genescorelist):
    w_lr1 = Standard_W1(ecc1)
    print(len(esspro_inppi1))
    print('result ')
    variable = globals()
    ACCMAX = 0
    BETAMAX = 0
    for i in range(1, 10):
        beta1 = i * 0.1
        beta1 = round(beta1, 1)
        #print('beta=', beta1)
        r1 = rwr_homo(w_lr1, beta1, biasR1, score_sublist1, score_complexlist1,genescorelist)
        variable['score2' + str(beta1)] = proscore(G1, r1, ID_NUM1)
        variable['score2_sorted' + str(beta1)] = sorted(variable['score2' + str(beta1)].items(),
                                                        key=lambda x: x[1], reverse=True)
        ACC, TP = result_analyse(G1, variable['score2_sorted' + str(beta1)], esspro_inppi1, cc1)
        if ACC>ACCMAX:
            ACCMAX = ACC
            BETAMAX = beta1
    print(NetName1+'best ACC and beta',ACCMAX,BETAMAX)
    score_list = []
    for each in variable['score2' + str(beta1)].keys():
        score_list.append(variable['score2' + str(beta1)][each])
    jackknifeCurve(cc1, G1, esspro_inppi1, BETAMAX)
    return variable['score2_sorted' + str(BETAMAX)], variable['score2' + str(beta1)], BETAMAX, score_list


def rwr_homo_fruitfly(w_lr1, beta1, biasR1, score_sublist1, score_complexlist1,genescorelist):
    score_sublist1 = np.array(score_sublist1)
    genescorelist = np.array(genescorelist)
    biasR1 = np.array(biasR1)
    #biasR1 = biasR1*score_sublist1
    #biasR1 = biasR1*score_sublist1*(score_complexlist1+genescorelist)
    biasR1 = biasR1 * score_complexlist1#fruitfly
    biasR1 = biasR1/max(biasR1)
    #biasR1 = score_sublist1 * score_complexlist1
    #biasR1 = biasR1 * score_complexlist1
    #biasR1 = np.array(score_complexlist1)
    r1 = biasR1
    l = 0.0000001
    change1 = 100000
    t = 0
    while (change1 > l):
        r_new1 = (1 - beta1) * np.dot(w_lr1, r1) + beta1 * biasR1
        change1_list = r_new1 - r1
        change1 = max(map(abs, change1_list))
        r1 = r_new1
        t = t + 1
        #print('1The number of iterations is：%s' % str(t), r_new1, change1)
    r1 = list(r1)
    return r1

def output_rwr_homo_fruitfly(NetName1, G1, numAllGene1, ID_NUM1, esspro_inppi1, biasR1, cc1, ecc1, score_sublist1,score_complexlist1, genescorelist):
    w_lr1 = Standard_W1(ecc1)
    print(len(esspro_inppi1))
    print('result')
    variable = globals()
    ACCMAX = 0
    BETAMAX = 0
    for i in range(1, 10):
        beta1 = i * 0.1
        beta1 = round(beta1, 1)
        #print('beta=', beta1)
        r1 = rwr_homo_fruitfly(w_lr1, beta1, biasR1, score_sublist1, score_complexlist1,genescorelist)
        variable['score2' + str(beta1)] = proscore(G1, r1, ID_NUM1)
        variable['score2_sorted' + str(beta1)] = sorted(variable['score2' + str(beta1)].items(),
                                                        key=lambda x: x[1], reverse=True)
        ACC, TP = result_analyse(G1, variable['score2_sorted' + str(beta1)], esspro_inppi1, cc1)
        if ACC>ACCMAX:
            ACCMAX = ACC
            BETAMAX = beta1
    print(NetName1+'best ACC and beta',ACCMAX,BETAMAX)
    score_list = []
    for each in variable['score2' + str(beta1)].keys():
        score_list.append(variable['score2' + str(beta1)][each])
    jackknifeCurve(cc1, G1, esspro_inppi1, BETAMAX)
    return variable['score2_sorted' + str(BETAMAX)], variable['score2' + str(beta1)], BETAMAX, score_list


def rwr_homo_123(beta1, w_lr1, biasR1, score_sublist1, score_complexlist1, beta11, biasR2, biasR3, w12, w13, score_sublist2, score_sublist3, score_complexlist2, score_complexlist3, scoregenelist1,scoregenelist2,scoregenelist3):
    biasR2 = np.array(biasR2) * (np.array(score_complexlist2)+scoregenelist2)
    #biasR2 =  np.array(biasR2) * np.array(score_sublist2) * np.array(score_complexlist2)
    biasR3 = np.array(biasR3) * np.array(score_sublist3) * (np.array(score_complexlist3)+np.array(scoregenelist3))
    bi1 = beta11*np.dot(w12, biasR2) + (1-beta11)*np.dot(w13, biasR3)
    #bi1 = bi1 / max(bi1)
    score_sublist1 = np.array(score_sublist1)
    biasR1 = np.array(biasR1)
    biasR1 = biasR1*score_sublist1 * (np.array(score_complexlist1)+np.array(scoregenelist1))

    r1 = biasR1
    l = 0.0000001
    change1 = 100000
    t = 0

    while (change1 > l):
        r_new1 = (1 - 0.4) * (beta1*np.dot(w_lr1, r1)+(1-beta1)*bi1) + 0.4 * biasR1
        #r_new1 = (1 - beta1) * np.dot(w_lr1, r1) + beta1 * biasR1

        change1_list = r_new1 - r1
        change1 = max(map(abs, change1_list))
        r1 = r_new1
        t = t + 1
        #print('1The number of iterations is：%s' % str(t), r_new1, change1)
    r1 = list(r1)
    return r1

def output_rwr_homo123(NetName1, G1, ID_NUM1, esspro_inppi1, biasR1, cc1, ecc1, score_sublist1, score_complexlist1,  biasR2, biasR3, w12, w13, score_sublist2, score_sublist3, score_complexlist2, score_complexlist3, scoregenelist1,scoregenelist2,scoregenelist3):
    w_lr1 = Standard_W1(ecc1)
    # W12 = Standard_W(w12)
    # W13 = Standard_W(w13)
    W12 = Standard_W1(w12)
    W13 = Standard_W1(w13)
    print(len(esspro_inppi1))
    print('result')
    variable = globals()
    ACCMAX = 0
    BETAMAX1 = 0
    BETAMAX11 = 0
    ACC1 = []
    for i in range(1, 10):
        beta1 = i * 0.1
        beta1 = round(beta1, 1)
        acc = []
        for j in range(1, 10):
            beta11 = j * 0.1
            beta11 = round(beta11, 1)
            #print('beta1=', beta1, 'beta11=', beta11)
            r1 = rwr_homo_123(beta1, w_lr1, biasR1, score_sublist1, score_complexlist1, beta11, biasR2, biasR3, W12, W13, score_sublist2, score_sublist3, score_complexlist2, score_complexlist3, scoregenelist1,scoregenelist2,scoregenelist3)
            variable['score2' + str(beta1)+str(beta11)] = proscore(G1, r1, ID_NUM1)
            variable['score2_sorted' + str(beta1)+str(beta11)] = sorted(variable['score2' + str(beta1)+str(beta11)].items(),
                                                            key=lambda x: x[1], reverse=True)
            ACC, TP = result_analyse(G1, variable['score2_sorted' + str(beta1)+str(beta11)], esspro_inppi1, cc1)
            acc.append(ACC)
            if ACC>ACCMAX:
                ACCMAX = ACC
                BETAMAX1 = beta1
                BETAMAX11 = beta11
        ACC1.append(acc)
    print(NetName1+'best ACC and beta',ACCMAX,BETAMAX1,BETAMAX11)
    return variable['score2_sorted' + str(BETAMAX1)+str(BETAMAX11)], variable['score2' + str(beta1)+str(beta11)], BETAMAX1, BETAMAX11, ACC1

def rwr_homo_213(beta1, w_lr2, biasR2, score_sublist2, score_complexlist2, beta11, biasR1,
                              biasR3, w21, w23, score_sublist1, score_sublist3, score_complexlist1,
                              score_complexlist3,score_genenodelist1,score_genenodelist2,score_genenodelist3):
    biasR1 = np.array(biasR1) * np.array(score_sublist1) * (np.array(score_complexlist1)+score_genenodelist1)
    biasR3 = np.array(biasR3) * np.array(score_sublist3) * (np.array(score_complexlist3)+score_genenodelist3)
    biasR2 = np.array(biasR2) * np.array(score_complexlist2)
    # biasR2 =  np.array(biasR2) * np.array(score_sublist2) * np.array(score_complexlist2)

    bi2 = beta11 * np.dot(w21, biasR1) + (1 - beta11) * np.dot(w23, biasR3)
    # bi2 = bi2 / max(bi2)

    r2 = biasR2
    l = 0.0000001
    change1 = 100000
    t = 0

    while (change1 > l):
        r_new2 = (1 - 0.9) * (beta1 * np.dot(w_lr2, r2) + (1 - beta1) * bi2) + 0.9 * biasR2
        # r_new1 = (1 - beta1) * np.dot(w_lr1, r1) + beta1 * biasR1

        change1_list = r_new2 - r2
        change1 = max(map(abs, change1_list))
        r2 = r_new2
        t = t + 1
        # print('1The number of iterations is：%s' % str(t), r_new1, change1)
    r2 = list(r2)
    return r2

def output_rwr_homo213(NetName2, G2, ID_NUM2, esspro_inppi2,biasR2, cc2, ecc2, score_sublist2,
                       score_complexlist2,biasR1, biasR3, w21, w23,score_sublist1, score_sublist3,
                       score_complexlist1, score_complexlist3,score_genenodelist1,score_genenodelist2,score_genenodelist3):
    w_lr2 = Standard_W1(ecc2)
    # W12 = Standard_W(w12)
    # W13 = Standard_W(w13)
    W21 = Standard_W1(w21)
    W23 = Standard_W1(w23)
    print(len(esspro_inppi2))
    print('result')
    variable = globals()
    ACCMAX = 0
    BETAMAX1 = 0
    BETAMAX11 = 0
    ACC1 = []
    for i in range(1, 10):
        beta1 = i * 0.1
        beta1 = round(beta1, 1)
        acc = []
        for j in range(1, 10):
            beta11 = j * 0.1
            beta11 = round(beta11, 1)
            #print('beta1=', beta1, 'beta11=', beta11)
            r2 = rwr_homo_213(beta1, w_lr2, biasR2, score_sublist2, score_complexlist2, beta11, biasR1,
                              biasR3, W21, W23, score_sublist1, score_sublist3, score_complexlist1,
                              score_complexlist3,score_genenodelist1,score_genenodelist2,score_genenodelist3)
            variable['score2' + str(beta1) + str(beta11)] = proscore(G2, r2, ID_NUM2)
            variable['score2_sorted' + str(beta1) + str(beta11)] = sorted(
                variable['score2' + str(beta1) + str(beta11)].items(),
                key=lambda x: x[1], reverse=True)
            ACC, TP = result_analyse(G2, variable['score2_sorted' + str(beta1) + str(beta11)], esspro_inppi2, cc2)
            acc.append(ACC)
            if ACC > ACCMAX:
                ACCMAX = ACC
                BETAMAX1 = beta1
                BETAMAX11 = beta11
        ACC1.append(acc)
    print(NetName2 + 'best ACC and beta', ACCMAX, BETAMAX1, BETAMAX11)
    return variable['score2_sorted' + str(BETAMAX1) + str(BETAMAX11)], variable[
            'score2' + str(beta1) + str(beta11)], BETAMAX1, BETAMAX11,ACC1

def rwr_homo_312(beta1, w_lr3, biasR3, score_sublist3, score_complexlist3, beta11, biasR1,
                              biasR2, w31, w32, score_sublist1, score_sublist2, score_complexlist1,
                              score_complexlist2,score_genenodelist1,score_genenodelist2,score_genenodelist3):
    biasR1 = np.array(biasR1) * np.array(score_sublist1) * (np.array(score_complexlist1)+score_genenodelist1)
    biasR3 = np.array(biasR3) * np.array(score_sublist3) * (np.array(score_complexlist3)+score_genenodelist3)
    biasR2 = np.array(biasR2) * (np.array(score_complexlist2)+score_genenodelist2)
    # biasR2 =  np.array(biasR2) * np.array(score_sublist2) * np.array(score_complexlist2)

    bi3 = beta11 * np.dot(w31, biasR1) + (1 - beta11) * np.dot(w32, biasR2)
    # bi3 = bi3 / max(bi3)

    r3 = biasR3
    l = 0.0000001
    change1 = 100000
    t = 0

    while (change1 > l):
        r_new3 = (1 - 0.9) * (beta1 * np.dot(w_lr3, r3) + (1 - beta1) * bi3) + 0.9 * biasR3
        change1_list = r_new3 - r3
        change1 = max(map(abs, change1_list))
        r3 = r_new3
        t = t + 1
        # print('1The number of iterations is：%s' % str(t), r_new1, change1)
    r3 = list(r3)
    return r3

def output_rwr_homo312(NetName3, G3, ID_NUM3, esspro_inppi3,biasR3, cc3, ecc3, score_sublist3,
                       score_complexlist3,biasR1, biasR2, w31, w32,score_sublist1, score_sublist2,
                       score_complexlist1, score_complexlist2,score_genenodelist1,score_genenodelist2,score_genenodelist3):
    w_lr3 = Standard_W1(ecc3)
    # W12 = Standard_W(w12)
    # W13 = Standard_W(w13)
    W31 = Standard_W1(w31)
    W32 = Standard_W1(w32)
    print(len(esspro_inppi3))
    print('result')
    variable = globals()
    ACCMAX = 0
    BETAMAX1 = 0
    BETAMAX11 = 0
    ACC1 = []
    for i in range(1, 10):
        beta1 = i * 0.1
        beta1 = round(beta1, 1)
        acc = []
        for j in range(1, 10):
            beta11 = j * 0.1
            beta11 = round(beta11, 1)
            #print('beta1=', beta1, 'beta11=', beta11)
            r3 = rwr_homo_312(beta1, w_lr3, biasR3, score_sublist3, score_complexlist3, beta11, biasR1,
                              biasR2, W31, W32, score_sublist1, score_sublist2, score_complexlist1,
                              score_complexlist2,score_genenodelist1,score_genenodelist2,score_genenodelist3)
            variable['score3' + str(beta1) + str(beta11)] = proscore(G3, r3, ID_NUM3)
            variable['score3_sorted' + str(beta1) + str(beta11)] = sorted(
                variable['score3' + str(beta1) + str(beta11)].items(),
                key=lambda x: x[1], reverse=True)
            ACC, TP = result_analyse(G3, variable['score3_sorted' + str(beta1) + str(beta11)], esspro_inppi3, cc3)
            acc.append(ACC)
            if ACC > ACCMAX:
                ACCMAX = ACC
                BETAMAX1 = beta1
                BETAMAX11 = beta11
        ACC1.append(acc)
    print(NetName3 + 'best ACC and beta', ACCMAX, BETAMAX1, BETAMAX11)
    return variable['score3_sorted' + str(BETAMAX1) + str(BETAMAX11)], variable[
        'score3' + str(beta1) + str(beta11)], BETAMAX1, BETAMAX11, ACC1


if __name__ == '__main__':
    netPath = r'D:\code\MLPR'
    NetName1 = 'DIP.txt'
    NetName2 = 'ppi_fruitfly.txt'
    NetName3 = 'humanppi.txt'
    homoName21 = 'fruitfly-yeast.txt'
    homoName31 = 'human-yeast.txt'
    homoName23 = 'fruitfly-human.txt'
    esspro_filename1 = 'essentialprotein.txt'
    esspro_filename2 = 'essentialprotein_fruitfly.txt'
    esspro_filename3 = 'esspro_human.txt'

    G1 = load_network(netPath, NetName1)
    ID_NUM1, NUM_ID1, numAllGene1 = mapping(G1)
    G2 = load_network(netPath, NetName2)
    ID_NUM2, NUM_ID2, numAllGene2 = mapping(G2)
    G3 = load_network(netPath, NetName3)
    ID_NUM3, NUM_ID3, numAllGene3 = mapping(G3)
    esspro_inppi1 = readEssentialProtein(G1, netPath, esspro_filename1)
    esspro_inppi2 = readEssentialProtein(G2, netPath, esspro_filename2)
    esspro_inppi3 = readEssentialProtein(G3, netPath, esspro_filename3)

    tupleHomoInfo21 = load_homo_info(netPath, homoName21)
    tupleHomoInfo31 = load_homo_info(netPath, homoName31)
    tupleHomoInfo23 = load_homo_info(netPath, homoName23)
    homobias1, homobias2, homobias3 = ini_homo(G1, G2, G3)
    homobias1, biasR1, homobias2, biasR2, homobias3, biasR3 = multi_homo_s(tupleHomoInfo21, tupleHomoInfo31,
                                                                           tupleHomoInfo23, homobias1, homobias2,
                                                                           homobias3, G1, G2, G3)
    homo_e(G1, homobias1)
    homo_e(G2, homobias2)
    homo_e(G3, homobias3)

    GOname1 = 'go_yeast.tsv'
    pro_go1 = GO_e(G1, GOname1)
    GOname2 = 'go_fruitfly.tsv'
    pro_go2 = GO_e(G2, GOname2)
    GOname3 = 'go_human.tsv'
    pro_go3 = GO_e(G3, GOname3)

    with open('humangenescore.pkl', 'rb') as f:
        loaded_dict_human = pickle.load(f)
    with open('fruitflygenescore.pkl', 'rb') as f:
        loaded_dict_fruitfly = pickle.load(f)
    with open('yeastgenescore.pkl', 'rb') as f:
        loaded_dict_yeast = pickle.load(f)

    cc1 = 0.25
    cc2 = 0.1
    cc3 = 0.35

    complexflename1 = 'complex_yeast.txt'
    complexscore_node1, score_complexlist1 = complex_node(G1, complexflename1)
    complexflename2 = 'complex_fruitfly.txt'
    complexscore_node2, score_complexlist2 = complex_node(G2, complexflename2)
    complexflename3 = 'complex_human.txt'
    complexscore_node3, score_complexlist3 = complex_node(G3, complexflename3)

    ecc1 = eccMatrix(numAllGene1, G1, ID_NUM1, loaded_dict_yeast)
    ecc2 = eccMatrix(numAllGene2, G2, ID_NUM2, loaded_dict_fruitfly)
    ecc3 = eccMatrix(numAllGene3, G3, ID_NUM3, loaded_dict_human)


    #subceller_11 = {'Nucleus': 0, 'Cytosol': 0, 'Cytoskeleton': 0, 'Endoplasmic reticulum': 0, 'Golgi apparatus': 0}
    subceller_11 = {'Nucleus': 0, 'Cytosol': 0, 'Cytoskeleton': 0, 'Peroxisome': 0,
                    'Vacuole': 0, 'Endoplasmic reticulum': 0, 'Golgi apparatus': 0,
                    'Plasma membrane': 0, 'Endosome': 0, 'Extracellular region': 0,
                    'Mitochondrion': 0}

    subcellername1 = 'subceller_yeast.tsv'
    subceller_selected1 = subceller_select(G1, subcellername1, subceller_11, esspro_inppi1)
    score_subceller1, score_sublist1, edge_sub1 = subceller_node(G1, subcellername1, subceller_selected1)

    subcellername2 = 'subceller_fruitfly.tsv'
    subceller_selected2 = subceller_select(G2, subcellername2, subceller_11, esspro_inppi2)
    score_subceller2, score_sublist2, edge_sub2 = subceller_node(G2, subcellername2, subceller_selected2)

    subcellername3 = 'subceller_human.tsv'
    subceller_selected3 = subceller_select(G3, subcellername3, subceller_11, esspro_inppi3)
    score_subceller3, score_sublist3, edge_sub3 = subceller_node(G3, subcellername3, subceller_selected3)

    score_genenode1, score_genenodelist1 = genescore_node(G1, loaded_dict_yeast)
    score_genenode2, score_genenodelist2 = genescore_node(G2, loaded_dict_fruitfly)
    score_genenode3, score_genenodelist3 = genescore_node(G3, loaded_dict_human)

    w21, w12 = initial_Wlayers1(G2, G1, tupleHomoInfo21, numAllGene2, numAllGene1, ID_NUM2, ID_NUM1, pro_go2, pro_go1)
    w31, w13 = initial_Wlayers1(G3, G1, tupleHomoInfo31, numAllGene3, numAllGene1, ID_NUM3, ID_NUM1, pro_go3, pro_go1)
    w23, w32 = initial_Wlayers1(G2, G3, tupleHomoInfo23, numAllGene2, numAllGene3, ID_NUM2, ID_NUM3, pro_go2, pro_go3)

    score1_sorted1, score_1, BETAMAX1, score_list1 = output_rwr_homo(NetName1, G1, numAllGene1, ID_NUM1, esspro_inppi1,
                                                                     biasR1, cc1, ecc1, score_sublist1,
                                                                     score_complexlist1, score_genenodelist1)
    score1_sorted1_multi, score1_multi, beta1max, BETAMAX11, ACC1 = output_rwr_homo123(NetName1, G1,
                    ID_NUM1,esspro_inppi1, biasR1, cc1, ecc1, score_sublist1, score_complexlist1, biasR2,
                    biasR3,w12, w13, score_sublist2, score_sublist3, score_complexlist2, score_complexlist3,
                    score_genenodelist1,score_genenodelist2,score_genenodelist3)
    with open('./result_experiment/yeast.pkl', 'wb') as f:
        pickle.dump(score1_sorted1_multi, f)

    score2_sorted2, score_2, BETAMAX2, score_list2 = output_rwr_homo_fruitfly(NetName2, G2, numAllGene2, ID_NUM2, esspro_inppi2,
                                                                     biasR2, cc2, ecc2, score_sublist2,
                                                                     score_complexlist2, score_genenodelist2)

    score2_sorted2_multi, score2_multi, beta2max, BETAMAX12, ACC2 = output_rwr_homo213(NetName2, G2,
                     ID_NUM2,esspro_inppi2,biasR2, cc2, ecc2,score_sublist2,score_complexlist2,biasR1,
                     biasR3, w21, w23,score_sublist1, score_sublist3,score_complexlist1,
                     score_complexlist3,score_genenodelist1,score_genenodelist2,score_genenodelist3)
    with open('./result_experiment/fruitfly.pkl', 'wb') as f:
        pickle.dump(score2_sorted2_multi, f)

    score3_sorted3, score_3, BETAMAX3, score_list3 = output_rwr_homo(NetName3, G3, numAllGene3, ID_NUM3, esspro_inppi3,
                                                                     biasR3, cc3, ecc3, score_sublist3,
                                                                     score_complexlist3, score_genenodelist3)
    score3_sorted3_multi, score3_multi, beta3max, BETAMAX13, ACC3 = output_rwr_homo312(NetName3, G3,
                     ID_NUM3,esspro_inppi3,biasR3, cc3, ecc3,score_sublist3,score_complexlist3,biasR1,
                     biasR2, w31, w32,score_sublist1, score_sublist2,score_complexlist1,score_complexlist2,
                     score_genenodelist1,score_genenodelist2,score_genenodelist3)
    with open('./result_experiment/human.pkl', 'wb') as f:
        pickle.dump(score3_sorted3_multi, f)

    print()