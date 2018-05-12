#!/usr/bin/python
#title           :SeqCluster.py
#description     :Implementation of Algorithm RaTa
#author          :Tanjil Ahmed
#date            :05/11/2018
#version         :1.0
#usage           :python SeqCluster.py <filename.fasta>
#==============================================================================

import sys
import numpy as np
import math
import pandas as pd

#################################BEGIN##########################################
#check if the input DNA sequences contain characters from {A,C,T,G,R,Y,M}
def file_format_check(fasta):
        string="ACTGRYM"
        for line_number,line in enumerate(fasta):
            if line_number%3==1:
                sequence=line.rstrip()
                if any(c not in string for c in sequence):
                    return 0
        return 1
#################################END############################################

#################################BEGIN##########################################
#calculate weight of the sequence based on frequency of the each base type
def seq_weight(sequence,w_A,w_C,w_T,w_G):
    gweight=0
    for base in sequence:
        if base=='A':
            gweight+=w_A
        if base=='C':
            gweight+=w_C
        if base=='T':
            gweight+=w_T
        if base=='G':
            gweight+=w_G
    return gweight
#################################END############################################


#################################BEGIN##########################################
#calculate weighted length of the cluster
def calculate_wl(cluster_gw,cluster_gl,c,wln):
    wt=0
    for j in range (0,wln+1):
        wt+=cluster_gw[c][j]*cluster_gl[c][j]
    return wt
#################################END############################################


#################################BEGIN##########################################
#Measure similarity between input sequence  and reference gene of the cluster
def sim(temp_w,temp_l,lref_cl,wref_cl):
    return abs(wref_cl*lref_cl-temp_w*temp_l)
#################################END############################################


#################################BEGIN##########################################
#Find the position of gene with closest gene; find the position of the sequence
#in a cluster whose length  is closest to the average sequence length of the
#cluster
def closest(list, Number):
    aux = []
    for gsq in list:
        gsq_l=len(gsq)
        aux.append(abs(Number-gsq_l))
    return aux.index(min(aux))
#################################END############################################



if __name__=="__main__":
    if len(sys.argv) == 1:
        f=open("nd2.fasta","r")
    else:
        f = open(sys.argv[1],"r")
    fasta = f.readlines()
    total_input_seq=0
#####################Start:File Format Check & Base Count#######################
    if file_format_check(fasta)==True:
        print("File contains appropriate character")
        counter ={'A':0, 'C':0,'T':0,'G':0}
        for line_number,line in enumerate(fasta):
            if line_number%3==0:
                fname=line.strip('>')
                fname=fname.strip('\n')
            if line_number%3==1:
                row=0
                sequence=line.rstrip()
                sequence=sequence.replace("R","A")
                sequence=sequence.replace("Y","T")
                sequence=sequence.replace("M","C")
                total_input_seq+=1
                for base in sequence:
                    counter[base]+=1
        #print(counter)
        N=counter['A']+counter['C']+counter['T']+counter['G']
        #print(N)
#######################END:File Format Check & Base Count#######################


#######################Start: Initialization, Algorithm Step-1-2 ###############
        w_A=counter['A']/N
        w_C=counter['C']/N
        w_T=counter['T']/N
        w_G=counter['G']/N

        print("Total Input sequence",total_input_seq)


        c=0
        n=0
        wln=0#index for seq weight and length
        seq_n=0
        w, h = total_input_seq, total_input_seq; # Max possible cluster for given file_format_check
        ref=[]
        cluster=[[] for i in range(w)]
        cluster_gw=[[] for i in range(w)]
        cluster_gl=[[] for i in range(w)]
        cluster_wl=[]
        cluster_th=[]
        w_cmax=[]
        l_cmax=[]
        cur_size=[]
        avg_length_c=[]
        thresh_p=0.55 # Important Parameter; Change in value result different output


        for line_number,line in enumerate(fasta):
            if line_number%3==0:
                fname=line.strip('>')
                fname=fname.strip('\n')
            if line_number%3==1:
                row=0
                sequence=line.rstrip()
                sequence=sequence.replace("R","A")
                sequence=sequence.replace("Y","T")
                sequence=sequence.replace("M","C")
                seq_n=seq_n+1
                L_seq=len(sequence)

                assigned_cluster=200
                min_dism=math.inf
#######################End: Initialization, Algorithm Step-1-2 #################

#######################Start: Algorithm Step-3-4 ###############################

                if seq_n==1:
                    ref.append(sequence)
                    cluster[c].append(sequence) #sequence cluster
                    cluster_gw[c].append(seq_weight(sequence,w_A,w_C,w_T,w_G)) # gene weight
                    cluster_gl[c].append(len(sequence)) # gene length
                    post=cluster_gl[c].index(max(cluster_gl[c]))
                    #print("post",post)
                    w_cmax.append(cluster_gw[c][post])
                    l_cmax.append(max(cluster_gl[c]))
                    #threshold
                    #print(w_cmax[c])
                    #print(l_cmax[c])
                    cluster_th.append(w_cmax[c]*l_cmax[c]*0.55)
                    #print(w_cmax[c]*l_cmax[c])
                    cur_size.append(np.count_nonzero(cluster_gl[c]))

                    cluster_wl.append(calculate_wl(cluster_gw,cluster_gl,c,wln))
                    total_cl_seq=len(cluster[c])
                    total_cl_seq_length=sum(cluster_gl[c])
                    avg_length_c.append(total_cl_seq_length/total_cl_seq)

#######################End: Algorithm Step-3-4 #################################
                else:
#######################Start: Algorithm Step-5.1 ###############################
                    for cl in range(0,c+1):
                        temp_w=seq_weight(sequence,w_A,w_C,w_T,w_G)
                        temp_l=len(sequence)
                        lref_cl=len(ref[cl])
                        wref_cl=seq_weight(ref[cl],w_A,w_C,w_T,w_G)
                        sim_measure=sim(temp_w,temp_l,lref_cl,wref_cl)
                        if sim_measure<cluster_th[cl] and sim_measure<min_dism:
                            min_dism=sim_measure
                            assigned_cluster=cl
########################End: Algorithm Step-5.1 ################################

########################Start: Algorithm Step-5.3 ##############################
                    if assigned_cluster==200:#new cluster create
                        c=c+1
                        wln=0########################check
                        ref.append(sequence)
                        cluster[c].append(sequence) #sequence cluster
                        cluster_gw[c].append(seq_weight(sequence,w_A,w_C,w_T,w_G)) # gene weight
                        cluster_gl[c].append(len(sequence)) # gene length
                        post=cluster_gl[c].index(max(cluster_gl[c]))
                        #print("post",post)
                        w_cmax.append(cluster_gw[c][post])
                        l_cmax.append(max(cluster_gl[c]))
                        cluster_th.append(w_cmax[c]*l_cmax[c]*thresh_p)
                        cur_size.append(np.count_nonzero(cluster_gl[c]))
                        cluster_wl.append(calculate_wl(cluster_gw,cluster_gl,c,wln))
                        total_cl_seq=len(cluster[c])
                        total_cl_seq_length=sum(cluster_gl[c])
                        avg_length_c.append(total_cl_seq_length/total_cl_seq)

########################End: Algorithm Step-5.3 ################################

                    else:

########################Start: Algorithm Step-5.2 ##############################

                        cluster[assigned_cluster].append(sequence) #sequence cluster
                        temp_w=seq_weight(sequence,w_A,w_C,w_T,w_G)
                        lg=len(sequence)
                        temp1=temp_w * lg
                        temp2=cluster_wl[assigned_cluster]+temp1
                        cur_size[assigned_cluster]=cur_size[assigned_cluster] +1
                        avg_length_c[assigned_cluster]=temp2/cur_size[assigned_cluster]
                        cluster_wl[assigned_cluster]=temp2
                        pos=closest(cluster[assigned_cluster],avg_length_c[assigned_cluster])
                        ref[assigned_cluster]=cluster[assigned_cluster][pos]

                        cluster_gw[assigned_cluster].append(seq_weight(sequence,w_A,w_C,w_T,w_G)) # gene weight
                        cluster_gl[assigned_cluster].append(len(sequence)) # gene length
                        post=cluster_gl[assigned_cluster].index(max(cluster_gl[assigned_cluster]))
                        w_cmax[assigned_cluster]=cluster_gw[assigned_cluster][post]
                        l_cmax[assigned_cluster]=max(cluster_gl[assigned_cluster])
                        cluster_th[assigned_cluster]=w_cmax[assigned_cluster]*l_cmax[assigned_cluster]*0.55

########################End: Algorithm Step-5.2 ################################

####################### Write output File ######################################
        df=pd.DataFrame()
        df[0]=[thresh_p,total_input_seq,c+1]
        df=df.T
        df.columns=['Cluster Threshold','Total Input sequence','Number of Generated Cluster']
        #fname="Threshold"+thresh_p+".txt"
        df.to_csv("0.55-10.txt",index=False, sep=' ', header=True) # change filename according to threshold value and number of input sequences


        #print(w_cmax)
        #print(l_cmax)
        for t in range(0,c+1):
            print("Cluster:",t+1,len(cluster[t]))
