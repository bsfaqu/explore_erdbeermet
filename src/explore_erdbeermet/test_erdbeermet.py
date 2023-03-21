from erdbeermet_pkg import simulation as sim
from erdbeermet_pkg import recognition as rec
from contextlib import redirect_stdout
import numpy as np
from decimal import Decimal as d
import math
import subprocess
from sys import argv
import sys
from itertools import permutations
import random
from solve import check_candidate,check_simple_6,rank_candidates, rev_a_checking_new, compare_candidates, test_deltas_new, visualize_splitstree, get_boxes
import matplotlib.pyplot as plt


def rev_a_val(D,a,x,y,z,u,v):
    return((-a*D[x,z]+a*D[z,y]+a*D[x,y]+D[u,z]-D[z,y]-D[u,y])/(2*a*D[x,y]+D[u,x]-D[u,y]-D[x,y]))

def init_scenario(filepath):
    scenario=sim.load(argv[1])
    V=[i for i in range(0,len(scenario.D))]
    # rank_candidates finds alpha values for every possible rev tuple
    candidates=rec.rank_candidates(scenario.D,V)
    fail_wp3 = False
    with open('output/logging.txt', 'w',encoding="utf-8") as f:
        with redirect_stdout(f):
            print("============================================== ALL")
            tree = rec.recognize(D=scenario.D, B=[], first_candidate_only=False, small_spike=False, print_info=True)
            #valid_divs, invalid_divs = classify_divergence(tree)
            print("============================================== WP3")
            tree3 = rec.recognize(D=scenario.D, B=[0,1,2,3], first_candidate_only=True, small_spike=False, print_info=True)
            fail_wp3 = True if tree3.root.valid_ways == 0 else False
    # if not fail_wp3:
    #     print("VALID SCENARIO")
    #     sys.exit()
    tree.visualize(save_as='output/vis_all.pdf')

    return candidates


def get_matrix(filepath):
    scenario=sim.load(filepath)
    return scenario.D

def init_scenario_sel(filepath,in_candidates):
    scenario=sim.load(filepath)
    V=[i for i in range(0,len(scenario.D))]
    # rank_candidates finds alpha values for every possible rev tuple
    candidates=rec.rank_candidates_selective(scenario.D,V,in_candidates)
    return candidates

def rnf_candidates_table(x,y,z,candidates,nodes):
    out_table=""

    out_table += str((x,y,z)) + "\t"
    for u in nodes:
        if(u in [x,y,z]):
            continue
        out_table += str(u) + "\t"

    out_table = out_table[0:-1]
    out_table += "\n"

    for u in nodes:
        if(u in [x,y,z]):
            continue
        out_table += str(u) + "\t\t"
        for v in nodes:
            if v in [x,y,z,u]:
                if (v == u):
                    out_table += "-\t"
                continue
            if candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))] < 0:
                    out_table += "< 0\t"
            elif candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))] > 1 :
                    out_table += "> 1\t"
            elif math.isnan(candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))]) :
                    out_table += "NaN\t"
            else:
                    out_table += str(round(candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))],2))+"\t"
        out_table += "\n"

    out_table += "\n\n******************************************** \n\n"
    out_table += str((x,y,z)) + "\t"
    for u in nodes:
        if(u in [x,y,z]):
            continue
        out_table += str(u) + "\t"

    out_table = out_table[0:-1]
    out_table += "\n"

    for u in nodes:
        if(u in [x,y,z]):
            continue
        out_table += str(u) + "\t\t"
        for v in nodes:
            if v in [x,y,z,u]:
                if (v == u):
                    out_table += "-\t"
                continue
            if candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))] < 0:
                    out_table += str(round(candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))],2))+"\t"
            elif candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))] > 1 :
                    out_table += str(round(candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))],2))+"\t"
            elif math.isnan(candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))]) :
                    out_table += str(round(candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))],2))+"\t"
            else:
                    out_table += str(round(candidates[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))],2))+"\t"
        out_table += "\n"
    out_table = out_table[0:-1]
    return out_table

def rnf_candidates(candidates):
    sort_cand={}
    agree_cand=[]
    agree_cons_arr=[]
    for k in candidates.keys():
        acc_alpha=0
        for k2 in candidates[k].keys():
            if(candidates[k][k2]>1 or candidates[k][k2]<0):
                acc_alpha+=abs(candidates[k][k2])
        sort_cand[k]=round(acc_alpha,3)

    # sort candidates by descending acc_alpha value
    sort_cand_out=sorted(sort_cand.items(),key=lambda x: x[1],reverse=True)

    # We count consensus of alphas per candidate
    for k in sort_cand_out:

        # We format the consens in a print string
        print_str=""
        print_str += str(k[0])+" - " + str(k[1]) + " || "
        # finding the number of witnesses that agree on an alpha is
        # a lot easier when the witnesses are sorted by their alpha
        sorted_c=sorted(candidates[k[0]].items(),key=lambda x: x[1])

        # saving the last occured alpha, and its support in the tuple
        last_comp=("notset",0)
        nan_counter=0
        inf_counter=0
        neg_inf_counter=0
        clean_sorted_c=[]
        agree=True
        agree_cons=-1
        for i in range(0,len(sorted_c)):
            if(math.isnan(sorted_c[i][1])):
                nan_counter+=1
            elif(math.inf==sorted_c[i][1]):
                #agree=False
                inf_counter+=1
            elif(-math.inf==sorted_c[i][1]):
                #agree=False
                neg_inf_counter+=1
            else:
                clean_sorted_c+=[sorted_c[i]]
        sorted_c=clean_sorted_c
        if nan_counter>0:
            print_str+="nan:"+str(nan_counter)+" "
        if inf_counter>0:
            print_str+="inf:"+str(inf_counter)+ " "
        if neg_inf_counter>0:
            print_str+="-inf:"+str(neg_inf_counter)+ " "
        for i in range(0,len(sorted_c)):
            if(last_comp[0]=="notset"):
                agree_cons=sorted_c[i][1]
                last_comp=(sorted_c[i][1],1)
            elif math.isclose(last_comp[0],sorted_c[i][1]):
                last_comp=(sorted_c[i][1],last_comp[1]+1)
            else:
                agree=False
                print_str+= str(round(last_comp[0],7)) + ":" + str(last_comp[1]) + " // "
                last_comp=(sorted_c[i][1],1)
        if(last_comp[0]!="notset"):
            print_str+= str(round(last_comp[0],7)) + ":" + str(last_comp[1])
        if(len(sorted_c)==0):
            print(print_str)
            continue
        else:
            if(agree and 0<=agree_cons<=1):
                agree_cand+=[k[0]]
                agree_cons_arr+=[agree_cons]
        print(print_str)
    return (agree_cand,agree_cons_arr)

candidates=init_scenario(argv[1])
agree_t=rnf_candidates(candidates)
agree_cand=agree_t[0]
agree_cand_alphas=agree_t[1]

for i in range(0,len(agree_cand)):
    print("["+str(i)+"] - " + str(agree_cand[i]))
scenario_string=""
with open(argv[1],"r") as f:
    scenario_string=f.read()
    #print(scenario_string)
spikes_comb=[(0,0),(1,1),(2,2),(3,3)]
hybr_comb=[(0,2),(0,3),(1,2),(1,3)]
curr_N=len(scenario_string.split("\n"))-1
#print(curr_N)
curr_D=[]
curr_V=[]

D_original = []

# subprocess.call("start firefox file:///C:/Users/Brujo/Documents/GitHub/explore_erdbeermet/output/vis_all.pdf &", shell=True,stdout=subprocess.DEVNULL)
# subprocess.call("firefox output/vis_all.pdf &", shell=True,stdout=subprocess.DEVNULL)
subprocess.call("google-chrome output/vis_all.pdf &", shell=True,stdout=subprocess.DEVNULL)


overview_str=""
while True:
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    dec_string=input("Do you wish to examine consens candidates or perform revR steps? cC/rR for choice.\n")
    #print(dec_string)
    if('c' in dec_string or 'C' in dec_string):
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        for i in range(0,len(agree_cand)):
            print("["+str(i)+"] - " + str(agree_cand[i]) + " " + str(agree_cand_alphas[i]))
        cand_string=input("Please Enter two number (i.e. 1,2) to test candidates against each other.\n")
        #print(cand_string)
        comp_cand=[agree_cand[int(x.strip())] for x in cand_string.split(",")]
        parents=[comp_cand[0][0],comp_cand[0][1],comp_cand[1][0],comp_cand[1][1]]
        print(comp_cand)
        for i in range(0,len(spikes_comb)):
            s=spikes_comb[i]
            update_scen_str=scenario_string
            update_scen_str+="("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(curr_N+1)+")"
            update_scen_str+=" 1.0; ["
            for j in range(0,curr_N+1):
                update_scen_str+="0.0,"
            update_scen_str+="1.0]\n"
            with open("rmet_tmp_scen","w+") as f:
                f.write(update_scen_str)
            overview_str+="-----------------------------------------------------"
            overview_str+="SPIKE - " + "("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(curr_N+1)+")"
            print("-----------------------------------------------------")
            print("SPIKE - " + "("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(curr_N+1)+")")
            rnf_candidates(init_scenario_sel("rmet_tmp_scen",comp_cand))
        for i in range(0,len(hybr_comb)):
            s=hybr_comb[i]
            update_scen_str=scenario_string
            update_scen_str+="("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(curr_N+1)+")"
            update_scen_str+=" 0.5; ["
            for j in range(0,curr_N+1):
                update_scen_str+="0.0,"
            update_scen_str+="1.0]\n"
            with open("rmet_tmp_scen","w+") as f:
                f.write(update_scen_str)
            overview_str+="-----------------------------------------------------"
            overview_str+="HYBRIDIZATION - " + "("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(curr_N+1)+")"
            print("-----------------------------------------------------")
            print("HYBRIDIZATION - " + "("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(curr_N+1)+")")
            rnf_candidates(init_scenario_sel("rmet_tmp_scen",comp_cand))
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        comb_string=input("Do you wish to examine a COMBINATION? y/Y for yes.\n")
        if('y' not in comb_string and 'Y' not in comb_string):
            continue
        while True:
            print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            comb_string=input("Please Enter the numbers of resolving events to combine (Zero-based, comma-separated). X for abort, a for all.\n")
            if("x" in comb_string or "X" in comb_string):
                break
            elif("a" in comb_string or "A" in comb_string):
                events=[x for x in range(0,8)]
            else:
                events=[int(x.strip()) for x in comb_string.split(",")]
            print(events)
            next_N=curr_N
            update_scen_str=scenario_string
            for i in events:
                if(i in [0,1,2,3]):
                    s=spikes_comb[i]
                    update_scen_str+="("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(next_N+1)+")"
                    update_scen_str+=" 1.0; ["
                    for j in range(0,next_N+1):
                        update_scen_str+="0.0,"
                    update_scen_str+="1.0]\n"
                if(i in [4,5,6,7]):
                    s=hybr_comb[i-4]
                    update_scen_str+="("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(next_N+1)+")"
                    update_scen_str+=" 0.5; ["
                    for j in range(0,next_N+1):
                        update_scen_str+="0.0,"
                    update_scen_str+="1.0]\n"
                next_N=next_N+1
            with open("rmet_tmp_scen","w+") as f:
                f.write(update_scen_str)
            print(update_scen_str)
            rnf_candidates(init_scenario_sel("rmet_tmp_scen",comp_cand))
            #print(overview_str)

    elif("r" in dec_string or "R" in dec_string):
        r_cand_string=input("Please Enter ONE number (i.e. 0) to remove a candidate pair.\n")
        print(r_cand_string)
        r_cand=agree_cand[int(r_cand_string.strip())]
        r_alpha=agree_cand_alphas[int(r_cand_string.strip())]
        with open("rmet_tmp_scen","w+") as f:
            f.write(scenario_string)
        scen=sim.load("rmet_tmp_scen",)
        V=[i for i in range(0,len(scen.D))]
        u=-1
        for j in range(0,len(V)):
                    if(j!=r_cand[0] and j!=r_cand[1] and j!=r_cand[2]):
                        u=j
                        break
        deltas=rec._compute_deltas(V,scen.D,r_alpha,r_cand[0],r_cand[1],r_cand[2],u)
        V_copy=[x for x in range(0,curr_N)]
        #V_copy=V.copy()
        #V_copy.remove(r_cand[2])
        print("WITNESS " + str(u))
        print(deltas)
        print(scen.D)

        D_original = scen.D.copy()

        input()
        D_copy=scen.D.copy()
        D_copy=rec._matrix_without_index(D_copy,V.index(r_cand[2]))
        print(len(D_copy))
        rec._update_matrix(V_copy,D_copy,r_cand[0],r_cand[1],deltas[2],deltas[3])
        t=rnf_candidates(rec.rank_candidates(D_copy,V_copy))
        agree_cand=t[0]
        agree_cand_alphas=t[1]
        curr_D=D_copy.copy()
        curr_V=V_copy.copy()
        curr_N=curr_N-1
        if(curr_N==3):
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print(str( (-0.546*curr_D[0,2]+0.546*curr_D[2,3]+curr_D[0,1]-curr_D[1,3]-0.454*curr_D[0,3])/(0.546*curr_D[2,3]+curr_D[1,2]-curr_D[1,3]-0.454*curr_D[2,3]) ))
                print(str())
                #print(str((-0.546*11.786738069999998+0.546*22.233512520000005+0.454*18.272028071093594+23.027000000000008-7.2362874800000085-22.233512520000005)/(2*0.546*18.272028071093594+17.506921930000004-7.2362874800000085-8.272028071093594)))
                print(str( ((-0.546*curr_D[0,2])+(0.546*curr_D[0,3])+(0.546*curr_D[2,3])+curr_D[0,1]-curr_D[0,3]-curr_D[1,3])/(2*0.546*curr_D[2,3]+curr_D[1,2]-curr_D[1,3]-curr_D[2,3]) ))
                print(str( (-0.546*curr_D[1,2]+0.546*curr_D[1,3]+0.546*curr_D[2,3]+curr_D[0,1]-curr_D[1,3]-curr_D[0,3])/(2*0.546*curr_D[2,3]+curr_D[0,2]-curr_D[0,3]-curr_D[2,3]) ))
                sys.exit()
        while True:
            # prnt_str=""
            # for i in range(0,len(curr_D)):
            #    for j in range(0,len(curr_D)):
            #        prnt_str+=str(curr_D[i,j]) +"\t"
            #    prnt_str+="\n"
            # print(prnt_str)
            print(curr_D)
            print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            dec_string=input("Do you wish to examine consens candidates, perform revR steps, examine a single candidate, or a group of candidates, or test 2 candidates on 6 leaves? cC/rR/sS/gG/tT for choice.\n")
            if('c' in dec_string or 'C' in dec_string):
                print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                for i in range(0,len(agree_cand)):
                    print("["+str(i)+"] - " + str(agree_cand[i]) + " " + str(agree_cand_alphas[i]))
                cand_string=input("Please Enter two number (i.e. 1,2) to test candidates against each other.\n")
                comp_cand=[agree_cand[int(x.strip())] for x in cand_string.split(",")]
                parents=[comp_cand[0][0],comp_cand[0][1],comp_cand[1][0],comp_cand[1][1]]
                next_N=curr_N+1
                for i in range(0,len(spikes_comb)):
                    s=spikes_comb[i]
                    V=[l for l in range(0,next_N+1)]
                    dt=[0 for l in range(0,next_N+1)]
                    #dt[-1]=1.0
                    dt[-1]=0.0
                    print(len(curr_D))
                    D_copy=rec.add_child(curr_D.copy(),parents[s[0]],parents[s[1]],next_N,0.5,dt)
                    print(len(D_copy))
                    print("-----------------------------------------------------")
                    print("SPIKE - " + "("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(curr_N)+")")
                    rnf_candidates(rec.rank_candidates_selective(D_copy,V,comp_cand))
                for i in range(0,len(hybr_comb)):
                    s=hybr_comb[i]
                    V=[l for l in range(0,next_N+1)]
                    dt=[0 for l in range(0,next_N+1)]
                    #dt[-1]=1.0
                    dt[-1]=0.0
                    print(len(curr_D))
                    D_copy=rec.add_child(curr_D.copy(),parents[s[0]],parents[s[1]],next_N,0.5,dt)
                    print(len(D_copy))
                    print("-----------------------------------------------------")
                    print("HYBRIDIZATION - " + "("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(next_N)+")")
                    rnf_candidates(rec.rank_candidates_selective(D_copy,V,comp_cand))
                print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                comb_string=input("Do you wish to examine a COMBINATION? y/Y for yes.\n")
                if('y' not in comb_string and 'Y' not in comb_string):
                    continue
                while True:
                    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                    comb_string=input("Please Enter the numbers of resolving events to combine (Zero-based, comma-separated). X or x for abort, a or A for all combinations.\n")
                    if("x" in comb_string or "X" in comb_string):
                        break
                    elif("a" in comb_string or "A" in comb_string):
                        events=[x for x in range(0,8)]
                    else:
                        events=[int(x.strip()) for x in comb_string.split(",")]
                    next_N=curr_N+1
                    D_copy=curr_D.copy()
                    for i in events:
                        if(i in [0,1,2,3]):
                            s=spikes_comb[i]
                            V=[l for l in range(0,next_N+1)]
                            dt=[0 for l in range(0,next_N+1)]
                            #dt[-1]=1.0
                            dt[-1]=0.0
                            D_copy=rec.add_child(D_copy.copy(),parents[s[0]],parents[s[1]],next_N,1,dt)
                            next_N=next_N+1
                        if(i in [4,5,6,7]):
                            s=hybr_comb[i-4]
                            V=[l for l in range(0,next_N+1)]
                            dt=[0 for l in range(0,next_N+1)]
                            #dt[-1]=1.0
                            dt[-1]=0.0
                            D_copy=rec.add_child(D_copy.copy(),parents[s[0]],parents[s[1]],next_N,0.5,dt)
                            next_N=next_N+1
                    rnf_candidates(rec.rank_candidates_selective(D_copy,V,comp_cand))
            elif("v" in dec_string or "V" in dec_string):
                print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                cand_string_4 = input("Please enter candidate indizes csv to examine box graph.")
                arr = [int(x) for x in cand_string_4.split(",")]
                visualize_splitstree(curr_D,"COMMANDSOUT.nex","4BOX.png",arr)
            elif "m" in dec_string or "M" in dec_string: # MAP OUT
                # cand_string_3 = input("Please enter children to examine as CSV!")
                #
                # cand_array = [int(x) for x in cand_string_3.split(",")]
                #
                # alphas = rec.rank_candidates(curr_D.copy(),curr_V.copy())

                cand_string_3 = input("Please enter children to examine as CSV!")

                cand_array = [(int(x.split(",")[0]),int(x.split(",")[1])) for x in cand_string_3.split(";")]

                alphas = rec.rank_candidates(curr_D.copy(),curr_V.copy())

                for t in cand_array:
                    z = t[0]
                    v = t[1]
                    for x in range(0,len(curr_V)):
                        for y in range(0,len(curr_V)):
                            if y <= x or z in [x,y] or v in [x,y]:
                                continue
                            # print("--------------------")
                            # print("(" + str(x) + ", " + str(y) + ": " + str(z) + ")  ")
                            invalid_counter = 0
                            negative_counter = 0
                            positive_counter = 0
                            for u in range(0,len(curr_V)):
                                if u not in [x,y,z,v]:
                                    if alphas[(x,y,z)][(min(u,v),max(u,v))]<0:
                                        # print(str(u) + " - " + str(v) + " <0")
                                        invalid_counter += 1
                                        negative_counter += 1
                                    if alphas[(x,y,z)][(min(u,v),max(u,v))]>1:
                                        invalid_counter += 1
                                        positive_counter += 1
                                        # print(str(u) + " - " + str(v) + " >1")
                            if invalid_counter < 1:
                                print("--------------------")
                                print("(" + str(x) + ", " + str(y) + ": " + str(z) + ")  ")
                                print("INVALID - " + str(invalid_counter))

                    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

                    x = t[1]
                    for y in range(0,len(curr_V)):
                        if z in [x,y]:
                            continue
                        # print("--------------------")
                        # print("(" + str(x) + ", " + str(y) + ": " + str(z) + ")  ")
                        invalid_counter = 0
                        negative_counter = 0
                        positive_counter = 0
                        for u in range(0,len(curr_V)):
                            for v in range(0,len(curr_V)):
                                if v<=u or u in [x,y,z] or v in [x,y,z]:
                                    continue
                                if alphas[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))]<0:
                                    # print(str(u) + " - " + str(v) + " <0")
                                    invalid_counter += 1
                                    negative_counter += 1
                                if alphas[(min(x,y),max(x,y),z)][(min(u,v),max(u,v))]>1:
                                    # print(str(u) + " - " + str(v) + " >1")
                                    invalid_counter += 1
                                    positive_counter += 1
                        if invalid_counter < 1:
                            print("--------------------")
                            print("(" + str(x) + ", " + str(y) + ": " + str(z) + ")  ")
                            print("INVALID - " + str(invalid_counter))
                        # print(">1 - "+ str(positive_counter))
                        # print("<0 - "+ str(negative_counter))


                                        # print("(" + str(x) + ", " + str(y) + ": " + str(z) + ")  " + str(u) + "-" + str(v) + ":  " + str(round(alphas[(x,y,z)][(min(u,v),max(u,v))],2)) )
                                         # print(str(u) + "-" + str(v) + ":  " + str(round(alphas[(x,y,z)][(min(u,v),max(u,v))],2)) )
                                    # print(alphas[(x,y,z)][(min(u,v),max(u,v))])


                # for t in cand_array:
                #     z = t[0]
                #     v = t[1]
                #     for x in range(0,len(curr_V)):
                #         for y in range(0,len(curr_V)):
                #             if y <= x or z in [x,y] or v in [x,y]:
                #                 continue
                #             print("------------------")
                #             print(str(x) + "," + str(y) + ":" + str(z))
                #             for u in range(0,len(curr_V)):
                #                 if u not in [x,y,z,v]:
                #                     print(str(u) + "-" + str(v) + ":  " + str(round(alphas[(x,y,z)][(min(u,v),max(u,v))],2)))
                #                     # print(alphas[(x,y,z)][(min(u,v),max(u,v))])



                # for z in cand_array:
                #     for x in range(0,len(curr_V)):
                #         for y in range(0,len(curr_V)):
                #             if y <= x and z in [x,y]:
                #                 continue
                #             print("------------------")
                #             print(str(x) + "," + str(y) + ":" + str(z))
                #             for u in range(0,len(curr_V)):
                #                 for v in range(u,len(curr_V)):
                #                     if u not in [x,y,z] and v not in [x,y,z] and u != v:
                #                         print(str(u) + "-" + str(v) + ":")
                #                         print(alphas[(x,y,z)][(u,v)])

            elif("e" in dec_string or "E" in dec_string):
                cand_string_3 = input("Please enter three numbers to examine chosen pairs alpha witnesses.")
                cand_tuple = [int(x) for x in cand_string_3.split(",")]
                # nodes_here = [x for x in range(0,curr_N)]
                cands = rec.rank_candidates(curr_D.copy(),curr_V.copy())
                out_table = rnf_candidates_table(cand_tuple[0],cand_tuple[1],cand_tuple[2],cands,curr_V)

                print("\n\n----------------------------------------------------\n")
                print(out_table)
                print("\n----------------------------------------------------\n\n")

            elif("t" in dec_string or "T" in dec_string):
                for i in range(0,len(agree_cand)):
                    print("["+str(i)+"] - " + str(agree_cand[i]) + " " + str(agree_cand_alphas[i]))
                print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                cand_string_6 = input("Please enter two numbers to examine candidate pairs separately.")
                c_tuple = agree_cand[int(cand_string_6.split(",")[0])]
                f_tuple = agree_cand[int(cand_string_6.split(",")[1])]

                c_alpha = agree_cand_alphas[int(cand_string_6.split(",")[0])]
                f_alpha = agree_cand_alphas[int(cand_string_6.split(",")[1])]

                print(c_tuple)
                print(c_alpha)
                print(f_tuple)
                print(f_alpha)

                check_simple_6(curr_D.copy(),c_tuple[0],c_tuple[1],c_tuple[2],f_tuple[0],f_tuple[1],f_tuple[2],c_alpha,f_alpha)

                alphas = rank_candidates(curr_D.copy(), curr_V.copy())
                rev_a_checking_new(curr_D.copy(),c_tuple[0],c_tuple[1],c_tuple[2],f_tuple[2],c_alpha,alphas)

                compare_candidates(curr_D.copy(),c_alpha,alphas,c_tuple[0],c_tuple[1],c_tuple[2],f_tuple[0],f_tuple[1],f_tuple[2])

                test_deltas_new(curr_D.copy(),f_alpha,c_alpha,f_tuple[0],f_tuple[1],f_tuple[2],c_tuple[0],c_tuple[1],c_tuple[2])



            elif("x" in dec_string or "X" in dec_string):
                curr_D = D_original.copy()
                curr_V = [i for i in range(0,len(curr_D))]
                curr_N = len(curr_D)-1
                t=rnf_candidates(rec.rank_candidates(curr_D,curr_V))
                agree_cand=t[0]
                agree_cand_alphas=t[1]

            elif "l" in dec_string or "L" in dec_string:
                print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

                tuple_data = []

                for z in curr_V:
                    z_min = -1
                    parent_tuple = []

                    for x in curr_V:
                        for y in curr_V:
                            if x == y or x == z or y == z:
                                continue
                            dz = rec._compute_delta_z(curr_D[x, y], curr_D[x, z], curr_D[y, z])

                            if dz < z_min or z_min < 0:
                                z_min = dz
                                parent_tuple = (x, y, z)
                    tuple_data += [[parent_tuple, z_min]]

                for t in tuple_data:
                    print(str(t[0]) + " --- " + str(t[1]))

                print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

            elif("d" in dec_string or "D" in dec_string):
                for i in range(0,len(agree_cand)):
                    print("["+str(i)+"] - " + str(agree_cand[i]) + " " + str(agree_cand_alphas[i]))
                print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                cand_string = input("Please enter the number of the candidate pair to examine.")
                choice = int(cand_string)

                chosen_child = agree_cand[choice][2]

                deltas_zs = []

                for x in range(0,len(curr_D)):
                    for y in range(0,len(curr_D)):

                        if x == y or x == chosen_child or y == chosen_child or x > y:
                            continue

                        dz = rec._compute_delta_z(curr_D[x,y],curr_D[x,chosen_child],curr_D[y,chosen_child])

                        deltas_zs += [(x,y,dz)]

                print("------------------------")

                deltas_zs = sorted(deltas_zs,key=lambda x: x[2])

                for t in deltas_zs:
                    print(t)

                time = [x for x in range(0,len(deltas_zs))]
                data = [x[2] for x in deltas_zs]

                print("------------------------")
                func_arr = []

                for u in range(0,len(curr_D)):
                    if u == agree_cand[choice][0] or u == agree_cand[choice][1] or u == chosen_child:
                        continue
                    t = get_boxes(curr_D,[agree_cand[choice][0],agree_cand[choice][1],agree_cand[choice][2],u],agree_cand[choice][0],agree_cand[choice][1],agree_cand[choice][2],u,agree_cand_alphas[choice])
                    func_arr += [t]

                func_arr = sorted(func_arr,key=lambda x: x[1])

                print("*******************")
                for t in func_arr:
                    print(t)
                print("*******************")

                time = [x[1] for x in func_arr]
                data = [x[0] for x in func_arr]

                fig, ax = plt.subplots()
                ax.plot(time,data)

                ax.set(xlabel="time",ylabel="delta_z",title="maybe we see something, maybe we dont")
                ax.grid()

                fig.savefig("PLOTTING_DELTAS.png")
                plt.show()




            elif('g' in dec_string or 'G' in dec_string):
                all_alphas = rec.rank_candidates(curr_D,curr_V)
                for i in range(0,len(agree_cand)):
                    print("["+str(i)+"] - " + str(agree_cand[i]) + " " + str(agree_cand_alphas[i]))
                print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                cand_string=input("Please enter candidates to check, comma separated.")
                comp_cand=[agree_cand[int(x.strip())] for x in cand_string.split(",")]
                if len(comp_cand)==1:
                    parents=[comp_cand[0][0],comp_cand[0][1]]
                    child=comp_cand[0][2]
                    us=[]
                    for j in range(0,len(curr_D)):
                        if(j!=parents[0] and j!=parents[1] and j!=child):
                            us+=[j]
                            u=j
                            #break
                    alpha=agree_cand_alphas[int(cand_string.strip())]
                    # print(us)
                    # print("^!")
                    # u_perms = list(permutations(us,2))
                    # print(u_perms)

                    deltas=rec._compute_deltas(curr_V,curr_D,alpha,parents[0],parents[1],child,us[0])
                    if (deltas[1] < 0 and not np.isclose(deltas[1],0.0) ) or (deltas[2] < 0 and not np.isclose(deltas[2],0.0) ) or (deltas[3] < 0 and not np.isclose(deltas[3],0.0) ) :
                        print("---------------------")
                        print("ALPHAS - " + str(parents[0]) + ", " + str(parents[1]) + ":" + str(child))
                        print("DELTAS SMALLER ZERO")
                        print(deltas)

                    # Delete child from matrix
                    curr_D2=curr_D.copy()
                    curr_V2=curr_V.copy()
                    print("-----------" + str(check_candidate(curr_D.copy(),parents[0],parents[1],child,alpha,all_alphas,print_info=True)) + "-----------")
                    # continue

                    curr_D2=curr_D.copy()
                    curr_V2=curr_V.copy()
                    curr_D2=rec._update_matrix_return(curr_V2,curr_D2,parents[0],parents[1],deltas[2],deltas[3])
                    curr_D2=rec._matrix_without_index(curr_D,curr_V.index(child))
                    curr_V2=[l for l in range(0,len(curr_D2))]

                    for j in range(0,len(curr_V)):
                        if(j!=parents[0] and j!=parents[1] and j!=child):
                            # print(j)
                            # print("J")

                            # deltas=rec._compute_deltas(curr_V,curr_D,alpha,parents[0],parents[1],child,j)
                            # if (deltas[1] < 0 and not np.isclose(deltas[1],0.0) ) or (deltas[2] < 0 and not np.isclose(deltas[2],0.0) ) or (deltas[3] < 0 and not np.isclose(deltas[3],0.0) ) :
                            #     print("---------------------")
                            #     print("ALPHAS - " + str(parents[0]) + ", " + str(parents[1]) + ":" + str(child))
                            #     print("DELTAS SMALLER ZERO")
                            #     print(deltas)
                            # # Delete child from matrix
                            # curr_D2=curr_D.copy()
                            # curr_V2=curr_V.copy()
                            # rec._update_matrix(curr_V2,curr_D2,parents[0],parents[1],deltas[2],deltas[3])
                            # curr_D2=rec._matrix_without_index(curr_D,curr_V.index(child))
                            # curr_V2=[l for l in range(0,len(curr_D2))]

                            print("---------------------")
                            print("ALPHAS - " + str(parents[0]) + ", " + str(parents[1]) + ":" + str(child))

                            # Check alphas for all legit witness combinations
                            # TODO: CHANGE THE COMBINATIONS BACK TO U; WERE ONLY DEPENDANT ON U.
                            for u in us:
                                # print(comb)
                                # if( parents[0] not in comb and parents[1] not in comb and child not in comb and j not in comb):
                                if u != parents[0] and u!= parents[1] and u!=child and u!=j:
                                    # w_0 = comb[0]
                                    # w_1 = comb[1]
                                    w_0 = u
                                    p_0 = parents[0]
                                    p_1 = parents[1]
                                    mock_child = j

                                    if p_0 > child: p_0-=1
                                    if p_1 > child: p_1-=1
                                    # if w_0 > child: w_0-=1
                                    # if w_1 > child: w_1-=1
                                    if u > child: w_0-=1
                                    if mock_child > child: mock_child-=1

                                    #print("CHECKING - " + str(parents[0]) + ", " + str(parents[1]) + ":" + str(mock_child) + "(WITNESS " + str(w_0) + " " + str(w_1) + " )" )
                                    print("CHECKING - " + str(parents[0]) + ", " + str(parents[1]) + ":" + str(j) + " ( WITNESS " + str(u) + " )" )
                                    check_alpha = rev_a_val(curr_D2,alpha,p_0,p_1,mock_child,w_0,w_0)
                                    chk_1=min(u,child)
                                    chk_2=max(u,child)
                                    print(str(check_alpha) + " <-> " + str(all_alphas[ (parents[0],parents[1],j) ][(chk_1,chk_2)]) + " - " + str(np.isclose(all_alphas[ (parents[0],parents[1],j) ][(chk_1,chk_2)],check_alpha)))

                        # for u in us:
                        #     if(j!=parents[0] and j!=parents[1] and j!=child and j!=u):
                        #         print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                        #         print("ALPHAS - " + str(parents[0]) + ", " + str(parents[1]) + ":" + str(child) + "(WITNESS " + str(u)+" )" )
                        #         deltas=rec._compute_deltas(curr_V,curr_D,alpha,parents[0],parents[1],child,u)
                        #         curr_D2=curr_D.copy()
                        #         curr_D2=rec._matrix_without_index(curr_D,curr_V.index(child))
                        #         curr_V2=[l for l in range(0,len(curr_D2))]
                        #         rec._update_matrix(curr_V2,curr_D2,parents[0],parents[1],deltas[2],deltas[3])
                        #         p0=parents[0]
                        #         p1=parents[1]
                        #         c=j
                        #         if(parents[0]>child):
                        #             p0-=1
                        #         if(parents[1]>child):
                        #             p1-=1
                        #         if(j>child):
                        #             c-=1
                        #         u_new=u
                        #         if u_new>child:
                        #             u_new-=u-1
                        #             #u-=1
                        #         #print("-------")
                        #         #print(u_new)
                        #         #print(u)
                        #         #print(type(u_new))
                        #         #print(type(u))
                        #         print("CHECKING - " + str(parents[0]) + ", " + str(parents[1]) + ":" + str(j) + "(WITNESS " + str(u_new) + " " + str(c) + " )" )
                        #         print(str(rev_a_val(curr_D2,alpha,p0,p1,c,u_new,child)))
                else:
                    parents=[comp_cand[0][0],comp_cand[0][1]]
                    children=[x[2] for x in comp_cand]
                    alphas=[]
                    for x in range(0,len(comp_cand)):
                        alphas+=[agree_cand_alphas[int(cand_string.split(",")[x].strip())]]
                    for x in range(0,len(children)):
                        for j in range(0,len(curr_V)):
                            if(j!=parents[0] and j!=parents[1] and j!=children[x]):
                                u=j
                                break
                        #print(curr_V)
                        #print(curr_D)
                        #print(alphas)
                        #print(parents)
                        #print(children)
                        #print(u)
                        deltas=rec._compute_deltas(curr_V,curr_D,alphas[x],parents[0],parents[1],children[x],u)
                        curr_D2=curr_D.copy()
                        curr_D2=rec._matrix_without_index(curr_D,curr_V.index(children[x]))
                        curr_V2=[l for l in range(0,len(curr_D2))]
                        rec._update_matrix(curr_V2  ,curr_D2,parents[0],parents[1],deltas[2],deltas[3])
                        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                        print("ALPHAS - " + str(parents[0]) + ", " + str(parents[1]) + ":" + str(children[x]))
                        for check_child in children:
                            if check_child==children[x]:
                                continue
                            for j in range(0,len(curr_V)):
                                if(j!=parents[0] and j!=parents[1] and j!=children[x] and j!=check_child):
                                    u=j
                                    break
                            # NEED TO FIT FOR REMOVED CHILD
                            p0=parents[0]
                            p1=parents[1]
                            c=check_child
                            if(parents[0]>children[x]):
                                p0-=1
                            if(parents[1]>children[x]):
                                p1-=1
                            if(check_child>children[x]):
                                c-=1
                            if u>children[x]:
                                u-=1
                            print("-------")
                            print("CHECKING - " + str(parents[0]) + ", " + str(parents[1]) + ":" + str(check_child))
                            #print("rev_a_val(D,"+str(alphas[x])+","+str(parents[0])+","+str(parents[1])+","+str(check_child)+","+str(u)+","+str(children[x]))
                            print(str(rev_a_val(curr_D2,alphas[x],p0,p1,c,u,children[x])))
            elif("r" in dec_string or "R" in dec_string):
                for i in range(0,len(agree_cand)):
                    print("["+str(i)+"] - " + str(agree_cand[i]))
                r_cand_string=input("Please Enter ONE number (i.e. 0) to remove a candidate pair.\n")
                r_cand=agree_cand[int(r_cand_string.strip())]
                r_alpha=agree_cand_alphas[int(r_cand_string.strip())]
                u=-1
                for j in range(0,len(curr_V)):
                    if(j!=r_cand[0] and j!=r_cand[1] and j!=r_cand[2]):
                        u=j
                        break
                deltas=rec._compute_deltas(curr_V,curr_D,r_alpha,r_cand[0],r_cand[1],r_cand[2],u)
                curr_D2=curr_D.copy()
                curr_D=rec._matrix_without_index(curr_D,curr_V.index(r_cand[2]))
                curr_V=[l for l in range(0,len(curr_D))]
                rec._update_matrix(curr_V,curr_D,r_cand[0],r_cand[1],deltas[2],deltas[3])
                t=rnf_candidates(rec.rank_candidates(curr_D,curr_V))
                agree_cand=t[0]
                agree_cand_alphas=t[1]
                print("Still Pseudometrik? - " + str(rec.is_pseudometric(curr_D)))
                if(len(curr_D)==4):
                    print("Valid 4? - " + str(rec.recognize4_matrix_only(curr_D)))
                curr_N=curr_N-1
                # if(curr_N==3 or True):
                #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 234")
                #     #print(str( (-0.546*curr_D[0,2]+0.546*curr_D[2,3]+curr_D[0,1]-curr_D[1,3]-0.454*curr_D[0,3])/(0.546*curr_D[2,3]+curr_D[1,2]-curr_D[1,3]-0.454*curr_D[2,3]) ))
                #     #print(str( ((-0.546*curr_D[0,2])+(0.546*curr_D[0,3])+(0.546*curr_D[2,3])+curr_D[0,1]-curr_D[0,3]-curr_D[1,3])/(2*0.546*curr_D[2,3]+curr_D[1,2]-curr_D[1,3]-curr_D[2,3]) ))
                #     #def rev_a_val(D,a,x,y,z,u,v):
                #     print(rev_a_val(curr_D,0.546,2,3,1,0,4))
                #    # print(str((-0.546*11.786738069999998+0.546*22.233512520000005+0.454*18.272028071093594+23.027000000000008-7.2362874800000085-22.233512520000005)/(2*0.546*18.272028071093594+17.506921930000004-7.2362874800000085-8.272028071093594)))
                #     #print(str( (-0.546*curr_D[1,2]+0.546*curr_D[1,3]+0.546*curr_D[2,3]+curr_D[0,1]-curr_D[1,3]-curr_D[0,3])/(2*0.546*curr_D[2,3]+curr_D[0,2]-curr_D[0,3]-curr_D[2,3]) ))
                #     print(rev_a_val(curr_D,0.546,2,3,0,1,4))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!r!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 231")
                #     print(rev_a_val(curr_D,0.7729,1,2,3,0,1))
                #     print(rev_a_val(curr_D,0.7729,1,2,0,3,1))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 230")
                #     print(rev_a_val(curr_D,0.2366,1,2,0,3,0))
                #     print(rev_a_val(curr_D,0.2366,1,2,3,0,0))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 012")
                #     print(rev_a_val(curr_D,0.1357,0,1,2,3,2))
                #     print(rev_a_val(curr_D,0.1357,0,1,3,2,2))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 013")
                #     print(rev_a_val(curr_D,0.1638,0,1,2,3,3))
                #     print(rev_a_val(curr_D,0.1638,0,1,3,2,3))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 014")
                #     print(rev_a_val(curr_D,0.9300,0,1,2,3,4))
                #     print(rev_a_val(curr_D,0.9300,0,1,3,2,4))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                #     #print("Alphas 345")
                #     #print(rev_a_val(curr_D,0.974,3,4,1,2,5))
                #     #print(rev_a_val(curr_D,0.974,3,4,1,2,5))
                #     print("Alphas 021")
                #     print(rev_a_val(curr_D,0.4117647,0,1,2,3,1))
                #     print(rev_a_val(curr_D,0.4117647,0,1,3,2,1))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 023")
                #     print(rev_a_val(curr_D,0.5,0,2,1,3,3))
                #     print(rev_a_val(curr_D,0.5,0,2,3,1,3))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 024")
                #     print(rev_a_val(curr_D,0.243,0,2,1,3,4))
                #     print(rev_a_val(curr_D,0.243,0,2,3,1,4))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 234")
                #     print(rev_a_val(curr_D,0.546,2,3,1,0,4))
                #     print(rev_a_val(curr_D,0.546,2,3,0,1,4))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!r!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 231")
                #     print(rev_a_val(curr_D,0.5841,1,2,3,0,1))
                #     print(rev_a_val(curr_D,0.5841,1,2,0,3,1))
                #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                #     print("Alphas 230")
                #     print(rev_a_val(curr_D,0.4885,1,2,0,3,0))
                #     print(rev_a_val(curr_D,0.4885,1,2,3,0,0))
                    #print(str(  (-0.772919*curr_D[0,2]+0.772919*curr_D[0,3]+0.772919*curr_D[2,3]+curr_D[0,4]-curr_D[0,3]-curr_D[3,4])/(2*0.772919*curr_D[2,3]+curr_D[2,4]-curr_D[3,4]-curr_D[2,3] ))
                    #print(str(  (-0.772919*curr_D[0,1]+0.772919*curr_D[0,2]+0.772919*curr_D[1,2]+curr_D[0,3]-curr_D[0,2]-curr_D[2,3])/(2*0.772919*curr_D[1,2]+curr_D[1,3]-curr_D[2,3]-curr_D[1,2] )))
                    #print(str( print(str( ((-0.772919*curr_D[2,4])+(0.772919*curr_D[3,4])+(0.772919*curr_D[2,3])+curr_D[0,4]-curr_D[3,4]-curr_D[0,3])/(2*0.772919*curr_D[2,3]+curr_D[0,2]-curr_D[0,3]-curr_D[2,3]) ))  ))
                    #print(str( print(str( ((-0.772919*curr_D[1,3])+(0.772919*curr_D[2,3])+(0.772919*curr_D[1,2])+curr_D[0,3]-curr_D[2,3]-curr_D[0,2])/(2*0.772919*curr_D[1,2]+curr_D[0,1]-curr_D[0,2]-curr_D[1,2]) ))  ))
                if(curr_N==3):
                    curr_N+=1
                    curr_D=curr_D2
                    curr_V+=[len(curr_V)]
                    t=rnf_candidates(rec.rank_candidates(curr_D,curr_V))
                    agree_cand=t[0]
                    agree_cand_alphas=t[1]
            elif('s' in dec_string or 'S' in dec_string):
                print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                for i in range(0,len(agree_cand)):
                    print("["+str(i)+"] - " + str(agree_cand[i]) + " " + str(agree_cand_alphas[i]))
                cand_string=input("Please Enter a number (i.e. 1) to test a candidate.\n")
                cand=int(cand_string.strip())
                comp_cand=agree_cand[cand]
                print(comp_cand)
                next_N=curr_N+1
                spks=[(0,0),(1,1)]
                parent=[comp_cand[0],comp_cand[1],comp_cand[2]]
                x=comp_cand[0]
                y=comp_cand[1]
                z=comp_cand[2]
                D_copy=curr_D.copy()
                a=agree_cand_alphas[cand]
                print(a)
                V=[i for i in range(0,len(D_copy))]
                u=-1
                for j in range(0,len(V)):
                    if(j!=x and j!=y and j!=z):
                        u=j
                        break
                del_z,xy,del_x,del_y=rec._compute_deltas(V,D_copy,a,x,y,z,u)
                V+=[len(D_copy)]
                dt=[0 for l in range(0,next_N+1)]
                #dt[-1]=1.0
                print("-----------------------------------------------------")
                print("HYBRIDIZATION - " + "("+str(parent[0])+", "+str(parent[1])+": "+str(next_N)+")")
                D_copy2=rec.add_child(D_copy.copy(),parent[0],parent[1],next_N,0.5,dt)
                rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand]))

                print("-----------------------------------------------------")
                print("HYBRIDIZATION - " + "("+str(parent[0])+", "+str(parent[2])+": "+str(next_N)+")")
                D_copy2=rec.add_child(D_copy.copy(),parent[0],parent[2],next_N,0.5,dt)
                rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand]))

                print("-----------------------------------------------------")
                print("HYBRIDIZATION - " + "("+str(parent[1])+", "+str(parent[2])+": "+str(next_N)+")")
                D_copy2=rec.add_child(D_copy.copy(),parent[1],parent[2],next_N,0.5,dt)
                rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand]))

                #for i in range(0,len(D_copy)):
                #    for j in range(0,len(D_copy)):
                #        if((i==x and j!=x) or (i!=x and j==x)):
                #            D_copy[i,j]=D_copy[i,j]-del_x
                print("-----------------------------------------------------")
                print("SPIKE - " + "("+str(parent[0])+", "+str(parent[0])+": "+str(next_N)+")")
                D_copy2=rec.add_child(D_copy.copy(),parent[0],parent[0],next_N,0.5,dt)
                rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand]))


                print("-----------------------------------------------------")
                print("SPIKE - " + "("+str(parent[1])+", "+str(parent[1])+": "+str(next_N)+")")
                D_copy2=rec.add_child(D_copy.copy(),parent[1],parent[1],next_N,0.5,dt)
                rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand]))

                print("-----------------------------------------------------")
                print("SPIKE - " + "("+str(parent[2])+", "+str(parent[2])+": "+str(next_N)+")")
                D_copy2=rec.add_child(D_copy.copy(),parent[2],parent[2],next_N,0.5,dt)
                rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand]))

                print("-----------------------------------------------------")
                print("SPIKE - " + "("+str(parent[0])+", "+str(parent[0])+": "+str(next_N)+")")
                print("SPIKE - " + "("+str(parent[1])+", "+str(parent[1])+": "+str(next_N)+")")
                #print("SPIKE - " + "("+str(parent[2])+", "+str(parent[2])+": "+str(next_N)+")")
                D_copy2=rec.add_child(D_copy.copy(),parent[0],parent[0],next_N,0.5,dt)
                V+=[len(V)]
                dt+=[0.0]
                D_copy2=rec.add_child(D_copy2.copy(),parent[1],parent[1],next_N+1,0.5,dt)
                V+=[len(V)]
                dt+=[0.0]
                D_copy2=rec.add_child(D_copy2.copy(),parent[0],parent[1],next_N+2,0.5,dt)
                #V+=[len(V)]
                #dt+=[0.0]
                #D_copy2=rec.add_child(D_copy2.copy(),parent[2],parent[2],next_N+2,0.5,dt)
                rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand]))
                #def _compute_deltas(V, D, alpha, x, y, z, u):
                #return delta_z, d_xy, delta_x, delta_y
                del_z1,xy1,del_x1,del_y1=rec._compute_deltas(V,D_copy2,a,x,y,z,next_N)
                del_z2,xy2,del_x2,del_y2=rec._compute_deltas(V,D_copy2,a,x,y,z,next_N+1)
                del_z3,xy3,del_x3,del_y3=rec._compute_deltas(V,D_copy2,a,x,y,z,next_N+2)
                dt=dt[0:-2]
                V=V[0:-2]
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("DELTA_X: "+ str(del_x))
                print("DELTA_Y: "+ str(del_y))
                print("DELTA_Z: "+ str(del_z))
                print("------------------")
                print("DELTA_X1: "+ str(del_x1))
                print("DELTA_Y1: "+ str(del_y1))
                print("DELTA_Z1: "+ str(del_z1))
                print("------------------")
                print("DELTA_X2: "+ str(del_x2))
                print("DELTA_Y2: "+ str(del_y2))
                print("DELTA_Z2: "+ str(del_z2))
                print("------------------")
                print("DELTA_X: "+ str(del_x3))
                print("DELTA_Y: "+ str(del_y3))
                print("DELTA_Z: "+ str(del_z3))
                ########## ADD DUMMY CHILD!!!!!!!!!!!!!!!
                #D_copy2=rec.add_child(D_copy.copy(),parent[2],parent[2],next_N,0.5,dt)

                #V+=[len(V)]
                #dt+=[0.0]
                #rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand]))
                #D_copy2=rec.add_child(D_copy2.copy(),parent[1],parent[1],next_N+1,0.5,dt)

                #V+=[len(V)]
                #dt+=[0.0]
                #for i in range(0,len(V)-1):
                #    if(i==parent[0] or i==parent[1] or i==comp_cand[2]):
                #        continue
                #    print("-----------------------------------------------------")
                #    print("HYBRIDIZATION - " + "("+str(parent[0])+", "+str(i)+": "+str(next_N)+")")
                #    D_copy2=rec.add_child(D_copy.copy(),parent[0],i,next_N,0.5,dt)
                #    rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand]))
                #    print("-----------------------------------------------------")
                #    print("HYBRIDIZATION - " + "("+str(parent[1])+", "+str(i)+": "+str(next_N)+")")
                #    D_copy2=rec.add_child(D_copy.copy(),parent[1],i,next_N,0.5,dt)
                #    rnf_candidates(rec._selective(D_copy2,V,[comp_cand]))
                dec_str2=input("Y/y for adding spikes?")
                if("y" in dec_str2 or "Y" in dec_str2):
                    while True:
                        D_copy2=D_copy.copy()
                        #D_copy2=rec.add_child(D_copy2,parent[0],parent[0],next_N,0.5,dt)
                        #D_copy2=rec.add_child(D_copy2,parent[1],parent[1],next_N+1,0.5,dt+[0.0])
                        for i in range(0,len(D_copy2)):
                            for j in range(0,len(D_copy2)):
                                if(j>i and i!=parent[0] and i!=parent[1] and i!=parent[2] and j!=parent[0] and j!=parent[1] and j!=parent[2]):
                                    al=agree_cand_alphas[cand]
                                    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                                    print(str(parent[0])+", "+str(parent[1]) +" : "+ str(parent[2]) + " - " + str(i) + " ; " + str(j))
                                    print("a:"+ str(D_copy2[parent[0],parent[1]]-del_x-del_y))
                                    print("b1: "+str((D_copy2[i,parent[0]]-del_x)-(D_copy2[i,parent[1]]-del_y)-(D_copy2[parent[0],parent[1]]-del_x-del_y)))
                                    print("b2: "+str((D_copy2[j,parent[0]]-del_x)-(D_copy2[j,parent[1]]-del_y)-(D_copy2[parent[0],parent[1]]-del_x-del_y)))
                                    print("c: "+ str(del_y))
                                    print("d: "+ str(del_x))
                                    print("ALPHA:" + str(agree_cand_alphas[cand]))
                                    print("bx: " + str(-2*(D_copy2[parent[0],parent[1]]-del_x-del_y)))
                                    print("by: " + "0")
                                    print("----------------------------")
                                    a=D_copy2[parent[0],parent[1]]-del_x-del_y
                                    b1=(D_copy2[i,parent[0]]-del_x)-(D_copy2[i,parent[1]]-del_y)-(D_copy2[parent[0],parent[1]]-del_x-del_y)
                                    b2=(D_copy2[j,parent[0]]-del_x)-(D_copy2[j,parent[1]]-del_y)-(D_copy2[parent[0],parent[1]]-del_x-del_y)
                                    bx=-2*(D_copy2[parent[0],parent[1]]-del_x-del_y)
                                    c=del_y
                                    d=del_x
                                    print("LOWER-BOUND-1: " + str( (b1-al*b1)/(2*(al*a-a-d)) ))
                                    print("UPPER-BOUND-1: " + str( (2*c-al*b1)/(2*(al*a+c)) ))
                                    print("LOWER-BOUND-2: " + str( (b2-al*b2)/(2*(al*a-a-d)) ))
                                    print("UPPER-BOUND-2: " + str( (2*c-al*b2)/(2*(al*a+c)) ))
                                   # print("LOWER-BOUND-X: " + str( (bx-al*b1)/(2*(al*a-a-d)) ))
                                    #print("UPPER-BOUND-X: " + str( (2*c-al*bx)/(2*(al*a+c)) )
                                    print("BOUND-X: " + str( (al*bx-bx)/(al*bx-bx+2*d) ))
                                    print("BOUND-Y: " + str( (c)/(a*al+c) ))
                        dec_str3=input("Please enter candidates to hybridize!")
                        p1=dec_str3.split(",")[0].strip()
                        p2=dec_str3.split(",")[1].strip()
                        a=float(dec_str3.split(",")[2].strip())
                        #p4=dec_str3.split(",")[3].strip()
                        print("-----------------------------------------------------")
                        print("HYBRIDIZATION - " + "("+str(p1)+", "+str(p2)+": "+str(next_N)+")")
                        #print("SPIKE - " + "("+p3+", "+str(p4)+": "+str(next_N+1)+")")
                        #D_copy2=rec.add_child(D_copy2.copy(),int(p1),int(p2),next_N+2,a,dt+[0.0,0.0])
                        print(next_N)
                        print(dt)
                        print(len(dt))
                        D_copy2=rec.add_child(D_copy2.copy(),int(p1),int(p2),next_N,a,dt)
                        #dt+=[0.0]
                        #D_copy2=rec.add_child(D_copy2.copy(),int(p3),int(p4),next_N+1,0.5,dt)
                        rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand]))
                        del_z1,xy1,del_x1,del_y1=rec._compute_deltas(V,D_copy2,a,x,y,z,next_N)
                        print("------------------")
                        print("DELTA_XH: "+ str(del_x1))
                        print("DELTA_YH: "+ str(del_y1))
                        print("DELTA_ZH: "+ str(del_z1))
                        print("------------------")
                        d_str=input("Ccontinue?")
                        uv=[]
                        if("c" in d_str or "C" in d_str):
                            continue
                        else:
                            break
                        #rnf_candidates(rec.rank_candidates_selective(D_copy2,V,[comp_cand,(comp_cand[0],comp_cand[1],next_N)]))
    else:
        break
