from erdbeermet import simulation as sim
from erdbeermet import recognition as rec
from contextlib import redirect_stdout
import numpy as np
from decimal import Decimal as d
import math
import subprocess
from sys import argv

def init_scenario(filepath):
    scenario=sim.load(argv[1])
    V=[i for i in range(0,len(scenario.D))]
    # rank_candidates finds alpha values for every possible rev tuple
    candidates=rec.rank_candidates(scenario.D,V)
    with open('output/logging.txt', 'w') as f:
        with redirect_stdout(f):
            print("============================================== ALL")
            tree = rec.recognize(D=scenario.D, B=[], first_candidate_only=False, small_spike=False, print_info=True)
            #valid_divs, invalid_divs = classify_divergence(tree)         
            print("============================================== WP3")
            tree3 = rec.recognize(D=scenario.D, B=[0,1,2,3], first_candidate_only=True, small_spike=False, print_info=True)    
            fail_wp3 = True if tree3.root.valid_ways == 0 else False  
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
  #  if(math.isclose(acc_alpha,0)):
  #      agree_cand[k]=candidates[k]

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
                agree=False
                inf_counter+=1
            elif(-math.inf==sorted_c[i][1]):
                agree=False
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

subprocess.call("firefox erdbeermet/output/vis_all.pdf &", shell=True,stdout=subprocess.DEVNULL)
overview_str=""
while True:
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    dec_string=input("Do you wish to examine consens candidates or perform revR steps? cC/rR for choice.\n")
    #print(dec_string)
    if('c' in dec_string or 'C' in dec_string):
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        for i in range(0,len(agree_cand)):
            print("["+str(i)+"] - " + str(agree_cand[i]))
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
        r_cand=agree_cand[int(r_cand_string.strip())]
        r_alpha=agree_cand_alphas[int(r_cand_string.strip())]
        with open("rmet_tmp_scen","w+") as f:
            f.write(scenario_string)
        scen=sim.load("rmet_tmp_scen",)
        V=[i for i in range(0,len(scen.D))]
        deltas=rec._compute_deltas(V,scen.D,r_alpha,r_cand[0],r_cand[1],r_cand[2],0)
        V_copy=[x for x in range(0,curr_N)]
        #V_copy=V.copy()
        #V_copy.remove(r_cand[2])
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
        while True:
            print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            dec_string=input("Do you wish to examine consens candidates or perform revR steps? cC/rR for choice.\n")
            if('c' in dec_string or 'C' in dec_string):
                print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                for i in range(0,len(agree_cand)):
                    print("["+str(i)+"] - " + str(agree_cand[i]))
                cand_string=input("Please Enter two number (i.e. 1,2) to test candidates against each other.\n")
                comp_cand=[agree_cand[int(x.strip())] for x in cand_string.split(",")]
                parents=[comp_cand[0][0],comp_cand[0][1],comp_cand[1][0],comp_cand[1][1]]
                for i in range(0,len(spikes_comb)):
                    s=spikes_comb[i]
                    V=[l for l in range(0,curr_N+2)]
                    dt=[0 for l in range(0,curr_N+2)]
                    dt[-1]=1.0
                    print(len(curr_D))
                    D_copy=rec.add_child(curr_D.copy(),parents[s[0]],parents[s[1]],curr_N+1,0.5,dt)
                    print(len(D_copy))
                    print("-----------------------------------------------------")
                    print("SPIKE - " + "("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(curr_N+1)+")")
                    rnf_candidates(rec.rank_candidates_selective(D_copy,V,comp_cand))
                for i in range(0,len(hybr_comb)):
                    s=hybr_comb[i]
                    V=[l for l in range(0,curr_N+2)]
                    dt=[0 for l in range(0,curr_N+2)]
                    dt[-1]=1.0
                    print(len(curr_D))
                    D_copy=rec.add_child(curr_D.copy(),parents[s[0]],parents[s[1]],curr_N+1,0.5,dt)
                    print(len(D_copy))
                    print("-----------------------------------------------------")
                    print("HYBRIDIZATION - " + "("+str(parents[s[0]])+", "+str(parents[s[1]])+": "+str(curr_N+1)+")")
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
                    next_N=curr_N
                    D_copy=curr_D.copy()
                    for i in events:
                        if(i in [0,1,2,3]):
                            s=spikes_comb[i]
                            V=[l for l in range(0,next_N+2)]
                            dt=[0 for l in range(0,next_N+2)]
                            dt[-1]=1.0
                            D_copy=rec.add_child(D_copy.copy(),parents[s[0]],parents[s[1]],next_N+1,1,dt)
                            next_N=next_N+1
                        if(i in [4,5,6,7]):
                            s=hybr_comb[i-4]
                            V=[l for l in range(0,next_N+2)]
                            dt=[0 for l in range(0,next_N+2)]
                            dt[-1]=1.0
                            D_copy=rec.add_child(D_copy.copy(),parents[s[0]],parents[s[1]],next_N+1,0.5,dt)
                            next_N=next_N+1
                    rnf_candidates(rec.rank_candidates_selective(D_copy,V,comp_cand))
            elif("r" in dec_string or "R" in dec_string):
                for i in range(0,len(agree_cand)):
                    print("["+str(i)+"] - " + str(agree_cand[i]))
                r_cand_string=input("Please Enter ONE number (i.e. 0) to remove a candidate pair.\n")
                r_cand=agree_cand[int(r_cand_string.strip())]
                r_alpha=agree_cand_alphas[int(r_cand_string.strip())]
                deltas=rec._compute_deltas(curr_V,curr_D,r_alpha,r_cand[0],r_cand[1],r_cand[2],0)
                curr_D=rec._matrix_without_index(curr_D,curr_V.index(r_cand[2]))
                curr_V=[l for l in range(0,len(curr_D))]
                rec._update_matrix(curr_V,curr_D,r_cand[0],r_cand[1],deltas[2],deltas[3])
                t=rnf_candidates(rec.rank_candidates(curr_D,curr_V))
                agree_cand=t[0]
                agree_cand_alphas=t[1]
                print("Still Pseudometrik? - " + str(rec.is_pseudometric(curr_D)))
                if(len(curr_D)==4):
                    print("Valid 4? - " + str(rec.recognize4_matrix_only(curr_D)))
    else:
        break     
