
# from erdbeermet_pkg import simulation as sim
from erdbeermet_pkg import recognition as rec
import random
import numpy as np
import copy
import sys
import time
import math
from itertools import permutations, combinations
from sys import argv

def rev_a_val(D,a,x,y,z,u):
    return((-a*D[x,z]+a*D[z,y]+a*D[x,y]+D[u,z]-D[z,y]-D[u,y])/(2*a*D[x,y]+D[u,x]-D[u,y]-D[x,y]))

def generate_random_scenario(n):
    # init output and alphabet
    scen_string = ""
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    already_added = []

    # add participating letters to first line
    for i in range(0,n):
        scen_string+=alphabet[i] + ","

    # remove last comma, add starting element (a)
    scen_string=scen_string[0:-1]
    scen_string += "\n"
    scen_string += alphabet[0]
    scen_string += "\n"

    # the first step is always direct branching of a (random deltas)
    # .0 is added at the end since otherwise float point calculations fail
    scen_string += "a,a,b;" + str(round(random.random(),5)) + "," + str(random.randint(0,10)) + ".0," + str(random.randint(0,10)) + ".0," + str(random.randint(0,10)) + ".0\n"
    already_added += ["a"]
    already_added += ["b"]

    # randomly draw parents, add new child letters,
    # parameters are generated randomly.
    for i in range(2,n):

        p_0 = already_added[random.randint(0,len(already_added)-1)]
        p_1 = copy.copy(p_0)

        while p_1 == p_0:
            p_1 = already_added[random.randint(0,len(already_added)-1)]

        c = alphabet[i]

        alpha = round(random.random(),5)
        del_0 = random.randint(0,10)
        del_1 = random.randint(0,10)
        del_2 = random.randint(0,10)

        scen_string += p_0 + "," + p_1 + "," + c + ";"+ str(alpha) + "," + str(del_0) + ".0," + str(del_1) + ".0," + str(del_2) + ".0\n"
        already_added += [c]
    # print(scen_string[0:-1])
    return scen_string[0:-1]

def generate_random_scenario(n):

    scen_string = ""
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    already_added = []

    for i in range(0,n):
        scen_string+=alphabet[i] + ","
    scen_string=scen_string[0:-1]
    scen_string += "\n"
    scen_string += alphabet[0]
    scen_string += "\n"

    scen_string += "a,a,b;" + str(round(random.random(),5)) + "," + str(random.randint(0,10)) + ".0," + str(random.randint(0,10)) + ".0," + str(random.randint(0,10)) + ".0\n"
    already_added += ["a"]
    already_added += ["b"]

    for i in range(2,n):
        p_0 = already_added[random.randint(0,len(already_added)-1)]
        p_1 = copy.copy(p_0)
        while p_1 == p_0:
            p_1 = already_added[random.randint(0,len(already_added)-1)]
        c = alphabet[i]
        alpha = round(random.random(),5)
        del_0 = round(round(random.random(),5)*random.randint(0,10),5)
        del_1 = round(round(random.random(),5)*random.randint(0,10),5)
        del_2 = round(round(random.random(),5)*random.randint(0,10),5)
        # del_0=0
        # del_1=0
        # scen_string += p_0 + "," + p_1 + "," + c + ";"+ str(alpha) + "," + str(del_0) + ".0," + str(del_1) + ".0," + str(del_2) + ".0\n"
        scen_string += p_0 + "," + p_1 + "," + c + ";"+ str(alpha) + "," + str(del_0) + "," + str(del_1) + "," + str(del_2) + "\n"

        already_added += [c]
    return scen_string


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
                # agree=False
                inf_counter+=1
            elif(-math.inf==sorted_c[i][1]):
                # agree=False
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
            if last_comp[0]=="notset":
                # if math.isnan(sorted_c[i][1]) or math.inf == sorted_c[i][1] or -math.inf==sorted_c[i][1] or sorted_c[i][1]==-0.0:
                #     continue
                agree_cons=sorted_c[i][1]
                last_comp=(sorted_c[i][1],1)
            # elif math.isnan(sorted_c[i][1]) or math.inf == sorted_c[i][1] or -math.inf==sorted_c[i][1] or sorted_c[i][1]==-0.0:
            #     continue
            elif np.isclose(last_comp[0],sorted_c[i][1]):
                last_comp=(sorted_c[i][1],last_comp[1]+1)
            else:
                agree=False
                print_str+= str(round(last_comp[0],7)) + ":" + str(last_comp[1]) + " // "
                last_comp=(sorted_c[i][1],1)
        if(last_comp[0]!="notset"):
            print_str+= str(round(last_comp[0],7)) + ":" + str(last_comp[1])
        if(len(sorted_c)==0):
            # print(print_str)
            continue
        else:
            if(agree and 0<=agree_cons<=1):
                agree_cand+=[k[0]]
                agree_cons_arr+=[agree_cons]
        # print(print_str)
    else:
        return (agree_cand,agree_cons_arr)

def make_distance_matrix(scenario):

    lines = scenario.split("\n")
    dist=[[0 for x in range(0,len(lines[0].split(",")))] for x in range(0,len(lines[0].split(",")))]
    nodes = [x for x in range(0,len(dist))]
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    seen = set()
    seen.add(0)

    for t in range(2,len(lines)-1):
        # Parse part of each insert line
        insert=lines[t].split(";")[0]
        params=lines[t].split(";")[1]

        # Parse parents and child instruction and transl
        lparent=nodes[alphabet.index(insert.split(",")[0])]
        rparent=nodes[alphabet.index(insert.split(",")[1])]
        child=nodes[alphabet.index(insert.split(",")[2])]

        seen.add(child)

        if(lparent not in seen or rparent not in seen):
            print(lparent)
            print(rparent)
            print(lines[t])
            print("Instruction in LINE " + str(t) + " contains parent nodes that do not exist!")
            sys.exit()
        alpha=np.float(params.split(",")[0])
        dl=np.float(params.split(",")[1])
        dr=np.float(params.split(",")[2])
        dc=np.float(params.split(",")[3])

        #R1
        for x in range(0,len(dist)):
            if(x==child):
                continue
            else:
                dist[x][child]=dist[x][lparent]*alpha + dist[x][rparent]*(1-alpha)
                dist[child][x]=dist[x][lparent]*alpha + dist[x][rparent]*(1-alpha)
        #R2
        for x in range(0,len(dist)):
           for y in range(x+1,len(dist)):
               if(y in seen and x in seen):
                   if(x==lparent or y==lparent):
                       dist[x][y]+=dl
                       dist[y][x]+=dl
                   if(x==rparent or y==rparent):
                       dist[x][y]+=dr
                       dist[y][x]+=dr
                   if(x==child or y==child):
                       dist[x][y]+=dc
                       dist[y][x]+=dc
    return dist

def _compute_alpha_xoff(V, D, x, y, z, u, v):

    x = V.index(x)
    y = V.index(y)
    z = V.index(z)
    u = V.index(u)
    v = V.index(v)

    numerator   = (D[u,z] + D[v,y]) - (D[v,z] + D[u,y])
    denominator = (D[u,x] + D[v,y]) - (D[v,x] + D[u,y])

    # print(f'({x},{y}:{z}),u={u},v={v}:   ({D[u,z]} + {D[v,y]}) - ({D[v,z]} + {D[u,y]})  // ({D[u,x]} + {D[v,y]}) - ({D[v,x]} + {D[u,y]}) = {numerator} // {denominator} = {numerator / denominator}')

    if not np.isclose(denominator, 0.0):#, rtol=glob_rtol, atol=glob_atol):
        if numerator/denominator==-0.0:
            return -0.0
        return numerator / denominator
    else:
        if numerator/denominator==math.inf:
            return math.inf
        elif numerator/denominator==-math.inf:
            return -math.inf
        return np.nan

def check_triangle(distance, x,y,z, competing_pairs,print_info=True):
    if print_info:
        print(str(x) + "," + str(y) + ":" + str(z) +" VS " + str(competing_pairs))
    p_0 = x
    p_1 = y

    if p_0 > z:
        p_0 -= 1
    if p_1 > z:
        p_1 -= 1

    # check if we hurt triangle condition of competing pairs
    for pair in competing_pairs:

        pair_0 = pair[0]
        pair_1 = pair[1]
        pair_c = pair[2]

        check_0_1 = pair_0 != z and pair_1 != z
        check_0_c = pair_0 != z and pair_c != z
        check_1_c = pair_1 != z and pair_c != z

        # print("--")
        # print(pair)
        # print(check_0_1)
        # print(check_0_c)
        # print(check_1_c)
        # print("--")

        if pair_0 > z:
            pair_0-=1
        if pair_1 > z:
            pair_1-=1
        if pair_c > z:
            pair_c-=1

        if check_0_1:
            # if distance[p_0,p_1] > distance[p_0,pair_0] + distance[pair_0,p_1]\
            # or  distance[p_0,p_1] > distance[p_0,pair_1] + distance[pair_1,p_1]:
            #     print("REVERSE CHECK FIRED------------------------------------------------->")
            #     return False

            # if print_info:
            #     print("TRIANGLE HURT " + str(pair_0) + ", " + str(pair_1) + " over " + str(p_0))
            #     print(str(distance[pair_0,pair_1]))
            #     print(str(distance[pair_0,p_0] + distance[p_0,pair_1]))
            # if print_info:
            #     print("TRIANGLE HURT " + str(pair_0) + ", " + str(pair_1) + " over " + str(p_1))
            #     print(str(distance[pair_0,pair_1]))
            #     print(str(distance[pair_0,p_1] + distance[p_1,pair_1]))
            if distance[pair_0,pair_1] > distance[pair_0,p_0] + distance[p_0,pair_1]\
            and not np.isclose(distance[pair_0,pair_1],distance[pair_0,p_0] + distance[p_0,pair_1]):
                if print_info:
                    print("TRIANGLE HURT " + str(pair_0) + ", " + str(pair_1) + " over " + str(p_0))
                    print(str(distance[pair_0,pair_1]))
                    print(str(distance[pair_0,p_0] + distance[p_0,pair_1]))
                return False

            if distance[pair_0,pair_1] > distance[pair_0,p_1] + distance[p_1,pair_1]\
            and not np.isclose(distance[pair_0,pair_1],distance[pair_0,p_1] + distance[p_1,pair_1]):
                if print_info:
                    print("TRIANGLE HURT " + str(pair_0) + ", " + str(pair_1) + " over " + str(p_1))
                    print(str(distance[pair_0,pair_1]))
                    print(str(distance[pair_0,p_1] + distance[p_1,pair_1]))
                return False

            for i in range(0,len(distance)):
                if i in [p_0,p_1,pair_0,pair_1,pair_c]:
                    continue
                else:
                    if distance[i,pair_0] > distance[i,p_0] + distance[p_0,pair_0] or\
                    distance[i,pair_0] > distance[i,p_1] + distance[p_1,pair_0] or\
                    distance[i,pair_1] > distance[i,p_0] + distance[p_0,pair_1] or\
                    distance[i,pair_1] > distance[i,p_1] + distance[p_1,pair_1]:
                        if print_info:
                            print("TRIANGLE HURT")
                        return False
        if check_0_c:
            # if distance[p_0,p_1] > distance[p_0,pair_0] + distance[pair_0,p_1]\
            # or  distance[p_0,p_1] > distance[p_0,pair_c] + distance[pair_c,p_1]:
            #     print("REVERSE CHECK FIRED------------------------------------------------->")
            #     return False
            # if print_info:
            #     print("TRIANGLE HURT " + str(pair_0) + ", " + str(pair_c) + " over " + str(p_0))
            #     print(str(distance[pair_0,pair_c]))
            #     print(str(distance[pair_0,p_0] + distance[p_0,pair_c]))
            # if print_info:
            #     print("TRIANGLE HURT " + str(pair_0) + ", " + str(pair_c) + " over " + str(p_1))
            #     print(str(distance[pair_0,pair_c]))
            #     print(str(distance[pair_0,p_1] + distance[p_1,pair_c]))
            if distance[pair_0,pair_c] > distance[pair_0,p_0] + distance[p_0,pair_c]\
            and not np.isclose(distance[pair_0,pair_c],distance[pair_0,p_0] + distance[p_0,pair_c]):
                if print_info:
                    print("TRIANGLE HURT " + str(pair_0) + ", " + str(pair_c) + " over " + str(p_0))
                    print(str(distance[pair_0,pair_c]))
                    print(str(distance[pair_0,p_0] + distance[p_0,pair_c]))
                return False
            if distance[pair_0,pair_c] > distance[pair_0,p_1] + distance[p_1,pair_c]\
            and not np.isclose(distance[pair_0,pair_c],distance[pair_0,p_1] + distance[p_1,pair_c]):
                if print_info:
                    print("TRIANGLE HURT " + str(pair_0) + ", " + str(pair_c) + " over " + str(p_1))
                    print(str(distance[pair_0,pair_c]))
                    print(str(distance[pair_0,p_1] + distance[p_1,pair_c]))
                return False

            for i in range(0,len(distance)):
                if i in [p_0,p_1,pair_0,pair_1,pair_c]:
                    continue
                else:
                    if distance[i,pair_0] > distance[i,p_0] + distance[p_0,pair_0] or\
                    distance[i,pair_0] > distance[i,p_1] + distance[p_1,pair_0] or\
                    distance[i,pair_c] > distance[i,p_0] + distance[p_0,pair_c] or\
                    distance[i,pair_c] > distance[i,p_1] + distance[p_1,pair_c]:
                        if print_info:
                            print("TRIANGLE HURT")
                        return False
        if check_1_c:
            # if distance[p_0,p_1] > distance[p_0,pair_1] + distance[pair_1,p_1]\
            # or  distance[p_0,p_1] > distance[p_0,pair_c] + distance[pair_c,p_1]:
            #     print("REVERSE CHECK FIRED------------------------------------------------->")
            #     return False
            # if print_info:
            #     print("TRIANGLE HURT " + str(pair_1) + ", " + str(pair_c) + " over " + str(p_0))
            #     print(str(distance[pair_1,pair_c]))
            #     print(str(distance[pair_1,p_0] + distance[p_0,pair_c]))
            # if print_info:
            #     print("TRIANGLE HURT " + str(pair_1) + ", " + str(pair_c) + " over " + str(p_1))
            #     print(str(distance[pair_1,pair_c]))
            #     print(str(distance[pair_1,p_1] + distance[p_1,pair_c]))
            if distance[pair_1,pair_c] > distance[pair_1,p_0] + distance[p_0,pair_c]\
            and not np.isclose(distance[pair_1,pair_c],distance[pair_1,p_0] + distance[p_0,pair_c]):
                if print_info:
                    print("TRIANGLE HURT " + str(pair_1) + ", " + str(pair_c) + " over " + str(p_0))
                    print(str(distance[pair_1,pair_c]))
                    print(str(distance[pair_1,p_0] + distance[p_0,pair_c]))
                return False
            if distance[pair_1,pair_c] > distance[pair_1,p_1] + distance[p_1,pair_c]\
            and not np.isclose(distance[pair_1,pair_c],distance[pair_1,p_1] + distance[p_1,pair_c]):
                if print_info:
                    print("TRIANGLE HURT " + str(pair_1) + ", " + str(pair_c) + " over " + str(p_1))
                    print(str(distance[pair_1,pair_c]))
                    print(str(distance[pair_1,p_1] + distance[p_1,pair_c]))
                return False
            for i in range(0,len(distance)):
                if i in [p_0,p_1,pair_0,pair_1,pair_c]:
                    continue
                else:
                    if distance[i,pair_1] > distance[i,p_0] + distance[p_0,pair_1] or\
                    distance[i,pair_1] > distance[i,p_1] + distance[p_1,pair_1] or\
                    distance[i,pair_c] > distance[i,p_0] + distance[p_0,pair_c] or\
                    distance[i,pair_c] > distance[i,p_1] + distance[p_1,pair_c]:
                        if print_info:
                            print("TRIANGLE HURT")
                        return False
    # print("TRIANGLE HOLDS!")
    return True


def check_candidate(distance,x,y,z,alpha,alphas,print_info=False,competing_pairs=[]):

    nodes = [i for i in range(0,len(distance)+1)]
    n = len(nodes)

    mock_children = []

    for c in range(0,n):
        if(c not in [x,y,z]):
            mock_children += [c]

    # print(mock_children[0])
    # print("^")

    deltas = rec._compute_deltas(nodes,distance,alpha,x,y,z,mock_children[0])
    if (deltas[1] < 0 and not np.isclose(deltas[1],0.0) ) or (deltas[2] < 0 and not np.isclose(deltas[2],0.0) ) or (deltas[3] < 0 and not np.isclose(deltas[3],0.0) ) :
        # print(deltas)
        return False

    distance = rec._update_matrix_return(nodes,distance,x,y,deltas[2],deltas[3])
    distance = rec._matrix_without_index(distance,nodes.index(z))
    #
    # print("JUST BEFORE CHECK TRIANGLE")
    # print(len(competing_pairs))
    # print(len(distance))
    if len(competing_pairs) != 0 and len(distance) > 4:
        # print("TRY CHECK TRIANGLE")
        check = check_triangle(distance,x,y,z,competing_pairs,print_info=print_info)
        if check:
            pass
        else:
            # print("TRIANGLE HURT")
            return False

    nodes = nodes[0:-1]
    # CHECK WHY WE REMOVE THIS MOCK CHILD
    mock_children = mock_children[0:-1]

    witnesses = []

    for node in nodes:
        if node not in [x,y,z]:
            witnesses += [node]

    p_0 = x
    p_1 = y

    if p_0 > z:
        p_0 -= 1
    if p_1 > z:
        p_1 -= 1

    for mock_child in mock_children:

        mock_child_adjusted = mock_child

        if mock_child > z:
            mock_child_adjusted -= 1

        for w_0 in witnesses:
            if w_0 != mock_child:
                w_0_adjusted = w_0
                if w_0_adjusted > z:
                    w_0_adjusted -= 1

                check_alpha = rev_a_val(distance,alpha,p_0,p_1,mock_child_adjusted,w_0_adjusted)

                chk_1=min(w_0,z)
                chk_2=max(w_0,z)

                # numerically print all the compare values
                # if print_info:
                #     print(str(check_alpha) + " <-> " + str(alphas[ (x,y,mock_child) ][(chk_1,chk_2)]))

                if np.isclose(alphas[ (x,y,mock_child) ][(chk_1,chk_2)],check_alpha):
                    continue
                else:
                    return False
    return True

def check_spikes(distance,x,y,z):
    distance_update = [[distance[0,0] for x in range(0,len(distance)+2)] for x in range(0,len(distance)+2)]
    distance_update = np.asarray(distance_update)

    for i in range(0,len(distance)):
        for j in range(0,len(distance)):
            distance_update[i,j] = distance[i,j]

    for j in range(0,len(distance)):
        distance_update[j,len(distance_update)-2] = distance_update[j,x]
        distance_update[len(distance_update)-2,j] = distance_update[j,x]

        distance_update[j,len(distance_update)-1] = distance_update[j,y]
        distance_update[len(distance_update)-1,j] = distance_update[j,y]

    distance_update[len(distance_update)-2,len(distance_update)-1] = distance_update[x, y]
    distance_update[len(distance_update)-1,len(distance_update)-2] = distance_update[x, y]

    # print(distance_update)

    nodes = [i for i in range(0,len(distance_update))]
    t = [(x,y,z)]

    candidates = rec.rank_candidates_selective(distance_update,nodes,t,print_info=False)
    k = list(candidates.keys())[0]
    found_one_valid = False
    # print("CHECK SPIKES " + str(x) + "," + str(y) + ":" + str(z))
    for a in candidates[k]:
        if candidates[k][a] < 1 or candidates[k][a] > 0 or np.isclose(candidates[k][a],0.0) or np.isclose(candidates[k][a],1.0):
            found_one_valid = True
        if not (math.isnan(candidates[k][a]) or math.inf == candidates[k][a] or -math.inf==candidates[k][a] or candidates[k][a]==-0.0):
            if (candidates[k][a] > 1 or candidates[k][a] < 0) and not (np.isclose(candidates[k][a],0.0) or np.isclose(candidates[k][a],1.0)):
                # print(False)
                return False
    return True
    # ranking = rnf_candidates(rec.rank_candidates_selective(distance_update,nodes,t),evaluate_spikes=True)

    # print(ranking)


def rank_candidates(D,V):
    candidates={}
    for x, y, z in permutations(V, 3):
        # considering x < y suffices
        if x > y:
            continue
        candidates[(x,y,z)]={}
        for u, v in combinations(V, 2):
            #if(x==1 and z==5):
            #    print(str(x)+","+str(x)+":"+str(z)+" - "+str(_compute_alpha(V, D, x, x, z, u, v)))
            if u in  (x, y, z) or v in (x, y, z):
                continue
            candidates[(x,y,z)][(u,v)] = _compute_alpha_xoff(V, D, x, y, z, u, v)
           # if(math.isnan(candidates[(x,y,z)][(u,v)])):
           #     candidates[(x,y,z)][(u,v)] = 1
    # for x in range(0,len(D)):
    #     for z in range(0,len(D)):
    #         if(x==z):
    #             continue
    #         else:
    #             candidates[(x,x,z)]={}
    #             for u, v in combinations(V, 2):
    #                 if u in (x,z) or v in (x,z):
    #                     continue
    #                 else:
    #                     candidates[(x,x,z)][(u,v)] = _compute_alpha_xoff(V, D, x, x, z, u, v)
    return candidates
def is_pseudometric(D, print_info=False, V=None,
                    return_info=False):
    """Check whether a given distance matrix is a pseudometric.

    Parameters
    ----------
    D : 2-dimensional numpy array
        Distance matrix
    rtol : float, optional
        Relative tolerance for equality. The default is 1e-05.
    atol : float, optional
        Absolute tolerance for equality. The default is 1e-08.
    print_info : bool, optional
        If True, print first encountered violation of the triangle inequality
        if any.
    V : list, optional
        List of items (used for info output).
    return_info : bool, optional
        If True, return an info string as a second return value. The default
        is False.

    Return
    ------
    bool or tuple of bool and str
        True if D is a pseudometric and optionally an info string.
    """

    N = D.shape[0]

    # check whether all entries are non-negative
    if not np.all(np.logical_or(np.isclose(D, 0.0),
                                D > 0.0)):
        return False if not return_info else (False, 'negative distances')

    # check whether all diagonal entries are zero
    if np.any(np.diagonal(D)):
        return False if not return_info else (False, 'non-zero diagonal')

    # check whether the matrix is symmetric
    if not np.allclose(D, D.T):
        return False if not return_info else (False, 'not symmetric')

    # check the triangle inequality
    for i in range(N-1):
        for j in range(i+1, N):
            minimum = np.min(D[i, :] + D[:, j])
            if minimum < D[i, j] and not np.isclose(minimum, D[i, j]):
                if print_info or return_info:
                    argmin = np.argmin(D[i, :] + D[:, j])
                    if not V:
                        info = f'triangle inequality violation: D[{i},'\
                               f'{j}]={D[i,j]} > {minimum} over {argmin}'
                    else:
                        info = f'triangle inequality violation: D[v{V[i]},'\
                               f'v{V[j]}]={D[i,j]} > {minimum} over v{V[argmin]}'
                        if print_info:
                            print(info)
                return False if not return_info else (False, info)

    return True if not return_info else (True, 'passed')

def solve_greedy(n,out,infile="",report_dead_ends=False,print_info=True):

    # # TODO: CHECK EXAMPLE_2, SPIKE TEST FAILS

    if(infile==""):
        scenario = generate_random_scenario(n)
    else:
        with open(infile,"r") as file:
            scenario = file.read()

    distance = make_distance_matrix(scenario)
    nodes = [x for x in range(0,len(distance))]

    distance = np.asarray(distance)
    swap_distance = distance.copy()

    while len(distance) > 4:
        if(len(distance) == 5):
            swap_distance = distance.copy()

        all_alphas = rank_candidates(distance,nodes)

        agree_cand, agree_cand_alphas = rnf_candidates(all_alphas)

        valid = []

        # print(agree_cand)

        for i in range(0,len(agree_cand)):
            candidate = agree_cand[i]

            competing_pairs = []
            for c in agree_cand:
                if c!=candidate:
                    competing_pairs += [c]

            # if len(distance) > 5 and print_info:
                # print(candidate)
                # print("VS")
                # print(competing_pairs)


            chk_cand = check_candidate(distance.copy(),candidate[0],candidate[1],candidate[2],
                            agree_cand_alphas[i], all_alphas, competing_pairs=competing_pairs, print_info=print_info)
            chk_spks = check_spikes(distance.copy(),candidate[0],candidate[1],candidate[2])
            if print_info:
                print("CANDIDATE - " + str(candidate))
                print("CHECK CANDIDATES - " + str(chk_cand))
                print("CHECK SPIKES - " + str(chk_spks))

            if chk_cand and chk_spks:
                valid += [i]
        # print(valid)
        if len(valid) == 0:
            if print_info:
                print("DEAD END! NO VALID CANDIDATES")
            if report_dead_ends:
                with open(out,"w+") as file:
                    file.write(scenario)
                return (False,"_NO_VALID_ALPHAS")
            else:
                return (True,"")
        else:
            valid_index = random.randint(0,len(valid)-1)

            t = copy.copy(agree_cand[valid[valid_index]])
            alpha = agree_cand_alphas[valid[valid_index]]
            witness = -1


            for node in nodes:
                if node != t[0] and node != t[1] and node != t[2]:
                    witness = node
                    break

            if print_info:
                print(str(len(distance)) + " LEFT")
                print("REMOVING " + str(t))
                print("WITNESS " + str(witness))

            if witness == -1:
                if print_info:
                    print("COULD NOT FIND WITNESS!")
                with open(out,"w+") as file:
                    file.write(scenario)
                return (False, "_NO_WITNESS")

            deltas = rec._compute_deltas(nodes.copy(),distance.copy(),alpha,t[0],t[1],t[2],witness)

            distance = rec._matrix_without_index(distance.copy(),nodes.copy().index(t[2]))

            nodes = nodes[0:-1]

            lparent = t[0]
            rparent = t[1]
            if lparent > t[2]:
                lparent -= 1
            if rparent > t[2]:
                rparent -= 1

            distance = rec._update_matrix_return(nodes.copy(),distance.copy(),lparent,rparent,deltas[2],deltas[3])

            still_metric = rec.is_pseudometric(distance.copy(),return_info=True)

            if(len(distance)==4):
                chk = rec.recognize4_matrix_only(distance.copy())
                if print_info:
                    print("METRIC ON 4? - " + str(chk))
                if(chk):
                    if print_info:
                        print("DONE!")
                else:
                    for i in range(0,len(valid)):

                        distance = swap_distance.copy()

                        nodes = [x for x in range(0,len(distance))]

                        t = agree_cand[valid[i]]
                        alpha = agree_cand_alphas[valid[i]]

                        for node in nodes:
                            if node != t[0] and node != t[1] and node != t[2]:
                                witness = node
                                break

                        if print_info:
                            print("REMOVING " + str(t))
                            print("WITNESS " + str(witness))

                        deltas = rec._compute_deltas(nodes.copy(),distance.copy(),alpha,t[0],t[1],t[2],witness)

                        distance = rec._matrix_without_index(distance.copy(),nodes.copy().index(t[2]))

                        nodes = nodes[0:-1]

                        lparent = t[0]
                        rparent = t[1]
                        if lparent > t[2]:
                            lparent -= 1
                        if rparent > t[2]:
                            rparent -= 1

                        distance = rec._update_matrix_return(nodes.copy(),distance.copy(),lparent,rparent,deltas[2],deltas[3])

                        # print(rec.recognize4_matrix_only(distance.copy()))

                        if rec.recognize4_matrix_only(distance.copy()):
                            if print_info:
                                print("STILL PSEUDOMETRIC")
                            return (True,"")
                    with open(out,"w+") as file:
                        file.write(scenario)
                    return (False, "_LAST_4")

            if(still_metric[0]):
                if print_info:
                    print("STILL PSEUDOMETRIC")
                    print()
                continue
            else:
                if print_info:
                    print("DEAD END!")
                    print(still_metric)
                with open(out,"w+") as file:
                    file.write(scenario)
                return (False, "_NO_PSEUDOMETRIC" + still_metric[1])
    return (True,"")

if __name__ == '__main__':

    base_string = "EXAMPLE_"
    invalid_counter = 0

    general_counter = 0

    path = argv[1]

    if path[-1] == "/":
        pass
    else:
        path += "/"

    while (True):
        print("RUN " + str(general_counter))
        general_counter += 1
        check = solve_greedy(10,path+base_string+str(invalid_counter)+".txt",print_info=False,report_dead_ends=False)
        if check[0]:
            pass
        else:
            invalid_counter += 1
        if(general_counter > 100000):
            break
