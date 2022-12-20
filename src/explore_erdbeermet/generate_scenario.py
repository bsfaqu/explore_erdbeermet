import random
import sys
import copy
import math

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

    scen_string += "a,a,b;" + str(round(random.random(),5)) + "," + str(round(round(random.random(),5)*random.randint(0,10),5)) + ".0," + str(round(round(random.random(),5)*random.randint(0,10),5)) + ".0," + str(round(round(random.random(),5)*random.randint(0,10),5)) + ".0\n"
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
        scen_string += p_0 + "," + p_1 + "," + c + ";"+ str(alpha) + "," + str(del_0) + ".0," + str(del_1) + ".0," + str(del_2) + ".0\n"
        already_added += [c]
    print(scen_string[0:-1])
    return scen_string

# print(sys.argv)
generate_random_scenario(int(sys.argv[1]))
