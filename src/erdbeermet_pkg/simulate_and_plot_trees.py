# -*- coding: utf-8 -*-

import os,sys
from contextlib import redirect_stdout

from erdbeermet.simulation import simulate, load
from erdbeermet.recognition import recognize, is_pseudometric
from erdbeermet.tools.Tree import Tree, TreeNode
from erdbeermet.visualize.BoxGraphVis import plot_box_graph


def classify_divergence(recognition_tree):
    
    valid_divs = []
    invalid_divs = [] 
    
    def classify_children(v):
        if v.Divergence:
            if v.valid_ways:
                valid_divs.append(v.Divergence)
            else:
                invalid_divs.append(v.Divergence)
                
        for c in v.children:
            classify_children(c)
                
    classify_children(recognition_tree.root)
    
    return valid_divs, invalid_divs


valid_mins = [None, None, None]
valid_maxs = [None, None, None]
valid_avrgs = [None, None, None]

invalid_mins = [None, None, None]
invalid_maxs = [None, None, None]
invalid_avrgs = [None, None, None]

prec_string=["0.01", "0.0001", "0.000001"]

j = 0
while True:

  j = j + 1
    
  if j % 100 == 0:
    print(j)
    for i in range(3):
        vavrg = valid_avrgs[i]
        if not vavrg:
            vavrg = -1
        ivavrg = invalid_avrgs[i]
        if not ivavrg:
            ivavrg = -1
        print("Precision", prec_string[i], "Valid", valid_mins[i], valid_maxs[i], "{:.4f}".format(vavrg), "Invalid", invalid_mins[i], invalid_maxs[i], "{:.4f}".format(ivavrg))    

  for i, folder, prec in [ (0, "res_prec2", 0.01), (1, "res_prec4", 0.0001), (2, "res_prec6", 0.000001) ]:

    scenario = simulate(8, delta_min=prec, alpha_tol=prec)
    scenario.write_history('sim_history')
        
    with open('logging.txt', 'w') as f:
        with redirect_stdout(f):

            print("============================================== ALL")
            tree = recognize(D=scenario.D, B=[], first_candidate_only=False, small_spike=False, print_info=True)
            valid_divs, invalid_divs = classify_divergence(tree)
            
            print("============================================== WP3")
            tree3 = recognize(D=scenario.D, B=[0,1,2,3], first_candidate_only=True, small_spike=False, print_info=True)    
            fail_wp3 = True if tree3.root.valid_ways == 0 else False            

            fname = folder+"/Example_"+str(j)
            if fail_wp3:
                fname = fname + "_WP3"
         

    if len(valid_divs) > 0:
        valid_min = min(valid_divs)
        valid_max = max(valid_divs)
        valid_avrg = sum(valid_divs)/float(len(valid_divs))

        fname = fname + "_va" + str(valid_min) + "-" + str(valid_max)

        if not valid_mins[i]:
            valid_mins[i] = valid_min
        else:
            valid_mins[i] = min(valid_min, valid_mins[i])

        if not valid_maxs[i]:
            valid_maxs[i] = valid_max
        else:
            valid_maxs[i] = max(valid_max, valid_maxs[i])   

        if not valid_avrgs[i]:
            valid_avrgs[i] = valid_avrg
        else:
            valid_avrgs[i] = (valid_avrg + valid_avrgs[i])/2.0

    if len(invalid_divs) > 0:
        invalid_min = min(invalid_divs)
        invalid_max = max(invalid_divs)
        invalid_avrg = sum(invalid_divs)/float(len(invalid_divs))

        fname = fname + "_iv" + str(invalid_min) + "-" + str(invalid_max)

        if not invalid_mins[i]:
            invalid_mins[i] = invalid_min
        else:
            invalid_mins[i] = min(invalid_min, invalid_mins[i])

        if not invalid_maxs[i]:
            invalid_maxs[i] = invalid_max
        else:
            invalid_maxs[i] = max(invalid_max, invalid_maxs[i])   

        if not invalid_avrgs[i]:
            invalid_avrgs[i] = invalid_avrg
        else:
            invalid_avrgs[i] = (invalid_avrg + invalid_avrgs[i])/2.0
      
    os.mkdir(fname)
    os.popen('cp sim_history '+fname)
    os.popen('mv logging.txt '+fname)  
    
    tree.visualize(save_as=fname+'/vis_all.pdf')
    tree3.visualize(save_as=fname+'/vis_wp3.pdf')

    with open(folder+'/divergence_log.txt', 'a') as f:
        if len(valid_divs) > 0:
           f.write("Valid"+str(j)+"\t"+ str(len(valid_divs)) +"\t"+str( valid_divs[int(len(valid_divs)/2)] )+"\t"+str(min(valid_divs))+"\t"+str(max(valid_divs))+"\t"+str(valid_divs)+"\n")
        else:
           f.write("Valid"+str(j)+"\t"+ str(len(valid_divs)) +"\t-1\t-1\t-1\t"+str(valid_divs)+"\n")
        if len(invalid_divs) > 0:
           f.write("Invalid"+str(j)+"\t"+ str(len(invalid_divs)) +"\t"+str( invalid_divs[int(len(invalid_divs)/2)] )+"\t"+str(min(invalid_divs))+"\t"+str(max(invalid_divs))+"\t"+str(invalid_divs)+"\n")
        else:
           f.write("Invalid"+str(j)+"\t"+ str(len(invalid_divs)) +"\t-1\t-1\t-1\t"+str(invalid_divs)+"\n")

