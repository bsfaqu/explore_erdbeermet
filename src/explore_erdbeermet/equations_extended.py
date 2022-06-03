from sympy import *
from sys import argv
import os
import pathlib
from texttable import Texttable
from itertools import combinations,permutations
import time
import sys

# Simply performs the alpha calculation
def equatify(a,b,c,d,e,f,g,h):
    return ( ((a+b)-(c+d)) / ((e+f)-(g+h)) )

# pretty print matrix with col and row headers
def ppmat(matrix,seen):
    col_width=[0 for x in range(0,len(matrix))]
    header = [""]+seen
    t = Texttable()
    t.add_row(header)
    for x in range(0,len(matrix)):
        for y in range(0,len(matrix)):
            col_width[y]=max(col_width[y],len(str(matrix[x][y])))
        ins_row=[seen[x]]+matrix[x]
        print(ins_row)
        t.add_row(ins_row)
    col_width=[1]+col_width
    t.set_cols_width(col_width)
    print(t.draw())

# calculate alpha and output some info
def get_alpha(matrix,seen,x,y,z,u,v,print_info):
    if(print_info):
        print("alpha: ("+seen[x]+","+seen[y]+":"+seen[z]+") - WITNESS "+ seen[u] +","+seen[v] )
        print()
    return( ( (matrix[u][z]+matrix[v][y])-( matrix[v][z]+matrix[u][y] ) )/( (matrix[u][x]+matrix[v][y] )-(matrix[v][x]+matrix[u][y]) ) )

# read scenario from argv1
with open(argv[1],"r") as f:
    lines=f.read().split("\n")

# init nodes, labels, visited nodes and distance matrix
nodes={}
labels=[]    
in_nodes=lines[0].split(",")
start_node=lines[1]
seen = []
seen+=[start_node]
iterations=1

# matrix stores console output symbols
matrix=[[0 for x in range(0,len(in_nodes))] for y in range(0,len(in_nodes))]
matrix_update=[[0 for x in range(0,len(in_nodes))] for y in range(0,len(in_nodes))]

# lmatrix stores latex code for pdf output
lmatrix=[[0 for x in range(0,len(in_nodes))] for y in range(0,len(in_nodes))]
lmatrix_update=[[0 for x in range(0,len(in_nodes))] for y in range(0,len(in_nodes))]

# preprocess nodes well encounter
for i in range(0,len(in_nodes)):
    nodes[in_nodes[i]]=i
    labels+=[in_nodes[i]]

# read scenario and compute distances
for t in range(2,len(lines)-1):
    #print(lines[t])
    iterations+=1
    # Parse part of each insert line
    insert=lines[t].split(";")[0]
    params=lines[t].split(";")[1]
    
    # Parse parents and child instruction
    lparent=insert.split(",")[0]
    rparent=insert.split(",")[1]
    child=insert.split(",")[2]
    
    # Append new node to seen
    seen += [child]
    
    # Initialize matrix with "base distance" d0xy
    # For the root once we reached iteration 4
    if(iterations==4):
        for x in range(0,len(seen)):
            for y in range(0,len(seen)):
                if x==y:
                    continue
                else:
                    # console output formulas
                    exec("d0"+seen[x]+seen[y]+" = symbols("+"\"d0"+seen[x]+seen[y]+"\")")
                    exec("d0"+seen[y]+seen[x]+" = symbols("+"\"d0"+seen[y]+seen[x]+"\")")
                    matrix[x][y]=eval("d0"+seen[y]+seen[x])
                    matrix[y][x]=eval("d0"+seen[y]+seen[x])
                    
                    # latex formulas
                    exec("ld0"+seen[x]+seen[y]+" = symbols("+"\"d^{\\scriptstyle0}_{"+seen[x]+""+seen[y]+"}\")")
                    exec("ld0"+seen[y]+seen[x]+" = symbols("+"\"d^{\\scriptstyle0}_{"+seen[y]+""+seen[x]+"}\")")
                    lmatrix[x][y]=eval("ld0"+seen[y]+seen[x]+"")
                    lmatrix[y][x]=eval("ld0"+seen[y]+seen[x]+"")
                    
    # Update the matrix when were beyond the root construction
    elif(iterations>4):
        # Save updated values in updated matrix
        matrix_update=matrix.copy()
        lmatrix_update=lmatrix.copy()
        
        ####### CONSOLE OUT #########
         # Declare symbolic alpha for corresponding R-Step
        exec("a"+str(iterations-4)+"=symbols("+"\"a"+str(iterations-4)+"\")")
        
        # Declare symbolic deltas for corresponding R-Step
        exec("del"+str(iterations-4)+"_"+lparent+"=symbols("+"\"del"+str(iterations-4)+"_"+lparent+"\")")
        exec("del"+str(iterations-4)+"_"+rparent+"=symbols("+"\"del"+str(iterations-4)+"_"+rparent+"\")")
        exec("del"+str(iterations-4)+"_"+child+"=symbols("+"\"del"+str(iterations-4)+"_"+child+"\")")
        
        
        ####### LATEX OUT #########
        # Declare symbolic alpha for corresponding R-Step
        exec("la"+str(iterations-4)+"=symbols("+"\"alpha_"+str(iterations-4)+"\")")
        
        # Declare symbolic deltas for corresponding R-Step
        exec("ldel"+str(iterations-4)+"_"+lparent+"=symbols("+"\"\\delta^"+str(iterations-4)+"_{"+lparent+"}\")")
        exec("ldel"+str(iterations-4)+"_"+rparent+"=symbols("+"\"\\delta^"+str(iterations-4)+"_{"+rparent+"}\")")
        exec("ldel"+str(iterations-4)+"_"+child+"=symbols("+"\"\\delta^"+str(iterations-4)+"_{"+child+"}\")")
        
        #R1 - hybridization
        for x in range(0,len(seen)):
            if(seen[x]==child):
                continue
            else:
                # update child distances
                ####### CONSOLE OUT #########
                matrix_update[x][nodes[child]]=matrix[x][nodes[lparent]]*eval("a"+str(iterations-4)) + matrix[x][nodes[rparent]]*(1-eval("a"+str(iterations-4)))
                matrix_update[nodes[child]][x]=matrix[x][nodes[lparent]]*eval("a"+str(iterations-4)) + matrix[x][nodes[rparent]]*(1-eval("a"+str(iterations-4)))
                
                ####### LATEX OUT #########
                lmatrix_update[x][nodes[child]]=lmatrix[x][nodes[lparent]]*eval("la"+str(iterations-4)) + lmatrix[x][nodes[rparent]]*(1-eval("la"+str(iterations-4)))
                lmatrix_update[nodes[child]][x]=lmatrix[x][nodes[lparent]]*eval("la"+str(iterations-4)) + lmatrix[x][nodes[rparent]]*(1-eval("la"+str(iterations-4)))
        
        #R2 - mutation
        for x in range(0,len(seen)):
            for y in range(x+1,len(seen)):
                if(x==y):
                    continue
                if(seen[x]==lparent or seen[y]==lparent):
                    ####### CONSOLE OUT #########
                    matrix_update[x][y]+=eval("del"+str(iterations-4)+"_"+lparent)
                    matrix_update[y][x]+=eval("del"+str(iterations-4)+"_"+lparent)
                    
                    ####### LATEX OUT #########
                    lmatrix_update[x][y]+=eval("ldel"+str(iterations-4)+"_"+lparent)
                    lmatrix_update[y][x]+=eval("ldel"+str(iterations-4)+"_"+lparent)
                if(seen[x]==rparent or seen[y]==rparent):
                    ####### CONSOLE OUT #########
                    matrix_update[x][y]+=eval("del"+str(iterations-4)+"_"+rparent)
                    matrix_update[y][x]+=eval("del"+str(iterations-4)+"_"+rparent)
                    
                    ####### LATEX OUT #########
                    lmatrix_update[x][y]+=eval("ldel"+str(iterations-4)+"_"+rparent)
                    lmatrix_update[y][x]+=eval("ldel"+str(iterations-4)+"_"+rparent)
                if(seen[x]==child or seen[y]==child):
                    ####### CONSOLE OUT #########
                    matrix_update[x][y]+=eval("del"+str(iterations-4)+"_"+child)
                    matrix_update[y][x]+=eval("del"+str(iterations-4)+"_"+child)
                    
                    ####### LATEX OUT #########
                    lmatrix_update[x][y]+=eval("ldel"+str(iterations-4)+"_"+child)
                    lmatrix_update[y][x]+=eval("ldel"+str(iterations-4)+"_"+child)
                    
        # Finalize the update
        matrix=matrix_update.copy()
        lmatrix=lmatrix_update.copy()
        
# Look at these pretty matrices
ppmat(matrix.copy(),seen)
ppmat(lmatrix.copy(),seen)

indices=[x for x in range(0,len(matrix))]
counter=0

# LaTeX preamble and table init
preamble="\\documentclass[12pt]{article}\\usepackage{amsmath,amsfonts}\n\\usepackage{longtable}\n\\usepackage{enumerate}\n\\begin{document}\n"
latex_str="\\renewcommand*{\\arraystretch}{2.3}\\begin{longtable}{l|c}\\hline\n"
init_printing()

# output all the alpha-calculation combinations
for x,y,z,u,v in permutations(indices,5):
    if y>x and v>u:
        counter+=1
        eq=get_alpha(matrix,seen,x,y,z,u,v,True).simplify()
        eql=get_alpha(lmatrix,seen,x,y,z,u,v,False).simplify()
        
        # Creates Latex table row, FORMAT HERE
        latex_str+="("+labels[x]+","+labels[y]+":"+labels[z]+") - "+ labels[u]+","+labels[v]+ "& {$\\displaystyle " + latex(eql.simplify())+" $}\\\\[0.4cm]\\hline \n"
        pprint(eq)
        print()
        print("-----------------------------")
    else:
        continue
latex_str+="\\end{longtable}\n"
latex_str=latex_str.replace("scriptstyle0", "scriptscriptstyle 0")
with open("alphas.tex","w+") as tfile:
    tfile.write(latex_str)
    
# Show Latex pdf
preview(latex_str,output="pdf",filename="alphas.pdf",preamble=preamble,euler=False)
