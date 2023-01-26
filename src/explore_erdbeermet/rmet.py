import os
import sys
from sys import argv
import subprocess
import pathlib
from pyvis.network import Network
import pyvis.options
import string
import json
import numpy
from erdbeermet_pkg import recognition as rec

importfile=str(pathlib.Path().resolve())+"/"+argv[2]


with open(argv[1],"r") as f:
    lines=f.read().split("\n")

net=Network(directed=True,height="900px",width="1200px")
options=pyvis.options
nodes={}
labels=[]
seen=set()
in_nodes=lines[0].split(",")
start_node=lines[1]

erdbeer_string=""

for i in range(0,len(in_nodes)):
    nodes[in_nodes[i]]=i
    labels+=[in_nodes[i]]
    net.add_node(i,label=str(i),color="#009000",font="22px futura black",labelHighlightBold=True,physics=False)

seen.add(nodes[start_node])
dist=[[0 for x in range(0,len(nodes))] for x in range(0,len(nodes))]

for t in range(2,len(lines)-1):
    print(lines[t])
    # Parse part of each insert line
    insert=lines[t].split(";")[0]
    params=lines[t].split(";")[1]

    # Parse parents and child instruction and transl
    lparent=nodes[insert.split(",")[0]]
    rparent=nodes[insert.split(",")[1]]
    child=nodes[insert.split(",")[2]]
    if(lparent==rparent):
        net.add_edge(lparent,child,physics=True,label=str(params.split(",")[0]),font="14px futura grey align",alignlabel="middle",shadow=True)
    else:
        net.add_edge(lparent,child,physics=True,label=params.split(",")[0],font="14px futura grey",alignfont="middle",shadow=True)
        net.add_edge(rparent,child,physics=True,label=str(1-float(params.split(",")[0])),font="14px futura grey",align="middle",shadow=True)
    erdbeer_string+= "("+str(lparent)+", "+str(rparent)+": "+ str(child)+") "
    seen.add(child)
    if(lparent not in seen or rparent not in seen):
        print("Instruction in LINE " + str(t) + " contains parent nodes that do not exist!")
        sys.exit()
    alpha=float(params.split(",")[0])
    dl=float(params.split(",")[1])
    dr=float(params.split(",")[2])
    dc=float(params.split(",")[3])
    erdbeer_string+=str(alpha)+"; ["
    for x in range(0,t):
        if x==lparent and x==rparent:
            erdbeer_string+=str(dl+dr)+","
            continue
        if x==rparent:
            erdbeer_string+=str(dr)+","
            continue
        if x==lparent:
            erdbeer_string+=str(dl)+","
            continue
        if x==child:
            erdbeer_string+=str(dc)+","
            continue
        else:
            erdbeer_string+="0.0,"
            continue
    erdbeer_string=erdbeer_string[0:-1]
    erdbeer_string+="]\n"
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
dist_str=""
for x in range(0,len(dist)):
    for y in range(0,len(dist)):
        if(y>=x):
            dist_str+=str(dist[x][y])+"\t"
        else:
            dist_str+="\t"
    dist_str+="\n"

print(dist_str)
last_e=len(dist)-1
my_np_arr=numpy.asarray(dist)
for x in range(0,len(dist)-1):
    for y in range(x+1,len(dist)-1):
        for z in range(y+1,len(dist)-1):
            for last_e in range(z+1,len(dist)-1):
                if(x<y<z<last_e):
                    print(str(x) + " - " + str(y) + " - " +str(z) + " - " + str(last_e) )
                    print("0\t"+str(dist[x][y])+"\t"+str(dist[x][z])+"\t"+str(dist[x][last_e]))
                    print("\t0\t"+str(dist[y][z])+"\t"+str(dist[y][last_e]))
                    print("\t\t0\t"+str(dist[z][last_e]))
                    print("\t\t\t0")
                    print(rec.recognize4_new(my_np_arr,x,y,z,last_e))
                    print("----------------")
                else:
                    continue
with open(argv[2],"w+") as outfile:
    outfile.write("#nexus\n")
    outfile.write("BEGIN Taxa;\n")
    outfile.write("DIMENSIONS ntax=" + str(len(dist))+";\n")
    outfile.write("TAXLABELS\n")
    for x in range(0,len(labels)):
        outfile.write("["+str(x+1)+"]"+"'"+str(labels[x])+"'\n")
    outfile.write(";\n")
    outfile.write("END; [Taxa]\n")
    outfile.write("BEGIN Distances;\n")
    outfile.write("DIMENSIONS ntax=" + str(len(dist))+";\n")
    outfile.write("FORMAT labels=no diagonal triangle=upper;\n")
    outfile.write("MATRIX\n")
    outfile.write(dist_str)
    outfile.write(";\n")
    outfile.write("END; [Distances]\n")
    outfile.write("BEGIN st_Assumptions;\nuptodate;\ndisttransform=SplitDecomposition;\nsplitstransform=EqualAngle;\nSplitsPostProcess filter=dimension value=4;\n exclude  no missing;\nautolayoutnodelabels;\nEND; [st_Assumptions]")
with open("splitstree_commands.nex","w+") as f:
    f.write("LOAD FILE="+importfile+";\n")
    f.write("UPDATE;\n")
    f.write("UPDATE;\n")
    f.write("EXPORTGRAPHICS format=PNG file="+str(pathlib.Path().resolve())+"/"+argv[3]+" REPLACE=yes;\n")
    f.write("QUIT\n")

with open(argv[4],"w+") as simfile:
    simfile.write(erdbeer_string)
net.toggle_physics(True)
net.toggle_stabilization(True)
net.show_buttons()
#print(str(options.Options()).replace("'",'"'))
net.set_options(
"""
{"interaction": {"hideEdgesOnDrag": false, "hideNodesOnDrag": false, "dragNodes": true}, "configure": {"enabled": false}, "physics": {"enabled": true, "stabilization": true}, "edges": {"font": {"align": "middle"}, "smooth": {"enabled": true, "type": "dynamic"}, "color": {"inherit": true}}}
"""
)
net.show("pynet.html")
subprocess.call("python src/explore_erdbeermet/test_erdbeermet.py " + argv[4], shell=True)
#subprocess.call("firefox src/erdbeermet/output/vis_all.pdf &", shell=True,stdout=subprocess.DEVNULL)

# Let this point to your splitstree installation
subprocess.call("/scratch/bruno/SplitsTree/splitstree4/./SplitsTreeCMD -c " + "splitstree_commands.nex" ,shell=True)
# subprocess.call("start firefox "+argv[3],shell=True,stdout=subprocess.DEVNULL)
