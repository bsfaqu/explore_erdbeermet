from sympy import *
from sys import argv
from texttable import Texttable
from itertools import permutations


# Simply performs the alpha calculation
def equatify(a, b, c, d, e, f, g, h):
    return (((a + b) - (c + d)) / ((e + f) - (g + h)))


# pretty print matrix with col and row headers
def ppmat(matrix, seen):
    col_width = [0 for x in range(0, len(matrix))]
    header = [""] + seen
    t = Texttable()
    t.add_row(header)
    for x in range(0, len(matrix)):
        for y in range(0, len(matrix)):
            col_width[y] = max(col_width[y], len(str(matrix[x][y])))
        ins_row = [seen[x]] + matrix[x]
        print(ins_row)
        t.add_row(ins_row)
    col_width = [1] + col_width
    t.set_cols_width(col_width)
    print(t.draw())


# calculate alpha and output some info
def get_alpha(matrix, seen, x, y, z, u, v, print_info):
    if (print_info):
        print("alpha: (" + seen[x] + "," + seen[y] + ":" + seen[z] + ") - WITNESS " + seen[u] + "," + seen[v])
        print()
    return (((matrix[u][z] + matrix[v][y]) - (matrix[v][z] + matrix[u][y])) / (
                (matrix[u][x] + matrix[v][y]) - (matrix[v][x] + matrix[u][y])))


# def _compute_delta_x(alpha, xz, d_xy, delta_z):
#     # print("--")
#     # print(alpha)
#     # print(xz)
#     # print(d_xy)
#     # print(delta_z)
#     # print("--")
#     return xz - (1-alpha) * d_xy - delta_z

def _compute_d_xy(alpha, xz, yz, ux, uy, uz, delta_z):
    return ((uz - alpha * ux - (1 - alpha) * uy
             - 2 * delta_z + alpha * xz + (1 - alpha) * yz)
            / (2 * alpha * (1 - alpha)))


# TODO: INCLUDE DXY CALCULATIONS!!!
def get_delta_x(matrix, x, y, z, alpha, delta_z):
    d_xy = _compute_d_xy(alpha, matrix[x][z], matrix[y][z], matrix[u][x], matrix[u][y], matrix[u][z], delta_z)
    return (matrix[x][z].simplify() - (1 - alpha) * d_xy.simplify() - delta_z)


def get_delta_y(matrix, x, y, z, alpha, delta_z):
    d_xy = _compute_d_xy(alpha, matrix[x][z], matrix[y][z], matrix[u][x], matrix[u][y], matrix[u][z], delta_z)
    return (matrix[y][z].simplify() - (alpha * d_xy.simplify()) - delta_z)


def get_delta_z(matrix, x, y, z):
    return (0.5 * (matrix[x][z].simplify() + matrix[y][z].simplify() - matrix[x][y].simplify()))


# def _compute_delta_y(alpha, yz, d_xy, delta_z):
#     # print("--")
#     # print(alpha)
#     # print(yz)
#     # print(d_xy)
#     # print(delta_z)
#     # print("--")
#     return yz - alpha * d_xy - delta_z


# def _compute_delta_z(xy, xz, yz):
#
#     return 0.5 * (xz + yz - xy)

# read scenario from argv1
with open(argv[1], "r") as f:
    lines = f.read().split("\n")

# init nodes, labels, visited nodes and distance matrix
nodes = {}
labels = []
in_nodes = lines[0].split(",")
start_node = lines[1]
seen = []
seen += [start_node]
iterations = 1

# matrix stores console output symbols
matrix = [[0 for x in range(0, len(in_nodes))] for y in range(0, len(in_nodes))]
matrix_update = [[0 for x in range(0, len(in_nodes))] for y in range(0, len(in_nodes))]

# lmatrix stores latex code for pdf output
lmatrix = [[0 for x in range(0, len(in_nodes))] for y in range(0, len(in_nodes))]
lmatrix_update = [[0 for x in range(0, len(in_nodes))] for y in range(0, len(in_nodes))]

lmatrix1_based = [[0 for x in range(0, len(in_nodes))] for y in range(0, len(in_nodes))]
lmatrix1_based_update = [[0 for x in range(0, len(in_nodes))] for y in range(0, len(in_nodes))]

# preprocess nodes well encounter
for i in range(0, len(in_nodes)):
    nodes[in_nodes[i]] = i
    labels += [in_nodes[i]]

# read scenario and compute distances
for t in range(2, len(lines) - 1):
    # print(lines[t])
    iterations += 1
    # Parse part of each insert line
    insert = lines[t].split(";")[0]
    params = lines[t].split(";")[1]

    # Parse parents and child instruction
    lparent = insert.split(",")[0]
    rparent = insert.split(",")[1]
    child = insert.split(",")[2]

    # Append new node to seen
    seen += [child]

    # Initialize matrix with "base distance" d0xy
    # For the root once we reached iteration 4
    if (iterations == 5):
        for x in range(0, len(seen)):
            for y in range(0, len(seen)):
                if x == y:
                    continue
                else:
                    # console output formulas
                    exec("d0" + seen[x] + seen[y] + " = symbols(" + "\"d0" + seen[x] + seen[y] + "\")")
                    exec("d0" + seen[y] + seen[x] + " = symbols(" + "\"d0" + seen[y] + seen[x] + "\")")
                    matrix[x][y] = eval("d0" + seen[y] + seen[x])
                    matrix[y][x] = eval("d0" + seen[y] + seen[x])

                    # latex formulas
                    exec("ld0" + seen[x] + seen[y] + " = symbols(" + "\"d^{\\scriptstyle0}_{" + seen[x] + "" + seen[
                        y] + "}\")")
                    exec("ld0" + seen[y] + seen[x] + " = symbols(" + "\"d^{\\scriptstyle0}_{" + seen[y] + "" + seen[
                        x] + "}\")")
                    lmatrix[x][y] = eval("ld0" + seen[y] + seen[x] + "")
                    lmatrix[y][x] = eval("ld0" + seen[y] + seen[x] + "")
    # Update the matrix when were beyond the root construction
    elif (iterations > 5):

        # Initialize 1-based latex matrix
        # if iterations ==5:
        #     for x in range(0,len(seen)):
        #         for y in range(0,len(seen)):
        #             if x==y:
        #                 continue
        #             else:
        #                 exec("ld1"+seen[x]+seen[y]+" = symbols("+"\"d^{\\scriptstyle1}_{"+seen[x]+""+seen[y]+"}\")")
        #                 exec("ld1"+seen[y]+seen[x]+" = symbols("+"\"d^{\\scriptstyle1}_{"+seen[y]+""+seen[x]+"}\")")
        #                 lmatrix1_based[x][y]=eval("ld1"+seen[y]+seen[x]+"")
        #                 lmatrix1_based[y][x]=eval("ld1"+seen[y]+seen[x]+"")
        # if iterations > 5:
        #     lmatrix1_based_update = lmatrix1_based.copy()
        # Save updated values in updated matrix

        matrix_update = matrix.copy()
        lmatrix_update = lmatrix.copy()

        ####### CONSOLE OUT #########
        # Declare symbolic alpha for corresponding R-Step
        exec("a" + str(iterations - 4) + "=symbols(" + "\"a" + str(iterations - 4) + "\")")

        # Declare symbolic deltas for corresponding R-Step
        exec("del" + str(iterations - 4) + "_" + lparent + "=symbols(" + "\"del" + str(
            iterations - 4) + "_" + lparent + "\")")
        exec("del" + str(iterations - 4) + "_" + rparent + "=symbols(" + "\"del" + str(
            iterations - 4) + "_" + rparent + "\")")
        exec("del" + str(iterations - 4) + "_" + child + "=symbols(" + "\"del" + str(
            iterations - 4) + "_" + child + "\")")

        ####### LATEX OUT #########
        # Declare symbolic alpha for corresponding R-Step
        exec("la" + str(iterations - 4) + "=symbols(" + "\"alpha_" + str(iterations - 4) + "\")")

        # Declare symbolic deltas for corresponding R-Step
        exec("ldel" + str(iterations - 4) + "_" + lparent + "=symbols(" + "\"\\delta^" + str(
            iterations - 4) + "_{" + lparent + "}\")")
        exec("ldel" + str(iterations - 4) + "_" + rparent + "=symbols(" + "\"\\delta^" + str(
            iterations - 4) + "_{" + rparent + "}\")")
        exec("ldel" + str(iterations - 4) + "_" + child + "=symbols(" + "\"\\delta^" + str(
            iterations - 4) + "_{" + child + "}\")")

        # R1 - hybridization
        for x in range(0, len(seen)):
            if (seen[x] == child):
                continue
            else:
                # update child distances
                ####### CONSOLE OUT #########
                matrix_update[x][nodes[child]] = matrix[x][nodes[lparent]] * eval("a" + str(iterations - 4)) + \
                                                 matrix[x][nodes[rparent]] * (1 - eval("a" + str(iterations - 4)))
                matrix_update[nodes[child]][x] = matrix[x][nodes[lparent]] * eval("a" + str(iterations - 4)) + \
                                                 matrix[x][nodes[rparent]] * (1 - eval("a" + str(iterations - 4)))

                ####### LATEX OUT #########
                lmatrix_update[x][nodes[child]] = lmatrix[x][nodes[lparent]] * eval("la" + str(iterations - 4)) + \
                                                  lmatrix[x][nodes[rparent]] * (1 - eval("la" + str(iterations - 4)))
                lmatrix_update[nodes[child]][x] = lmatrix[x][nodes[lparent]] * eval("la" + str(iterations - 4)) + \
                                                  lmatrix[x][nodes[rparent]] * (1 - eval("la" + str(iterations - 4)))

                # if iterations > 5:
                #     lmatrix1_based_update[x][nodes[child]]=lmatrix1_based[x][nodes[lparent]]*eval("la"+str(iterations-5)) + lmatrix1_based[x][nodes[rparent]]*(1-eval("la"+str(iterations-5)))
                #     lmatrix1_based_update[nodes[child]][x]=lmatrix1_based[x][nodes[lparent]]*eval("la"+str(iterations-5)) + lmatrix1_based[x][nodes[rparent]]*(1-eval("la"+str(iterations-5)))

        # R2 - mutation
        for x in range(0, len(seen)):
            for y in range(x + 1, len(seen)):
                if (x == y):
                    continue
                if (seen[x] == lparent or seen[y] == lparent):
                    ####### CONSOLE OUT #########
                    matrix_update[x][y] += eval("del" + str(iterations - 4) + "_" + lparent)
                    matrix_update[y][x] += eval("del" + str(iterations - 4) + "_" + lparent)

                    ####### LATEX OUT #########
                    lmatrix_update[x][y] += eval("ldel" + str(iterations - 4) + "_" + lparent)
                    lmatrix_update[y][x] += eval("ldel" + str(iterations - 4) + "_" + lparent)

                if (seen[x] == rparent or seen[y] == rparent):
                    ####### CONSOLE OUT #########
                    matrix_update[x][y] += eval("del" + str(iterations - 4) + "_" + rparent)
                    matrix_update[y][x] += eval("del" + str(iterations - 4) + "_" + rparent)

                    ####### LATEX OUT #########
                    lmatrix_update[x][y] += eval("ldel" + str(iterations - 4) + "_" + rparent)
                    lmatrix_update[y][x] += eval("ldel" + str(iterations - 4) + "_" + rparent)
                if (seen[x] == child or seen[y] == child):
                    ####### CONSOLE OUT #########
                    matrix_update[x][y] += eval("del" + str(iterations - 4) + "_" + child)
                    matrix_update[y][x] += eval("del" + str(iterations - 4) + "_" + child)

                    ####### LATEX OUT #########
                    lmatrix_update[x][y] += eval("ldel" + str(iterations - 4) + "_" + child)
                    lmatrix_update[y][x] += eval("ldel" + str(iterations - 4) + "_" + child)

        # Finalize the update
        matrix = matrix_update.copy()
        lmatrix = lmatrix_update.copy()

# Look at these pretty matrices
ppmat(matrix.copy(), seen)
ppmat(lmatrix.copy(), seen)

indices = [x for x in range(0, len(matrix))]
counter = 0

# LaTeX preamble and table init
preamble = "\\documentclass[12pt]{article}\n\\usepackage[a3paper]{geometry}\n\\usepackage{amsmath,amsfonts}\n\\usepackage{longtable}\n\\usepackage{enumerate}\n\\begin{document}\n"
latex_str = "\\renewcommand*{\\arraystretch}{2.3}\n\\begin{longtable}{l|l}\\hline\n"
init_printing()

if False:
    # noinspection PyUnreachableCode
    for x, y, z, u in permutations(indices, 4):
        if y > x:
            if (x == 0 and y == 1 and z == 2):
                alpha = get_alpha(lmatrix, seen, x, y, z, 1, 4, True).simplify()
            elif x == 3 and y == 4 and z == 5:
                alpha = get_alpha(lmatrix, seen, x, y, z, 0, 2, True).simplify()
            else:
                alpha = symbols("\\alpha_{" + str(x) + ";" + str(y) + ";" + str(z) + "}")
            # deltaz
            delta_z = get_delta_z(lmatrix, x, y, z).simplify().evalf(2)
            # deltax
            delta_x = get_delta_x(lmatrix, x, y, z, alpha, delta_z).simplify().evalf(2)
            # deltay
            delta_y = get_delta_y(lmatrix, x, y, z, alpha, delta_z).simplify().evalf(2)

            latex_str += "$\delta^{" + str(labels[x]) + "}_{" + str(labels[x]) + "," + str(labels[y]) + "," + str(
                labels[z]) + ";" + str(labels[u]) + "}$ "
            latex_str += "& {$\\displaystyle " + latex(delta_x.simplify()) + " $}\\\\[0.4cm]\\hline \n"

            latex_str += "$\delta^{" + str(labels[y]) + "}_{" + str(labels[x]) + "," + str(labels[y]) + "," + str(
                labels[z]) + ";" + str(labels[u]) + "}$ "
            latex_str += "& {$\\displaystyle " + latex(delta_y.simplify()) + " $}\\\\[0.4cm]\\hline \n"

            latex_str += "$\delta^{" + str(labels[z]) + "}_{" + str(labels[x]) + "," + str(labels[y]) + "," + str(
                labels[z]) + ";" + str(labels[u]) + "}$ "
            latex_str += "& {$\\displaystyle " + latex(delta_z.simplify()) + " $}\\\\[0.4cm]\\hline \n"

            sum = delta_x + delta_y + delta_z
            sum.simplify().evalf()
            latex_str += "SUM & {$\\displaystyle " + latex(sum.simplify()) + " $}\\\\[0.4cm]\\hline \n"

# output all the alpha-calculation combinations
for x, y, z, u, v in permutations(indices, 5):
    if y > x and v > u or (x, y in [0, 2, 3, 4] and z == 3 and v > u):
        counter += 1
        eq = get_alpha(matrix, seen, x, y, z, u, v, True).simplify()
        eql = get_alpha(lmatrix, seen, x, y, z, u, v, False).simplify()

        # Creates Latex table row, FORMAT HERE
        latex_str += "(" + labels[x] + "," + labels[y] + ":" + labels[z] + ") - " + labels[u] + "," + labels[
            v] + "& {$\\displaystyle " + latex(eql.simplify()) + " $}\\\\[0.4cm]\\hline \n"
        pprint(eq)
        print()
        print("-----------------------------")
    else:
        continue

latex_str += "\\end{longtable}\n"
latex_str = latex_str.replace("scriptstyle0", "scriptscriptstyle 0")
with open(argv[2], "w+") as tfile:
    tfile.write(preamble + latex_str + "\\end{document}")

# Show Latex pdf
preview(latex_str, output="pdf", filename="alphas.pdf", preamble=preamble, euler=False)
