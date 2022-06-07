# -*- coding: utf-8 -*-

from itertools import combinations, permutations
import numpy as np

from erdbeermet.tools.Tree import Tree, TreeNode

from collections import defaultdict

import networkx as nx

from decimal import Decimal as d

import random, sys

import math



__author__ = 'David Schaller'


glob_rtol=1e-09
glob_atol=1e-10


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
    if not np.all(np.logical_or(np.isclose(D, 0.0, rtol=glob_rtol, atol=glob_atol),
                                D > 0.0)):
        return False if not return_info else (False, 'negative distances')
    
    # check whether all diagonal entries are zero
    if np.any(np.diagonal(D)):
        return False if not return_info else (False, 'non-zero diagonal')
    
    # check whether the matrix is symmetric
    if not np.allclose(D, D.T, rtol=glob_rtol, atol=glob_atol):
        return False if not return_info else (False, 'not symmetric')
    
    # check the triangle inequality
    for i in range(N-1):
        for j in range(i+1, N):
            minimum = np.min(D[i, :] + D[:, j])
            if minimum < D[i, j] and not np.isclose(minimum, D[i, j],
                                                    rtol=glob_rtol, atol=glob_atol):
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


def distance_sums_matrix(D, x, y, z, u):
    
    xy_zu = D[x,y] + D[z,u]
    xz_yu = D[x,z] + D[y,u]
    xu_yz = D[x,u] + D[y,z]
    
    return xy_zu, xz_yu, xu_yz


def restrict_matrix(D, indices):
    
    if min(indices) < 0 or max(indices) >= D.shape[0]:
        raise IndexError("List contains index that is out of range!")
    
    n = len(indices)
    D_new = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            D_new[i, j] = D[indices[i], indices[j]]
    
    return D_new


def _recognize4_parent_xy(D, x, y, z, u):
    
    left = D[x,y] * (D[x,y] + 2 * D[z,u] - D[x,z] - D[y,u] - D[x,u] - D[y,z])
    right = (D[x,z] - D[y,z]) * (D[y,u] - D[x,u])
    
    return np.isclose(left, right, rtol=glob_rtol, atol=glob_atol) or left < right

def _recognize4_xy_zu(D, x, y, z, u):
    
    return (_recognize4_parent_xy(D, x, y, z, u) or 
            _recognize4_parent_xy(D, z, u, x, y))
    

def recognize4_new(D, x, y, z, u):
    
    if not is_pseudometric(restrict_matrix(D, [x, y, z, u])):
        return False
    
    dsums = distance_sums_matrix(D, x, y, z, u)
    
    if dsums[0] == max(dsums):
        return _recognize4_xy_zu(D, x, y, z, u)
    elif dsums[1] == max(dsums):
        return _recognize4_xy_zu(D, x, z, y, u)
    else:
        return _recognize4_xy_zu(D, x, u, y, z)
    
    
def recognize4_matrix_only(D):
    
    return recognize4_new(D, 0, 1, 2, 3)


def _compute_delta_x(alpha, xz, d_xy, delta_z):
    
    return xz - (1-alpha) * d_xy - delta_z


def _compute_delta_y(alpha, yz, d_xy, delta_z):
    
    return yz - alpha * d_xy - delta_z


def _compute_delta_z(xy, xz, yz):
    
    return 0.5 * (xz + yz - xy)


def _compute_d_xy(alpha, xz, yz, ux, uy, uz, delta_z):
    
    return (   (uz - alpha * ux - (1-alpha) * uy 
                - 2 * delta_z + alpha * xz + (1-alpha) * yz)
            / (2 * alpha * (1-alpha))   )

  
def _close_to_equal(a):
    
    if np.isclose(a, 0.0, rtol=glob_rtol, atol=glob_atol):
        return 0.0
    elif np.isclose(a, 1.0, rtol=glob_rtol, atol=glob_atol):
        return 1.0
    else:
        return a
    

def _non_negative(a):
    
    return np.isclose(a, 0.0, rtol=glob_rtol, atol=glob_atol) or a > 0.0


def _all_non_negative(a):
    
    for val in a:
        if not _non_negative(val):
            return False
        
    return True

   
def _compute_alpha(V, D, x, y, z, u, v):
    
    x = V.index(x)
    y = V.index(y)
    z = V.index(z)
    u = V.index(u)
    v = V.index(v)
    
    numerator   = (D[u,z] + D[v,y]) - (D[v,z] + D[u,y])
    denominator = (D[u,x] + D[v,y]) - (D[v,x] + D[u,y])
    
    print(f'({x},{y}:{z}),u={u},v={v}:   ({D[u,z]} + {D[v,y]}) - ({D[v,z]} + {D[u,y]})  // ({D[u,x]} + {D[v,y]}) - ({D[v,x]} + {D[u,y]}) = {numerator} // {denominator} = {numerator / denominator}')
    
    if not np.isclose(denominator, 0.0, rtol=glob_rtol, atol=glob_atol):
        return numerator / denominator
    else:
        return np.nan
def _compute_alpha_xoff(V, D, x, y, z, u, v):
    
    x = V.index(x)
    y = V.index(y)
    z = V.index(z)
    u = V.index(u)
    v = V.index(v)
    
    numerator   = (D[u,z] + D[v,y]) - (D[v,z] + D[u,y])
    denominator = (D[u,x] + D[v,y]) - (D[v,x] + D[u,y])
    
    print(f'({x},{y}:{z}),u={u},v={v}:   ({D[u,z]} + {D[v,y]}) - ({D[v,z]} + {D[u,y]})  // ({D[u,x]} + {D[v,y]}) - ({D[v,x]} + {D[u,y]}) = {numerator} // {denominator} = {numerator / denominator}')
    
    if not np.isclose(denominator, 0.0, rtol=glob_rtol, atol=glob_atol):
        if numerator/denominator==-0.0:
            #print("RETTTT 0.0")
            return 0.0
        return numerator / denominator
    else:
        if numerator/denominator==math.inf:
            return math.inf
        elif numerator/denominator==-math.inf:
            return -math.inf
        return np.nan

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
    for x in range(0,len(D)):
        for z in range(0,len(D)):
            if(x==z):
                continue
            else:
                candidates[(x,x,z)]={}
                for u, v in combinations(V, 2):
                    if u in (x,z) or v in (x,z):
                        continue
                    else:
                        candidates[(x,x,z)][(u,v)] = _compute_alpha_xoff(V, D, x, x, z, u, v)
    return candidates
def rank_candidates_selective(D,V,candidates_in):
    candidates={}
    for t in candidates_in:
        x=t[0]
        y=t[1]
        z=t[2]        
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
#    for x in range(0,len(D)):
#        for z in range(0,len(D)):
#            if(x==z):
#                continue
#            else:
#                candidates[(x,x,z)]={}
#                for u, v in combinations(V, 2):
#                    if u in (x,z) or v in (x,z):
#                        continue
#                    else:
#                        candidates[(x,x,z)][(u,v)] = _compute_alpha_xoff(V, D, x, x, z, u, v)
    return candidates
def _find_candidates(D, V, print_info, B):
    
    candidates = []
    n = len(V)
    
    if print_info: print(f'-----> n = {n}, V = {V} ---> Candidates')
    
    for x, y, z in permutations(V, 3):
        
        # Test
        if B and z in B:
            continue

        # considering x < y suffices
        if x > y:
            continue
        
        alpha = np.zeros(( (n-3) * (n-4) // 2 ,))
        
        pos = 0
        u_witness = None
        for u, v in combinations(V, 2):
            if u in  (x, y, z) or v in (x, y, z):
                continue
            
            alpha[pos] = _compute_alpha(V, D, x, y, z, u, v)
            
            if not u_witness and not np.isnan(alpha[pos]):
                u_witness = u
                
            pos += 1
        
        nan_mask = np.isnan(alpha)
        
        if not np.any(nan_mask) and np.allclose(alpha, alpha[0], rtol=glob_rtol, atol=glob_atol):
                        
            alpha[0] = _close_to_equal(alpha[0])
            
            if alpha[0] >= 0.0 and alpha[0] <= 1.0:
            
                str_alphas = [ "{0:0.15f}".format(a).split('.')[1] for a in alpha ]
                #print(str_alphas, file=sys.stderr)

                divergence_position = None
                for i in range(len(str_alphas[0])):
                    c = str_alphas[0][i]
                    for j in range(1, len(str_alphas)):
                        if str_alphas[j][i] != c:
                            divergence_position = i
                            break
                    if divergence_position != None:
                        break            
            
                candidates.append((x, y, z, u_witness, alpha[0], divergence_position))
                deltas = _compute_deltas(V, D, alpha[0], x, y, z, u_witness)
                
                if print_info: 
                    print(f'({x}, {y}: {z})', end='   ')
                    print('alpha=',  ["{0:0.15f}".format(i) for i in alpha] , end=' ,')
                    print('δx = {:.15f}, δy = {:.15f}, '\
                          'δz = {:.15f}, dxy = {:.15f}'.format(deltas[2],
                                                             deltas[3],
                                                             deltas[0],
                                                             deltas[1]))
            
        elif not np.all(nan_mask):
                    
            ref_alpha = alpha[ np.argmin(nan_mask) ]
            masked_alpha = np.ma.array(alpha, mask=nan_mask)
            
            if np.ma.allclose(masked_alpha, ref_alpha, masked_equal=True, rtol=glob_rtol, atol=glob_atol):
            
                ref_alpha = _close_to_equal(ref_alpha)
                
                if ref_alpha >= 0.0 and ref_alpha <= 1.0:
                    candidates.append((x, y, z, u_witness, ref_alpha, None))
                    
        else:

            # choose an arbitrary alpha (e.g. 0.5) and witness u (?)
            ref_alpha, u_witness = 0.5, None
            for u in V:
                if u not in (x, y, z):
                    u_witness = u
                    break
            candidates.append((x, y, z, u_witness, ref_alpha, None))
            
    return candidates


def _compute_deltas(V, D, alpha, x, y, z, u):
    
    x = V.index(x)
    y = V.index(y)
    z = V.index(z)
    u = V.index(u)
    
    delta_z = _compute_delta_z(D[x,y], D[x,z], D[y,z])
    
    # handle alpha in {0, 1}
    if alpha == 0.0 or alpha == 1.0:
        return delta_z, D[x,y], 0.0, 0.0
    
    d_xy = _compute_d_xy(alpha, D[x,z], D[y,z], D[u,x], D[u,y], D[u,z], delta_z)
    delta_x = _compute_delta_x(alpha, D[x,z], d_xy, delta_z)
    delta_y = _compute_delta_y(alpha, D[y,z], d_xy, delta_z)
    
    return delta_z, d_xy, delta_x, delta_y


def _update_matrix(V, D, x, y, delta_x, delta_y):
    
    x = V.index(x)
    y = V.index(y)
    
    if delta_x:             # if not 0.0
        D[:, x] -= delta_x
        D[x, :] -= delta_x
        D[x, x] = 0.0
    
    if delta_y:             # if not 0.0
        D[:, y] -= delta_y
        D[y, :] -= delta_y
        D[y, y] = 0.0
        
        
def _matrix_without_index(D, index):
    
    n = D.shape[0]
    
    if index < 0 or index >= n:
        raise IndexError(f"Index {index} is out of range!")
    
    D_new = np.zeros((n-1, n-1))
    
    indices = [i for i in range(n) if i != index]
    
    for i in range(n-1):
        for j in range(n-1):
            D_new[i, j] = D[indices[i], indices[j]]
    
    return D_new


def _finalize_tree(recognition_tree):
    
    def _sort_children(v):
        v.children.sort(key=lambda c: c.R_step)
        for c in v.children:
            _sort_children(c)
    
    for v in recognition_tree.postorder():
        if v.valid_ways and v.parent:
            v.parent.valid_ways += v.valid_ways
            
    recognition_tree.valid_ways = recognition_tree.root.valid_ways
    recognition_tree.successes = recognition_tree.root.valid_ways
            
    _sort_children(recognition_tree.root)
    
def add_child(D_in,x,y,z,alpha,delta):
    #print("------------")
    #print(z)
    #print(len(delta))
    #print(type(D_in))
    D=np.zeros((len(D_in)+1, len(D_in)+1))
    #print(len(D_in))
    #print(len(D))
    for i in range(0,len(D_in)):
        for j in range(0,len(D_in)):
            D[i, j]=D_in[i, j]
            D[j, i]=D_in[j, i]
    if (x == y or
        (x is None) or (y is None) or
        alpha == 1.0 or alpha == 0.0):
        
        if x is None or alpha == 0.0:
            x = y
            
        D[x, z] = 0.0
        D[z, x] = 0.0
        
        for u in range(z):
            D[u, z] = D[x, z]
            D[z, u] = D[x, z]
    # recombination event      
    else:
        for u in range(z):
            if u != x and u != y:
                
                d = alpha * D[x, u] + (1 - alpha) * D[y, u]
                D[u, z] = d
                D[z, u] = d
                
        D[z, x] = (1 - alpha) * D[x, y]
        D[x, z] = (1 - alpha) * D[x, y]
        D[z, y] = alpha * D[x, y]
        D[y, z] = alpha * D[x, y]
    
    # distance increment, i.e., independent evolution after event
    if len(delta) != z + 1:
        raise RuntimeError(f'invalid length of delta array for z={z}')
                
    for p in range(z):
        for q in range(p+1, z+1):
            D[p, q] += delta[p] + delta[q]
            D[q, p] = D[p, q] 
    return (D)
def recognize(D, first_candidate_only=False, print_info=False,  B=[], small_spike=False):
    """Recognition of type R matrices.
    
    Parameters
    ----------
    D : 2-dimensional numpy array
        A distance matrix.
    first_candidate_only : bool, optional
        If True, only consider the first found candidate for a merge event.
        The default is False.
    print_info : bool, True
        If True, print the recognition history. The default is False.
    
    Returns
    -------
    Tree
        The recognition tree.
    
    See also
    --------
    tools.Tree
    """

    n = D.shape[0]
    V = [i for i in range(n)]
    
    recognition_tree = Tree(TreeNode(n, V, D=D))
    stack = []
    
    # trivial failure if not a pseudometric
    if not is_pseudometric(D):
        if print_info: print('no pseudometric')
        recognition_tree.root.info = 'no pseudometric'
    
    # every pseudometric is additve and thus also an R matrix
    elif n <= 3:
        if print_info: print(print(f'SUCCESS on {V}'))
        recognition_tree.root.valid_ways = 1
    
    # otherwise start the recognition algorithm
    else:
        stack.append(recognition_tree.root)
    
    
    while stack:
        
        parent = stack.pop()
        V, D = parent.V, parent.D
        n = len(V)
  
        if print_info: print(D)
        
        if n > 4:

            # cadidates: (x, y ,z , u_witness, alpha)

            # sort list so that A < B if A(delta_y) < B(delta'_y) && A(y) == B(y)    
            # => "oldest" candidate at bottom of list, "youngest" candidate at top of list
        
            candidates = _find_candidates(D, V, print_info, B)
            
            found_valid = False

            #spikelength
            if small_spike and n > 5:
                 
                new_candidates = []
                for x, y, z, u_witness, alpha, div in candidates:
                    V_copy = V.copy()
                    V_copy.remove(z)
                                    
                    deltas = _compute_deltas(V, D, alpha, x, y, z, u_witness)
                    
                    if not _all_non_negative(deltas):
                        continue
                    
                    D_copy = _matrix_without_index(D, V.index(z))
                    _update_matrix(V_copy, D_copy, x, y, deltas[2], deltas[3])
                    
                    still_metric, metric_info = is_pseudometric(D_copy,
                                                                return_info=True,
                                                                V=V_copy)
                    if not still_metric:
                        continue
                        
                    new_candidates.append( (x, y, z, u_witness, alpha, div) )
                    
                candidates = new_candidates

                edges = []
                delta = []
                new_candidates = []

                for i, c1 in enumerate(candidates):
                    x1 = c1[0]
                    y1 = c1[1]
                    z1 = c1[2]
                    u1 = c1[3]
                    alpha1 = c1[4]
                    delta_z1, d_xy1, delta_x1, delta_y1 = _compute_deltas(V, D, alpha1, x1, y1, z1, u1)

                    if print_info:
                        print("Spikelength Node",i, x1, delta_x1, y1, delta_y1, z1, delta_z1)

                    for j, c2 in enumerate(candidates):
                        x2 = c2[0]
                        y2 = c2[1]
                        z2 = c2[2]
                        u2 = c2[3]
                        alpha2 = c2[4]
                        delta_z2, d_xy2, delta_x2, delta_y2 = _compute_deltas(V, D, alpha2, x2, y2, z2, u2)

                        if i >= j:
                            continue
                        elif not set([x1, y1, z1]) & set([x2, y2, z2]):
                            continue
                        else:
                            for a, delta_a in zip([x1, y1, z1], [delta_x1, delta_y1, delta_z1]):
                                for b, delta_b in zip([x2, y2, z2], [delta_x2, delta_y2, delta_z2]):
                                    if a == b and not np.isclose(delta_a, delta_b, rtol=glob_rtol, atol=glob_atol):
                                        if delta_a < delta_b:
                                            edges.append((i, j))
                                            delta.append((delta_x1,delta_x2))
                                        else:
                                            edges.append((j, i))
                                            delta.append((delta_x2,delta_x1))

                if print_info:
                    print("CycleGraph")
                    for e in edges:
                        print(e)
                    
                graph = nx.DiGraph(edges)
                for i in range(len(candidates)):
                    if not graph.has_node(i):
                        graph.add_node(i)
                for node, degree in graph.in_degree():
                    if degree == 0:
                        new_candidates.append(candidates[node])

                if new_candidates:
                    candidates = []
                    candidates.append(random.choice(new_candidates))
                else:
                    for cycle in nx.simple_cycles(graph):
                        if len(cycle) >= 2:
                            print("ERROR: cycle detected. Aborting.")
                            parent.info = "cycle detected"
                            _finalize_tree(recognition_tree)    
                            return recognition_tree

            if print_info: 
                print(f'-----> n = {n}, V = {V} ---> R-steps actually carried out')
            for x, y, z, u_witness, alpha, div in candidates:
                V_copy = V.copy()
                V_copy.remove(z)
                
                child = TreeNode(n-1, V_copy, R_step=(x, y, z, alpha), Divergence=div)
                parent.add_child(child)
                
                deltas = _compute_deltas(V, D, alpha, x, y, z, u_witness)
                
                if print_info:
                    print('({}, {}: {}) alpha={:.5f} | u={}'.format(x, y, z, alpha, u_witness),
                          end='   ')
                    print('δx = {:.3f}, δy = {:.3f}, '\
                          'δz = {:.3f}, dxy = {:.3f}'.format(deltas[2],
                                                             deltas[3],
                                                             deltas[0],
                                                             deltas[1]))
                
                if not _all_non_negative(deltas):
                    if print_info: print('         |___ negative δ/dxy')
                    child.info = 'negative delta/dxy'
                    child.Divergence=None
                    continue
                
                D_copy = _matrix_without_index(D, V.index(z))
                _update_matrix(V_copy, D_copy, x, y, deltas[2], deltas[3])
                child.D = D_copy
                
                still_metric, metric_info = is_pseudometric(D_copy,
                                                            return_info=True,
                                                            V=V_copy)
                
                if not still_metric:
                    if print_info: print( '         |___ no pseudometric')
                    if print_info: print(f'         |___ {metric_info}')
                    child.info = 'no pseudometric'
                    child.Divergence=None
                    continue


                
                found_valid = True
                if print_info: print(f'         |___ STACKED {V_copy}')
                stack.append(child)
                                
                # for n = 5 always check all candidates
                if first_candidate_only and n > 5:
                    break
                
            if not candidates or not found_valid:
                parent.info = 'no candidate'
                
        else:
            if print_info: print(f'-----> n = {n} R-map test')
            if recognize4_matrix_only(D):
                if print_info: print(f'SUCCESS on {V}')
                parent.valid_ways = 1
            else:
                if print_info: print(f'NO R-MAP on {V}')
                parent.info = 'spikes too short'
    
    _finalize_tree(recognition_tree)    
    return recognition_tree
