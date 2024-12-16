#!/usr/bin/env sage

import sys

def simplify_polynomial(terms):
    total_length = 0
    for i in range(len(terms)):
        if terms[i] != 0:
            total_length = i + 1
    return terms[:total_length]

def multiply_polynomials(term1, term2):
    result = {}
    for i1,t1 in enumerate(term1):
        for i2,t2 in enumerate(term2):
            i_res = i1 +  i2
            if i_res not in result:
                result[i_res] = 0
            result[i_res] = result[i_res] + (t1 * t2)
    result_terms = [0] * (max(result)+1)
    for i,k in result.items():
        result_terms[i] = k
    return result_terms

def get_dimension_of_lattice(terms):
    return terms[1]

def clear_denominators(terms):
    result = terms
    dim = get_dimension_of_lattice(terms)
    for i in range(dim):
        result = multiply_polynomials(result, [1, -1])
    return result

def get_smallest(terms):
    for i in range(1,len(terms)):
        if terms[i] != 0:
            return i
    return None

def long_division(term, divisor):
    term = simplify_polynomial(term)
    divisor = simplify_polynomial(divisor)
    
    if len(term) < len(divisor):
        return [],term

    working = [0]*(len(term) - len(divisor) + 1)

    while len(term) >= len(divisor):
        term_leading = term[-1]
        divisor_leading = divisor[-1]

        multiplier = term_leading / divisor_leading
        working[len(term) - len(divisor)] = multiplier

        for i in range(len(divisor)):
            term[len(term) - i - 1] = term[len(term) - i - 1] - ( divisor[len(divisor) - i - 1] * multiplier )
        term = simplify_polynomial(term)

    working = simplify_polynomial(working)

    return working, term


def compute_quantum(decomp):
    result = [1]
    for d in decomp:
        result = multiply_polynomials(result, [1] + [0]*(d-1) + [-1] )
        result,ignore = long_division(result, [1,-1])
        assert(ignore == [])
    return result

def quantum_poly_decomp(terms):
    cleared = clear_denominators(terms)
    #
    dim = get_dimension_of_lattice(terms)
    #
    working = cleared
    decomp = []
    for i in range(dim):
        m = get_smallest(working)
        if m is None: return None
        decomp = decomp + [m]
        working, sanity = long_division(working, [1] + [0]*(m-1) + [-1] )
        if not (sanity == []):
            return None
    if not (working == [1]):
        return None
    if not (terms == compute_quantum(decomp)):
        return None
    return decomp


def get_ranks(poset):
    result = [0]*(poset.rank()+1)
    for h in poset:
        result[h.rank()] = result[h.rank()] + 1
    return result


def iterate_over_coordinates(max_vec):
    tmp = [0]*len(max_vec)
    
    while tmp != max_vec:
        yield tuple(tmp)
        error = True
        for i in range(len(max_vec)):
            if tmp[i] >= max_vec[i] - 1:
                tmp[i] = 0
            else:
                tmp[i] = tmp[i] + 1
                error = False
                break
        if error:
            break

def lattice(sizes):
    edges = []
    
    for i in range(len(sizes)):
        tmp_sizes = list(sizes)
        tmp_sizes[i] = tmp_sizes[i] - 1
        tmp_sizes = tuple(tmp_sizes)
        
        for source in iterate_over_coordinates(tmp_sizes):
            target = list(source)
            target[i] = target[i] + 1
            target = tuple(target)
            
            edges.append([source, target])
    
    return DiGraph(edges)

def is_lattice(poset):
    h = poset.hasse_diagram()
    
    sizes = quantum_poly_decomp(get_ranks(poset))
    
    if sizes is None:
        return False
    
    l = lattice(sizes)
    
    sub = h.subgraph_search(l)
    
    if sub is not None:
        return (sizes, l, h, sub)

    return False

def print_usage():
    print(
"""Usage: generate-sublattice.sage [OPTIONS] TYPE NUM

Attempts to find a spanning cubical lattice in a given group

Arguments:

    TYPE    one of A, B, C, D, E, F, G, or H
    NUM     this is the rank of the finite Coxeter (sub)group

Options:

    --help              print these usage instructions
    --draw-graph        creates a graph of the embedding if one is found
    --affine=UPPER      use the affine Coxeter group. In this case, you
                         need to give an upper bound for the Bruhat graph

Example:
    
    sage generate-sublattice.sage D 5

    sage generate-sublattice.sage --draw-graph --affine=1021020 A 2"""
        ,file=sys.stderr)
    sys.exit(1r)

def print_err(message):
    print(message, file=sys.stderr)
    sys.exit(1r)

def numeric_check(value):
    # note that we cannot use .is_numeric here since this will accept
    #  numeric values which are not a part of ASCII
    for v in value:
        if str(v) not in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
            return False
    return True

import datetime
def print_status_line(file, info):
    datestring = '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
    file.write("[" + datestring + "]: " + info + "\n")
    file.flush()

def open_file_or_die(fname):
    try: return open(fname, "w")
    except:
        print_err("Failed to open file '" + fname + "'")

def coxeter_or_die(ctype):
    try:
        W = CoxeterGroup(ctype)
        if W is not None: return W
    except: pass
    print_err("Cannot construct Coxeter group as specified")


def coxeter3_or_die(W, q):
    try:
        from coxeter3_sage import Coxeter3
        return Coxeter3(W, q)
    except:
        print_err("Coxeter3 might not be installed")

def element_or_die(W, element_desc):
    gens = W.gens()

    e = W.one()
    for s in element_desc:
        i = int(s)
        if len(gens) <= i: print_err("Cannot construct specified element")
        e = e * gens[i]

    return e

def file_prefix(ctype, max_element):
    if len(ctype) == 2:
        return ctype[0] + str(ctype[1])
    elif len(ctype) == 3:
        return ctype[0] + str(ctype[1]) + "_tilde_" + "".join([str(s) for s in max_element.reduced_word()])
    else:
        print_err("An unknown error occured in 'file_prefix'")


def compute_position(coordinate,zz=1):
    if len(coordinate) == 2:
        return tuple(coordinate)
    elif len(coordinate) == 3:
        z,y,x = tuple(coordinate)
        
        #SQRT2 = 1.41421356237310
        
        return (
            (x-y)*(zz-1),
            (x+y)*(zz-1)+0.9*z,
        )
    else:
        raise "err"

def write3d_graph(graph, cert, pdf_fname, sizes):
    zz = 1
    if len(sizes) == 3: zz = sizes[2]

    size_multi = sizes[1]
    if len(sizes) == 3: size_multi = size_multi + sizes[2]

    position = {}
    for sn, ln in cert.items():
        position[sn] = compute_position(tuple(ln), zz)
    G = graph.plot(
        vertex_labels=False,
        pos=position,
        figsize=(2 * size_multi, 2 * size_multi))
    G.save_image(pdf_fname)

def write_graph(graph, cert, pdf_fname, sizes, rank):
    if sizes is not None and len(sizes) <= 3 and cert is not None:
        write3d_graph(graph,cert,pdf_fname, sizes)
    else:
        G = graph.plot(vertex_labels=False)
        G.save_image(pdf_fname, figsize=(2*rank, 2*rank))

#############################################################
# read the command line arguments
#############################################################


inp_type = None
inp_num = None
inp_draw_graph = False
inp_affine = None

for line in sys.argv[1:]:
    if line == "--help": print_usage()
    elif line == "--draw-graph": inp_draw_graph = True
    elif line.startswith("--affine="):
        inp_affine = line[len("--affine="):]
        if not numeric_check(inp_affine):
            print_err("The UPPER element must be an integer sequence")
    elif line in ["A", "B", "C", "D", "E", "F", "G", "H"]:
        if inp_type is not None:
            print_err("Repeated Coxter type: only enter the Coxeter type once")
        inp_type = line
    elif numeric_check(line):
        inp_num = int(line)
    else:
        print_err("Unrecognised commandline option '"+line+"'")

if inp_type is None:
    print_err("You must specify a Coxeter type")
if inp_num is None:
    print_err("You must specify a rank for the finite Coxeter (sub)group")

#############################################################
#############################################################

ctype = None
if inp_affine is None:
    ctype = [inp_type, inp_num]
else:
    ctype = [inp_type, inp_num, 1]

W = coxeter_or_die(ctype)

R.<q> =  LaurentPolynomialRing(ZZ)

cox = coxeter3_or_die(W, q)

if (inp_affine) is None:
    e = W.long_element()
else:
    e = element_or_die(W, inp_affine)

#############################################################
#############################################################

with open_file_or_die(file_prefix(ctype,e) + "_status.txt") as status_file:
    cox_matrix = W.coxeter_matrix()

    with open(file_prefix(ctype,e) + "_matrix_desc.txt", "w") as f:
        f.write(str(cox_matrix))

    print_status_line(status_file, "output coxeter matrix")

    print_status_line(status_file, "beginning to construct Hasse diagram")

    bg = W.bruhat_graph(W.one(), e)
    pos = W.bruhat_interval_poset(W.one(), e)
    hasse = pos.hasse_diagram()

    print_status_line(status_file, "finished constructing Hasse diagram")

    print_status_line(status_file, "checking if there is a spanning cubical lattice")
    tmp = is_lattice(pos)
    
    if tmp is False:
        print_status_line(status_file, "!! there is no spanning cubical lattice")
        
        if inp_draw_graph:
          print_status_line(status_file, "beginnng to write graph")
          write_graph(hasse, None, file_prefix(ctype,e)+".pdf",None, pos.rank())
          print_status_line(status_file, "finished writing graph")

        print_status_line(status_file, "bye")
        exit()

    # there is a lattice let's get it
    sizes, l, h, lattice_subgraph = tmp

    print_status_line(status_file, "calculating an embedding of the lattice")
    a,cert = lattice_subgraph.is_isomorphic(l, certificate=True)
    print_status_line(status_file, "finished calculating an embedding of the lattice")

    invert_cert = { v: k for k,v in cert.items()}

    with open(file_prefix(ctype, e) + "_map.txt", "w") as f:
        for r in range(pos.rank()+1):
            f.write("\n")
            f.write("Vertices of rank "+str(r)+"\n")
            for k,v in invert_cert.items():
                if sum(k) != r:
                    continue
                line = str(k) + " ~> " + str(v.element.reduced_word())
                f.write(line+"\n")

    if inp_draw_graph:
        print_status_line(status_file, "beginnng to write graph")
        write_graph(hasse, cert, file_prefix(ctype,e)+".pdf",sizes, pos.rank())
        print_status_line(status_file, "finished writing graph")

    print_status_line(status_file, "Finished writing an embedding of the subgraph")
    print_status_line(status_file, "bye")


