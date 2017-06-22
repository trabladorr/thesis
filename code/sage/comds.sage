load("xorcount.sage")

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import functools, time, collections

#inspired by https://stackoverflow.com/a/12605800/
#similar to itertools.product, but only calls a generator's next() when required, while caching previous generated elements in a list
def cache_gen_product(caching=True,*iterables):

    class Cache_gen:
        def __init__(self,gen):
            if isinstance(gen,collections.Iterator):
                self.gen = gen
            else:
                self.gen = iter(gen)
            self.data = []
            self.i = -1

        def __iter__(self): #restart iteration
            self.i = -1
            return self

        def next(self):
            self.i += 1
            # assert(self.i <= len(self.data))
            if self.i == len(self.data):
                self.data.append(next(self.gen)) #potentially raises StopIteration to signify end 
            return self.data[self.i]

    def expand(vals, prev_vals=None):
        if prev_vals is None:
            for val in vals:
                yield (val,)
        else:
            for prev_val in prev_vals:
                for val in vals:
                    yield prev_val + (val,)
            

    n = len(iterables)
    assert(n>0)

    ret = None

    for i,iterable in enumerate(iterables):
        if caching:
            ret = expand(Cache_gen(iterable), ret)
        else:
            ret = expand(iterable, ret)

    return ret

#generator for all polynomial ring elements up to degree n
def gen_poly_ring_elements(poly_ring,n):
    field_ring = poly_ring.base_ring()
    generators = [iter(field_ring) for _ in range(n)]
    elems = cache_gen_product(field_ring.degree() == 1,*generators)
    return itertools.imap(lambda l:poly_ring(list(reversed(l))),elems)

# tests for gen_poly_ring_elements
for F_test in (GF(3),GF(2^2, 'a')):
    for i in range(1,5):
        polys = set(gen_poly_ring_elements(F_test['x'],i))
        # print "All polynomials in ",F_test,"with degree up to",i,":",len(polys)
        assert(len(polys) == F_test.order()**i)

#transpose a polynomial representing a circulant matrix
def poly_transpose(p,x,n):
    if p == 0:
        return 0
    coeff_list = p.list()+([0]*(n-len(p.list())))
    constant = coeff_list.pop(0)
    ret = sum ([e*x**(i+1) for i,e in enumerate(reversed(coeff_list))])
    return ret + (constant*x**0)

#calculate the reciprocal polynomial of p
def rev_factor(p,x):
    coeff_list = p.list()
    constant = coeff_list[0]
    ret = sum ([e*x**(i) for i,e in enumerate(reversed(coeff_list))])
    return ret/constant

#helper function for printing elements of extension fields
def q_to_list_hex(g,n):
    return [field_to_int(e) for e in g.list()+([0]*(n-len(g.list())))]

#convert a polynomial representation of a circulant matrix to the matrix itself
def q_to_mat(q):
    return circ(q.list()[::-1])

#convert a circulant matrix to its polynomial representation
def mat_to_q(mat,x):
    return sum([e*x**(i) for i,e in enumerate(mat[0])])

#convert a circulant matrix first row coefficients to its polynomial representation
def coeff_to_q(coeff,q_x):
    return sum([e*q_x**(i) for i,e in enumerate(coeff)])

#run a number of checks to assert that we found the ideal multiplicative units correctly
def assert_Eis_correct(Eis,q_x,paired_ideal_indices,n):
    #assert that the sum of the multiplicative units, Ei, is 1
    sum_ideals = 0
    for Ei in Eis:
        sum_ideals += Ei
    assert(sum_ideals == 1)

    #assert that multiplying all Eis by all Eis produces a diagonal matrix of them
    for Ei in Eis:
        for Ej in Eis:
            assert(((Ei == Ej) and Ei*Ej == Ei) or ((Ei != Ej) and Ei*Ej == 0))

    #assert that the paired Eis pair up nicely
    for k,(i,j) in enumerate(paired_ideal_indices):
        assert(Eis[i] == poly_transpose(Eis[j],q_x,n))
        assert(Eis[j] == poly_transpose(Eis[i],q_x,n))
        assert(Eis[i] * poly_transpose(Eis[j],q_x,n) == Eis[i])
        assert(poly_transpose(Eis[i],q_x,n) * Eis[j] == Eis[j])

def gen_ideal_orth_elements_genopt(poly_ring,q_ring,q_x,polys,n,unpaired_ideal_indices,paired_ideal_indices,len_unpaired_orth,len_paired_orth,factor_deg,verbose=False,asserts=True):
    field_ring = poly_ring.base_ring()
    Eis = [None for _ in polys]
    q_polys = [q_ring(p) for p in polys]

    
    #optimization attempt to find the multiplicative units without iterating through all ideal elements
    for i in range(len(q_polys)):
        candidate = mat_to_q(q_to_mat(q_polys[i]).echelon_form(),q_x)
        for j in iter(field_ring):
            if (j*candidate)*q_polys[i] == q_polys[i]:
                # print candidate*j
                Eis[i] = candidate*j
                break

    #find the multiplicative units by iterating through all ideal elements
    if not all([Ei is not None for Ei in Eis]):
        if verbose:
            print "Iterating over",field_ring.order()**max(factor_deg),"elements to find ideal multiplicative units...\n"
        for i,p in enumerate(polys):
            if Eis[i] is not None:
                continue
            j= 0
            for e in gen_poly_ring_elements(poly_ring,factor_deg[i]):
                j += 1
                ei = q_ring(p*e)
                if ei != 0 and ei*ei == ei:
                    Eis[i] = ei
                    break

    #check that the multiplicative units are correct
    if asserts:
        assert_Eis_correct(Eis,q_x,paired_ideal_indices,n)

    # generator for unpaired ideal orthogonal elements
    # iterates through all ideal elements until all orthogonal elements are found
    def gen_unpaired_orth_elements(Ei,poly,factor_deg,poly_ring,q_ring,len_unpaired_orth,n,asserts=True):
        elems = gen_poly_ring_elements(poly_ring,factor_deg)
        i = 0
        for e in elems:
            ei = q_ring(poly*e)
            if ei*poly_transpose(ei,q_x,n) == Ei:
                yield ei
                i += 1
                if (i == len_unpaired_orth):
                    return
        assert(False)

    #optimization function, attempts to find the corresponding elements in paired ideals
    def solve_AXeqB_circ(qA,qB,n,q_ring):
        q_x = q_ring.gen()

        mat = q_to_mat(qA).T
        sol = mat.solve_right(vector(qB.list()))
        return coeff_to_q(sol,q_x)

    # generator for paired (i,j) ideal orthogonal elements,
    def gen_paired_orth_elements(Ei,Ej,poly,poly_star,factor_deg,poly_ring,q_ring,n,asserts=True):

        q_x = q_ring.gen()


        # iterates through all elements of the ith ideal
        for e in gen_poly_ring_elements(poly_ring,factor_deg):
            ei = q_ring(poly*e)
            if ei == 0: 
                continue

            # assert that doesn't always work
            # if asserts:
            #     assert (solve_AXeqB_circ(q_ring(e),ei,n,q_ring) == q_ring(poly))


            #optimization attempt to avoid iterating through all elements of the jth ideal
            e2 = solve_AXeqB_circ(ei,Ei,n,q_ring)
            ej = e2*Ej


            if ei*poly_transpose(ej,q_x,n) != Ei:
                #fall back to iterating through all elements of the jth ideal to find the corresponding element
                for e2 in gen_poly_ring_elements(poly_ring,n):
                    ej = q_ring(e2*poly_star)
                    if ei*poly_transpose(ej,q_x,n) == Ei:
                        break

            #check that we did indeed find the pair (too paranoid)
            # if asserts:
            #     assert(ei*poly_transpose(ej,q_x,n) == Ei)

            yield ei + ej

    #initialize generators for all ideals
    unpaired_orth = [gen_unpaired_orth_elements(Eis[i],polys[i],factor_deg[i],poly_ring,q_ring,len_unpaired_orth[unpaired_ideal_indices.index(i)],n) for i in unpaired_ideal_indices]
    
    paired_orth = [gen_paired_orth_elements(Eis[i],Eis[j],polys[i],polys[j],factor_deg[i],poly_ring,q_ring,n) for (i,j) in paired_ideal_indices]

    return Eis,unpaired_orth,paired_orth



def find_COMDS_matrices(n,r=None,field=None,mod=None,verbose=False,find_all=False,asserts=True,update_every = 10000):

    #initialize galois field
    if field is not None:
        field_ring = field
        field_x = field.gen()
    elif mod is not None and r is not None:
        field_ring,field_x = make_field(r,mod)
    elif r is not None:
        field_ring,field_x = GF(2^r).objgen()
    else:
        print "No Field or binary extension field exponent given"
        return
    mod = field_ring.modulus().list()

    #assert field order is coprime to n
    if(gcd(field_ring.order(),n)>1):
        print "q=",field_ring.order()," and n=",n," are not coprime: gcd=",gcd(field_ring.order(),n)
        return

    #initialize XOR count generator, only applies for binary extension fields
    if (field_ring.characteristic() == 2):
        r = field_ring.degree()
        mult_xor_counts = gen_xor_count(r,mod,True)
        if verbose is not None and verbose:
            print "Generated XOR counts for",field_ring
        if verbose is not None:
            print "\n\nChecking",field_ring,", mod:",field_ring.modulus()
    else:
        mult_xor_counts = {i:int(0) for i in range(int(field_ring.order()))}

    #create polynomial ring
    poly_ring, poly_x = field_ring['p'].objgen()
    poly_mod = poly_x^n - 1

    #calculate the factors of x^n-1
    factors = poly_mod.factor()

    if asserts:
        for _,fc in factors:
            assert(fc == 1)

    factors = [f for f,_ in factors]

    #create polynomial quotient ring, each element representing a circulant matrix
    q_ring, q_x = poly_ring.quo(poly_mod).objgen()

    if verbose is not None and verbose:
        print "There are",len(q_ring.base_ring())**n,"circulant",n,"x",n," matrices in",field_ring,"\n"

    #create the elements generating the ideals
    ideal_gens = [q_ring(poly_mod/f) for f in factors]
    factor_deg = [f.degree() for f in factors]

    #find paired and unpaired ideals, defined by indices
    pair_factors = [rev_factor(f,poly_x) for f in factors]

    unpaired_ideal_indices = []
    paired_ideal_indices = []

    for i,(f,fi) in enumerate(zip(factors,pair_factors)):
        if f == fi: 
            unpaired_ideal_indices.append(i)
        elif fi in factors:
            ifi = factors.index(fi)
            if (ifi,i) not in paired_ideal_indices:
                paired_ideal_indices.append((i,ifi))
        else:
            print "DBG: ideal pairing problem"
            return

    #calculate the number of orthogonal elements within ideals
    len_unpaired_orth = []
    for i in unpaired_ideal_indices:
        if factors[i].degree() == 1 and factors[i].list()[0] == field_ring(-1): #x-1 factor
            if field_ring.order()%2 == 0:
                len_unpaired_orth.append(1)
            else:
                len_unpaired_orth.append(2)
        elif factors[i].degree() == 1: #x+1 factor
            if field_ring.order()%2 == 1:
                len_unpaired_orth.append(2)
            else:
                len_unpaired_orth.append(1)
        else:
            len_unpaired_orth.append(1+int(factors[i].degree()/2))

    len_paired_orth = [field_ring.order()**factor_deg[i]-1 for i,_ in paired_ideal_indices]


    #find E_i for every ideal, and create generators for the orthogonal elements
    Eis,unpaired_orth,paired_orth = gen_ideal_orth_elements_genopt(poly_ring,q_ring,q_x,ideal_gens,n,unpaired_ideal_indices,paired_ideal_indices,len_unpaired_orth,len_paired_orth,factor_deg,verbose is not None and verbose,asserts)


    #print stats
    if verbose is not None and verbose:
        print poly_mod, " has %d irreducible factors forming ideals: \n"%(len(factors))
        for i,(f,d,gen) in enumerate(zip(factors,factor_deg,ideal_gens)):
            print "\t",i,": Ideal formed by:",poly_mod,"/",f,", factor degree: ",d
            print "\t Multiplicative unit:",q_to_list_hex(Eis[i],n)
            if i in unpaired_ideal_indices:
                print "\t Unpaired Ideal, gen:",list_field_to_int(gen)
                print "\t Ideal has",len_unpaired_orth[unpaired_ideal_indices.index(i)],"Orthogonal Elements"
            else:
                for (o,p) in paired_ideal_indices:
                    if o == i or p == i:
                        break
                print "\t Paired Ideal ",(o,p),", gen:",list_field_to_int(gen)
                print "\t Ideal pair has",len_paired_orth[paired_ideal_indices.index((o,p))],"Orthogonal Elements"
            print ""
        
        
    #create orthogonal polynomial generator
    orth_gen = cache_gen_product(True,*(unpaired_orth+paired_orth))

    n_mds = 0
    n_orth = 0
    n_mat = 0
    min_xor = None
    min_xor_mds = None

    if verbose is not None and verbose:
        print "Iterating over all",reduce(lambda k,l:k*l,len_unpaired_orth+len_paired_orth),"orthogonal matrices :\n"

    #iterate over all circulant orthogonal matrices
    for o in orth_gen:

        o = functools.reduce(lambda z,y: z+y, o)

        coeff_list = o.list()[::-1]

        n_mat += 1
        if verbose and n_mat%update_every==0:
            print "Checked",n_mat,"orthogonal circulant matrices..."

        #check the matrix XOR count
        if (field_ring.characteristic() == 2):
            table_xor_count = circmat_xor_count(mult_xor_counts,n,r,coeff_list)
        else:
            table_xor_count = 0
        if not find_all and min_xor is not None and min_xor <= table_xor_count:
            continue

        #create matrix
        m = circ(coeff_list)

        #assert it is orthogonal
        if asserts:
            assert(m*m.T == matrix.identity(field_ring,n))
        n_orth += 1

        #check for the MDS property
        if mds(m):
            if min_xor is None or min_xor > table_xor_count:
                min_xor = table_xor_count
                min_xor_mds = coeff_list
            n_mds += 1
            if verbose is not None and verbose:
                print "MDS! : circ(",coeff_list,"), XOR count:",table_xor_count
        else :
            pass

    #print search results
    if verbose is not None:
        print "Checked",n_mat," circulant",n,"x",n," matrices for XOR count in",field_ring,", mod:",field_to_hex_str(GF(2)['x'](mod)),"\n"
        print "Checked",n_orth,"orthogonal circulant",n,"x",n," matrices for MDS property in",field_ring,", mod:",field_to_hex_str(GF(2)['x'](mod)),"\n"
        print "Found",n_mds,"MDS orthogonal circulant",n,"x",n," matrices in",field_ring,", mod:",field_to_hex_str(GF(2)['x'](mod)),"\n"
        if min_xor_mds is not None:
            print "circ(",list_field_to_hex_str(min_xor_mds),"): MDS with lowest XOR count : ",min_xor

    #assert we checked all circulant orthogonal matrices
    if asserts:
        assert(reduce(lambda k,l:k*l,len_unpaired_orth+len_paired_orth)==n_mat)

    #calculate the standard deviation of the field
    if (field_ring.characteristic() == 2):
        std = np.std(np.array(gen_xor_count(r,mod,False).values()))
    else:
        std = None

    return (n_mat,n_orth,n_mds,min_xor,std,min_xor_mds)


# Run a number of searches for COMDS matrices, checking that it all works correctly


start = time.time()
assert(find_COMDS_matrices(5,r=4,verbose=None,find_all=True)[:-1] == (225,225,100,150,2.817356917396161))
assert(find_COMDS_matrices(3,field=GF(7),verbose=None,find_all=True)[:-1] == (12,12,6,0,None))
assert(find_COMDS_matrices(5,field=GF(9),verbose=None,find_all=True)[:-1] == (8,8,0,None,None))
assert(find_COMDS_matrices(7,field=GF(3),verbose=None,find_all=True)[:-1] == (8,8,0,None,None))
assert(find_COMDS_matrices(4,field=GF(7),verbose=None,find_all=True)[:-1] == (8,8,0,None,None))
assert(find_COMDS_matrices(9,field=GF(5),verbose=None,find_all=True)[:-1] == (16,16,0,None,None))
assert(find_COMDS_matrices(4,field=GF(3),verbose=None,find_all=True)[:-1] == (8,8,0,None,None))
assert(find_COMDS_matrices(7,field=GF(2^2),verbose=None,find_all=True)[:-1] == (63,63,0,None,0.5))
assert(find_COMDS_matrices(8,field=GF(3),verbose=None,find_all=True)[:-1] == (64,64,0,None,None))
end = time.time()
print "Tests took %0.2f sec(s)"%(end-start)

# Searches for 5x5 GF(2^r) matrices

# find_COMDS_matrices(5,r=4,verbose=True)
# find_COMDS_matrices(5,r=8,verbose=True)
# find_COMDS_matrices(5,r=12,verbose=True)
# find_COMDS_matrices(5,r=16,verbose=True)
# find_COMDS_matrices(5,field=GF(2^24, 'a'),verbose=True,asserts=False)



#run nxn COMDS search for all modulos in GF(2^r)
def find_COMDS_matrices_GF2r(n,r,verbose=True):
    res = []
    poly_ring, poly_x = GF(2)['p'].objgen()
    for mod in gen_poly_ring_elements(poly_ring,r):
        mod += poly_x^r
        if mod.is_irreducible():
            data = find_COMDS_matrices(n,r,mod=mod.list(),verbose=verbose,find_all=False)
            print "\nSearched for %dx%d COMDS matrices in GF(2^%d)/%s:"%(n,n,r,str(mod))
            print data
            res.append(data)
    return res

#find_COMDS_matrices_GF2r(5,4,verbose=None)

#run nxn COMDS search for all modulos in GF(2^r), plot minxor against field xor stddev
def plot_COMDS_GF2r_minxor_xorstddev(n,r,formats,verbose=True):
    data = [(res[3],res[4]) for res in find_COMDS_matrices_GF2r(n,r,verbose=verbose)]
    data.sort(key=lambda l: l[0])
    plt.scatter(*zip(*data))
    plt.title("Min XORCount %dx%d COMDS matrix in GF(2_%d), all modulos"%(n,n,r))
    plt.xlabel("Min XORCount COMDS in Field")
    plt.ylabel("Field XORCount Standard Deviation")

    for form in formats:
        plt.savefig("results/minxor_stddev_scatter_COMDS_%dx%d_GF2_%d.%s"%(n,n,r,form))

# plot_COMDS_GF2r_minxor_xorstddev(5,4,["png"],verbose=None)
# plot_COMDS_GF2r_minxor_xorstddev(5,6,["png"],verbose=None)
# plot_COMDS_GF2r_minxor_xorstddev(5,8,["png"],verbose=None)



