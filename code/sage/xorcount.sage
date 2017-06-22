import copy, itertools
from collections import Counter

load("tools.sage")

#convert a list of sets into a binary matrix
def sets_to_binary_matrix(sets,r):
    ret = [[0 for _ in range(r)] for _ in sets]
    for i,s in enumerate(sets):
        for bit in s:
            ret[i][bit] = 1
    return matrix(GF(2),ret)

#convert a binary matrix into list of sets
def matrix_to_sets(mat):
    ret = [set() for _ in range(mat.nrows())]
    for ir,row in enumerate(mat):
        for ic,bit in enumerate(row):
            if bit:
                ret[ir].add(ic)
    return ret

#given the bits to be XORed to perform a multiplication, 
#find the the bits to be XORed to perform the inverse of that multiplication
def calculate_inverse_xor_bits(xors,r):
    mat = sets_to_binary_matrix(xors,r)
    if (mat.rank() < r):
        return xors
    reverse_xors = matrix_to_sets(mat.T.inverse().T)

    return reverse_xors

#find the bits that need to be XORed to perform a multiplication with the given element in GF(2^r)/[mod_bits] 
def xor_bits(elem,r,mod_bits):
    res = [set() for _ in range(r)]
    overflow = [set() for _ in range(r-1)]

    b = [i for i in range(r)]

    bits = to_bool_list(elem,r)

    for i in range(r):
        if bits[i]:
            for j in range(r):
                if j<i:
                    overflow[r-1-i+j].add(b[j])
                else:
                    res[j-i].add(b[j])


    mod_bits = list(reversed(mod_bits))

    #reduce
    for i in range(r-1):
        reduction = copy.copy(overflow[i])
        if len(overflow[i])>0:
            for k in range(r+1):
                if mod_bits[k] and i+k<r-1:
                    overflow[i+k] = overflow[i+k].symmetric_difference(reduction)
                elif mod_bits[k] and i+k>=r-1:
                    res[i+k-r+1] = res[i+k-r+1].symmetric_difference(reduction)
    return res


#calculate the XOR count for multiplication by elem, potentially first converted into the bits to be XORed
def xor_count(elem,r,mod_bits=None):
    if not isinstance(elem,list):
        assert(mod_bits is not None)
        elem = xor_bits(elem,r,mod_bits)
    cnt = sum([len(elem[i])-1 if len(elem[i])>1 else 0 for i in range(r)])
    return cnt


#calculate (or create a generator) the XOR bits for all elements in GF(2^r)/[mod] 
def gen_xor_bits(r,mod,big_field_opt=False):

    _,field_x = make_field(r,mod)
    
    if not big_field_opt:
        mult_xor_counts = {}
        for e in range(2^r):
            elem = int_to_poly(e,field_x)
            mult_xor_counts[e] = xor_bits(elem,r,mod)
        return mult_xor_counts

    class delay_dict(dict):
        def __init__(self,r,mod,field_x):
            self.__dict__ = {}
            self.r = r
            self.mod = mod
            self.field_x = field_x
            
        def __getitem__(self, key):
            if key not in self.__dict__:
                elem = int_to_poly(key,self.field_x)
                self.__dict__[key] = xor_bits(elem,self.r,self.mod)
            return self.__dict__[key]
        def values(self):
            return None

    return delay_dict(r,mod,field_x)

#calculate (or create a generator) the XOR count for all elements in GF(2^r)/[mod] 
def gen_xor_count(r,mod,big_field_opt=False):

    _,field_x = make_field(r,mod)
    
    if not big_field_opt:
        mult_xor_counts = {}
        for e in range(2^r):
            elem = int_to_poly(e,field_x)
            mult_xor_counts[e] = xor_count(elem,r,mod)
        return mult_xor_counts

    class delay_dict(dict):
        def __init__(self,r,mod,field_x):
            self.__dict__ = {}
            self.r = r
            self.mod = mod
            self.field_x = field_x
            
        def __getitem__(self, key):
            if key not in self.__dict__:
                elem = int_to_poly(key,self.field_x)
                self.__dict__[key] = xor_count(elem,self.r,self.mod)
            return self.__dict__[key]
        def values(self):
            return None

    return delay_dict(r,mod,field_x)

#sum up XOR counts of elemeents in a list, mult_xor_counts being the output of gen_xor_count
def list_xor_count(mult_xor_counts,l):
    return sum([mult_xor_counts[field_to_int(e)] for e in l])

#returns XOR count of a circulant matrix (or any matrix with rows being permutations of each other)
def circmat_xor_count(mult_xor_counts,n,r,l):
    return n * (list_xor_count(mult_xor_counts,l)+(n-l.count(0)-1)*r)


# #test XOR count generation against table in Sim, Khoo (2015) for GF(2^8)
def assert_XORcount_GF256_correct():
    F,xx = GF(2)['xx'].objgen()
    xors = [[0x11d,0x177,0x1f3,0x169,0x1bd,0x1e7,0x12b,0x1d7],[0,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[2,3,5,5,3,5,5,3,5],[3,11,13,11,11,11,11,11,11],[4,6,10,7,6,9,7,6,9],[5,14,14,13,10,15,13,14,15],[6,13,15,14,15,14,14,15,14],[7,21,19,18,19,22,18,23,18],[8,9,14,9,10,12,9,9,13],[9,17,16,15,12,16,15,11,17],[10,16,17,16,11,17,16,16,16],[11,24,19,20,13,19,20,18,22],[12,15,16,18,18,17,18,19,16],[13,23,22,22,24,19,22,21,22],[14,22,19,19,21,22,19,24,19],[15,30,25,21,27,26,21,26,27],[16,12,15,11,15,15,11,13,16],[17,12,19,17,19,21,13,21,18],[18,21,18,18,14,20,16,12,19],[19,21,22,22,18,24,20,20,19],[20,18,21,16,11,18,18,17,17],[21,18,21,20,19,22,22,25,21],[22,27,20,21,12,19,21,18,24],[23,27,20,23,20,25,27,26,26],[24,17,17,22,21,19,22,24,17],[25,17,19,26,27,21,26,26,19],[26,22,26,23,26,20,21,23,22],[27,22,28,25,32,20,27,25,26],[28,23,19,21,23,22,23,24,20],[29,23,25,23,25,22,29,26,20],[30,32,24,20,26,27,20,25,29],[31,32,30,20,28,29,28,27,31],[32,16,16,13,19,18,16,17,19],[33,14,22,13,19,18,22,15,23],[34,13,21,20,22,23,13,24,20],[35,11,27,22,22,25,21,22,22],[36,26,20,22,17,21,17,13,22],[37,24,22,24,21,23,21,11,24],[38,23,25,23,18,26,20,22,19],[39,21,27,27,22,26,26,20,19],[40,21,24,16,11,24,23,18,18],[41,19,28,18,13,28,27,22,18],[42,18,23,21,20,21,22,25,21],[43,16,27,25,22,27,28,29,23],[44,27,20,23,11,19,22,18,27],[45,25,28,27,9,25,24,22,29],[46,28,19,22,22,24,27,27,26],[47,26,27,28,20,28,31,31,30],[48,18,19,26,24,21,23,26,19],[49,24,25,28,28,19,27,24,25],[50,17,18,27,27,22,28,29,18],[51,23,24,31,31,22,30,27,22],[52,24,27,25,28,22,20,24,22],[53,30,29,29,28,22,26,22,30],[54,19,30,24,33,19,27,25,29],[55,25,32,30,33,17,31,23,35],[56,23,19,23,24,23,24,25,22],[57,29,19,27,22,25,30,29,28],[58,22,28,22,25,20,31,24,19],[59,28,28,28,23,24,35,28,27],[60,33,23,20,26,28,19,27,31],[61,39,27,26,28,32,27,31,35],[62,32,28,17,25,29,28,24,28],[63,38,32,25,27,31,34,28,34],[64,19,17,16,24,20,20,22,22],[65,15,15,22,30,20,18,24,18],[66,16,26,13,21,19,25,17,25],[67,12,24,21,27,17,25,19,23],[68,17,23,23,24,25,13,26,23],[69,13,25,27,26,23,9,28,21],[70,10,28,22,23,28,24,23,22],[71,6,30,28,25,28,22,25,22],[72,30,23,23,20,22,21,15,25],[73,26,27,27,24,26,21,23,25],[74,27,22,26,23,25,20,10,26],[75,23,26,32,27,27,22,18,24],[76,24,25,24,18,27,20,23,20],[77,20,25,26,26,29,18,31,18],[78,21,28,29,19,26,25,20,17],[79,17,28,33,27,30,25,28,13],[80,25,26,21,11,25,25,19,24],[81,29,24,25,13,27,25,21,26],[82,20,29,16,16,32,30,26,17],[83,24,27,22,18,32,28,28,21],[84,19,24,22,21,20,22,25,21],[85,23,26,24,27,20,24,27,21],[86,14,27,23,24,27,29,30,22],[87,18,29,27,30,29,29,32,24],[88,28,20,26,11,19,24,18,27],[89,32,20,28,19,25,22,26,29],[90,23,29,27,6,26,23,21,30],[91,27,29,31,14,30,19,29,30],[92,26,18,21,23,26,27,28,26],[93,30,14,21,27,30,27,36,30],[94,25,27,28,20,25,28,29,29],[95,29,23,30,24,31,26,37,31],[96,19,21,27,27,24,24,25,21],[97,25,21,29,25,24,24,29,21],[98,26,26,30,28,19,29,24,28],[99,32,26,30,26,21,27,28,30],[100,17,17,28,27,23,31,31,20],[101,23,21,32,29,25,29,35,22],[102,24,26,33,30,22,30,28,19],[103,30,30,35,32,22,26,32,23],[104,26,29,24,29,24,19,26,22],[105,32,31,28,33,20,21,24,26],[106,29,28,29,28,23,26,21,31],[107,35,30,31,32,21,26,19,33],[108,16,29,23,35,19,28,28,27],[109,22,27,29,35,17,28,26,29],[110,23,32,30,32,14,29,21,36],[111,29,30,34,32,10,27,19,36],[112,23,22,26,24,25,25,26,23],[113,21,18,30,26,27,23,30,21],[114,32,17,27,21,24,30,29,28],[115,30,13,29,23,28,30,33,28],[116,21,26,21,26,18,32,22,18],[117,19,26,27,24,22,32,26,14],[118,26,29,28,21,25,35,27,27],[119,24,29,32,19,27,37,31,25],[120,34,22,21,26,29,18,29,32],[121,32,24,27,26,27,14,27,30],[122,39,27,24,29,28,25,32,35],[123,37,29,28,29,28,23,30,31],[124,32,26,14,22,30,27,21,25],[125,30,24,22,26,30,25,19,25],[126,37,31,23,27,29,32,26,32],[127,35,29,29,31,27,32,24,30],[128,21,23,20,29,22,23,26,26],[129,29,25,16,29,28,25,20,26],[130,16,14,23,32,25,20,27,19],[131,24,16,21,32,29,24,21,17],[132,17,27,13,23,21,28,20,27],[133,25,25,7,27,25,28,14,29],[134,12,26,22,28,16,27,19,24],[135,20,24,18,32,22,29,13,24],[136,20,25,27,25,28,16,27,27],[137,28,21,25,27,30,20,27,27],[138,15,26,28,26,23,7,28,22],[139,23,22,28,28,23,13,28,24],[140,12,29,22,25,31,27,25,22],[141,20,29,18,23,31,29,25,20],[142,3,30,29,24,26,20,24,21],[143,11,30,27,22,28,24,24,21],[144,35,26,25,22,23,22,17,28],[145,35,24,23,26,27,26,11,34],[146,28,27,26,25,26,23,26,27],[147,28,25,26,29,28,25,20,31],[148,27,22,28,26,28,19,9,25],[149,27,16,24,26,30,25,3,29],[150,24,27,31,27,27,26,20,24],[151,24,21,29,27,31,30,14,26],[152,26,28,26,18,29,21,24,21],[153,26,24,26,16,29,23,24,23],[154,19,23,25,27,28,16,29,18],[155,19,19,27,25,26,16,29,22],[156,22,28,31,16,26,24,20,20],[157,22,28,29,18,24,28,20,24],[158,15,27,32,27,29,25,27,9],[159,15,27,32,29,29,27,27,15],[160,29,29,23,12,26,27,19,25],[161,27,29,23,20,24,31,23,29],[162,30,24,28,11,29,28,20,30],[163,28,24,26,19,29,30,24,32],[164,21,31,14,18,31,30,27,16],[165,19,27,16,22,31,32,31,22],[166,26,26,25,19,34,29,30,21],[167,24,22,25,23,32,29,34,25],[168,20,25,24,22,22,22,28,24],[169,18,23,22,24,24,28,26,28],[170,25,26,23,27,17,25,25,19],[171,23,24,19,29,21,29,23,25],[172,12,27,21,26,27,27,32,21],[173,10,29,21,32,31,31,30,23],[174,17,28,26,29,30,28,31,24],[175,15,30,24,35,32,30,29,28],[176,29,20,26,11,19,26,18,27],[177,35,20,24,15,15,28,22,29],[178,32,21,29,22,26,23,27,30],[179,38,21,25,26,24,27,31,30],[180,21,30,27,3,26,25,20,30],[181,27,26,27,11,24,29,24,30],[182,24,27,32,12,29,16,27,29],[183,30,23,30,20,25,22,31,27],[184,24,20,21,25,27,27,29,26],[185,30,14,17,31,27,27,27,24],[186,31,11,18,26,30,26,38,31],[187,37,5,12,32,32,28,36,31],[188,24,26,28,19,22,28,27,27],[189,30,24,26,21,24,30,25,27],[190,27,21,27,22,29,21,34,28],[191,33,19,23,24,29,25,32,30],[192,20,26,28,29,28,25,24,22],[193,16,26,26,27,30,25,24,22],[194,25,21,31,26,25,26,29,21],[195,21,21,27,24,25,24,29,23],[196,28,26,29,29,19,28,24,31],[197,24,30,25,31,19,26,24,29],[198,33,25,30,24,20,27,31,30],[199,29,29,24,26,22,23,31,30],[200,17,16,29,27,24,34,33,21],[201,13,22,29,31,30,32,27,21],[202,26,21,34,30,25,29,34,22],[203,22,27,32,34,29,25,28,20],[204,25,28,32,29,23,31,29,16],[205,21,30,30,29,27,27,23,18],[206,30,29,35,34,20,24,32,25],[207,26,31,31,34,26,18,26,25],[208,28,31,23,30,23,18,27,22],[209,32,31,23,32,27,12,27,16],[210,31,32,28,31,20,23,24,27],[211,35,32,26,33,22,19,24,23],[212,28,27,30,28,24,25,21,31],[213,32,31,28,26,26,21,21,27],[214,35,28,29,31,21,24,16,32],[215,39,32,25,29,25,22,16,30],[216,13,29,22,36,19,29,30,25],[217,17,31,24,36,27,25,24,23],[218,20,24,29,35,16,28,27,28],[219,24,26,29,35,22,26,21,24],[220,21,33,31,32,16,30,20,36],[221,25,31,31,36,22,28,14,32],[222,28,28,32,29,5,23,15,35],[223,32,26,30,33,13,23,9,29],[224,24,24,29,24,26,25,27,25],[225,30,26,27,30,28,23,25,21],[226,19,19,30,29,27,22,32,20],[227,25,21,30,35,31,22,30,18],[228,32,18,28,20,23,30,29,28],[229,38,24,28,22,27,26,27,22],[230,31,9,27,23,28,29,32,27],[231,37,15,29,25,30,27,30,23],[232,21,24,20,28,16,32,20,18],[233,27,28,16,32,14,28,24,14],[234,16,25,27,23,21,31,25,11],[235,22,29,25,27,21,29,29,5],[236,25,30,25,18,25,35,26,27],[237,31,30,23,26,25,29,30,25]]
    mods = [int_to_poly(xor,xx).list() for xor in xors[0]]
    xorcounts = [gen_xor_count(8,mod,big_field_opt=False) for mod in mods]
    for xor in xors[1:]:
        val = xor[0]
        for i,(mod,pval) in enumerate(zip(mods,xor[1:])):
            assert(xorcounts[i][val]==pval)

assert_XORcount_GF256_correct()

#print XOR count table for GF(2^r)/mod for every mod in mods
#mods expected as integer representation of the binary polynomial
def print_XOR_count_table(r,mods):
    F,xx = GF(2)['xx'].objgen()
    xor_counts = [gen_xor_count(r,int_to_poly(mod,xx).list()) for mod in mods]

    print "mod:\t",
    for mod in mods:
            print "0x%x\t"%(mod),
    print

    for e in range(2^r):
        print e,":\t",
        for xor_count in xor_counts:
            print xor_count[e],"\t",
        print

# print_XOR_count_table(r=4,mods=[0x13,0x19,0x1f])
# print_XOR_count_table(r=8,mods=[0x11d,0x177,0x1f3,0x169,0x1bd,0x1e7,0x12b,0x1d7])

#calculate xor bits for multiplying by init*applied, both provided as xor bits
def combine_factors(init_xors,applied_xors):
    # print(init_xors,applied_xors)
    ret = [set() for _ in range(len(init_xors))]

    for i,bit_xors in enumerate(applied_xors):
        for bit in bit_xors:
            ret[i] = ret[i].symmetric_difference(init_xors[bit])
    return ret

#multiplying by a ratio means multiplying by x/y = x*y^(-1)
#calculate the XOR bits and count for all ratios in GF(2^r)/mod
def calculate_ratio_XOR_table_and_count(r,mod):

    field_ring,field_x = make_field(r,mod)
    mult_xors = {}
    div_xors = {}
    for e in range(2^r):
        elem = int_to_poly(e,field_x)
        # print hex(e),count
        mult_xors[e] = xor_bits(elem,r,mod)
        div_xors[e] = calculate_inverse_xor_bits(mult_xors[e],r)

    ratio = {}
    ratio_xor_counts = {}

    for multi in range(1,2^r):
        for divi in range(1,2^r):
            ratio[(multi,divi)] = combine_factors(mult_xors[multi],div_xors[divi])
            ratio_xor_counts[(multi,divi)] = xor_count(ratio[(multi,divi)],r)

    return (ratio,ratio_xor_counts)

#print the xor counts for ratios in GF(2^r)/mod, where ratio = x/y and y<div_print_limit
def print_XOR_RATIO_count_table(r,mod,div_print_limit):


    ratio,ratio_xor_counts = calculate_ratio_XOR_table_and_count(r,mod)

    for multi in range(0,2^r):
        if multi == 0:
            print(" DIV "+" ".join(["%02x"%(divi) for divi in range(0,div_print_limit)]))
            print("MULT")
            print("%02x:   -- "%(multi)+" ".join(["%02d"%(0) for _ in range(1,div_print_limit)]))
        else:
            print("%02x:   -- "%(multi)+" ".join(["%02d"%(ratio_xor_counts[(multi,divi)]) for divi in range(1,div_print_limit)]))
    

# print_XOR_RATIO_count_table(8,[1,0,1,1,1,0,0,0,1],div_print_limit=20)