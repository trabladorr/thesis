import binascii, copy, subprocess, random, os

load("xorcount.sage")

def x2(r,mod,num):
    field_ring,field_x = make_field(r,mod)
    v = int_to_poly(num,field_x)
    x2 = int_to_poly(0x02,field_x)
    return field_to_hex_str(v*x2)
assert(x2(8,[1,0,1,0,0,1,1,0,1],0xf0)=='0x85')

def mat_vec_mult(r,mod,mat_hex,vector_hex):
    field_ring,field_x = make_field(r,mod)
    mat = matrix(field_ring,list_int_to_poly(mat_hex,field_x))
    v = matrix(list_int_to_poly(vector_hex,field_x)).T
    # print_matrix_hex(mat)
    # print_matrix_hex(v)
    # print_matrix_hex(mat*v)
    return mat*v


def mix_columns_multiplication(vector_hex):
    m = [[2,3,1,1],[1,2,3,1],[1,1,2,3],[3,1,1,2]]
    #irr poly: x8 + x4 + x3 + x + 1, 
    z = mat_vec_mult(8,[1,1,0,1,1,0,0,0,1],m,vector_hex)
    return list_field_to_hex_str(z.T[0])
assert(mix_columns_multiplication([0xdb,0x13,0x53,0x45])==['0x8e','0x4d','0xa1','0xbc'])

def inv_mix_columns_multiplication(vector_hex):
    m = [[14,11,13,9],[9,14,11,13],[13,9,14,11],[11,13,9,14]]
    z = mat_vec_mult(8,[1,1,0,1,1,0,0,0,1],m,vector_hex)
    return list_field_to_hex_str(z.T[0])
assert(inv_mix_columns_multiplication([0x8e,0x4d,0xa1,0xbc])==['0xdb','0x13','0x53','0x45'])

def h1_multiplication(vector_hex):
    #Anubis Hadamard matrix
    m = [[0x01,0x02,0x04,0x06],[0x02,0x01,0x06,0x04],[0x04,0x06,0x01,0x02],[0x06,0x04,0x02,0x01]]
    #irr poly: 0x11d, 100011101, x8 + x4 + x3 + x2 + x0
    z = mat_vec_mult(8,[1,0,1,1,1,0,0,0,1],m,vector_hex)
    return list_field_to_hex_str(z.T[0])

def h2_multiplication(vector_hex):
    #Hadamard matrix
    m = [[0x01,0x02,0xb0,0xb2],[0x02,0x01,0xb2,0xb0],[0xb0,0xb2,0x01,0x02],[0xb2,0xb0,0x02,0x01]]
    #irr poly: 0x165, 101100101, x8 + x6 + x3 + x2 + x0
    z = mat_vec_mult(8,[1,0,1,0,0,1,1,0,1],m,vector_hex)
    return list_field_to_hex_str(z.T[0])


assert((h1_multiplication([0x01,0x02,0x03,0x04]))==['0x11', '0x1a', '0x3', '0xc'])
assert((h1_multiplication([0x10,0x20,0x30,0x40]))==['0xd', '0xbd', '0x30', '0xc0'])
assert((h2_multiplication([0x01,0x02,0x03,0x04]))==['0xb2', '0xb9', '0xba', '0xb5'])
assert((h2_multiplication([0x10,0x20,0x30,0x40]))==['0x8', '0xb8', '0x88', '0x78'])

def print_bytearray(input,bytes_per_row=16):

    for i in range(0,len(input),bytes_per_row):
        if i+bytes_per_row>=len(input):
            endi = len(input)
        else:
            endi = i+bytes_per_row
        print ' '.join(format(x, '02x') for x in input[i:endi])

def var_hex_int(values):
    ret = 0
    mult = int(256^len(values))
    for v in values:
        mult /= int(256)
        if mult == 0:
            ret += int(hex(v),16)
        else:
            ret += int(hex(v),16)*mult
    return ret

def chunk(iterable,chunksize=1):
    return zip(*[iter(iterable)]*chunksize)

#checks that the c implementation for correctnesss by using mult_io, and random vectors
def test_c_mult(n_tests,verbose):

    vector_len = 4
    mat = {0:[[0x01,0x02,0x04,0x06],[0x02,0x01,0x06,0x04],[0x04,0x06,0x01,0x02],[0x06,0x04,0x02,0x01]],\
        1:[[0x01,0x02,0xb0,0xb2],[0x02,0x01,0xb2,0xb0],[0xb0,0xb2,0x01,0x02],[0xb2,0xb0,0x02,0x01]]}
    

    for gf2_order,n_vectors,table_type,mult_type,impl_type,mod in [\
        (8,4,0,0,0,[1,0,1,1,1,0,0,0,1]),\
        (8,4,0,0,1,[1,0,1,1,1,0,0,0,1]),\
        (8,16,0,0,0,[1,0,1,1,1,0,0,0,1]),\
        (8,16,0,1,0,[1,0,1,1,1,0,0,0,1]),\
        (8,16,1,0,0,[1,0,1,0,0,1,1,0,1]),\
        (8,16,1,1,0,[1,0,1,0,0,1,1,0,1]),\
        (16,2,0,0,0,GF(2)['x'].irreducible_element(16).list()),\
        (16,2,0,0,1,GF(2)['x'].irreducible_element(16).list()),\
        (8,128,0,2,0,[1,0,1,1,1,0,0,0,1]),\
        (8,128,0,2,2,[1,0,1,1,1,0,0,0,1]),\
        (8,128,0,2,3,[1,0,1,1,1,0,0,0,1])]:

        field_ring,field_x = make_field(gf2_order,mod)
        
        elem_data_len = int(gf2_order/8)
        vector_data_len = vector_len*elem_data_len
        random_data_len = n_tests*n_vectors*vector_data_len

        random_data = bytearray(os.urandom(random_data_len))

        if verbose:
            print("Generated %d bytes, in %d vectors:"%(random_data_len,n_vectors))
            print_bytearray(random_data,random_data_len/n_vectors)


        vector_hex = [var_hex_int(i) for i in chunk(random_data,elem_data_len)]
        vectors = [vector_hex[i:i+vector_len] for i in range(0,random_data_len,vector_len)]


        if verbose:
            print("./mult_io %d %d %d %d %d %d"%(gf2_order,n_vectors,n_tests,table_type,mult_type,impl_type))
        process=subprocess.Popen(['../c/mult_io',str(gf2_order),str(n_vectors),str(n_tests),str(table_type),str(mult_type),str(impl_type)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdoutdata,stderrdata=process.communicate(input=random_data)

        
        if len(stdoutdata) != random_data_len and verbose:
            print("C Debug:\n")
            print(str(stdoutdata[random_data_len:]))

        correct = True
        for i,v in zip(range(0,random_data_len,vector_data_len),vectors):

            z = mat_vec_mult(gf2_order,mod,mat[table_type],v).T[0]
            vector_output = bytearray(stdoutdata[i:i+vector_data_len])
            recvd = [var_hex_int(vec) for vec in chunk(vector_output,elem_data_len)]
            if (verbose):
                print(list_field_to_hex_str(list_int_to_poly(v,field_x)),list_field_to_hex_str(z),list_field_to_hex_str(list_int_to_poly(recvd,field_x)))
            correct = correct and list_field_to_int(z)==recvd

        if len(stderrdata) > 0:
            print("Error:\n",str(stderrdata))

        assert(correct)


test_c_mult(5,False)

#print c intrinsics code for GF(2^8) multiplication by mat in bitslice arrangement
#if intermediates_opt is set to True, attempt to find all common intermediates, and print them as temporary variables
def generate_bitslice_xormult_code(r,mod,mat,intermediates_opt=True):
    bits = r*len(mat)

    xor_bits_dict = gen_xor_bits(r,mod,big_field_opt=True)

    xors = [[] for i in range(bits)]
    for ir,row in enumerate(mat):
        for ic,mult in enumerate(row):
            for i in range(r):
                xors[ir*r+i].append(xor_bits_dict[(mult)][i])

    for i in range(bits):
        xors[i] = [list(j) for j in xors[i]]
        xors[i] = [item+b*r for b,sublist in enumerate(xors[i]) for item in sublist]

    # print without intermediates
    if not intermediates_opt:
        for i in range(bits):
            print("output[%d] = %s;"%(i," ^ ".join(["input[%d]"%(xor) for xor in xors[i]])))
        return

    #intermediates optimisation
    pairs = [[] for i in range(bits)]
    pair_count = Counter()
    for i in range(bits):
        pairs[i] = set(itertools.combinations(xors[i], 2))
        pair_count.update(pairs[i])


    marked_bits = set()
    tmps = []
    xors_saved = 0

    while True:
        
        [(pair,count)] = pair_count.most_common(1)
        pair_count[pair] = 0
        if count < 2:
            break
        tmps.append(pair)
        xors_saved += count - 1

        pair_count = Counter()
        for i in range(bits):
            if pair in pairs[i]:
                xors[i].remove(pair[0])
                xors[i].remove(pair[1])
                xors[i].append("tmp%d"%(len(tmps)-1))
                pairs[i] = set(itertools.combinations(xors[i], 2))
            pair_count.update(pairs[i])

    def var_to_str(x):
        if isinstance(x,str):
            return x
        elif isinstance(x,sage.rings.integer.Integer):
            return "input[%d]"%(x)
        print type(x)
        assert(False)

    # print with intermediates
    for i,t in enumerate(tmps):
        print("__m128i tmp%d = %s;"%(i," ^ ".join([var_to_str(xor) for xor in t])))
    for i in range(bits):
        print("output[%d] = %s;"%(i," ^ ".join([var_to_str(xor) for xor in xors[i]])))


# mat = [[0x01,0x02,0x04,0x06],[0x02,0x01,0x06,0x04],[0x04,0x06,0x01,0x02],[0x06,0x04,0x02,0x01]] #anubis
# generate_bitslice_xormult_code(8,[1,0,1,1,1,0,0,0,1],mat)


#macro generation 
# for r in [range(i*4,i*4+4) for i in range(8)]:
#     print("sse_trans_slice_x4_macro("+",".join(["m128i_input%d"%(i) for i in r])+");\\")


#print multiplication shuffles for anubis, gf(16)
# field_ring,field_x = make_field(16,GF(2)['x'].irreducible_element(16).list())
# for s,z in zip(['x2','x4','x6'],map(lambda l:int_to_poly(l,field_x),[0x02,0x04,0x06])):
#   for upper,upperstr in zip([False,True],['lower','upper']):
#       for rng,rngstr in zip([range(0,16),range(0,256,16),range(0,4096,256),range(0,65536,4096)],("cdef", "89ab", "4567", "0123")):
#           print "const __m128i mul_"+s+"_"+upperstr+"_"+rngstr+" = _mm_setr_epi8(",
#           if upper:
#               print ', '.join(map(lambda l:field_to_hex_str(z*int_to_poly(l,field_x),llim=8),rng)),
#           else:
#               print ', '.join(map(lambda l:field_to_hex_str(z*int_to_poly(l,field_x),ulim=8),rng)),
#           print ");"

#print multiplication shuffle, anubis, 6div4
# field_ring,field_x = make_field(8,[1,0,1,1,1,0,0,0,1])
# for z in map(lambda l:int_to_poly(l,field_x),[0x06]):
#   print "const __m128i mul_"+field_to_hex_str(z)[1:]+"div4_low = _mm_setr_epi8(",
#   print ', '.join(map(lambda l:field_to_hex_str(z*int_to_poly(l,field_x)/int_to_poly(0x04, field_x)),range(0,16))),
#   print ");"
#   print "const __m128i mul_"+field_to_hex_str(z)[1:]+"div4_high = _mm_setr_epi8(",
#   print ', '.join(map(lambda l:field_to_hex_str(z*int_to_poly(l,field_x)/int_to_poly(0x04, field_x)),range(0,256,16))),
#   print ");"

# # print multiplication shuffles for anubis
# field_ring,field_x = make_field(8,[1,0,1,1,1,0,0,0,1])
# for z in map(lambda l:int_to_poly(l,field_x),[0x02,0x04,0x06]):
#   print "const __m128i mul_"+field_to_hex_str(z)[1:]+"_low = _mm_setr_epi8(",
#   print ', '.join(map(lambda l:field_to_hex_str(z*int_to_poly(l,field_x)),range(0,16))),
#   print ");"
# for z in map(lambda l:int_to_poly(l,field_x),[0x02,0x04,0x06]):
#   print "const __m128i mul_"+field_to_hex_str(z)[1:]+"_high = _mm_setr_epi8(",
#   print ', '.join(map(lambda l:field_to_hex_str(z*int_to_poly(l,field_x)),range(0,256,16))),
#   print ");"

# # print multiplication shuffles for h2
# field_ring,field_x = make_field(8,[1,0,1,0,0,1,1,0,1])
# for z in map(lambda l:int_to_poly(l,field_x),[0x02,0xb0,0xb2]):
#   print "const __m128i mul_"+field_to_hex_str(z)[1:]+"_low = _mm_setr_epi8(",
#   print ', '.join(map(lambda l:field_to_hex_str(z*int_to_poly(l,field_x)),range(0,16))),
#   print ");"
# for z in map(lambda l:int_to_poly(l,field_x),[0x02,0xb0,0xb2]):
#   print "const __m128i mul_"+field_to_hex_str(z)[1:]+"_high = _mm_setr_epi8(",
#   print ', '.join(map(lambda l:field_to_hex_str(z*int_to_poly(l,field_x)),range(0,256,16))),
#   print ");"


#print additional shuffles
# for i in range(8):
#   for j in range(4):
#       print ', '.join(["%d%d%d"%(b,j,i) for b in range(8)])
# print
# for i in range(8):
#   for j in range(4):
#       print ', '.join(["%d%d%d"%(b,j,i) for b in range(8)])

# for z in range(0,4):
#   print "(",
#   for i in range(0,4):
#       for j in range (0,4):
#           print  "0x%x,"%(i+z*4),
#   print ");"


#print direct rearrange shuffles
# for i in range(0,16,4):
#   for j in map(lambda l:int(l),[1,0,3,2]):
#       print  "%s,"%(hex(i+j)),
# print
# for i in range(0,16,4):
#   for j in map(lambda l:int(l),[2,3,0,1]):
#       print  "%s,"%(hex(i+j)),
# print
# for i in range(0,16,4):
#   for j in map(lambda l:int(l),[3,2,1,0]):
#       print  "%s,"%(hex(i+j)),
# print
