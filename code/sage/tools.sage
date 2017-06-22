

#converts an integer into a binary polynomial
def int_to_poly(num,poly_x):
    binary_string = bin(num)[2:]
    return sum([poly_x^i for i,j in enumerate(binary_string[::-1]) if j!='0'])

#int_to_poly applied on lists recursively
def list_int_to_poly(m,field_x):
    if hasattr(m, "__getitem__") and not isinstance(m, basestring):
        return [list_int_to_poly(i,field_x) for i in m]
    return int_to_poly(m,field_x)

#converts a binary polynomial into an integer
def field_to_int(elem,llim=None,ulim=None):
    if hasattr(elem,'list'):
        return int(sum([int(j)*2**i for i,j in enumerate(elem.list())]))
    if isinstance(elem,int):
        return elem
    if hasattr(elem,'_int_repr'):
        return int(elem._int_repr())
    if llim is not None:
        return int(sum([int(j)*2**i for i,j in enumerate(elem.polynomial().list()[llim:])]))
    if ulim is not None:
        return int(sum([int(j)*2**i for i,j in enumerate(elem.polynomial().list()[:ulim])]))
    return int(sum([int(j)*2**i for i,j in enumerate(elem.polynomial().list())]))

#field_to_hex applied on lists recursively
def list_field_to_int(m):
    if hasattr(m, "__getitem__"):
        return [list_field_to_int(i) for i in m]
    return field_to_int(m)

#converts a binary polynomial into a hexadecimal representation string
def field_to_hex_str(elem,llim=None,ulim=None):
    return hex(field_to_int(elem,llim=llim,ulim=ulim))

#field_to_hex_str applied on lists recursively
def list_field_to_hex_str(m):
    if hasattr(m, "__getitem__"):
        return [list_field_to_hex_str(i) for i in m]
    return field_to_hex_str(m)


#converts a galois field element into a polynomial coefficient list of set length r
def to_bool_list(elem,r):
    if elem == 0:
        coeff = []
    else:
        coeff= elem.polynomial().list()
    return map(lambda a, b: b if a is None else a,coeff,[0]*r)

#print matrix with binary polynomials as hexadecimal representation strings
def print_matrix_hex(mat):
    res = ""
    for i in mat:
        res += str(list_field_to_hex_str(i))+'\n'
    print res

#create a binary galois field using the given modulus and field degree
def make_field(r,mod=None):
    poly_ring, poly_x = GF(2**r)['a'].objgen()
    if mod is None:
        return GF(2**r, name='o').objgen()
    return GF(2**r, name='o', modulus=poly_ring(mod)).objgen()

#random vector in the defined field and size, as list of hexadecimal strings 
def random_vector(vector_size,field_ring,field_x):
    return [hex(random.randrange(0,field_ring.order()))[2:] for vi in range(vector_size)]

#turn coefficients into a circulant matrix
#code provided by Elmar Wolfgang Tischhauser
def circ(row):
    field = row[0].parent()
    k = len(row)
    M = matrix(field, k)
    for i in range(k):
        for j in range(k):
            M[i,j] = row[(i - j - 1) % k]
    return M

#check a matrix for the MDS property
#code provided by Elmar Wolfgang Tischhauser
def mds(M):
    assert M.is_square()
    for k in range(M.nrows()+1):
        if 0 in M.minors(k):
            return False
    return True