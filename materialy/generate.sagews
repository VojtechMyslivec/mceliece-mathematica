# A little help from Sage to generate irreducible polynomial over GF(2^m)
# and import them to Mathematica (in form of list of coef. i.e. {{1,1,0},{0,0,0},{0,1,1}} and so on).

# This file was created with Tomas Kalvoda and Miro Hroncok assistance
# Thanks to both of you
# ps: This code is written in English to make Miro happy ;]

# ===================================================================
# will print the list in {1,2,...,n} format (without a newline)
def print_item( list_to_print ):
    # double {{ and }} to print actually a { and } respecively
    sys.stdout.write( '{{{}}}'.format(','.join(str(i) for i in list_to_print)) )

# transform the polynomial coefficient coef into a list in "normal" order and pad it with zeros to total length of m
def norm_coef( coef, m ):
    pole = coef.polynomial().coefficients( sparse=False )
    pole += [0] * (m - len(pole))

    return (reversed(pole))

# ===================================================================
# this will define the desired field and ring of polynomials over it
m = 8
T = GF(2^m, 'x')
R = PolynomialRing(T, 'y')
x = T.gen()
y = R.gen()

# -------------------------------------------------------------------
# the irreducible polynomial of the field
g = T.modulus()
pole = reversed(g.coefficients( sparse=False ))
g
print_item(pole)

# -------------------------------------------------------------------
# random irreducible polynomial from the ring
p = R.irreducible_element(3)
print(p)

sys.stdout.write('{ ')
for coef in list(reversed(p.coefficients( sparse=False )))[:-1]:
    pole = norm_coef( coef, m )

    print_item(pole)
    sys.stdout.write(', ')
pole = norm_coef( p.coefficients()[0], m )
print_item(pole)
sys.stdout.write(' }')




