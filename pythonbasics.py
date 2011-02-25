
# dictionaries
x = dict()
x.keys()
x.values()
x["key"] = "value"
x["key2"] = "value2"

# lists
y = list()
firstelement = y[0]
lastelement = y[-1]
subsequence = y[2:4]

# sets
z = set()
z2 = set()
z3 = z - z2
z4 = z + z2
z5 = z.intersection(z2)

y = [1,2,3,4,5]
y = range(1,6)
y.reverse()

# function definition
def function(arg1, arg2):
    print(arg1)
    print("this is arg2: %s" % arg2)
    print("this is arg1 %s and arg2 %s" % (arg1, arg2))
    return arg1

# conditionals
if "this string" in y:
    print("the list y contains the string this string")
else:
    print("did not find this stringin y")

# loops
for myelem in y:
    print(myelem)

# list comprehensions provide a shorthand for common looping/filtering idioms

         [f(i) for i in sequence if pred(i)]

#  is the same as ....

	results = list()
	for i in sequence:
	    if pred(i) is True:
	        result.append(f(i))

# modules are objects. Their functions etc. are members.
import random
random.choice(sdf) # <--- choice is a function provided by the module random

# objects provide functions and attributes as members
d = dict()
d.keys()    # <--- keys() is a method (i.e., member function) of d


