for i in range(0,10):
    print i
    
def func(x):
    return x*x+2

for i in range(0,10):
    print func(i)
    print "%d\t%d" % (i, func(i))