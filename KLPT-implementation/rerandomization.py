from sage.all import Zmod, ZZ
import random

def KLPT_input_rerandomization(I,order,Nc):
    w = order.random_element()
    nt = w.reduced_characteristic_polynomial().discriminant()
    while Nc.divides(nt) or (Zmod(Nc)(nt).is_square()):
        w = order.random_element()
        nt = w.reduced_characteristic_polynomial().discriminant()
    inter = I.intersection(order*1)
    notin = Nc*I.right_order()
    x = inter.random_element()
    while(x in notin):
        x = inter.random_element()
    r = random.randrange(0,Nc+1)
    s = 1
    if(r==Nc):
        r,s=1,0
    g = (r+s*w)*x
    J = I*(g.conjugate()/(I.norm()))
    return(J, ZZ(J.norm()),g)



    
