import cagl
ctx = cagl.fmpq_mpoly_ctx(2)
p = cagl.fmpq_mpoly(ctx)
p.set_str_pretty("x^2+y^2", ("x","y") )
q = cagl.fmpq_mpoly(ctx)
q.set_str_pretty("x+y", ("x","y") )
r = q+p
print(r.get_str_pretty(("x", "y")))
r = q*p
print(r.get_str_pretty(("x", "y")))
r = r.der(0)
r = r.integ(0)
print(r.get_str_pretty(("x", "y")))
r.gen(0)
print(r.get_str_pretty(("x", "y")))

p.set_str_pretty("x", ("x","y") )
q.set_str_pretty("y", ("x","y") )

s = cagl.fmpq_mpoly_matrix(2,2,ctx)
t = cagl.fmpq_mpoly_matrix(2,2,ctx)
s[0,0] = p
s[1,1] = q
t[0,1] = q
t[1,0] = p
v = s @ t

print([ [ v[i,j].get_str_pretty(("x","y")) for j in range(2) ] for i in range(2) ])
