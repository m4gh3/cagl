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
