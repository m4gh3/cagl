import cagl
import numpy as np

ctx = cagl.fmpq_mpoly_ctx(6)

gen_names = [ f"w{i}" for i in range(4) ] + [ "λ0", "λ1" ]

w = [ cagl.fmpq_mpoly_matrix(2, 1, ctx ) for i in range(2) ]

for i in range(2):
    w[i].gens_fill(i<<1)

x_train = cagl.fmpq_mpoly_matrix(4, 2, ctx )
x_train.set_from_np( np.array( [
                        [ 1., 1. ],
                        [ 1., -1. ],
                        [ -1., 1. ],
                        [ -1., -1. ]
                    ]), 16 )

y_train = cagl.fmpq_mpoly_matrix(4, 1, ctx )
y_train.set_from_np(np.array( [ [1.], [-1.], [-1.], [1.] ] ), 16 )

out = (x_train @ w[0])*(x_train @ w[1])
err = (y_train - out).sqfrob()

lagr = err + cagl.fmpq_mpoly(ctx, " λ0*w0^2 + λ0*w1^2 - λ0  +  λ1*w2^2 + λ1*w3^2 - λ1", gen_names )

print("lagrangian for the mul problem:")
print(f"\t{lagr.get_str_pretty(gen_names)}")

#compute the gradient of the lagrangian:
print("---------------------------------")
print("lagrangian gradient:")
for i in range(6):
    print(f"\t{lagr.der(i).get_str_pretty(gen_names)}")
print(cagl.solve_from_gens([lagr.der(i) for i in range(6)]))
