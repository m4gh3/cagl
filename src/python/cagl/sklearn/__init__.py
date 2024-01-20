''' This file is part of cagl.
 *
 * cagl is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * cagl is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cagl.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Lorenzo Magherini (m4gh3) '''


from sklearn.base import BaseEstimator, TransformerMixin, ClassifierMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from sklearn.metrics import accuracy_score
import numpy as np
import cagl
import time



class HalfSphereCompact(BaseEstimator, TransformerMixin ):
    def __init__(self, alpha=1 ):
        self.alpha = alpha
    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None ):
        scale = np.linalg.norm(X, axis=1, keepdims=True )**2 + self.alpha**2
        extra_alpha = self.alpha*np.ones(X.shape[:-1]+(1,), dtype=X.dtype )
        X = np.concatenate((X, extra_alpha), axis=1 )
        return X/np.sqrt(scale)

    def inverse_transform(self, X ):
        #the last components reveal about the norm:
        #scale = self.alpha/X[:,-1:] 
        return ((self.alpha/X[:,-1:])*X)[:,:-1]



class PXOR(BaseEstimator, ClassifierMixin ):
        

    def __init__(self, dval=0 ):
        super().__init__()
        self.dval = dval

    def fit(self, X, y ):

        ndims = X.shape[1]
        nsamples = X.shape[0]
        
        ctx = cagl.fmpq_mpoly_ctx(2*ndims+2)
        gen_names = [ f"w{i}" for i in range(2*ndims) ] + [ "λ0", "λ1" ]

        x_train = cagl.fmpq_mpoly_matrix(nsamples, ndims, ctx )
        x_train.set_from_np(X.astype(np.double), 8, 8 )

        w = [ cagl.fmpq_mpoly_matrix(ndims, 1, ctx ) for i in range(2) ]
        
        for i in range(2):
            w[i].gens_fill(i*ndims)

        out = (x_train @ w[0])*(x_train @ w[1])

        y_train = cagl.fmpq_mpoly_matrix(1, nsamples, ctx )
        y_train_np = (2*y-1).astype(np.double)
        y_train.set_from_np(y_train_np, 8, 8 )

        err = (y_train @ out)[0,0]

        cons = cagl.fmpq_mpoly(ctx, "0" )

        for i in range(ndims):
            cons += cagl.fmpq_mpoly(ctx, f"λ0*w{i}^2", gen_names )
        for i in range(ndims, 2*ndims ):
            cons += cagl.fmpq_mpoly(ctx, f"λ1*w{i}^2", gen_names )

        cons += cagl.fmpq_mpoly(ctx, "- λ0 - λ1", gen_names )

        lagr = err + cons

        a = cagl.solve_from_gens([lagr.der(i) for i in range(2*ndims+2)])

        time.sleep(0.5)

        eval_ = float('inf')
       
        print(a)

        for i in range(a.shape[0]):
        
            w0 = a[i,0:ndims]
            w1 = a[i,ndims:2*ndims]

            new_eval = - (y_train_np @ ((X @ w0)*(X @ w1)))
       
            print(i, ': ', new_eval )
            print(w0, w1 )
            if new_eval < eval_:
                
                eval_ = new_eval
                self.w0 = w0
                self.w1 = w1
        
        self._is_fitted = True


    def __sklearn_is_fitted__(self):
        return hasattr(self, "_is_fitted") and self._is_fitted


    def predict(self, X ):
        
        check_is_fitted(self)

        X = check_array(X)
        pred = 1*((X @ self.w0)*(X @ self.w1) > self.dval)

        return pred


    def score(self, X, y ):

        check_is_fitted(self)

        X, y = check_X_y(X, y )
        pred = self.predict(X)
        acc  = accuracy_score(y, pred )

        return acc

    def predict_val(self, X ):
        return (X @ self.w0)*(X @ self.w1)
