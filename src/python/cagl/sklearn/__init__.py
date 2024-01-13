from sklearn.base import BaseEstimator, TransformerMixin
import numpy as np

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
