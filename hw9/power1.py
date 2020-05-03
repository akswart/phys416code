
import numpy as np

def power1(A,X,eps,max1):
# function to compute the maximum eigenvalue (lambda) of a Matrix A and
# its corresponding eigen vector (eigenvector)
# X is some base vector row matrix, input
# eps is the tolerance you want the eigenvalue to be computed to
# max1 is the maximum number of iteractions allowed
    lamda=0.0
    cnt=1
    xnew=X
    xold=0*X
    err=np.linalg.norm(X)
    # xnew = xnew/max(abs(xnew))
    xnew = xnew/np.linalg.norm(xnew)
    
    while err>eps and cnt < max1:
        xold=xnew
        # ck = max(abs(xnew))
        ck = np.linalg.norm(xnew)
        xnew=np.dot(A,xnew)/ck
        cnt = cnt+1
        err = np.linalg.norm(xold-xnew)
    
    if (cnt >=max1):
        print('error in power2, max number of iterations exceeded')
    
    eigenvector = xnew/np.linalg.norm(xnew)
    lamda = ck

    return lamda,eigenvector

M=np.array([[2., -1., -1.], [-1., 2., -1], [1., -1., 2]])
x=np.array([[3.], [2.], [1.]])

eigenvalue,eigenvector = np.linalg.eig(M)
# get max eigenvalue
emax = np.max(eigenvalue)
emax_index = np.argmax(eigenvalue)
evmax = eigenvector[:,emax_index]


eigenvalue1,eigenvector1=power1(M,x,1.0e-3,20)

print('max eigenvalue from numpy eig =',emax)
print('corresponding eigenvector =',evmax)
print('max eigenvalue from powermethod =',eigenvalue1)
print('corresponding eigenvector =',eigenvector1)

