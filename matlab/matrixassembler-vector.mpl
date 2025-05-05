# Matrices assembling for Gauss-Bonnet black holes's quasi-normal modes computation
# Non-extreme case, spin s = 1

MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  N::integer, # N : number of Tchebyshev modes
  L::numeric, # L : angular momentum, L >= 1
  n::numeric, # n : number of compactified dimensions
  ah::numeric,# ah: physical string coupling parameter
  p::string   # p : string of characters containing the path where we save the assembled matrices
  )
  local f::function, fp::function, rhoh::numeric, alpha::numeric, kappa::numeric, lambda::numeric, a:: numeric, fpp::function, fppp::function, UT::function, KT::function, KTp::function, KTpp::function, qT::function, P::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string:

  with(LinearAlgebra): # To work with matrices

  Digits := d: # Setup multiprecision arithmetics

  # Determination of the event horizon:
  rhoh := fsolve(rhoh^(n - 1) + ah*rhoh^(n - 3) - 2 = 0, rhoh = 0 .. 3):
  printf("Found BH's event horizon rhoh = %f\n", rhoh):

  # Physical parameters:
  kappa := 8*ah/rhoh^(n + 1):
  lambda := (L - 1)*(L + n):
  a := (2*sqrt(1 + kappa)*(1 + sqrt(1 + kappa)))/(2*(n - 1)*(1 + sqrt(1 + kappa)) + kappa*(n - 3)):

  # Other functions appearing in the model:
  f := proc(y) option remember;
    1 - 4*(1 + sqrt(1 + kappa))*(1 - y)^(n - 1)/(2^(n + 1)*(1 + sqrt(1 + kappa*((1 - y)/2)^(n + 1)))):
  end proc:

  fp := proc(y) option remember;
    4*(1 + sqrt(1 + kappa))*(1 - y)^(n - 1)*(n - 1)/((1 - y)*2^(n + 1)*(1 + sqrt(1 + kappa*(1/2 - y/2)^(n + 1)))) - (1 + sqrt(1 + kappa))*(1 - y)^(n - 1)*kappa*(1/2 - y/2)^(n + 1)*(n + 1)/(2^(n + 1)*(1 + sqrt(1 + kappa*(1/2 - y/2)^(n + 1)))^2*sqrt(1 + kappa*(1/2 - y/2)^(n + 1))*(1/2 - y/2)):
  end proc:

  fpp := proc(y) option remember;
    -28*((-kappa^3*4^n*(-1 + y)^3*(n - 3)*(n - 5)*(1/4 - y/4)^(3*n)/14 + (-1 + y)^2*2^n*kappa^2*(n^2 - 46/7*n + 43/7)*(1/4 - y/4)^(2*n) - 4*(-1 + y)*(n^2 - 31/7*n + 22/7)*kappa*(1/4 - y/4)^n + (32*2^(-n)*(n - 1)*(-2 + n))/7)*sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n) - (3*kappa^3*4^n*(-1 + y)^3*(n - 1)*(n - 7)*(1/4 - y/4)^(3*n))/7 + (26*(-1 + y)^2*(n^2 - 74/13*n + 57/13)*2^n*kappa^2*(1/4 - y/4)^(2*n))/7 - 72*(n^2 - 37/9*n + 26/9)*(-1 + y)*kappa*(1/4 - y/4)^n/7 + (64*2^(-n)*(n - 1)*(-2 + n))/7)*(1 - y)^(n - 3)*(1 + sqrt(1 + kappa))/(sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n)*(2 + sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n))^3*(-2 + kappa*(-1 + y)*(1/2 - y/2)^n)^2):
  end proc:

  fppp := proc(y) option remember;
    -3072*(1 - y)^(n - 4)*((kappa^5*16^n*(-1 + y)^5*(n - 1)*(n^2 - 14*n + 57)*(1/4 - y/4)^(5*n)/384 - kappa^4*8^n*(-1 + y)^4*(n^3 - 18*n^2 + 65*n - 44)*(1/4 - y/4)^(4*n)/32 + 17*(n^3 - 291/17*n^2 + 823/17*n - 501/17)*kappa^3*(-1 + y)^3*4^n*(1/4 - y/4)^(3*n)/96 - (7*kappa^2*2^n*(-1 + y)^2*(n^3 - 183/14*n^2 + 212/7*n - 243/14)*(1/4 - y/4)^(2*n))/12 + (kappa*(n^2 - 8*n + 10)*(-1 + y)*(1/4 - y/4)^n - (2*2^(-n)*(n - 2)*(n - 3))/3)*(n - 1))*sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n) + ((n^3 - 33/2*n^2 + 68*n - 117/2)*kappa^5*(-1 + y)^5*16^n*(1/4 - y/4)^(5*n))/64 - 49*(n^3 - 879/49*n^2 + 2879/49*n - 279/7)*8^n*kappa^4*(-1 + y)^4*(1/4 - y/4)^(4*n)/384 + 19*(n^3 - 300/19*n^2 + 799/19*n - 482/19)*kappa^3*(-1 + y)^3*4^n*(1/4 - y/4)^(3*n)/32 - kappa^6*32^n*(-1 + y)^6*(n - 5)*(n - 7)*(n - 3)*(1/4 - y/4)^(6*n)/1536 - (13*kappa^2*2^n*(-1 + y)^2*(n^3 - 12*n^2 + 27*n - 200/13)*(1/4 - y/4)^(2*n))/8 + 7*(kappa*(-1 + y)*(n - 11/7)*(n - 6)*(1/4 - y/4)^n - (4*2^(-n)*(n - 2)*(n - 3))/7)*(n - 1)/3)*(1 + sqrt(1 + kappa))/(sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n)*(2 + sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n))^4*(-2 + kappa*(-1 + y)*(1/2 - y/2)^n)^4):
  end proc:

  KV := proc(y) option remember;
    (1 - y)^(n/2)/(2^((n - 2)/2)*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)))):
  end proc:

  qV := proc(y) option remember;
    lambda*f(y)*(1 - y)^2*(n*(n - 1)*(2*rhoh)^(n + 1) + 2*ah*(n - 3)*(1 - y)^(n + 1))/(4*(n - 1)*(n*(2*rhoh)^(n + 1) + 4*ah*(1 - y)^(n + 1))):
  end proc:

  KVp := proc(y) option remember;
    -(1 - y)^(n/2)*n/(2*(1 - y)*2^(n/2 - 1)*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)))) - (1 - y)^(n/2)*(-2*ah*(1 - y)*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)) + ah*(1 - y)^2*(-fp(y)*(n - 4) - (1 - y)*fpp(y))/((n - 1)*(n - 2)))/(2*2^(n/2 - 1)*(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)))^(3/2)):
  end proc:

  KVpp := proc(y) option remember;
    (1 - y)^(n/2)*n^2/(4*(1 - y)^2*2^(n/2 - 1)*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)))) - (1 - y)^(n/2)*n/(2*(1 - y)^2*2^(n/2 - 1)*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)))) + (1 - y)^(n/2)*n*(-2*ah*(1 - y)*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)) + ah*(1 - y)^2*(-fp(y)*(n - 4) - (1 - y)*fpp(y))/((n - 1)*(n - 2)))/(4*(1 - y)*2^(n/2 - 1)*(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)))^(3/2)) + (1 - y)^(n/2)*n*(-2*ah*(1 - y)*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)) + ah*(1 - y)^2*(-fp(y)*(n - 4) - (1 - y)*fpp(y))/((n - 1)*(n - 2)))/(4*(1 - y)*2^(n/2 - 1)*(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)))^(3/2)) - (1 - y)^(n/2)*(2*ah*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)) - 2*ah*(1 - y)*(-fp(y)*(n - 4) - (1 - y)*fpp(y))/((n - 1)*(n - 2)) - 2*ah*(1 - y)*(-fp(y)*(n - 4) - (1 - y)*fpp(y))/((n - 1)*(n - 2)) + ah*(1 - y)^2*(-fpp(y)*(n - 4) + fpp(y) - (1 - y)*fppp(y))/((n - 1)*(n - 2)))/(2*2^(n/2 - 1)*(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)))^(3/2)) + (3*(1 - y)^(n/2)*(-2*ah*(1 - y)*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)) + ah*(1 - y)^2*(-fp(y)*(n - 4) - (1 - y)*fpp(y))/((n - 1)*(n - 2)))*(-2*ah*(1 - y)*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)) + ah*(1 - y)^2*(-fp(y)*(n - 4) - (1 - y)*fpp(y))/((n - 1)*(n - 2))))/(4*2^(n/2 - 1)*(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/((n - 1)*(n - 2)))^(5/2)):
  end proc:

  UV := proc(y) option remember;
    qV(y) + (1 - y)^3*f(y)/(4*KV(y))*((1 - y)*fp(y)*KVp(y) + f(y)*((1 - y)*KVpp(y) - 2*KVp(y))):
  end proc:

  # Definition of the 2nd order ODE coefficients:
  L00 := proc(y) option remember;
    -4*UV(y)/(1 - y)^2:
  end proc:

  L01 := proc(y) option remember;
    (1 - y)*f(y)*((1 - y)*fp(y) - 2*f(y)):
  end proc:

  L02 := proc(y) option remember;
    (1 - y)^2*f(y)^2:
  end proc:

  L10 := proc(y) option remember;
    rhoh*f(y)*(1 - y)*(a*f(y)*(3 + y) - (1 + y)*(a*(1 - y)^2 - 2*(1 + y))*fp(y)/(1 - y))/(1 + y)^2:
  end proc:

  L11 := proc(y) option remember;
    2*rhoh*f(y)^2*(2*(1 + y) - a*(1 - y)^2)/(1 + y):
  end proc:

  L12 := proc(y) option remember;
    0:
  end proc:

  L20 := proc(y) option remember;
    rhoh^2*(4*(1 + y)^2 - f(y)^2*(a*(1 - y)^2 - 2*(1 + y))^2)/(-y^2 + 1)^2:
  end proc:

  L21 := proc(y) option remember;
    0:
  end proc:

  L22 := proc(y) option remember;
    0:
  end proc:

  # Precompute collocation points
  xi := [seq(evalf(cos((2.0*i-1.0)*Pi/(2.0*N))), i=1..N)]:

  # Precompute Chebyshev polynomial expressions and derivatives
  T_expr := [seq(ChebyshevT(j, x), j=0..N-1)]:
  DT_expr := [seq(diff(T_expr[j+1], x), j=0..N-1)]:
  D2T_expr := [seq(diff(T_expr[j+1], x$2), j=0..N-1)]:

  # Precompute matrices T_mat, DT_mat, D2T_mat
  T_mat := Matrix(N, N, (i,j) -> evalf(subs(x=xi[i], T_expr[j]))):
  DT_mat := Matrix(N, N, (i,j) -> evalf(subs(x=xi[i], DT_expr[j]))):
  D2T_mat := Matrix(N, N, (i,j) -> evalf(subs(x=xi[i], D2T_expr[j]))):

  # Precompute coefficient vectors
  L00_vec := Vector(N, i -> evalf(L00(xi[i]))):
  L01_vec := Vector(N, i -> evalf(L01(xi[i]))):
  L02_vec := Vector(N, i -> evalf(L02(xi[i]))):
  L10_vec := Vector(N, i -> evalf(L10(xi[i]))):
  L11_vec := Vector(N, i -> evalf(L11(xi[i]))):
  L12_vec := Vector(N, i -> evalf(L12(xi[i]))): # Will be zero
  L20_vec := Vector(N, i -> evalf(L20(xi[i]))):
  L21_vec := Vector(N, i -> evalf(L21(xi[i]))): # Will be zero
  L22_vec := Vector(N, i -> evalf(L22(xi[i]))): # Will be zero

  # Assemble matrices efficiently
  M0 := Matrix(N, N, (i,j) -> L00_vec[i]*T_mat[i,j] + L01_vec[i]*DT_mat[i,j] + L02_vec[i]*D2T_mat[i,j]):
  M1 := Matrix(N, N, (i,j) -> L10_vec[i]*T_mat[i,j] + L11_vec[i]*DT_mat[i,j] + L12_vec[i]*D2T_mat[i,j]):
  M2 := Matrix(N, N, (i,j) -> L20_vec[i]*T_mat[i,j] + L21_vec[i]*DT_mat[i,j] + L22_vec[i]*D2T_mat[i,j]):

  # We finally export the data from Maple and save in files:
  path := cat(p, "/data/"):
  nstr := convert(N, string):

  ExportMatrix(cat(path, "M0_", nstr, ".mat"), M0, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M1_", nstr, ".mat"), M1, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M2_", nstr, ".mat"), M2, target=MATLAB, mode=ascii):

end proc: