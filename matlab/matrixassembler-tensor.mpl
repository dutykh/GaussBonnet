# Matrices assembling for Gauss-Bonnet black holes's quasi-normal modes computation
# Non-extreme case, spin s = 2
MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  N::integer, # N : number of Tchebyshev modes
  L::numeric, # L : angular momentum, L >= 2
  n::numeric, # n : number of compactified dimensions
  ah::numeric,# ah: physical string coupling parameter
  p::string   # p : string of characters containing the path where we save the assembled matrices
  )
  local f::function, fp::function, rhoh::numeric, alpha::numeric, kappa::numeric, lambda::numeric, a:: numeric, fpp::function, fppp::function, UT::function, KT::function, KTp::function, KTpp::function, qT::function, P::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string:
  with(LinearAlgebra):
  Digits := d:
  # Determination of the event horizon:
  rhoh := fsolve(rhoh^(n - 1) + ah*rhoh^(n - 3) - 2 = 0, rhoh = 0 .. 3):
  printf("Found BH's event horizon rhoh = %f\n", rhoh):
  # Physical parameters:
  kappa := 8*ah/rhoh^(n + 1):
  lambda := L*(L + n - 1):
  a := (2*sqrt(1 + kappa)*(1 + sqrt(1 + kappa)))/(2*(n - 1)*(1 + sqrt(1 + kappa)) + kappa*(n - 3)):
  # Other parameters appearing in the model:
  f := y -> 1 - 4*(1 + sqrt(1 + kappa))*(1 - y)^(n - 1)/(2^(n + 1)*(1 + sqrt(1 + kappa*((1 - y)/2)^(n + 1)))):
  fp := y -> 4*(1 + sqrt(1 + kappa))*(1 - y)^(n - 1)*(n - 1)/((1 - y)*2^(n + 1)*(1 + sqrt(1 + kappa*(1/2 - y/2)^(n + 1)))) - (1 + sqrt(1 + kappa))*(1 - y)^(n - 1)*kappa*(1/2 - y/2)^(n + 1)*(n + 1)/(2^(n + 1)*(1 + sqrt(1 + kappa*(1/2 - y/2)^(n + 1)))^2*sqrt(1 + kappa*(1/2 - y/2)^(n + 1))*(1/2 - y/2)):
  fpp := y -> -28*((-kappa^3*4^n*(-1 + y)^3*(n - 3)*(n - 5)*(1/4 - y/4)^(3*n)/14 + (-1 + y)^2*2^n*kappa^2*(n^2 - 46/7*n + 43/7)*(1/4 - y/4)^(2*n) - 4*(-1 + y)*(n^2 - 31/7*n + 22/7)*kappa*(1/4 - y/4)^n + (32*2^(-n)*(n - 1)*(-2 + n))/7)*sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n) - (3*kappa^3*4^n*(-1 + y)^3*(n - 1)*(n - 7)*(1/4 - y/4)^(3*n))/7 + (26*(-1 + y)^2*(n^2 - 74/13*n + 57/13)*2^n*kappa^2*(1/4 - y/4)^(2*n))/7 - 72*(n^2 - 37/9*n + 26/9)*(-1 + y)*kappa*(1/4 - y/4)^n/7 + (64*2^(-n)*(n - 1)*(-2 + n))/7)*(1 - y)^(n - 3)*(1 + sqrt(1 + kappa))/(sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n)*(2 + sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n))^3*(-2 + kappa*(-1 + y)*(1/2 - y/2)^n)^2):
  fppp := y -> -3072*(1 - y)^(n - 4)*((kappa^5*16^n*(-1 + y)^5*(n - 1)*(n^2 - 14*n + 57)*(1/4 - y/4)^(5*n)/384 - kappa^4*8^n*(-1 + y)^4*(n^3 - 18*n^2 + 65*n - 44)*(1/4 - y/4)^(4*n)/32 + 17*(n^3 - 291/17*n^2 + 823/17*n - 501/17)*kappa^3*(-1 + y)^3*4^n*(1/4 - y/4)^(3*n)/96 - (7*kappa^2*2^n*(-1 + y)^2*(n^3 - 183/14*n^2 + 212/7*n - 243/14)*(1/4 - y/4)^(2*n))/12 + (kappa*(n^2 - 8*n + 10)*(-1 + y)*(1/4 - y/4)^n - (2*2^(-n)*(n - 2)*(n - 3))/3)*(n - 1))*sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n) + ((n^3 - 33/2*n^2 + 68*n - 117/2)*kappa^5*(-1 + y)^5*16^n*(1/4 - y/4)^(5*n))/64 - 49*(n^3 - 879/49*n^2 + 2879/49*n - 279/7)*8^n*kappa^4*(-1 + y)^4*(1/4 - y/4)^(4*n)/384 + 19*(n^3 - 300/19*n^2 + 799/19*n - 482/19)*kappa^3*(-1 + y)^3*4^n*(1/4 - y/4)^(3*n)/32 - kappa^6*32^n*(-1 + y)^6*(n - 5)*(n - 7)*(n - 3)*(1/4 - y/4)^(6*n)/1536 - (13*kappa^2*2^n*(-1 + y)^2*(n^3 - 12*n^2 + 27*n - 200/13)*(1/4 - y/4)^(2*n))/8 + 7*(kappa*(-1 + y)*(n - 11/7)*(n - 6)*(1/4 - y/4)^n - (4*2^(-n)*(n - 2)*(n - 3))/7)*(n - 1)/3)*(1 + sqrt(1 + kappa))/(sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n)*(2 + sqrt(4 - 2*kappa*(-1 + y)*(1/2 - y/2)^n))^4*(-2 + kappa*(-1 + y)*(1/2 - y/2)^n)^4):
  KT := y -> 2^(n/2 - 1)*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1))/(1 - y)^(n/2):
  qT := y -> lambda*(1 - y)^2*f(y)*(4*(n - 1)*(n - 2)*rhoh^2 - ah*(1 - y)^3*(fpp(y)*(1 - y) - 2*fp(y)) + ah*(n - 3)*(1 - y)^2*((n - 4)*(1 - f(y)) - 2*(1 - y)*fp(y)))/(4*(n - 2)*(4*(n - 1)*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y)))):
  KTp := y -> 2^(n/2 - 1)*(-2*ah*(1 - y)*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1) + ah*(1 - y)^2*(-(n - 3)*fp(y) + fp(y) - (1 - y)*fpp(y))/(n - 1))/(2*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1))*(1 - y)^(n/2)) + 2^(n/2 - 1)*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1))*n/(2*(1 - y)^(n/2)*(1 - y)):
  KTpp := y -> -2^(n/2 - 1)*(-2*ah*(1 - y)*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1) + ah*(1 - y)^2*(-(n - 3)*fp(y) + fp(y) - (1 - y)*fpp(y))/(n - 1))^2/(4*(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1))^(3/2)*(1 - y)^(n/2)) + 2^(n/2 - 1)*(-2*ah*(1 - y)*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1) + ah*(1 - y)^2*(-(n - 3)*fp(y) + fp(y) - (1 - y)*fpp(y))/(n - 1))*n/(2*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1))*(1 - y)^(n/2)*(1 - y)) + 2^(n/2 - 1)*(2*ah*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1) - 4*ah*(1 - y)*(-(n - 3)*fp(y) + fp(y) - (1 - y)*fpp(y))/(n - 1) + ah*(1 - y)^2*(-(n - 3)*fpp(y) + 2*fpp(y) - (1 - y)*fppp(y))/(n - 1))/(2*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1))*(1 - y)^(n/2)) + 2^(n/2 - 1)*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1))*n^2/(4*(1 - y)^(n/2)*(1 - y)^2) + 2^(n/2 - 1)*sqrt(4*rhoh^2 + ah*(1 - y)^2*((n - 3)*(1 - f(y)) - (1 - y)*fp(y))/(n - 1))*n/(2*(1 - y)^(n/2)*(1 - y)^2):
  UT := y -> qT(y) + (1 - y)^3*f(y)/(4*KT(y))*((1 - y)*fp(y)*KTp(y) + f(y)*((1 - y)*KTpp(y) - 2*KTp(y))):
  # Definition of the 2nd order ODE coefficients:
  L00 := y -> -4*UT(y)/(1 - y)^2:
  L01 := y -> (1 - y)*f(y)*((1 - y)*fp(y) - 2*f(y)):
  L02 := y -> (1 - y)^2*f(y)^2:
  L10 := y -> rhoh*f(y)*(1 - y)*(a*f(y)*(3 + y) - (1 + y)*(a*(1 - y)^2 - 2*(1 + y))*fp(y)/(1 - y))/(1 + y)^2:
  L11 := y -> 2*rhoh*f(y)^2*(2*(1 + y) - a*(1 - y)^2)/(1 + y):
  L12 := y -> 0:
  L20 := y -> rhoh^2*(4*(1 + y)^2 - f(y)^2*(a*(1 - y)^2 - 2*(1 + y))^2)/(-y^2 + 1)^2:
  L21 := y -> 0:
  L22 := y -> 0:
  P := y -> simplify(add(a[j]*ChebyshevT(j, y), j=0..N-1)):
  M0 := Matrix(N):
  M1 := Matrix(N):
  M2 := Matrix(N):
  for i from 1 to N do
    xi := cos((2.0*i-1.0)*Pi/(2.0*N)); # Chebyshev roots collocation points
    expr0 := evalf(L00(xi)*P(xi) + L01(xi)*subs(x=xi, diff(P(x),x)) + L02(xi)*subs(x=xi, diff(P(x),x$2))):
    expr1 := evalf(L10(xi)*P(xi) + L11(xi)*subs(x=xi, diff(P(x),x)) + L12(xi)*subs(x=xi, diff(P(x),x$2))):
    expr2 := evalf(L20(xi)*P(xi) + L21(xi)*subs(x=xi, diff(P(x),x)) + L22(xi)*subs(x=xi, diff(P(x),x$2))):
    for j from 1 to N do
      M0[i,j] := coeff(expr0, a[j-1]):
      M1[i,j] := coeff(expr1, a[j-1]):
      M2[i,j] := coeff(expr2, a[j-1]):
    end do:
  end do:
  # We finally export the data from Maple and save in files:
  path := cat(p, "/data/"):
  nstr := convert(N, string):
  ExportMatrix(cat(path, "M0_", nstr, ".mat"), M0, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M1_", nstr, ".mat"), M1, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M2_", nstr, ".mat"), M2, target=MATLAB, mode=ascii):
end proc: