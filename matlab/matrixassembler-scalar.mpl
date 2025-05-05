# Matrices assembling for dirty black-holes' quasi-normal modes computation
# Non-extreme case
MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  N::integer, # N : number of Tchebyshev modes
  L::numeric, # L : angular momentum, L >= 0
  n::numeric, # n : number of compactified dimensions
  ah::numeric,# ah: physical string coupling parameter
  p::string   # p : string of characters containing the path where we save the assembled matrices
  )
  local f::function, g::function, rhoh::numeric, alpha::numeric, fp::function, gp::function, V::function, F::function, Fp::function, P::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string:
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
  V := y -> (1 - y)^2/16*f(y)*(n*(n - 2)*f(y) + 2*n*fp(y)*(1 - y) + 4*lambda):
  # Definition of the 2nd order ODE coefficients:
  L00 := y -> -4*V(y)/(1 - y)^2:
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
  nstr := convert(N, string);
  ExportMatrix(cat(path, "M0_", nstr, ".mat"), M0, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M1_", nstr, ".mat"), M1, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M2_", nstr, ".mat"), M2, target=MATLAB, mode=ascii):
end proc: